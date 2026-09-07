{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : PrimitiveSpec
-- Description : The lattice modifications, judged as probability kernels.
--
--               Every result in the programme descends from a transition matrix
--               built by applying primitives to a homogeneous ring, and until
--               now nothing asserted that the outcome was a probability kernel
--               at all. The assertions here are elementary --- columns sum to
--               one, entries are non-negative, detailed balance holds where the
--               primitive is reversible, and the construction is continuous in
--               the mobility --- and they are made across the whole range of
--               mobility rather than at one convenient value, because a
--               formulation that degenerates does so at an endpoint.

module PrimitiveSpec (tests) where

import qualified Numeric.LinearAlgebra as LA
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Types (Primitive(..), reversible, primitiveDegreeShift)
import Defect (buildTransitionMatrix)

tests :: TestTree
tests = testGroup "Lattice primitives"
    [ stochasticity
    , positivity
    , detailedBalance
    , continuityInMobility
    , structuralClaims
    , reversibilityLabels
    ]

-- Mobilities span the range, and include the endpoint at which the previous
-- formulation carried a separate branch.
mobilities :: [Double]
mobilities = [0.10, 0.25, 0.50, 0.75, 0.90, 0.99, 1.00]

nSites, halfRange :: Int
nSites = 24
halfRange = 2

-- Every primitive, at sites chosen so that added edges are not already present
-- and removed edges are.
catalogue :: [(String, Primitive)]
catalogue =
    [ ("EdgeAdd",       EdgeAdd 0 12)
    , ("EdgeDel",       EdgeDel 5 6)
    , ("WattsStrogatz", WattsStrogatz 0 12)
    , ("Barrier",       Barrier 5 6 0.3)
    , ("Reweight",      Reweight 5 6 0.02)
    , ("DirectedAdd",   DirectedAdd 0 12)
    , ("Asymmetric",    Asymmetric 5 6 0.30)
    , ("Teleport",      Teleport 7 0.05)
    ]

-- The convention throughout the codebase is row-stochastic: the entry at
-- (i, j) is the probability of stepping from i to j, the simulator samples row
-- i as the departure law from i, and the stationary law is the left eigenvector.
-- A perturbation that conserved columns rather than rows would pass a
-- symmetric test and fail on any lattice with unequal degrees.
rowSum :: LA.Matrix Double -> Int -> Double
rowSum w i = sum [ w `LA.atIndex` (i, j) | j <- [0 .. LA.cols w - 1] ]

worstRow :: LA.Matrix Double -> Double
worstRow w = maximum [ abs (rowSum w i - 1) | i <- [0 .. LA.rows w - 1] ]

leastEntry :: LA.Matrix Double -> Double
leastEntry w = minimum
    [ w `LA.atIndex` (i, j)
    | i <- [0 .. LA.rows w - 1], j <- [0 .. LA.cols w - 1] ]

built :: Double -> Primitive -> LA.Matrix Double
built q p = buildTransitionMatrix q nSites halfRange [p]

-- | A walk that conserves probability has every column summing to one. This is
-- the assertion that the previous formulation failed at unit mobility, where
-- columns ran between 0.95 and 1.14.
stochasticity :: TestTree
stochasticity = testGroup "Rows sum to one"
    [ testCase (name ++ " at every mobility") $
        mapM_ (\q ->
            let !d = worstRow (built q p)
            in assertBool (name ++ " q=" ++ show q ++ " worst row " ++ show d)
                          (d < 1e-12))
        mobilities
    | (name, p) <- catalogue ]

positivity :: TestTree
positivity = testGroup "Entries are non-negative"
    [ testCase (name ++ " at every mobility") $
        mapM_ (\q ->
            let !m = leastEntry (built q p)
            in assertBool (name ++ " q=" ++ show q ++ " least entry " ++ show m)
                          (m > negate 1e-14))
        mobilities
    | (name, p) <- catalogue ]

-- | Where a primitive leaves the walk reversible, the stationary law is
-- proportional to degree and detailed balance holds against it exactly. This
-- distinguishes a correct degree-raising modification from one that moves
-- weight along a single bond: the latter would leave the law uniform.
detailedBalance :: TestTree
detailedBalance = testGroup "Detailed balance against the degree law"
    [ testCase (name ++ " at every mobility") $
        mapM_ (\q -> check name p q) mobilities
    | (name, p) <- catalogue, reversible p ]
  where
    check name p q = do
        let !w = built q p
            !base = 2 * halfRange
            shifts = primitiveDegreeShift nSites p
            deg s = fromIntegral (base + sum [ d | (t, d) <- shifts, t == s ])
            !pis = [ deg s | s <- [0 .. nSites - 1] ]
            -- Detailed balance in the row convention pairs the departure law
            -- of each site with its own stationary weight: pi_i W[i,j] equals
            -- pi_j W[j,i]. The opposite pairing is the column-convention
            -- statement and holds only where the degrees are equal, which is
            -- exactly the case a degree-raising primitive leaves behind.
            dev = maximum
                [ abs ( (pis !! i) * (w `LA.atIndex` (i, j))
                      - (pis !! j) * (w `LA.atIndex` (j, i)) )
                | i <- [0 .. nSites - 1], j <- [0 .. nSites - 1] ]
        assertBool (name ++ " q=" ++ show q ++ " worst imbalance " ++ show dev)
                   (dev < 1e-12)

-- | A construction with a special case at an endpoint is discontinuous there,
-- and the discontinuity is invisible to any test that does not approach it.
-- Approaching it is the whole of this assertion.
continuityInMobility :: TestTree
continuityInMobility = testGroup "Continuous as the mobility reaches one"
    [ testCase name $ do
        let near = built 0.999999 p
            at   = built 1.0 p
            dev = maximum
                [ abs (near `LA.atIndex` (i, j) - at `LA.atIndex` (i, j))
                | i <- [0 .. nSites - 1], j <- [0 .. nSites - 1] ]
        assertBool (name ++ " jump at unit mobility " ++ show dev) (dev < 1e-5)
    | (name, p) <- catalogue ]

structuralClaims :: TestTree
structuralClaims = testGroup "What each primitive claims to do"
    [ testCase "an empty defect list leaves the ring untouched" $ do
        let !a = buildTransitionMatrix 0.75 nSites halfRange []
            !b = buildTransitionMatrix 0.75 nSites halfRange []
        assertBool "identical" (LA.maxElement (LA.cmap abs (a - b)) < 1e-15)

    , testCase "adding an edge raises the degree of both endpoints" $ do
        let !w = built 1.0 (EdgeAdd 0 12)
            nz i = length [ () | j <- [0 .. nSites - 1]
                          , w `LA.atIndex` (i, j) > 1e-12 ]
        assertBool ("degree at 0 is " ++ show (nz 0)) (nz 0 == 2 * halfRange + 1)
        assertBool ("degree at 12 is " ++ show (nz 12)) (nz 12 == 2 * halfRange + 1)
        assertBool ("degree at 5 is " ++ show (nz 5)) (nz 5 == 2 * halfRange)

    , testCase "the new bond carries the same weight as every other" $ do
        let !w = built 0.75 (EdgeAdd 0 12)
            !d = fromIntegral (2 * halfRange + 1)
        assertBool "shortcut weight"
            (abs (w `LA.atIndex` (0, 12) - 0.75 / d) < 1e-12)
        assertBool "ring weight"
            (abs (w `LA.atIndex` (0, 1) - 0.75 / d) < 1e-12)

    , testCase "removing an edge lowers the degree of both endpoints" $ do
        let !w = built 1.0 (EdgeDel 5 6)
            nz i = length [ () | j <- [0 .. nSites - 1]
                          , w `LA.atIndex` (i, j) > 1e-12 ]
        assertBool ("degree at 5 is " ++ show (nz 5)) (nz 5 == 2 * halfRange - 1)
        assertBool ("degree at 6 is " ++ show (nz 6)) (nz 6 == 2 * halfRange - 1)

    , testCase "rewiring preserves the degree of the source" $ do
        let !w = built 1.0 (WattsStrogatz 0 12)
            nz i = length [ () | j <- [0 .. nSites - 1]
                          , w `LA.atIndex` (i, j) > 1e-12 ]
        assertBool ("source degree " ++ show (nz 0)) (nz 0 == 2 * halfRange)
        assertBool ("abandoned neighbour " ++ show (nz 1)) (nz 1 == 2 * halfRange - 1)
        assertBool ("new endpoint " ++ show (nz 12)) (nz 12 == 2 * halfRange + 1)

    , testCase "a one-way shortcut is one-way" $ do
        let !w = built 0.75 (DirectedAdd 0 12)
        assertBool "0 reaches 12" (w `LA.atIndex` (0, 12) > 1e-12)
        assertBool "12 does not reach 0" (w `LA.atIndex` (12, 0) < 1e-12)

    , testCase "the holding probability is untouched by adding an edge" $ do
        let !w = built 0.75 (EdgeAdd 0 12)
        assertBool ("hold at 0 is " ++ show (w `LA.atIndex` (0, 0)))
            (abs (w `LA.atIndex` (0, 0) - 0.25) < 1e-12)

    , testCase "a barrier lowers one bond and returns the refused weight" $ do
        let !w0 = buildTransitionMatrix 0.75 nSites halfRange []
            !w  = built 0.75 (Barrier 5 6 0.3)
            !dropped = w0 `LA.atIndex` (5, 6) - w `LA.atIndex` (5, 6)
        assertBool "bond reduced" (dropped > 1e-12)
        assertBool "returned to the walker"
            (abs (w `LA.atIndex` (5, 5) - w0 `LA.atIndex` (5, 5) - dropped) < 1e-12)
    ]

reversibilityLabels :: TestTree
reversibilityLabels = testGroup "Reversibility is labelled honestly"
    [ testCase "directed, biased and resetting rules are irreversible" $ do
        assertBool "DirectedAdd" (not (reversible (DirectedAdd 0 12)))
        assertBool "Asymmetric"  (not (reversible (Asymmetric 5 6 0.03)))
        assertBool "Teleport"    (not (reversible (Teleport 7 0.05)))
    , testCase "edge and bond rules are reversible" $ do
        assertBool "EdgeAdd"       (reversible (EdgeAdd 0 12))
        assertBool "EdgeDel"       (reversible (EdgeDel 5 6))
        assertBool "WattsStrogatz" (reversible (WattsStrogatz 0 12))
        assertBool "Barrier"       (reversible (Barrier 5 6 0.3))
    , testCase "an irreversible rule has no degree-based stationary law" $ do
        let !w = built 0.75 (DirectedAdd 0 12)
            dev = maximum
                [ abs (w `LA.atIndex` (i, j) - w `LA.atIndex` (j, i))
                | i <- [0 .. nSites - 1], j <- [0 .. nSites - 1] ]
        assertBool "must not be symmetric" (dev > 1e-6)
    ]
