{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : CellsSpec
-- Description : All six cells against an independent solve of the pair chain.
--               The dispatcher chooses a different route for each cell --- a
--               closed form, a direct evaluation, or a contact renewal --- so
--               agreeing with one reference in one cell says nothing about the
--               others, and each is checked separately.

module CellsSpec (tests) where

import Data.List (sortBy)
import Data.Ord (comparing, Down(..))
import qualified Numeric.LinearAlgebra as LA
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Cells
import Types (Primitive(..))

tests :: TestTree
tests = testGroup "The six cells"
    [ matrixConventions
    , closedFormAvailability
    , allSixAgainstOracle
    , weightsAndMeans
    , defectsInEveryCell
    ]

-- | Sum of a column of a transition matrix, which is the total probability
-- leaving the site it indexes. Hoisted to the top level and given a signature
-- because an inferred type here would be polymorphic in the element type, and
-- the constraint that carries is not one a local binding may hold.
rowSum :: LA.Matrix Double -> Int -> Double
rowSum w i = sum [ w `LA.atIndex` (i, j) | j <- [0 .. LA.cols w - 1] ]

conserves :: LA.Matrix Double -> Bool
conserves w = all (\i -> within 1e-13 1 (rowSum w i)) [0 .. LA.rows w - 1]

within :: Double -> Double -> Double -> Bool
within tol want got
    | abs want < 1e-14 = abs got <= tol
    | otherwise        = abs (got - want) / abs want <= tol

-- | Mean time to encounter with absorption, and the probability of encountering
-- at all, obtained by solving the pair chain outright. Nothing here is shared
-- with the dispatcher: no spectrum, no generating function, no contour.
pairOracle :: LA.Matrix Double -> LA.Matrix Double -> Double -> (Int, Int)
           -> (Double, Double)
pairOracle wA wB !rho (!a0, !b0) =
    let !ns = LA.rows wA
        !n2 = ns * ns
        idx !a !b = a * ns + b
        !entries =
            [ ((idx a b, idx a' b'), pa * pb * damp a' b')
            | a <- [0 .. ns - 1], b <- [0 .. ns - 1]
            , a' <- [0 .. ns - 1], b' <- [0 .. ns - 1]
            , let !pa = wA `LA.atIndex` (a', a)
            , let !pb = wB `LA.atIndex` (b', b)
            , pa /= 0, pb /= 0 ]
        damp !a' !b' = if a' == b' then 1 - rho else 1
        !srcList =
            [ ( idx a b
              , rho * sum [ (wA `LA.atIndex` (c, a)) * (wB `LA.atIndex` (c, b))
                          | c <- [0 .. ns - 1] ] )
            | a <- [0 .. ns - 1], b <- [0 .. ns - 1] ]
        !kern = LA.accum (LA.konst 0 (n2, n2)) (+) entries
        !src  = LA.accum (LA.konst 0 n2) (+) srcList
        !sys  = LA.ident n2 - kern
        !phi  = sys LA.<\> src
        !gvec = sys LA.<\> (src + kern LA.#> phi)
        !p = phi `LA.atIndex` idx a0 b0
        !g = gvec `LA.atIndex` idx a0 b0
    in if p <= 0 then (0, 0) else (p, g / p)

-- | Reflection returns refused weight to the walker; absorption does not.
matrixConventions :: TestTree
matrixConventions = testGroup "Transition matrices"
    [ testCase "periodic and reflecting conserve, absorbing does not" $ do
        assertBool "1D periodic"    (conserves (cellMatrix D1 Periodic 12 1 0.7))
        assertBool "1D reflecting"  (conserves (cellMatrix D1 Reflecting 12 1 0.7))
        assertBool "2D periodic"    (conserves (cellMatrix D2 Periodic 5 1 0.7))
        assertBool "2D reflecting"  (conserves (cellMatrix D2 Reflecting 5 1 0.7))
        let wa = cellMatrix D1 Absorbing 12 1 0.7
        assertBool "1D absorbing must leak at the edge" (rowSum wa 0 < 1 - 1e-9)
        assertBool "1D absorbing conserves in the interior"
            (within 1e-13 1 (rowSum wa 6))
    , testCase "state count is linear in one dimension and quadratic in two" $ do
        assertBool "1D" (cellStates D1 16 == 16)
        assertBool "2D" (cellStates D2 16 == 256)
    ]

-- | The reflecting and absorbing bases diagonalise the operator at unit range
-- and not beyond, and the module reports that rather than returning a formula
-- that does not hold.
closedFormAvailability :: TestTree
closedFormAvailability = testGroup "Closed-form availability"
    [ testCase "periodic at every range" $
        assertBool "" (closedSpectrumExists Periodic 1
                       && closedSpectrumExists Periodic 4)
    -- A mirror wall is diagonalised by the cosine basis at every range, so the
    -- reflecting cell keeps its closed form throughout. The absorbing cell does
    -- not: the odd extension leaves a residue one site beyond the boundary,
    -- which a step of range two reaches and a step of unit range cannot.
    , testCase "reflecting at every range, absorbing at unit range only" $ do
        assertBool "k=1 reflecting" (closedSpectrumExists Reflecting 1)
        assertBool "k=3 reflecting" (closedSpectrumExists Reflecting 3)
        assertBool "k=1 absorbing"  (closedSpectrumExists Absorbing 1)
        assertBool "k=2 absorbing"  (not (closedSpectrumExists Absorbing 2))

    , testCase "the cosine composition reproduces the reflecting spectrum at every range" $
        mapM_ (\k -> do
            let !n = 10
                !q = 0.75
                !w = cellMatrix D1 Reflecting n k q
                !got = reverse (LA.toList (LA.eigenvaluesSH (LA.trustSym
                          (LA.scale 0.5 (w + LA.tr w)))))
                want = [ (1 - q) + q * sum [ cos (pi * fromIntegral (l * m)
                                                 / fromIntegral n)
                                           | m <- [1 .. k] ] / fromIntegral k
                       | l <- [0 .. n - 1] ]
                dev = maximum (zipWith (\a b -> abs (a - b))
                                 (sortBy (comparing Down) got)
                                 (sortBy (comparing Down) want))
            assertBool ("k=" ++ show k ++ " deviation " ++ show dev) (dev < 1e-12))
        [1, 2, 3]

    , testCase "a mirror wall does not inflate the holding probability at the edge" $ do
        let !w = cellMatrix D1 Reflecting 10 2 0.75
        assertBool ("edge holding " ++ show (w `LA.atIndex` (0, 0)))
            (w `LA.atIndex` (0, 0) < 0.5)
    ]

allSixAgainstOracle :: TestTree
allSixAgainstOracle = testGroup "Against the pair chain"
    [ testCase "1D periodic"    $ check D1 Periodic   16 (1, 0) (9, 0)  1.0
    , testCase "1D reflecting"  $ check D1 Reflecting 14 (1, 0) (9, 0)  1.0
    , testCase "1D absorbing"   $ check D1 Absorbing  12 (2, 0) (8, 0)  1.0
    , testCase "2D periodic"    $ check D2 Periodic    5 (0, 0) (2, 2)  1.0
    , testCase "2D reflecting"  $ check D2 Reflecting  5 (0, 0) (2, 2)  1.0
    , testCase "2D absorbing"   $ check D2 Absorbing   4 (0, 0) (2, 2)  1.0
    , testCase "1D periodic, imperfect absorption"
                                $ check D1 Periodic   14 (1, 0) (8, 0)  0.4
    , testCase "2D reflecting, imperfect absorption"
                                $ check D2 Reflecting  4 (0, 0) (2, 1)  0.35
    ]
  where
    check dim bc size a b rho = do
        let !q1 = 0.75
            !q2 = 0.55
            !wA = cellMatrix dim bc size 1 q1
            !wB = cellMatrix dim bc size 1 q2
            !ia = cellIndex dim size a
            !ib = cellIndex dim size b
            (!pO, !tO) = pairOracle wA wB rho (ia, ib)
            !r = encounterCell dim bc size 1 q1 q2 rho (a, b) [] 0
        assertBool ("weight: oracle " ++ show pO
                    ++ " cell " ++ show (crSplittingWeight r))
            (within 1e-7 pO (crSplittingWeight r))
        assertBool ("mean: oracle " ++ show tO
                    ++ " cell " ++ show (crMean r)
                    ++ " via " ++ crProvenance r)
            (within 1e-6 tO (crMean r))

weightsAndMeans :: TestTree
weightsAndMeans = testGroup "Reporting convention"
    [ testCase "measure-preserving cells carry unit weight" $ do
        let ws = [ crSplittingWeight (encounterCell d b s 1 0.75 0.6 1.0
                                        ((1, 0), (5, 2)) [] 0)
                 | (d, b, s) <- [ (D1, Periodic, 12), (D1, Reflecting, 12)
                                , (D2, Periodic, 4),  (D2, Reflecting, 4) ] ]
        assertBool ("weights " ++ show ws) (all (within 1e-12 1) ws)
    , testCase "absorbing cells carry a weight strictly inside the unit interval" $ do
        let ws = [ crSplittingWeight (encounterCell d Absorbing s 1 0.75 0.6 1.0
                                        ((1, 0), (5, 2)) [] 0)
                 | (d, s) <- [(D1, 12), (D2, 4)] ]
        assertBool ("weights " ++ show ws) (all (\p -> p > 0 && p < 1) ws)
    , testCase "every cell returns a positive finite mean" $ do
        let ms = [ crMean (encounterCell d b s 1 0.75 0.6 1.0 ((1, 0), (5, 2)) [] 0)
                 | d <- [D1, D2]
                 , b <- [Periodic, Reflecting, Absorbing]
                 , let s = if d == D1 then 12 else 4 ]
        assertBool ("means " ++ show ms)
            (all (\m -> m > 0 && not (isNaN m) && not (isInfinite m)) ms)
    ]

-- | A local modification is a perturbation of a homogeneous neighbourhood, so
-- the same primitive should serve whatever boundary encloses it. Adding an edge
-- must leave every cell a probability kernel and must shorten the encounter,
-- since it supplies a route that was not there before.
defectsInEveryCell :: TestTree
defectsInEveryCell = testGroup "A defect in every cell"
    -- A measure-preserving cell conserves everywhere. An absorbing cell
    -- conserves in the interior and loses exactly the attempts that leave, so
    -- its rows are bounded above by one rather than equal to it; that deficit
    -- is the boundary doing its work and not an error to be tuned away.
    [ testCase "adding an edge leaves the measure-preserving cells stochastic" $
        mapM_ (\(d, b, s, u, v) -> do
            let !w = cellMatrixWith d b s 1 0.75 [EdgeAdd u v]
                dev = maximum [ abs (rowSum w i - 1) | i <- [0 .. LA.rows w - 1] ]
            assertBool (show (d, b) ++ " worst row " ++ show dev) (dev < 1e-12))
        [ (D1, Periodic, 12, 0, 6), (D1, Reflecting, 12, 0, 6)
        , (D2, Periodic, 4, 0, 10), (D2, Reflecting, 4, 0, 10) ]

    , testCase "adding an edge leaves the absorbing cells substochastic" $
        mapM_ (\(d, b, s, u, v, interior) -> do
            let !w = cellMatrixWith d b s 1 0.75 [EdgeAdd u v]
                worst = maximum [ rowSum w i | i <- [0 .. LA.rows w - 1] ]
                least = minimum [ rowSum w i | i <- [0 .. LA.rows w - 1] ]
            assertBool (show (d, b) ++ " no row may exceed one: " ++ show worst)
                (worst < 1 + 1e-12)
            assertBool (show (d, b) ++ " some row must leak: " ++ show least)
                (least < 1 - 1e-9)
            assertBool (show (d, b) ++ " the interior must conserve")
                (abs (rowSum w interior - 1) < 1e-12))
        [ (D1, Absorbing, 12, 0, 6, 3)
        , (D2, Absorbing, 5, 0, 12, 2 + 5 * 2) ]

    , testCase "an added edge shortens the encounter" $
        mapM_ (\(d, b, s, u, v, a, c) -> do
            let bare = crMean (encounterCell d b s 1 0.75 0.6 1.0 (a, c) [] 0)
                withE = crMean (encounterCell d b s 1 0.75 0.6 1.0 (a, c)
                                    [EdgeAdd u v] 0)
            assertBool (show (d, b) ++ " bare " ++ show bare
                        ++ " with edge " ++ show withE) (withE < bare))
        [ (D1, Periodic, 16, 0, 8, (1, 0), (7, 0))
        , (D1, Reflecting, 16, 0, 8, (1, 0), (7, 0)) ]
    ]
