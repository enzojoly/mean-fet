{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : AbsorbingSpec
-- Description : The absorbing cells judged against a pair chain that carries
--               the loss explicitly. Where the domain destroys walkers the
--               encounter law no longer sums to one, so there are two numbers
--               to establish rather than one, and both are obtained here by
--               solving the joint system outright.

module AbsorbingSpec (tests) where

import Data.List (sortBy)
import Data.Ord (comparing, Down(..))
import qualified Numeric.LinearAlgebra as LA
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Absorbing

tests :: TestTree
tests = testGroup "Absorbing boundaries"
    [ spectrumOneDimension
    , massIsLost
    , againstPairOracle
    , twoDimensions
    , schemaPair
    ]

within :: Double -> Double -> Double -> Bool
within tol want got
    | abs want < 1e-14 = abs got <= tol
    | otherwise        = abs (got - want) / abs want <= tol

-- | Splitting weight and conditional mean from the joint chain, with the loss
-- carried explicitly. Writing p for the probability of meeting before either
-- walker is lost, one step either produces the encounter outright or leaves a
-- surviving pair from which the same question is asked, so p solves a linear
-- system whose kernel is the joint step damped on the diagonal. The first
-- moment restricted to the meeting event solves the same system against a
-- source that carries the accumulated survival, and their ratio is the
-- conditional mean.
pairAbsorbing :: Int -> LA.Matrix Double -> LA.Matrix Double -> Double
              -> (Int, Int) -> (Double, Double)
pairAbsorbing !ns wA wB !rho (!a0, !b0) =
    let !n2 = ns * ns
        idx !a !b = a * ns + b

        !entries =
            [ ((idx a b, idx a' b'), pa * pb * damp a' b')
            | a <- [0 .. ns - 1], b <- [0 .. ns - 1]
            , a' <- [0 .. ns - 1], b' <- [0 .. ns - 1]
            , let !pa = wA `LA.atIndex` (a', a)
            , let !pb = wB `LA.atIndex` (b', b)
            , pa /= 0, pb /= 0
            ]
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

-- | At unit range the killed operator is diagonalised by the sine basis. At
-- greater range it is not, and the module reports that rather than returning a
-- formula that does not hold.
spectrumOneDimension :: TestTree
spectrumOneDimension = testGroup "Killed spectrum"
    [ testCase "k=1 closed form matches the killed matrix" $ do
        let !n = 12
            !q = 0.7
            !w = killedRing1D n 1 q
            !dense = LA.toList (LA.eigenvaluesSH (LA.trustSym w))
            !formula = maybe [] id (absorbingSpectrum1D n 1 q)
            desc = sortBy (comparing Down)
            !dev = maximum (zipWith (\a b -> abs (a - b))
                              (desc dense) (desc formula))
        assertBool ("max deviation " ++ show dev) (dev < 1e-12)
    , testCase "k=2 admits no closed form and says so" $
        assertBool "expected Nothing" (absorbingSpectrum1D 12 2 0.7 == Nothing)
    ]


-- | The domain destroys walkers, so the columns of the transition matrix do
-- not sum to one and the deficit is exactly the loss.
massIsLost :: TestTree
massIsLost = testGroup "Substochasticity"
    [ testCase "interior columns conserve, boundary columns leak" $ do
        let !n = 10
            !q = 0.8
            !w = killedRing1D n 1 q
            colSum j = sum [ w `LA.atIndex` (i, j) | i <- [0 .. n - 1] ]
        assertBool "interior must conserve" (within 1e-14 1 (colSum 5))
        assertBool ("left edge must leak: " ++ show (colSum 0))
            (colSum 0 < 1 - 1e-9)
        assertBool ("right edge must leak: " ++ show (colSum (n - 1)))
            (colSum (n - 1) < 1 - 1e-9)
    , testCase "a box leaks on every edge" $ do
        let !l = 5
            !w = killedBox2D l 1 0.75
            colSum j = sum [ w `LA.atIndex` (i, j) | i <- [0 .. l * l - 1] ]
        assertBool "centre conserves" (within 1e-14 1 (colSum (2 + l * 2)))
        assertBool "corner leaks twice over" (colSum 0 < 1 - 0.3)
    ]

againstPairOracle :: TestTree
againstPairOracle = testGroup "One dimension against the pair chain"
    [ testCase "N=10 q=0.75/0.6 rho=1" $ check 10 0.75 0.60 1.0 (2, 7)
    , testCase "N=8 q=0.9/0.9 rho=0.5"  $ check 8  0.90 0.90 0.5 (1, 6)
    , testCase "N=9 q=0.5/0.85 rho=1"   $ check 9  0.50 0.85 1.0 (0, 5)
    ]
  where
    check n q1 q2 rho starts = do
        let !wA = killedRing1D n 1 q1
            !wB = killedRing1D n 1 q2
            (!pO, !tO) = pairAbsorbing n wA wB rho starts
            (!pE, !tE) = encounterAbsorbing wA wB rho starts
        assertBool ("splitting weight: oracle " ++ show pO ++ " engine " ++ show pE)
            (within 1e-8 pO pE)
        assertBool ("conditional mean: oracle " ++ show tO ++ " engine " ++ show tE)
            (within 1e-6 tO tE)

twoDimensions :: TestTree
twoDimensions = testGroup "Two dimensions against the pair chain"
    [ testCase "L=4 q=0.75/0.75 rho=1" $ do
        let !l = 4
            !w = killedBox2D l 1 0.75
            (!pO, !tO) = pairAbsorbing (l * l) w w 1.0 (0, 2 + l * 2)
            (!pE, !tE) = encounterAbsorbing w w 1.0 (0, 2 + l * 2)
        assertBool ("weight: oracle " ++ show pO ++ " engine " ++ show pE)
            (within 1e-8 pO pE)
        assertBool ("mean: oracle " ++ show tO ++ " engine " ++ show tE)
            (within 1e-6 tO tE)
    , testCase "L=4 heterogeneous mobilities" $ do
        let !l = 4
            !wA = killedBox2D l 1 0.9
            !wB = killedBox2D l 1 0.4
            (!pO, !tO) = pairAbsorbing (l * l) wA wB 1.0 (0, 3 + l * 1)
            (!pE, !tE) = encounterAbsorbing wA wB 1.0 (0, 3 + l * 1)
        assertBool ("weight: oracle " ++ show pO ++ " engine " ++ show pE)
            (within 1e-8 pO pE)
        assertBool ("mean: oracle " ++ show tO ++ " engine " ++ show tE)
            (within 1e-6 tO tE)
    , testCase "the splitting weight is a probability strictly below one" $ do
        let !w = killedBox2D 4 1 0.75
            (!p, _) = encounterAbsorbing w w 1.0 (0, 2 + 4 * 2)
        assertBool ("weight " ++ show p) (p > 0 && p < 1)
    ]

-- | The two numbers a domain admits are reported together, and a
-- measure-preserving domain is the weight-one member of the same family rather
-- than a separate case. A consumer that reads the mean therefore reads the
-- conditional mean everywhere, and obtains the unconditional one exactly when
-- the weight is one, without having to know which kind of domain produced it.
schemaPair :: TestTree
schemaPair = testGroup "Reporting convention"
    [ testCase "an absorbing domain carries a weight strictly inside (0,1)" $ do
        let !w = killedRing1D 10 1 0.75
            (!p, !t) = encounterAbsorbing w w 1.0 (2, 7)
        assertBool ("weight " ++ show p) (p > 0 && p < 1)
        assertBool ("conditional mean " ++ show t) (t > 0)
    , testCase "raising the loss lowers the weight and shortens the conditional mean" $ do
        let (!pWide, !tWide) =
                let !w = killedRing1D 16 1 0.6 in encounterAbsorbing w w 1.0 (3, 12)
            (!pNarrow, !tNarrow) =
                let !w = killedRing1D 8 1 0.6 in encounterAbsorbing w w 1.0 (1, 6)
        assertBool ("weights " ++ show (pWide, pNarrow)) (pNarrow < pWide)
        assertBool ("means " ++ show (tWide, tNarrow)) (tNarrow < tWide)
    , testCase "the conditional mean is finite where the unconditional is not" $ do
        let !w = killedRing1D 12 1 0.5
            (!p, !t) = encounterAbsorbing w w 1.0 (1, 10)
        assertBool "loss must be genuine" (p < 1)
        assertBool ("conditional mean must be finite: " ++ show t)
            (not (isNaN t) && not (isInfinite t) && t > 0)
    ]
