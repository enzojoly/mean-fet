{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : GridSpec
-- Description : Correctness tests for the square-lattice first-encounter solver.
--               The periodic torus is checked through exact symmetries, the
--               closed-form mean against the inverted PMF, and simulation. The
--               reflecting box is checked through walker-relabel symmetry, the
--               generating-function-derivative mean against a converged PMF sum,
--               simulation, PMF normalisation, and separation monotonicity.
--               The production sweep configuration (L=16 reflecting box, starts
--               (0,0)-(8,8)) is checked at three reference cells spanning the
--               high-mobility corner and the large-mean extreme. The iterative
--               pair-chain PMF is cross-checked against the GF inversion and the
--               closed-form means on both boundary conditions. Finally a set of
--               anchors certifies the second-round accelerations: the Kronecker
--               composition of the reflecting eigensystem against a dense
--               diagonalisation, the k=1 sine-cell spectrum, the gambler's-ruin
--               fundamental-matrix identity, the geometric-law inverter, and
--               seed-protocol invariance of the simulator under batching.

module GridSpec (tests) where

import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)
import Data.Complex (Complex(..))
import Data.List (sort)
import qualified Data.Vector.Unboxed as V
import qualified Numeric.LinearAlgebra as LA

import Distribution (invertPGF, invertPGFWith, nextPow2)
import Lattice
    ( encounterPMFTorus, encounterMeanTorus, simulateTorus
    , encounterPMFReflect, encounterMeanReflect, simulateReflect
    , encounterPMFReflectIter, encounterPMFTorusIter
    , reflectMatrix, reflectEigensystem, latEigQ, latIndex
    )

truncatedMean :: V.Vector Double -> Double
truncatedMean = V.sum . V.imap (\i p -> fromIntegral i * p)

tests :: TestTree
tests = testGroup "Grid (square lattice)"
    [ symmetryTests
    , consistencyTests
    , simulationTests
    , physicsTests
    , reflectingTests
    , sweepConfigTests
    , iterativePMFTests
    , anchorTests
    ]

symmetryTests :: TestTree
symmetryTests = testGroup "Torus exact symmetries"
    [ testCase "walker relabel: mean(q1,q2) = mean(q2,q1)" $ do
        let !m1 = encounterMeanTorus 8 1 0.4 0.8 1.0 (0, 0) (4, 2)
            !m2 = encounterMeanTorus 8 1 0.8 0.4 1.0 (0, 0) (4, 2)
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 1e-9 * (1 + abs m1))

    , testCase "axis swap: separation (dx,dy) = (dy,dx)" $ do
        let !m1 = encounterMeanTorus 8 1 0.75 0.5 1.0 (0, 0) (4, 2)
            !m2 = encounterMeanTorus 8 1 0.75 0.5 1.0 (0, 0) (2, 4)
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 1e-9 * (1 + abs m1))

    , testCase "start swap: mean(a,b) = mean(b,a)" $ do
        let !m1 = encounterMeanTorus 8 1 0.6 0.6 1.0 (1, 2) (5, 6)
            !m2 = encounterMeanTorus 8 1 0.6 0.6 1.0 (5, 6) (1, 2)
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 1e-9 * (1 + abs m1))
    ]

consistencyTests :: TestTree
consistencyTests = testGroup "Torus closed-form mean vs inverted PMF"
    [ testCase "homogeneous: closed-form mean matches PMF sum" $ do
        let !p  = encounterPMFTorus 6 1 0.75 0.75 1.0 (0, 0) (3, 3) 4000
            !mt = truncatedMean p
            !me = encounterMeanTorus 6 1 0.75 0.75 1.0 (0, 0) (3, 3)
        assertBool ("closed=" ++ show me ++ " pmfsum=" ++ show mt)
            (abs (me - mt) < 0.02 * me + 0.5)

    , testCase "heterogeneous: closed-form mean matches PMF sum" $ do
        let !p  = encounterPMFTorus 6 1 0.5 0.9 1.0 (0, 0) (2, 4) 4000
            !mt = truncatedMean p
            !me = encounterMeanTorus 6 1 0.5 0.9 1.0 (0, 0) (2, 4)
        assertBool ("closed=" ++ show me ++ " pmfsum=" ++ show mt)
            (abs (me - mt) < 0.02 * me + 0.5)

    , testCase "PMF normalises to one (rho=1)" $ do
        let !p    = encounterPMFTorus 6 1 0.75 0.75 1.0 (0, 0) (3, 3) 4000
            !mass = V.sum p
        assertBool ("mass=" ++ show mass) (mass > 0.98 && mass < 1.02)
    ]

simulationTests :: TestTree
simulationTests = testGroup "Torus cross-validation versus simulation"
    [ testCase "homogeneous L=6: closed-form mean within 12% of sim" $ do
        let !me = encounterMeanTorus 6 1 0.75 0.75 1.0 (0, 0) (2, 2)
            (!ms, _, _) = simulateTorus 6 1 0.75 0.75 1.0 (0, 0) (2, 2) 60000 42 4000
            !relErr = abs (me - ms) / max 1 ms
        assertBool ("closed=" ++ show me ++ " sim=" ++ show ms
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.12)

    , testCase "heterogeneous L=6: closed-form mean within 12% of sim" $ do
        let !me = encounterMeanTorus 6 1 0.5 0.9 1.0 (0, 0) (2, 2)
            (!ms, _, _) = simulateTorus 6 1 0.5 0.9 1.0 (0, 0) (2, 2) 60000 42 4000
            !relErr = abs (me - ms) / max 1 ms
        assertBool ("closed=" ++ show me ++ " sim=" ++ show ms
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.12)
    ]

physicsTests :: TestTree
physicsTests = testGroup "Torus physics"
    [ testCase "rho < 1 increases the mean" $ do
        let !mFull = encounterMeanTorus 8 1 0.75 0.75 1.0 (0, 0) (4, 4)
            !mPart = encounterMeanTorus 8 1 0.75 0.75 0.5 (0, 0) (4, 4)
        assertBool ("partial=" ++ show mPart ++ " full=" ++ show mFull)
            (mPart > mFull)

    , testCase "farther even separation gives larger mean" $ do
        let !mNear = encounterMeanTorus 12 1 0.75 0.75 1.0 (0, 0) (2, 2)
            !mFar  = encounterMeanTorus 12 1 0.75 0.75 1.0 (0, 0) (6, 6)
        assertBool ("near=" ++ show mNear ++ " far=" ++ show mFar)
            (mFar > mNear)
    ]

reflectingTests :: TestTree
reflectingTests = testGroup "Reflecting box"
    [ testCase "walker relabel: mean(q1,q2,a,b) = mean(q2,q1,b,a)" $ do
        let !m1 = encounterMeanReflect 6 1 0.4 0.8 1.0 (0, 0) (2, 1)
            !m2 = encounterMeanReflect 6 1 0.8 0.4 1.0 (2, 1) (0, 0)
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 1e-6 * (1 + abs m1))

    , testCase "start swap (homogeneous): mean(a,b) = mean(b,a)" $ do
        let !m1 = encounterMeanReflect 6 1 0.7 0.7 1.0 (0, 1) (4, 3)
            !m2 = encounterMeanReflect 6 1 0.7 0.7 1.0 (4, 3) (0, 1)
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 1e-6 * (1 + abs m1))

    , testCase "homogeneous: GF-derivative mean matches PMF sum" $ do
        let !p  = encounterPMFReflect 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !mt = truncatedMean p
            !me = encounterMeanReflect 5 1 0.7 0.7 1.0 (0, 0) (2, 2)
        assertBool ("gf=" ++ show me ++ " pmfsum=" ++ show mt)
            (abs (me - mt) < 0.03 * me + 0.5)

    , testCase "heterogeneous: GF-derivative mean matches PMF sum" $ do
        let !p  = encounterPMFReflect 5 1 0.5 0.9 1.0 (0, 0) (2, 1) 500
            !mt = truncatedMean p
            !me = encounterMeanReflect 5 1 0.5 0.9 1.0 (0, 0) (2, 1)
        assertBool ("gf=" ++ show me ++ " pmfsum=" ++ show mt)
            (abs (me - mt) < 0.03 * me + 0.5)

    , testCase "PMF normalises to one (rho=1)" $ do
        let !p    = encounterPMFReflect 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !mass = V.sum p
        assertBool ("mass=" ++ show mass) (mass > 0.98 && mass < 1.02)

    , testCase "L=5: GF-derivative mean within 12% of sim" $ do
        let !me = encounterMeanReflect 5 1 0.7 0.7 1.0 (0, 0) (2, 2)
            (!ms, _, _) = simulateReflect 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 60000 42 4000
            !relErr = abs (me - ms) / max 1 ms
        assertBool ("gf=" ++ show me ++ " sim=" ++ show ms
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.12)

    , testCase "heterogeneous L=5: GF-derivative mean within 12% of sim" $ do
        let !me = encounterMeanReflect 5 1 0.5 0.9 1.0 (0, 0) (2, 1)
            (!ms, _, _) = simulateReflect 5 1 0.5 0.9 1.0 (0, 0) (2, 1) 60000 42 4000
            !relErr = abs (me - ms) / max 1 ms
        assertBool ("gf=" ++ show me ++ " sim=" ++ show ms
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.12)

    , testCase "farther separation gives larger mean" $ do
        let !mNear = encounterMeanReflect 6 1 0.7 0.7 1.0 (0, 0) (1, 1)
            !mFar  = encounterMeanReflect 6 1 0.7 0.7 1.0 (0, 0) (5, 5)
        assertBool ("near=" ++ show mNear ++ " far=" ++ show mFar)
            (mFar > mNear)
    ]

sweepConfigTests :: TestTree
sweepConfigTests = testGroup "Sweep config (L=16 box, starts (0,0)-(8,8))"
    [ testCase ("q1=" ++ show q1 ++ " q2=" ++ show q2
                ++ " finite, positive, ~" ++ show ref) $ do
          let !m = encounterMeanReflect 16 1 q1 q2 1.0 (0, 0) (8, 8)
          assertBool ("non-finite mean: " ++ show m)
              (not (isNaN m) && not (isInfinite m))
          assertBool ("non-positive mean: " ++ show m) (m > 0)
          assertBool ("m=" ++ show m ++ " expected ~" ++ show ref ++ " (off by >3%)")
              (abs (m - ref) <= 0.03 * ref)
    | (q1, q2, ref) <- referenceCells ]
  where
    referenceCells :: [(Double, Double, Double)]
    referenceCells =
        [ (0.95, 0.95, 593.0)
        , (0.95, 0.80, 592.1)
        , (0.05, 0.20, 3361.0)
        ]

iterativePMFTests :: TestTree
iterativePMFTests = testGroup "Iterative pair-chain PMF (power iteration)"
    [ testCase "reflecting: matches GF-inversion PMF (homogeneous L=5)" $ do
        let !pg  = encounterPMFReflect     5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !pit = encounterPMFReflectIter 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !len = min (V.length pg) (V.length pit)
            !mx  = V.maximum (V.zipWith (\x y -> abs (x - y))
                              (V.take len pg) (V.take len pit))
        assertBool ("max|GF - iter| = " ++ show mx) (mx < 1e-3)

    , testCase "reflecting: matches GF-inversion PMF (heterogeneous L=5)" $ do
        let !pg  = encounterPMFReflect     5 1 0.5 0.9 1.0 (0, 0) (2, 1) 500
            !pit = encounterPMFReflectIter 5 1 0.5 0.9 1.0 (0, 0) (2, 1) 500
            !len = min (V.length pg) (V.length pit)
            !mx  = V.maximum (V.zipWith (\x y -> abs (x - y))
                              (V.take len pg) (V.take len pit))
        assertBool ("max|GF - iter| = " ++ show mx) (mx < 1e-3)

    , testCase "reflecting: iterative mean matches closed-form mean (L=5)" $ do
        let !pit = encounterPMFReflectIter 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !mt  = truncatedMean pit
            !me  = encounterMeanReflect    5 1 0.7 0.7 1.0 (0, 0) (2, 2)
        assertBool ("iter mean=" ++ show mt ++ " closed=" ++ show me)
            (abs (mt - me) < 0.02 * me + 0.5)

    , testCase "reflecting: iterative PMF normalises to one (rho=1, L=5)" $ do
        let !pit  = encounterPMFReflectIter 5 1 0.7 0.7 1.0 (0, 0) (2, 2) 500
            !mass = V.sum pit
        assertBool ("mass=" ++ show mass) (mass > 0.98 && mass < 1.02)

    , testCase "torus: iterative mean matches exact closed-form mean (L=6)" $ do
        let !pit = encounterPMFTorusIter 6 1 0.6 0.6 1.0 (0, 0) (3, 3) 1000
            !mt  = truncatedMean pit
            !me  = encounterMeanTorus    6 1 0.6 0.6 1.0 (0, 0) (3, 3)
        assertBool ("iter mean=" ++ show mt ++ " closed=" ++ show me)
            (abs (mt - me) < 0.02 * me + 0.5)

    , testCase "torus: iterative matches GF-inversion PMF (L=6)" $ do
        let !pg  = encounterPMFTorus     6 1 0.6 0.6 1.0 (0, 0) (3, 3) 1000
            !pit = encounterPMFTorusIter 6 1 0.6 0.6 1.0 (0, 0) (3, 3) 1000
            !len = min (V.length pg) (V.length pit)
            !mx  = V.maximum (V.zipWith (\x y -> abs (x - y))
                              (V.take len pg) (V.take len pit))
        assertBool ("max|GF - iter| = " ++ show mx) (mx < 1e-3)
    ]

anchorTests :: TestTree
anchorTests = testGroup "Second-round anchors"
    [ testCase "Kronecker reflecting spectrum matches dense diagonalisation (L=5,k=1)" $ do
        let l = 5; k = 1; q = 0.7
            (!vals, _) = reflectEigensystem l k q
            !dense = LA.toList (LA.eigenvaluesSH (LA.trustSym (reflectMatrix l k q)))
            !a = sort vals
            !b = sort dense
            !mx = maximum (zipWith (\x y -> abs (x - y)) a b)
        assertBool ("max eigenvalue gap = " ++ show mx) (mx < 1e-9)

    , testCase "Kronecker reflecting spectrum matches dense (heterogeneous L=6,k=2)" $ do
        let l = 6; k = 2; q = 0.55
            (!vals, _) = reflectEigensystem l k q
            !dense = LA.toList (LA.eigenvaluesSH (LA.trustSym (reflectMatrix l k q)))
            !mx = maximum (zipWith (\x y -> abs (x - y)) (sort vals) (sort dense))
        assertBool ("max eigenvalue gap = " ++ show mx) (mx < 1e-9)

    , testCase "k=1 DCT spectrum: reflectMatrix eigenvalues equal cosine composition" $ do
        let l = 8; q = 0.6
            !dense = sort (LA.toList (LA.eigenvaluesSH (LA.trustSym (reflectMatrix l 1 q))))
            !formula = sort
                [ (1 - q) + (q / 2) * (cos (pi * fromIntegral mx / fromIntegral l)
                                     + cos (pi * fromIntegral my / fromIntegral l))
                | mx <- [0 .. l - 1], my <- [0 .. l - 1] ]
            !mxdev = maximum (zipWith (\x y -> abs (x - y)) dense formula)
        assertBool ("max deviation = " ++ show mxdev) (mxdev < 1e-9)

    , testCase "latEigQ tensor factorisation matches reflect-free torus stencil" $ do
        let l = 6; k = 1; q = 0.7
            !viaLat = latEigQ q l k (2, 3)
            !direct = (1 - q) + q * 0.5
                * (cos (2 * pi * 2 / fromIntegral l) + cos (2 * pi * 3 / fromIntegral l))
        assertBool ("latEigQ=" ++ show viaLat ++ " direct=" ++ show direct)
            (abs (viaLat - direct) < 1e-12)

    , testCase "gambler's-ruin: fundamental-matrix MFPT equals (x+1)(N-x)/q (k=1)" $ do
        let n = 10; q = 0.7 :: Double
            band :: LA.Matrix Double
            band = LA.build (n, n) $ \i j ->
                let x = round i :: Int; x' = round j :: Int
                in if x == x' then 0
                   else if abs (x - x') == 1 then q / 2
                   else 0
            interior = LA.ident n - band
            ones = LA.konst 1 n :: LA.Vector Double
            solved = LA.toList (interior LA.<\> ones)
            expected = [ fromIntegral ((x + 1) * (n - x)) / q | x <- [0 .. n - 1] ]
            mx = maximum (zipWith (\a b -> abs (a - b)) solved expected)
        assertBool ("max MFPT gap = " ++ show mx) (mx < 1e-6)

    , testCase "geometric inverter: recovers p(1-p)^(t-1) to budget" $ do
        let p = 0.3 :: Double
            gf z = ((p :+ 0) * z) / (1 - ((1 - p) :+ 0) * z)
            tmax = 400
            got = invertPGFWith tmax (nextPow2 (2 * tmax + 1)) 14.0 gf
            want t = if t == 0 then 0 else p * (1 - p) ** fromIntegral (t - 1)
            mx = V.maximum (V.imap (\t g -> abs (g - want t)) got)
        assertBool ("max sup-norm = " ++ show mx) (mx < 1e-11)

    , testCase "seed-protocol invariance: torus sim mean stable across walker counts" $ do
        let (!m1, _, _) = simulateTorus 6 1 0.75 0.75 1.0 (0, 0) (2, 2) 40000 7 3000
            (!m2, _, _) = simulateTorus 6 1 0.75 0.75 1.0 (0, 0) (2, 2) 80000 7 3000
        assertBool ("m1=" ++ show m1 ++ " m2=" ++ show m2)
            (abs (m1 - m2) < 0.05 * m1 + 1.0)
    ]
