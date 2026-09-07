{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : ConsistencySpec
-- Description : Cross-validation between independent computation routes.

module ConsistencySpec (tests) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC

import Data.Complex (Complex(..))
import qualified Data.Vector.Unboxed as V

import Types
import Ring
import Defect
import Passage
import Distribution
import Simulate
import Estimate (correctedFromPMF)

lazyTransition :: N -> K -> Pos -> Pos -> Matrix
lazyTransition n k u v = buildTransitionMatrix (qRing k) n k [EdgeAdd u v]

tests :: TestTree
tests = testGroup "Consistency"
    [ eigenvalueSanity
    , eigensystemVsDefect
    , pmfVsSimulation
    , fftVsNaive
    ]

eigenvalueSanity :: TestTree
eigenvalueSanity = testGroup "Eigenvalue Sanity"
    [ QC.testProperty "lazy: all |eigenvalue| <= 1" $
        forAll (chooseInt (6, 25)) $ \n ->
            let k = 2; u = 0; v = n `div` 2
            in ringDist n u v > k ==>
                all (\e -> abs e <= 1 + 1e-6) (exactEigenvalues $ lazyTransition n k u v)

    , QC.testProperty "lazy: exactly one eigenvalue = 1" $
        forAll (chooseInt (6, 25)) $ \n ->
            let k = 2; u = 0; v = n `div` 2
            in ringDist n u v > k ==>
                length (filter (\e -> abs (e - 1) < 1e-6) (exactEigenvalues $ lazyTransition n k u v)) == 1

    , testCase "N=36 non-lazy: no eigenvalue > 1" $ do
        let w = buildTransitionMatrix 1.0 36 2 [EdgeAdd 0 18]
            eigs = exactEigenvalues w
        assertBool "bounded" $ all (\e -> e <= 1 + 1e-6) eigs
    ]

eigensystemVsDefect :: TestTree
eigensystemVsDefect = testGroup "Eigensystem vs Defect GF"
    [ QC.testProperty "lazy: eigensystem GF ~= defect scalar GF" $
        forAll (chooseInt (8, 30)) $ \n ->
            let k = 2; q = qRing k; u = 0; v = n `div` 2; src = 1; tgt = v - 1
            in ringDist n u v > k && src /= tgt ==>
                let defs = primitivesToDefects q n k [EdgeAdd u v]
                    w = lazyTransition n k u v
                    (eigs, evecs) = exactEigensystem w
                    pi_ = stationaryDist w
                    piR = (pi_ !! src) / (pi_ !! tgt)
                    eigGF = mkEigensystemGF eigs evecs src tgt 1.0 piR
                    defGF z = firstPassageGF q n k defs src tgt z
                    z = 0.8 :+ 0.1
                    diff = abs (realPart (eigGF z) - realPart (defGF z))
                in diff < 0.01
    ]
  where realPart (r :+ _) = r

pmfVsSimulation :: TestTree
pmfVsSimulation = testGroup "PMF vs Simulation"
    [ testCase "lazy N=50: GF MFPT ~= simulation MFPT" $ do
        let n = 50; k = 2; q = qRing k; src = 1; tgt = 24
            defs = primitivesToDefects q n k [EdgeAdd 0 25]
            gf z = firstPassageGF q n k defs src tgt z
            pmf = invertPGF 1000 gf
            -- Corrected for the tail beyond the horizon, as the simulated mean
            -- now is. A raw sum over a finite horizon omits the slowest
            -- arrivals, and comparing it against a corrected mean compares two
            -- different quantities.
            gfMFPT = correctedFromPMF pmf
            sr = simulate $ defaultSimConfig
                { simN = n, simK = k, simQ = q
                , simSrc = src, simTgt = tgt
                , simDefects = [EdgeAdd 0 25]
                , simWalkers = 20000, simSeed = 42, simMaxT = 1000 }
        assertBool ("GF=" ++ show gfMFPT ++ " sim=" ++ show (srMeanFPT sr))
            (abs (gfMFPT - srMeanFPT sr) < max (5 * srStdErr sr) (0.1 * gfMFPT))
    ]

fftVsNaive :: TestTree
fftVsNaive = testGroup "FFT vs Naive Inversion"
    [ testCase "lazy N=30: FFT matches naive pointwise" $ do
        let n = 30; k = 2; q = qRing k; src = 1; tgt = 14
            defs = primitivesToDefects q n k [EdgeAdd 0 15]
            gf z = firstPassageGF q n k defs src tgt z
            fftPmf = invertPGF 200 gf
            naivePmf = invertPGFNaive 200 512 14.0 gf
        assertBool "pointwise close" $
            V.all (\(a, b) -> abs (a - b) < 1e-4) (V.zip fftPmf naivePmf)

    , testCase "normalisation preserved: FFT sum ~= 1" $ do
        let n = 20; k = 2; q = qRing k; src = 1; tgt = 9
            gf z = pureRingFPGF q n k src tgt z
            pmf = invertPGF 500 gf
            total = V.sum pmf
        assertBool ("sum=" ++ show total) $ total > 0.90 && total < 1.01
    ]
