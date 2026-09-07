{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : ReferenceSpec
-- Description : Assertions against values derived outside this codebase, and
--               the structural invariants the encounter machinery must satisfy
--               whatever route computes it. Every target here was obtained by
--               an independent method (symbolic limit over exact rationals,
--               sixty-digit arithmetic, exact linear solve of the pair chain,
--               direct time iteration to negligible residual mass, and Monte
--               Carlo), so a failure indicts the implementation rather than
--               the target.

module ReferenceSpec (tests) where

import Data.Complex (Complex(..))
import qualified Data.Vector.Unboxed as V
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist)
import Distribution (Modality(..), modality, noiseFloor)
import Encounter
    ( encounterMeanRingTwoQ
    , encounterMeanTwoQ
    , encounterGFTwoQ
    , encounterPMFRingTwoQ
    )
import Lattice (encounterMeanTorus)

tests :: TestTree
tests = testGroup "Reference values and invariants"
    [ closedFormTargets
    , rhoForm
    , torusTarget
    , swapSymmetry
    , frozenPartner
    , kacInvariant
    , meanSelfConsistency
    , modalityBehaviour
    ]

within :: Double -> Double -> Double -> Bool
within tol want got
    | want == 0 = abs got <= tol
    | otherwise = abs (got - want) / abs want <= tol

report :: String -> Double -> Double -> String
report label want got =
    label ++ ": want " ++ show want ++ ", got " ++ show got
    ++ ", relative " ++ show (abs (got - want) / abs want)

-- | Mean encounter time on the homogeneous ring. Derived in closed form and
-- confirmed against a sixty-digit evaluation of the generating-function
-- derivative, an exact solve of the relative chain, and a direct iteration
-- carried to residual mass below 1e-12.
closedFormTargets :: TestTree
closedFormTargets = testGroup "Homogeneous ring mean, closed form"
    [ testCase "N=96 K=8 q=0.75/0.75 rho=1 d0=46" $ do
        let !got = encounterMeanRingTwoQ 0.75 0.75 96 4 1.0 1 47
            !want = 284.171216051941
        assertBool (report "ring K=8" want got) (within 1e-10 want got)
    , testCase "N=96 K=12 q=0.75/0.75 rho=1 d0=46" $ do
        let !got = encounterMeanRingTwoQ 0.75 0.75 96 6 1.0 1 47
            !want = 187.359286866158
        assertBool (report "ring K=12" want got) (within 1e-10 want got)
    , testCase "N=96 K=8 q=0.75/0.25 rho=1 d0=46" $ do
        let !got = encounterMeanRingTwoQ 0.75 0.25 96 4 1.0 1 47
            !want = 388.862994054417
        assertBool (report "ring asymmetric" want got) (within 1e-10 want got)
    , testCase "N=96 K=8 q=0.90/0.10 rho=1 d0=46" $ do
        let !got = encounterMeanRingTwoQ 0.90 0.10 96 4 1.0 1 47
            !want = 376.082544078114
        assertBool (report "ring strongly asymmetric" want got)
            (within 1e-10 want got)
    ]

-- | The imperfect-absorption form of the encounter generating function,
-- against a brute-force iteration of the pair chain. This is the one assertion
-- that distinguishes the first-arrival law from the encounter law: the two
-- coincide only at rho = 1, and a form that omits rho is out by tens of per
-- cent here.
rhoForm :: TestTree
rhoForm = testGroup "Imperfect absorption"
    [ testCase "ring N=7 q=0.8/0.35 rho=0.45 starts (1,5) z=0.8" $ do
        let !wA = buildTransitionMatrix 0.80 7 1 []
            !wB = buildTransitionMatrix 0.35 7 1 []
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !piV = stationaryDist wA
            !got = realPartOf
                (encounterGFTwoQ eA vA eB vB piV 0.45 (1, 5) (0.8 :+ 0))
            !want = 0.111591469348980
        assertBool (report "rho form" want got) (within 1e-10 want got)
    ]
  where
    realPartOf (x :+ _) = x

torusTarget :: TestTree
torusTarget = testGroup "Torus mean, state-count Kac term"
    [ testCase "L=4 q=0.7/0.5 rho=0.6 d0=(2,1)" $ do
        let !got = encounterMeanTorus 4 1 0.7 0.5 0.6 (0, 0) (2, 1)
            !want = 33.408115911
        assertBool (report "torus" want got) (within 1e-8 want got)
    ]

-- | Relabelling the walkers exchanges both the mobilities and the starting
-- sites. Exchanging the mobilities alone is not a symmetry once the geometry
-- distinguishes the two starts.
swapSymmetry :: TestTree
swapSymmetry = testGroup "Walker relabelling"
    [ testCase "T(qA,qB;a,b) = T(qB,qA;b,a) on the ring" $ do
        let !l = encounterMeanRingTwoQ 0.9 0.3 48 2 1.0 5 29
            !r = encounterMeanRingTwoQ 0.3 0.9 48 2 1.0 29 5
        assertBool (report "swap" l r) (within 1e-12 l r)
    , testCase "T(qA,qB;a,b) = T(qB,qA;b,a) on the torus" $ do
        let !l = encounterMeanTorus 8 1 0.8 0.25 1.0 (1, 2) (5, 6)
            !r = encounterMeanTorus 8 1 0.25 0.8 1.0 (5, 6) (1, 2)
        assertBool (report "swap torus" l r) (within 1e-12 l r)
    ]

-- | With one walker frozen the encounter reduces to first passage to a fixed
-- target, so the mean must agree with the single-walker result.
frozenPartner :: TestTree
frozenPartner = testGroup "Frozen partner reduces to first passage"
    [ testCase "qB = 0, N=96 K=8 d0=46" $ do
        let !got = encounterMeanRingTwoQ 0.75 0.0 96 4 1.0 1 47
            !want = 488.7693534925
        assertBool (report "frozen partner" want got) (within 1e-9 want got)
    ]

-- | Each failed absorption costs one mean return time, and on a uniform
-- stationary law that cost is the number of states. The difference is
-- therefore exactly |Lambda| (1 - rho) / rho, independent of mobility,
-- connectivity and separation.
kacInvariant :: TestTree
kacInvariant = testGroup "Kac return cost"
    [ testCase "ring N=96, rho=0.5" $ check 96 0.5
    , testCase "ring N=96, rho=0.2" $ check 96 0.2
    , testCase "ring N=48, rho=0.8" $ check 48 0.8
    , testCase "torus L=8, rho=0.5" $ do
        let !d = encounterMeanTorus 8 1 0.7 0.4 0.5 (0, 0) (3, 4)
                 - encounterMeanTorus 8 1 0.7 0.4 1.0 (0, 0) (3, 4)
            !want = 64 * (1 - 0.5) / 0.5
        assertBool (report "kac torus" want d) (within 1e-9 want d)
    ]
  where
    check bigN rho = do
        let !d = encounterMeanRingTwoQ 0.75 0.4 bigN 2 rho 1 (bigN `div` 2)
                 - encounterMeanRingTwoQ 0.75 0.4 bigN 2 1.0 1 (bigN `div` 2)
            !want = fromIntegral bigN * (1 - rho) / rho
        assertBool (report "kac" want d) (within 1e-9 want d)

-- | The mean taken as a limit of the generating function and the mean summed
-- from the recovered distribution describe the same random variable, so they
-- must agree once the untabulated tail is accounted for. This is the check
-- that catches an error in either the limit or the inversion without needing
-- an external reference at all.
meanSelfConsistency :: TestTree
meanSelfConsistency = testGroup "Limit against summed distribution"
    [ testCase "ring N=48 K=4 q=0.75/0.6 rho=1" $ do
        let !pmf = encounterPMFRingTwoQ 0.75 0.6 48 2 1.0 1 25 4000
            !mass = V.sum pmf
            !summed = V.sum (V.imap (\t p -> fromIntegral t * p) pmf)
            !corrected = summed / mass
            !limitMean = encounterMeanRingTwoQ 0.75 0.6 48 2 1.0 1 25
        assertBool ("mass too low to compare: " ++ show mass) (mass > 0.999)
        assertBool (report "self-consistency" limitMean corrected)
            (within 5e-3 limitMean corrected)
    ]

-- | The modality of a discrete distribution is the number of sign changes of
-- its first difference. The cases below are the ones a naive local-maximum
-- test gets wrong.
modalityBehaviour :: TestTree
modalityBehaviour = testGroup "Modality"
    [ testCase "rise then monotone decay has one mode" $ do
        let !p = V.fromList (0 : [ exp (negate (fromIntegral i) / 8) | i <- [0 :: Int .. 60] ])
            !m = modality p
        assertBool ("peaks = " ++ show (mdPeaks m)) (length (mdPeaks m) == 1)
    , testCase "flat apex is still one mode" $ do
        let raw = [ exp (negate ((fromIntegral i - 20) ** 2) / 50) | i <- [0 :: Int .. 60] ]
            !p = V.fromList (0 : take 20 raw ++ [raw !! 20] ++ drop 20 raw)
            !m = modality p
        assertBool ("peaks = " ++ show (mdPeaks m)) (length (mdPeaks m) == 1)
    , testCase "tail ripple below the floor is not a mode" $ do
        let base = [ exp (negate (fromIntegral i) / 12) | i <- [0 :: Int .. 200] ]
            bump i x = if i == 150 then x * (1 + 1e-9) else x
            !p = V.fromList (0 : zipWith bump [0 :: Int ..] base)
            !m = modality p
        assertBool ("peaks = " ++ show (mdPeaks m)) (length (mdPeaks m) == 1)
    , testCase "two separated humps give two modes and a valley" $ do
        let g c w i = exp (negate ((fromIntegral i - c) ** 2) / w)
            !p = V.fromList
                (0 : [ g 15 40 i + 0.7 * g 120 300 i | i <- [0 :: Int .. 250] ])
            !m = modality p
        assertBool ("peaks = " ++ show (mdPeaks m)) (length (mdPeaks m) == 2)
        assertBool "valley must be located" (mdValley m /= Nothing)
        assertBool ("w2 = " ++ show (mdW2 m)) (mdW2 m > 0)
    , testCase "resolution floor tracks the accuracy parameter" $ do
        let !f14 = noiseFloor 14.0 1.0
            !f20 = noiseFloor 20.0 1.0
        assertBool "floor must be positive" (f14 > 0 && f20 > 0)
        assertBool "floor must fall as accuracy rises" (f20 < f14)
    ]
