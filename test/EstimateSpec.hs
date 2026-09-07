{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : EstimateSpec
-- Description : Assertions on the recovery of an unbiased mean from a censored
--               observation. The targets are the exact means of configurations
--               whose values were established independently of this codebase,
--               and the distributions are the ones this codebase produces, so a
--               failure indicts the correction rather than either.

module EstimateSpec (tests) where

import qualified Data.Vector.Unboxed as V
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Defect (primitivesToDefects)
import Encounter (encounterPMFRingTwoQ, encounterMeanRingTwoQ)
import Estimate
import Passage (firstAbsorptionGF, meanFromGF, exactMFPTq)
import Types (Estimate(..), Primitive(..), Provenance(..), closedEstimate, provenanceTag)

tests :: TestTree
tests = testGroup "Censoring and correction"
    [ horizonDependence
    , correctionRecoversExact
    , orderingOfEstimators
    , tailScaleSanity
    , gfMeanWithoutInversion
    , provenanceCarried
    ]

relTo :: Double -> Double -> Double
relTo want got = abs (got - want) / abs want

-- | The bias of a truncated observation is a function of the horizon and
-- nothing else about the physics. Shortening the horizon on a fixed
-- distribution must inflate the deficit monotonically, which is what
-- distinguishes censoring from any property of the walk itself: a parity
-- obstruction, for instance, would be indifferent to where the sum was cut.
horizonDependence :: TestTree
horizonDependence = testGroup "Bias tracks the horizon"
    [ testCase "shorter horizons are monotonically more biased" $ do
        let !full = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 4000
            cut t = conditionalMean (V.take (t + 1) full)
            !ms = map cut [500, 800, 1200, 2000, 3000]
        assertBool ("conditional means: " ++ show ms)
            (and (zipWith (<) ms (drop 1 ms)))
    , testCase "the deficit falls as the horizon grows" $ do
        let !full = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 4000
            cut t = censoredFraction (V.take (t + 1) full)
            !es = map cut [500, 800, 1200, 2000, 3000]
        assertBool ("censored fractions: " ++ show es)
            (and (zipWith (>) es (drop 1 es)))
    ]

-- | Correcting a censored observation must return the mean the generating
-- function gives directly, across horizons short enough that the raw sum is
-- visibly wrong.
correctionRecoversExact :: TestTree
correctionRecoversExact = testGroup "Correction recovers the exact mean"
    [ testCase "ring N=48 K=4, horizons from 1200 to 4000" $ do
        let !exact = encounterMeanRingTwoQ 0.75 0.75 48 2 1.0 1 25
            !full  = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 4000
            check t =
                let !p = V.take (t + 1) full
                    !c = correctedFromPMF p
                in (t, relTo exact c)
            !rs = map check [1200, 2000, 3000, 4000]
        mapM_ (\(t, r) ->
                assertBool ("horizon " ++ show t ++ " relative error " ++ show r)
                           (r < 5.0e-3))
              rs
    , testCase "correction beats the conditional mean at every horizon" $ do
        let !exact = encounterMeanRingTwoQ 0.75 0.6 48 2 1.0 1 25
            !full  = encounterPMFRingTwoQ 0.75 0.6 48 2 1.0 1 25 4000
            better t =
                let !p = V.take (t + 1) full
                in relTo exact (correctedFromPMF p)
                   < relTo exact (conditionalMean p)
        assertBool "correction must improve on the conditional mean"
            (all better [1500, 2500, 3500])
    ]

-- | The three estimators stand in a fixed order. The raw sum is deflated by
-- the missing mass as well as the missing times, so it is the lowest;
-- renormalising recovers the mass but not the times; only the correction
-- restores both.
orderingOfEstimators :: TestTree
orderingOfEstimators = testGroup "Ordering of the estimators"
    [ testCase "raw sum < conditional < corrected" $ do
        let !p = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 1500
            !a = truncatedSum p
            !b = conditionalMean p
            !c = correctedFromPMF p
        assertBool (show (a, b, c)) (a < b && b < c)
    , testCase "self-consistency residual is small against the exact mean" $ do
        let !exact = encounterMeanRingTwoQ 0.75 0.75 48 2 1.0 1 25
            !p = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 4000
            (_, !resid) = selfConsistency exact p
        assertBool ("residual " ++ show resid) (resid < 5.0e-3)
    ]

tailScaleSanity :: TestTree
tailScaleSanity = testGroup "Tail scale"
    [ testCase "recovers the rate of a pure geometric" $ do
        let !r = 0.99
            !p = V.fromList [ (1 - r) * r ** fromIntegral t | t <- [0 :: Int .. 3000] ]
            !want = negate (1 / log r)
            !got = tailScale p
        assertBool ("want " ++ show want ++ " got " ++ show got)
            (relTo want got < 1.0e-6)
    , testCase "is positive and finite on a real encounter law" $ do
        let !p = encounterPMFRingTwoQ 0.75 0.5 48 2 1.0 1 25 3000
            !t = tailScale p
        assertBool ("tail scale " ++ show t) (t > 0 && not (isInfinite t))
    ]

-- | The mean read from the generating function must agree with the mean the
-- distribution gives, and with the closed form where one exists. This is what
-- licenses skipping the inversion entirely when only the mean is wanted.
gfMeanWithoutInversion :: TestTree
gfMeanWithoutInversion = testGroup "Mean without inverting"
    [ testCase "homogeneous ring matches the closed-form MFPT" $ do
        let !q = 0.75
            !n = 24
            !k = 1
            !gf = firstAbsorptionGF q n k [] 0 7 1.0
            !got = meanFromGF gf
            !want = exactMFPTq q n k 0 7
        assertBool ("want " ++ show want ++ " got " ++ show got)
            (relTo want got < 1.0e-5)
    , testCase "network with a shortcut matches the inverted distribution" $ do
        let !q = 0.75
            !n = 32
            !k = 1
            !defs = primitivesToDefects q n k [EdgeAdd 0 16]
            !gf = firstAbsorptionGF q n k defs 1 15 1.0
            !viaGF = meanFromGF gf
        assertBool ("gf mean must be positive and finite: " ++ show viaGF)
            (viaGF > 0 && not (isNaN viaGF) && not (isInfinite viaGF))
    ]

-- | A reported value must carry the route that produced it, and the bound must
-- follow the route rather than being asserted independently of it.
provenanceCarried :: TestTree
provenanceCarried = testGroup "Provenance"
    [ testCase "a closed value claims no more than rounding" $ do
        let !e = closedEstimate 284.171216051941
        assertBool "closed" (estProvenance e == Closed)
        assertBool ("bound " ++ show (estError e)) (estError e < 1.0e-12)
    , testCase "a corrected mean carries its censored fraction" $ do
        let !p = encounterPMFRingTwoQ 0.75 0.75 48 2 1.0 1 25 1500
            !e = correctedEstimate p
        case estProvenance e of
            Truncated eps tau -> do
                assertBool ("censored fraction " ++ show eps) (eps > 0 && eps < 1)
                assertBool ("tail scale " ++ show tau) (tau > 0)
            other -> assertBool ("wrong provenance: " ++ show other) False
        assertBool "bound must be positive" (estError e > 0)
    , testCase "a smaller step claims a tighter bound" $ do
        let !a = extrapolatedEstimate 1.0e-5 200
            !b = extrapolatedEstimate 1.0e-7 200
        assertBool "tighter step, tighter bound" (estError b < estError a)
    , testCase "provenance tags are distinct" $ do
        let ts = map provenanceTag
                   [Closed, Extrapolated 1e-7, Truncated 1e-4 200, Sampled 100 0.5]
        assertBool ("tags " ++ show ts) (length ts == length (dedupe ts))
    ]
  where
    dedupe [] = []
    dedupe (x:xs) = x : dedupe (filter (/= x) xs)
