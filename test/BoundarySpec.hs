{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : BoundarySpec
-- Description : Extreme values and parameter-boundary stability.

module BoundarySpec (tests) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC

import qualified Data.Vector.Unboxed as V

import Types
import Ring
import Defect
import Passage
import Distribution
import Simulate

lazyTransition :: N -> K -> Pos -> Pos -> Matrix
lazyTransition n k u v = buildTransitionMatrix (qRing k) n k [EdgeAdd u v]

tests :: TestTree
tests = testGroup "Boundary"
    [ extremeValues
    , parameterBoundaries
    , simulationBoundaries
    ]

extremeValues :: TestTree
extremeValues = testGroup "Extreme Values"
    [ testCase "N=50 eigenvalues stable" $ do
        let eigs = exactEigenvalues $ lazyTransition 50 2 0 25
        assertBool "count" $ length eigs == 50
        assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) eigs

    , testCase "N=100 eigenvalues stable" $ do
        let eigs = exactEigenvalues $ lazyTransition 100 2 0 50
        assertBool "count" $ length eigs == 100
        assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) eigs

    , QC.testProperty "boundary configs are stable" $
        forAll (chooseInt (8, 40)) $ \n ->
            let k = 2; u = 0; v = n `div` 2
            in ringDist n u v > k ==>
                let eigs = exactEigenvalues $ lazyTransition n k u v
                in length eigs == n && all (\e -> abs e <= 1 + 1e-6) eigs
    ]

parameterBoundaries :: TestTree
parameterBoundaries = testGroup "Parameter Boundaries"
    [ testGroup "N"
        [ testCase "N=4 (minimum)" $
            length (exactEigenvalues $ buildTransitionMatrix (qRing 1) 4 1 []) @?= 4
        , testCase "N=6" $
            length (exactEigenvalues $ lazyTransition 6 2 0 3) @?= 6
        , testCase "N=8" $
            length (exactEigenvalues $ lazyTransition 8 2 0 4) @?= 8
        , testCase "N=50" $
            length (exactEigenvalues $ lazyTransition 50 2 0 25) @?= 50
        , testCase "N=100" $
            length (exactEigenvalues $ lazyTransition 100 2 0 50) @?= 100
        , testCase "N=200" $
            length (exactEigenvalues $ lazyTransition 200 2 0 100) @?= 200
        ]

    , testGroup "K"
        [ testCase "K=2 (standard)" $ do
            let eigs = exactEigenvalues $ lazyTransition 12 2 0 6
            assertBool "has 1" $ any (\e -> abs (e - 1) < 1e-6) eigs

        , testCase "K=3" $ do
            let eigs = exactEigenvalues $ buildTransitionMatrix (qRing 3) 18 3 [EdgeAdd 0 9]
            assertBool "has 1" $ any (\e -> abs (e - 1) < 1e-6) eigs
        ]

    , testGroup "Shortcut Distance"
        [ testCase "d=N/2 (antipodal)" $ do
            let eigs = exactEigenvalues $ lazyTransition 20 2 0 10
            assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) eigs

        , testCase "d=N/2-1 (near-anti.)" $ do
            let eigs = exactEigenvalues $ lazyTransition 20 2 0 9
            assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) eigs
        ]
    ]

simulationBoundaries :: TestTree
simulationBoundaries = testGroup "Simulation Boundaries"
    [ testCase "maxT=1: most walkers survive" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 10, simK = 2, simQ = 0.75
                , simSrc = 0, simTgt = 5, simDefects = []
                , simWalkers = 100, simSeed = 1, simMaxT = 1 }
        assertBool "survived > absorbed" $ srSurvived sr >= srAbsorbed sr

    , testCase "adjacent target: fast absorption" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 10, simK = 2, simQ = 0.75
                , simSrc = 0, simTgt = 1, simDefects = []
                , simWalkers = 1000, simSeed = 1, simMaxT = 500 }
        assertBool "most absorbed" $ srAbsorbed sr > 900

    , testCase "encounter on small ring: finite meeting time" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 8, simK = 1, simQ = 0.75
                , simSrc = 0, simSrcB = 4, simDefects = []
                , simWalkers = 100, simSeed = 42, simMaxT = 2000
                , simMode = Simulate.Encounter }
        assertBool "some absorbed" $ srAbsorbed sr > 0
    ]
