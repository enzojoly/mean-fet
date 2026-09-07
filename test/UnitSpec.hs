-- |
-- Module      : UnitSpec
-- Description : Spot checks against known answers for ring, defect, passage,
--               and simulation modules.

module UnitSpec (tests) where

import Test.Tasty
import Test.Tasty.HUnit

import Data.Complex (magnitude)
import qualified Data.Vector.Unboxed as V

import Types
import Ring
import Defect
import Passage
import Simulate

ringEigenvalue :: N -> K -> Int -> R
ringEigenvalue n k = ringEigQ (qRing k) n k

ringEigenvalues :: N -> K -> [R]
ringEigenvalues n k = [ringEigQ (qRing k) n k j | j <- [0..n-1]]

lazyRingMatrix :: N -> K -> Matrix
lazyRingMatrix n k = ringMatrix (qRing k) n k

lazyTransition :: N -> K -> Pos -> Pos -> Matrix
lazyTransition n k u v = buildTransitionMatrix (qRing k) n k [EdgeAdd u v]

tests :: TestTree
tests = testGroup "Unit"
    [ ringTests
    , matrixTests
    , transitionTests
    , eigenvalueTests
    , matrixAccessTests
    , simulationTests
    ]

ringTests :: TestTree
ringTests = testGroup "Ring"
    [ testCase "lam0 = 1 (stationary)" $
        ringEigenvalue 10 2 0 @?= 1.0

    , testCase "lam0 = 1 for various N" $ do
        ringEigenvalue 8 2 0 @?= 1.0
        ringEigenvalue 12 2 0 @?= 1.0
        ringEigenvalue 20 2 0 @?= 1.0

    , testCase "N eigenvalues returned" $
        length (ringEigenvalues 12 2) @?= 12

    , testCase "eigenvector has N components" $
        length (ringEigenvector 10 3) @?= 10

    , testCase "eigenvector normalisation" $ do
        let vec = ringEigenvector 8 0
            norm = sqrt $ sum [magnitude c ** 2 | c <- vec]
        assertBool "normalised to 1" (abs (norm - 1.0) < 1e-10)

    , testCase "K=2 lazy eigenvalue formula" $ do
        let n = 12; j = 3
            lamNL = 0.5 * (cos (2 * pi * fromIntegral j / fromIntegral n)
                         + cos (4 * pi * fromIntegral j / fromIntegral n))
            expected = 0.5 + 0.5 * lamNL
        assertBool "matches lazy formula" $
            abs (ringEigenvalue n 2 j - expected) < 1e-10

    , testCase "non-lazy eigenvalue formula" $ do
        let n = 12; j = 3
            expected = 0.5 * (cos (2 * pi * fromIntegral j / fromIntegral n)
                            + cos (4 * pi * fromIntegral j / fromIntegral n))
        assertBool "matches cosine sum" $
            abs (ringEigNonLazy n 2 j - expected) < 1e-10

    , testCase "K=3 explicit formula" $ do
        let n = 12; j = 2; selfLoop = 1 / 3 :: Double
            cosTerm = (cos (2*pi*2/12) + cos (4*pi*2/12) + cos (6*pi*2/12)) / 3
            expected = selfLoop + (1 - selfLoop) * cosTerm
        assertBool "match" $ abs (ringEigenvalue n 3 j - expected) < 1e-10
    ]

matrixTests :: TestTree
matrixTests = testGroup "Ring Matrix"
    [ testCase "size matches N" $
        matrixSize (lazyRingMatrix 12 2) @?= 12

    , testCase "diagonal = 1-q = 0.5 for K=2" $ do
        let m = lazyRingMatrix 10 2
        assertBool "W[0,0]=0.5" $ abs (matrixGet m 0 0 - 0.5) < 1e-10
        assertBool "W[5,5]=0.5" $ abs (matrixGet m 5 5 - 0.5) < 1e-10

    , testCase "neighbour weight = q/(2K) = 0.125 for K=2" $ do
        let m = lazyRingMatrix 10 2
        assertBool "W[0,1]=0.125" $ abs (matrixGet m 0 1 - 0.125) < 1e-10
        assertBool "W[0,9]=0.125" $ abs (matrixGet m 0 9 - 0.125) < 1e-10

    , testCase "non-neighbour = 0" $
        assertBool "W[0,5]=0" $ abs (matrixGet (lazyRingMatrix 10 2) 0 5) < 1e-10

    , testCase "all entries non-negative" $
        assertBool "non-neg" $ all (>= -1e-10) $ concat $ matrixToLists $ lazyRingMatrix 12 2

    , testCase "rows sum to 1" $
        assertBool "rows=1" $ all (\s -> abs (s-1) < 1e-10) $
            map sum $ matrixToLists $ lazyRingMatrix 12 2

    , testCase "columns sum to 1 (doubly stochastic)" $ do
        let m = lazyRingMatrix 12 2
        let colSums = [sum [matrixGet m i j | i <- [0..11]] | j <- [0..11]]
        assertBool "cols=1" $ all (\s -> abs (s-1) < 1e-10) colSums
    ]

transitionTests :: TestTree
transitionTests = testGroup "Transition Matrix"
    [ testCase "stochastic: rows sum to 1" $
        assertBool "rows=1" $ all (\s -> abs (s-1) < 1e-10) $
            map sum $ matrixToLists $ lazyTransition 18 2 0 9

    , testCase "non-negative" $
        assertBool "non-neg" $ all (>= -1e-10) $
            concat $ matrixToLists $ lazyTransition 18 2 0 9

    -- An added edge raises the degree of its endpoints and leaves the holding
    -- probability alone: the mobility is shared over one more destination, and
    -- the walker is no less likely to stay than it was.
    , testCase "shortcut endpoints keep their holding probability" $ do
        let w = lazyTransition 18 2 0 9
        let q = qRing 2
        assertBool "W[0,0] = 1-q" $ abs (matrixGet w 0 0 - (1 - q)) < 1e-10
        assertBool "W[9,9] = 1-q" $ abs (matrixGet w 9 9 - (1 - q)) < 1e-10

    , testCase "shortcut entry present" $
        assertBool "W[0,9]>0" $ matrixGet (lazyTransition 18 2 0 9) 0 9 > 0.01
    ]

eigenvalueTests :: TestTree
eigenvalueTests = testGroup "Eigenvalues"
    [ testCase "count = N" $
        length (exactEigenvalues $ lazyTransition 12 2 0 6) @?= 12

    , testCase "has lam=1" $
        assertBool "has 1" $ any (\e -> abs (e-1) < 1e-6) $
            exactEigenvalues $ lazyTransition 12 2 0 6

    , testCase "all |lam| <= 1" $
        assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-8) $
            exactEigenvalues $ lazyTransition 18 2 0 9

    , testCase "all finite" $
        assertBool "finite" $ all (\e -> not (isNaN e) && not (isInfinite e)) $
            exactEigenvalues $ lazyTransition 18 2 0 9

    , testCase "large N=100" $ do
        let eigs = exactEigenvalues $ lazyTransition 100 2 0 50
        assertEqual "count" 100 (length eigs)
        assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) eigs
    ]

matrixAccessTests :: TestTree
matrixAccessTests = testGroup "Matrix Access"
    [ testCase "matrixGet retrieves correct values" $ do
        let m = matrixFromLists [[1,2],[3,4]]
        matrixGet m 0 0 @?= 1
        matrixGet m 0 1 @?= 2
        matrixGet m 1 0 @?= 3
        matrixGet m 1 1 @?= 4

    , testCase "matrixSize correct" $
        matrixSize (matrixFromLists [[1,2,3],[4,5,6],[7,8,9]]) @?= 3

    , testCase "matrixToLists/matrixFromLists round-trip" $ do
        let orig = [[1,2],[3,4]]
        matrixToLists (matrixFromLists orig) @?= orig
    ]

simulationTests :: TestTree
simulationTests = testGroup "Simulation"
    [ testCase "deterministic: same seed gives same result" $ do
        let cfg = defaultSimConfig { simWalkers = 1000, simSeed = 99, simMaxT = 500 }
            r1 = simulate cfg
            r2 = simulate cfg
        assertEqual "absorbed" (srAbsorbed r1) (srAbsorbed r2)
        assertBool "mean equal" $ abs (srMeanFPT r1 - srMeanFPT r2) < 1e-12

    , testCase "all walkers accounted for" $ do
        let w = 2000
            sr = simulate $ defaultSimConfig { simWalkers = w, simSeed = 42, simMaxT = 500 }
        assertEqual "total" w (srAbsorbed sr + srSurvived sr)

    , testCase "MFPT positive when walkers absorbed" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 10, simK = 2, simQ = 0.75
                , simSrc = 0, simTgt = 5
                , simWalkers = 1000, simSeed = 7, simMaxT = 500 }
        assertBool "absorbed > 0" $ srAbsorbed sr > 0
        assertBool "mean > 0" $ srMeanFPT sr > 0

    , testCase "shortcut reduces MFPT vs pure ring" $ do
        let n = 20; k = 2; q = qRing k; src = 0; tgt = 10
            srRing = simulate $ defaultSimConfig
                { simN = n, simK = k, simQ = q
                , simSrc = src, simTgt = tgt
                , simDefects = []
                , simWalkers = 5000, simSeed = 88, simMaxT = 2000 }
            srSC = simulate $ defaultSimConfig
                { simN = n, simK = k, simQ = q
                , simSrc = src, simTgt = tgt
                , simDefects = [EdgeAdd 3 13]
                , simWalkers = 5000, simSeed = 88, simMaxT = 2000 }
        assertBool ("sc=" ++ show (srMeanFPT srSC) ++ " < ring=" ++ show (srMeanFPT srRing))
            (srMeanFPT srSC < srMeanFPT srRing)

    , testCase "encounter mode: walkers eventually meet" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 10, simK = 2, simQ = 1.0
                , simSrc = 0, simSrcB = 5
                , simDefects = []
                , simWalkers = 2000, simSeed = 77, simMaxT = 500
                , simMode = Simulate.Encounter }
        assertBool "some absorbed" $ srAbsorbed sr > 0
        assertBool "mean > 0" $ srMeanFPT sr > 0

    , testCase "encounter: same start, t>0 guard prevents t=0 absorption" $ do
        let sr = simulate $ defaultSimConfig
                { simN = 10, simK = 2, simQ = 1.0
                , simSrc = 3, simSrcB = 3
                , simDefects = []
                , simWalkers = 100, simSeed = 1, simMaxT = 500
                , simMode = Simulate.Encounter }
        assertBool "should NOT absorb at t=0 (guard t>0)" $ srAbsorbed sr > 0
    ]
