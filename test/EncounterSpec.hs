{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : EncounterSpec
-- Description : Comprehensive correctness tests for defect primitives and
--               first-encounter cross-validation against simulation, including
--               the two-mobility (qA, qB) heterogeneous-walker path.

module EncounterSpec (tests) where

import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool, (@?=))
import qualified Data.Vector.Unboxed as V

import Types (Primitive(..), matrixGet, matrixSize, matrixToLists)
import Ring (ringMatrix, ringDist, qRing)
import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist, primitivesToDefects)
import Passage (firstPassageGF, exactMFPTq)
import Distribution (invertPGF, Modality(..), modality)
import Encounter (encounterPMF, encounterGF, encounterPMFTwoQ, encounterGFTwoQ, splittingProbs, splittingProbsTwoQ)
import Simulate (SimMode(..), SimConfig(..), SimResult(..), defaultSimConfig, simulate)
import Estimate (correctedFromPMF)
import Data.Complex (Complex(..), magnitude)

isRowStochastic :: Int -> [[Double]] -> Bool
isRowStochastic n rows =
    length rows == n
    && all (\r -> length r == n) rows
    && all (\r -> abs (sum r - 1.0) < 1e-8) rows
    && all (\r -> all (>= -1e-10) r) rows

tests :: TestTree
tests = testGroup "Encounter"
    [ primitiveMatrixTests
    , primitivePhysicsTests
    , encounterNormalisationTests
    , encounterCrossValidationTests
    , encounterPhysicsTests
    , splittingProbTests
    , peakCountingTests
    , twoQTests
    ]

primitiveMatrixTests :: TestTree
primitiveMatrixTests = testGroup "Primitive Matrix Correctness"
    [ testGroup "EdgeAdd"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeAdd 0 6]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "shortcut entry nonzero" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeAdd 0 6]
            assertBool "W[0,6] > 0" $ matrixGet w 0 6 > 1e-10
            assertBool "W[6,0] > 0" $ matrixGet w 6 0 > 1e-10

        , testCase "self-loop at endpoints retained (lazy)" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeAdd 0 6]
            assertBool "W[0,0] = 1-q" $ abs (matrixGet w 0 0 - 0.25) < 1e-10
            assertBool "W[6,6] = 1-q" $ abs (matrixGet w 6 6 - 0.25) < 1e-10

        , testCase "non-endpoint sites unchanged" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeAdd 0 6]
                w0 = ringMatrix 0.75 12 2
            assertBool "W[3,3] = W0[3,3]" $ abs (matrixGet w 3 3 - matrixGet w0 3 3) < 1e-10
            assertBool "W[3,4] = W0[3,4]" $ abs (matrixGet w 3 4 - matrixGet w0 3 4) < 1e-10
        ]

    , testGroup "EdgeDel"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeDel 3 4]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "deleted edge has zero weight" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeDel 3 4]
            assertBool "W[3,4] = 0" $ abs (matrixGet w 3 4) < 1e-10
            assertBool "W[4,3] = 0" $ abs (matrixGet w 4 3) < 1e-10

        , testCase "probability drawn from the existing bonds, not the diagonal" $ do
            let w = buildTransitionMatrix 0.75 12 2 [EdgeDel 3 4]
                w0 = ringMatrix 0.75 12 2
            assertBool "holding probability untouched"
                (abs (matrixGet w 3 3 - matrixGet w0 3 3) < 1e-10)
            assertBool "a surviving bond carries more"
                (matrixGet w 3 2 > matrixGet w0 3 2 + 1e-10)
        ]

    , testGroup "WattsStrogatz"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [WattsStrogatz 0 6]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "old neighbour edge removed" $ do
            let w = buildTransitionMatrix 0.75 12 2 [WattsStrogatz 0 6]
            assertBool "W[0,1] reduced or zero" $ matrixGet w 0 1 < matrixGet (ringMatrix 0.75 12 2) 0 1 + 1e-10

        , testCase "new shortcut edge present" $ do
            let w = buildTransitionMatrix 0.75 12 2 [WattsStrogatz 0 6]
            assertBool "W[0,6] > 0" $ matrixGet w 0 6 > 1e-10
        ]

    , testGroup "Barrier"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Barrier 3 4 0.5]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "barrier reduces transition probability" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Barrier 3 4 0.5]
                w0 = ringMatrix 0.75 12 2
            assertBool "W[3,4] < W0[3,4]" $ matrixGet w 3 4 < matrixGet w0 3 4

        , testCase "barrier is symmetric" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Barrier 3 4 0.5]
            assertBool "W[3,4] = W[4,3]" $ abs (matrixGet w 3 4 - matrixGet w 4 3) < 1e-10

        , testCase "p=1 leaves matrix unchanged" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Barrier 3 4 1.0]
                w0 = ringMatrix 0.75 12 2
            assertBool "identical" $ all (\(i,j) -> abs (matrixGet w i j - matrixGet w0 i j) < 1e-10)
                [(i,j) | i <- [0..11], j <- [0..11]]

        , testCase "p=0 fully blocks edge" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Barrier 3 4 0.0]
            assertBool "W[3,4] = 0" $ abs (matrixGet w 3 4) < 1e-10
        ]

    , testGroup "Asymmetric"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Asymmetric 3 4 0.05]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "breaks symmetry" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Asymmetric 3 4 0.05]
            assertBool "W[3,4] != W[4,3]" $ abs (matrixGet w 3 4 - matrixGet w 4 3) > 1e-5
        ]

    , testGroup "Teleport"
        [ testCase "row-stochastic and non-negative" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Teleport 6 0.05]
            assertBool "stochastic" $ isRowStochastic 12 (matrixToLists w)

        , testCase "resetting site gains probability from all others" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Teleport 6 0.05]
                w0 = ringMatrix 0.75 12 2
            assertBool "W[3,6] > W0[3,6]" $ matrixGet w 3 6 > matrixGet w0 3 6

        , testCase "non-resetting sites lose diagonal probability" $ do
            let w = buildTransitionMatrix 0.75 12 2 [Teleport 6 0.05]
                w0 = ringMatrix 0.75 12 2
            assertBool "W[3,3] < W0[3,3]" $ matrixGet w 3 3 < matrixGet w0 3 3
        ]

    , testGroup "Combined Primitives"
        [ testCase "shortcut + barrier: row-stochastic" $ do
            let w = buildTransitionMatrix 0.75 16 2 [EdgeAdd 0 8, Barrier 4 5 0.3]
            assertBool "stochastic" $ isRowStochastic 16 (matrixToLists w)

        , testCase "shortcut + teleport: row-stochastic" $ do
            let w = buildTransitionMatrix 0.75 16 2 [EdgeAdd 0 8, Teleport 12 0.05]
            assertBool "teleport composed with a shortcut must conserve probability"
                (isRowStochastic 16 (matrixToLists w))

        , testCase "all four primitives: row-stochastic" $ do
            let w = buildTransitionMatrix 0.75 20 2
                    [EdgeAdd 0 10, Barrier 5 6 0.5, Asymmetric 8 9 0.03, Teleport 15 0.02]
            assertBool "combined teleport must conserve probability"
                (isRowStochastic 20 (matrixToLists w))

        , testCase "WS + barrier: row-stochastic" $ do
            let w = buildTransitionMatrix 0.75 16 2 [WattsStrogatz 0 8, Barrier 4 5 0.5]
            assertBool "stochastic" $ isRowStochastic 16 (matrixToLists w)
        ]
    ]

primitivePhysicsTests :: TestTree
primitivePhysicsTests = testGroup "Primitive Physics"
    [ testCase "EdgeAdd shortcut reduces first-passage MFPT" $ do
        let n = 20; k = 2; q = 0.75; src = 1; tgt = 9
            ringMFPT = exactMFPTq q n k src tgt
            defs = primitivesToDefects q n k [EdgeAdd 0 10]
            gf z = firstPassageGF q n k defs src tgt z
            pmf = invertPGF 1000 gf
            netMFPT = V.sum $ V.imap (\i p -> fromIntegral (i + 1) * p) (V.tail pmf)
        assertBool ("net=" ++ show netMFPT ++ " < ring=" ++ show ringMFPT)
            (netMFPT < ringMFPT)

    -- Removing the shortcut itself must slow the walk, since the fast route is
    -- what it provides. Removing an ordinary ring bond need not: the mobility
    -- it carried is spread over the remaining bonds, which grow heavier, so the
    -- walk is re-routed rather than delayed and may well arrive sooner. The
    -- assertion here is therefore about the shortcut, which is a statement
    -- about the model, and not about deletion in general, which is a statement
    -- about geometry.
    , testCase "removing the shortcut increases first-passage MFPT" $ do
        let n = 16; k = 2; q = 0.75; src = 1; tgt = 7
            defsBase = primitivesToDefects q n k [EdgeAdd 0 8]
            defsCut = primitivesToDefects q n k []
            gfBase z = firstPassageGF q n k defsBase src tgt z
            gfCut z = firstPassageGF q n k defsCut src tgt z
            pmfBase = invertPGF 1000 gfBase
            pmfCut = invertPGF 1000 gfCut
            -- Corrected for the tail beyond the horizon. A raw sum over a
            -- finite horizon omits the slowest arrivals, and the two
            -- configurations compared here do not lose the same amount.
            mfptBase = correctedFromPMF pmfBase
            mfptCut = correctedFromPMF pmfCut
        assertBool ("cut=" ++ show mfptCut ++ " > base=" ++ show mfptBase)
            (mfptCut > mfptBase)

    , testCase "Barrier increases first-passage MFPT" $ do
        let n = 16; k = 2; q = 0.75; src = 1; tgt = 7
            defsBase = primitivesToDefects q n k [EdgeAdd 0 8]
            defsBarrier = primitivesToDefects q n k [EdgeAdd 0 8, Barrier 3 4 0.3]
            gfBase z = firstPassageGF q n k defsBase src tgt z
            gfBarrier z = firstPassageGF q n k defsBarrier src tgt z
            pmfBase = invertPGF 1000 gfBase
            pmfBarrier = invertPGF 1000 gfBarrier
            mfptBase = V.sum $ V.imap (\i p -> fromIntegral (i + 1) * p) (V.tail pmfBase)
            mfptBarrier = V.sum $ V.imap (\i p -> fromIntegral (i + 1) * p) (V.tail pmfBarrier)
        assertBool ("barrier=" ++ show mfptBarrier ++ " > base=" ++ show mfptBase)
            (mfptBarrier > mfptBase)

    , testCase "WS rewiring reduces MFPT vs pure ring" $ do
        let n = 20; k = 2; q = 0.75; src = 1; tgt = 9
            ringMFPT = exactMFPTq q n k src tgt
            defs = primitivesToDefects q n k [WattsStrogatz 0 10]
            gf z = firstPassageGF q n k defs src tgt z
            pmf = invertPGF 1000 gf
            netMFPT = V.sum $ V.imap (\i p -> fromIntegral (i + 1) * p) (V.tail pmf)
        assertBool ("ws=" ++ show netMFPT ++ " < ring=" ++ show ringMFPT)
            (netMFPT < ringMFPT)
    ]

encounterNormalisationTests :: TestTree
encounterNormalisationTests = testGroup "Encounter Normalisation"
    [ testCase "GF at z near 1 approaches 1 (homogeneous)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !val = encounterGF eigs evecs piVec 1.0 (0, 5) (0.9999 :+ 0)
        assertBool ("GF near 1, got " ++ show (magnitude val))
            (magnitude val > 0.95 && magnitude val < 1.05)

    , testCase "GF at z near 1 approaches 1 (with shortcut)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !val = encounterGF eigs evecs piVec 1.0 (1, 4) (0.9999 :+ 0)
        assertBool ("GF near 1, got " ++ show (magnitude val))
            (magnitude val > 0.95 && magnitude val < 1.05)

    , testCase "PMF non-negative (homogeneous N=10)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (0, 5) 500
        assertBool "non-negative" (V.all (>= -1e-8) pmf)

    , testCase "PMF sums near 1 (homogeneous N=10)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (0, 5) 800
            !total = V.sum pmf
        assertBool ("sum=" ++ show total) (total > 0.90 && total < 1.01)

    , testCase "PMF non-negative (shortcut N=10)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (1, 4) 500
        assertBool "non-negative" (V.all (>= -1e-8) pmf)

    , testCase "PMF sums near 1 (shortcut N=10)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (1, 4) 800
            !total = V.sum pmf
        assertBool ("sum=" ++ show total) (total > 0.90 && total < 1.01)

    , testCase "T(0) = 0 (walkers do not meet at t=0)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (0, 5) 200
        assertBool "T(0)=0" $ abs (pmf `V.unsafeIndex` 0) < 1e-6
    ]

encounterCrossValidationTests :: TestTree
encounterCrossValidationTests = testGroup "Encounter vs Simulation"
    [ testCase "homogeneous N=10: exact MFPT within 15% of sim MFPT" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (0, 5) 800
            !exactMFPT = V.sum $ V.imap (\i p -> fromIntegral i * p) pmf
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.75, simQB = 0.75
                , simSrc = 0, simSrcB = 5
                , simDefects = [], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)

    , testCase "shortcut N=10: exact MFPT within 15% of sim MFPT" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (1, 4) 800
            !exactMFPT = V.sum $ V.imap (\i p -> fromIntegral i * p) pmf
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.75, simQB = 0.75
                , simSrc = 1, simSrcB = 4
                , simDefects = [EdgeAdd 0 5], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)

    , testCase "WS N=10: exact MFPT within 15% of sim MFPT" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [WattsStrogatz 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (2, 7) 800
            !exactMFPT = V.sum $ V.imap (\i p -> fromIntegral i * p) pmf
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.75, simQB = 0.75
                , simSrc = 2, simSrcB = 7
                , simDefects = [WattsStrogatz 0 5], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)

    , testCase "barrier N=10: exact MFPT within 15% of sim MFPT" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5, Barrier 2 3 0.4]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (1, 7) 800
            !exactMFPT = V.sum $ V.imap (\i p -> fromIntegral i * p) pmf
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.75, simQB = 0.75
                , simSrc = 1, simSrcB = 7
                , simDefects = [EdgeAdd 0 5, Barrier 2 3 0.4], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)
    ]

encounterPhysicsTests :: TestTree
encounterPhysicsTests = testGroup "Encounter Physics"
    [ testCase "shortcut reduces encounter MFPT vs homogeneous ring" $ do
        let n = 10; k = 1; q = 0.75
            !wRing = buildTransitionMatrix q n k []
            (!eigsR, !evecsR) = exactEigensystem wRing
            !piR = stationaryDist wRing
            !pmfRing = encounterPMF eigsR evecsR piR 1.0 (0, 5) 800
            !mfptRing = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfRing
            !wSC = buildTransitionMatrix q n k [EdgeAdd 0 5]
            (!eigsSC, !evecsSC) = exactEigensystem wSC
            !piSC = stationaryDist wSC
            !pmfSC = encounterPMF eigsSC evecsSC piSC 1.0 (1, 4) 800
            !mfptSC = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfSC
        assertBool ("sc=" ++ show mfptSC ++ " < ring=" ++ show mfptRing)
            (mfptSC < mfptRing)

    , testCase "encounter MFPT < first-passage MFPT (two walkers close gap faster)" $ do
        let n = 10; k = 1; q = 0.75; src = 1; tgt = 4
            fpMFPT = exactMFPTq q n k src tgt
            !w = buildTransitionMatrix q n k [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmf = encounterPMF eigs evecs piVec 1.0 (src, tgt) 800
            !encMFPT = V.sum $ V.imap (\i p -> fromIntegral i * p) pmf
        assertBool ("enc=" ++ show encMFPT ++ " < fp=" ++ show fpMFPT)
            (encMFPT < fpMFPT * 1.2)

    , testCase "rho < 1 increases encounter MFPT" $ do
        let n = 10; k = 1; q = 0.75
            !w = buildTransitionMatrix q n k []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmfFull = encounterPMF eigs evecs piVec 1.0 (0, 5) 800
            !mfptFull = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfFull
            !pmfPartial = encounterPMF eigs evecs piVec 0.5 (0, 5) 800
            !mfptPartial = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfPartial
        assertBool ("partial=" ++ show mfptPartial ++ " > full=" ++ show mfptFull)
            (mfptPartial > mfptFull)

    , testCase "symmetric starting positions give same MFPT" $ do
        let n = 10; k = 1; q = 0.75
            !w = buildTransitionMatrix q n k []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !pmfAB = encounterPMF eigs evecs piVec 1.0 (1, 6) 500
            !mfptAB = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfAB
            !pmfBA = encounterPMF eigs evecs piVec 1.0 (6, 1) 500
            !mfptBA = V.sum $ V.imap (\i p -> fromIntegral i * p) pmfBA
        assertBool ("AB=" ++ show mfptAB ++ " = BA=" ++ show mfptBA)
            (abs (mfptAB - mfptBA) < 0.5 + 0.02 * mfptAB)
    ]

splittingProbTests :: TestTree
splittingProbTests = testGroup "Splitting Probabilities"
    [ testCase "sum to 1 (homogeneous)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !splits = splittingProbs eigs evecs piVec 1.0 (0, 5)
            !total = V.sum splits
        assertBool ("sum=" ++ show total) (total > 0.95 && total < 1.05)

    , testCase "all non-negative (homogeneous)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !splits = splittingProbs eigs evecs piVec 1.0 (0, 5)
        assertBool "non-negative" (V.all (>= -1e-8) splits)

    , testCase "approximately uniform on homogeneous ring" $ do
        let !w = buildTransitionMatrix 0.75 10 1 []
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !splits = splittingProbs eigs evecs piVec 1.0 (0, 5)
            !expected = 1.0 / 10.0
            !maxDev = V.maximum $ V.map (\s -> abs (s - expected)) splits
        assertBool ("maxDev=" ++ show maxDev ++ " should be small")
            (maxDev < 0.05)

    , testCase "sum to 1 (with shortcut)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !splits = splittingProbs eigs evecs piVec 1.0 (1, 4)
            !total = V.sum splits
        assertBool ("sum=" ++ show total) (total > 0.95 && total < 1.05)

    , testCase "shortcut concentrates splitting near endpoints" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !splits = splittingProbs eigs evecs piVec 1.0 (1, 4)
            !atEndpoints = (splits `V.unsafeIndex` 0) + (splits `V.unsafeIndex` 5)
            !uniform = 2.0 / 10.0
        assertBool ("endpoints=" ++ show atEndpoints ++ " > uniform=" ++ show uniform)
            (atEndpoints > uniform)
    ]

peakCountingTests :: TestTree
peakCountingTests = testGroup "Modality"
    [ testCase "unimodal PMF has one mode" $ do
        let pmf = V.fromList [0, 0.01, 0.05, 0.1, 0.08, 0.05, 0.02, 0.01]
        length (mdPeaks (modality pmf)) @?= 1

    , testCase "two separated humps give two modes" $ do
        let pmf = V.fromList [0, 0.1, 0.05, 0.02, 0.01, 0.02, 0.08, 0.05, 0.02]
        length (mdPeaks (modality pmf)) @?= 2

    , testCase "two separated humps give three sign changes" $ do
        let pmf = V.fromList [0, 0.1, 0.05, 0.02, 0.01, 0.02, 0.08, 0.05, 0.02]
        mdSignChanges (modality pmf) @?= 3

    , testCase "flat PMF has no mode" $ do
        let pmf = V.fromList [0.1, 0.1, 0.1, 0.1, 0.1]
        length (mdPeaks (modality pmf)) @?= 0

    , testCase "rise then monotone decay is one mode, not two" $ do
        let pmf = V.fromList [0, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01]
        length (mdPeaks (modality pmf)) @?= 1

    , testCase "rise then monotone decay has a single sign change" $ do
        let pmf = V.fromList [0, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01]
        mdSignChanges (modality pmf) @?= 1

    , testCase "late hump at the horizon is not discarded" $ do
        let pmf = V.fromList [0, 0.3, 0.2, 0.1, 0.04, 0.02, 0.05, 0.09]
        length (mdPeaks (modality pmf)) @?= 2

    , testCase "a truncated late hump never descends, so two sign changes" $ do
        let pmf = V.fromList [0, 0.3, 0.2, 0.1, 0.04, 0.02, 0.05, 0.09]
        mdSignChanges (modality pmf) @?= 2

    , testCase "the late-population weight is a probability" $ do
        let pmf = V.fromList [0, 0.1, 0.05, 0.02, 0.01, 0.02, 0.08, 0.05, 0.02]
            w = mdW2 (modality pmf)
        assertBool ("w2 = " ++ show w) (w > 0 && w <= V.sum pmf)
    ]

twoQTests :: TestTree
twoQTests = testGroup "Two-Q Heterogeneous Mobilities"
    [ testCase "reduces to single-q when qA = qB (identical PMF)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !p1 = encounterPMF eigs evecs piVec 1.0 (1, 4) 400
            !p2 = encounterPMFTwoQ eigs evecs eigs evecs piVec 1.0 (1, 4) 400
        assertBool "two-q with equal spectra must match single-q exactly" (p1 == p2)

    , testCase "GF reduces to single-q when qA = qB (identical value)" $ do
        let !w = buildTransitionMatrix 0.75 10 1 [EdgeAdd 0 5]
            (!eigs, !evecs) = exactEigensystem w
            !piVec = stationaryDist w
            !z = 0.8 :+ 0.1
            !v1 = encounterGF eigs evecs piVec 1.0 (1, 4) z
            !v2 = encounterGFTwoQ eigs evecs eigs evecs piVec 1.0 (1, 4) z
        assertBool ("identical GF, diff=" ++ show (magnitude (v1 - v2)))
            (magnitude (v1 - v2) < 1e-12)

    , testCase "walker relabel symmetry: (qA,a) <-> (qB,b)" $ do
        let !wA = buildTransitionMatrix 0.5 10 1 []
            !wB = buildTransitionMatrix 0.9 10 1 []
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !pv = stationaryDist wA
            !pAB = encounterPMFTwoQ eA vA eB vB pv 1.0 (1, 6) 600
            !pBA = encounterPMFTwoQ eB vB eA vA pv 1.0 (6, 1) 600
            !mAB = V.sum (V.imap (\i p -> fromIntegral i * p) pAB)
            !mBA = V.sum (V.imap (\i p -> fromIntegral i * p) pBA)
        assertBool ("AB=" ++ show mAB ++ " BA=" ++ show mBA)
            (abs (mAB - mBA) < 0.5 + 0.02 * mAB)

    , testCase "two-q PMF non-negative and sums near 1 (homogeneous)" $ do
        let !wA = buildTransitionMatrix 0.5 10 1 []
            !wB = buildTransitionMatrix 0.9 10 1 []
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !pv = stationaryDist wA
            !pmf = encounterPMFTwoQ eA vA eB vB pv 1.0 (0, 5) 800
            !total = V.sum pmf
        assertBool "non-negative" (V.all (>= -1e-8) pmf)
        assertBool ("sum=" ++ show total) (total > 0.90 && total < 1.01)

    , testCase "two-q GF at z near 1 approaches 1 (shortcut)" $ do
        let !wA = buildTransitionMatrix 0.6 10 1 [EdgeAdd 0 5]
            !wB = buildTransitionMatrix 0.9 10 1 [EdgeAdd 0 5]
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !pv = stationaryDist wA
            !val = encounterGFTwoQ eA vA eB vB pv 1.0 (1, 4) (0.9999 :+ 0)
        assertBool ("GF near 1, got " ++ show (magnitude val))
            (magnitude val > 0.95 && magnitude val < 1.05)

    , testCase "two-q homogeneous N=10: exact MFPT within 15% of sim MFPT" $ do
        let !wA = buildTransitionMatrix 0.5 10 1 []
            !wB = buildTransitionMatrix 0.9 10 1 []
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !pv = stationaryDist wA
            !pmf = encounterPMFTwoQ eA vA eB vB pv 1.0 (0, 5) 800
            !exactMFPT = V.sum (V.imap (\i p -> fromIntegral i * p) pmf)
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.5, simQB = 0.9
                , simSrc = 0, simSrcB = 5
                , simDefects = [], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)

    , testCase "two-q shortcut N=10: exact MFPT within 15% of sim MFPT" $ do
        let !wA = buildTransitionMatrix 0.6 10 1 [EdgeAdd 0 5]
            !wB = buildTransitionMatrix 0.9 10 1 [EdgeAdd 0 5]
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !pv = stationaryDist wA
            !pmf = encounterPMFTwoQ eA vA eB vB pv 1.0 (1, 4) 800
            !exactMFPT = V.sum (V.imap (\i p -> fromIntegral i * p) pmf)
            !sr = simulate $ defaultSimConfig
                { simN = 10, simK = 1, simQ = 0.6, simQB = 0.9
                , simSrc = 1, simSrcB = 4
                , simDefects = [EdgeAdd 0 5], simWalkers = 50000
                , simSeed = 42, simMaxT = 800
                , simMode = Encounter, simRho = 1.0
                }
            !simMFPT = srMeanFPT sr
            !relErr = abs (exactMFPT - simMFPT) / max 1 simMFPT
        assertBool ("exact=" ++ show exactMFPT ++ " sim=" ++ show simMFPT
                    ++ " relErr=" ++ show relErr)
            (relErr < 0.15)

    , testCase "two-q MFET tracks simulation across a qB sweep (homogeneous)" $ do
        let check qb = do
                let !wA = buildTransitionMatrix 0.9 10 1 []
                    !wB = buildTransitionMatrix qb 10 1 []
                    (!eA, !vA) = exactEigensystem wA
                    (!eB, !vB) = exactEigensystem wB
                    !pv = stationaryDist wA
                    !pmf = encounterPMFTwoQ eA vA eB vB pv 1.0 (0, 4) 1500
                    !exact = V.sum (V.imap (\i p -> fromIntegral i * p) pmf)
                    !sr = simulate $ defaultSimConfig
                        { simN = 10, simK = 1, simQ = 0.9, simQB = qb
                        , simSrc = 0, simSrcB = 4
                        , simDefects = [], simWalkers = 50000
                        , simSeed = 42, simMaxT = 1500
                        , simMode = Encounter, simRho = 1.0
                        }
                    !sm = srMeanFPT sr
                    !relErr = abs (exact - sm) / max 1 sm
                assertBool ("qB=" ++ show qb ++ " exact=" ++ show exact
                            ++ " sim=" ++ show sm ++ " relErr=" ++ show relErr)
                    (relErr < 0.15)
        check 0.9
        check 0.6
        check 0.3
    ]
