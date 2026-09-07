{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : RegressionSpec
-- Description : Pinned reference values that must not drift across refactors.

module RegressionSpec (tests) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC

import Data.List (sort)
import qualified Data.Vector.Unboxed as V

import Types
import Ring
import Defect
import Passage
import Distribution

ringEigenvalue :: N -> K -> Int -> R
ringEigenvalue n k = ringEigQ (qRing k) n k

lazyTransition :: N -> K -> Pos -> Pos -> Matrix
lazyTransition n k u v = buildTransitionMatrix (qRing k) n k [EdgeAdd u v]

shortestPath :: N -> K -> Pos -> Pos -> Pos -> Pos -> Int
shortestPath bigN k src tgt scA scB = minimum [direct, viaAB, viaBA]
  where
    rd a b = ringDist bigN a b
    stepsOn d = (d + k - 1) `div` k
    direct = stepsOn (rd src tgt)
    viaAB  = stepsOn (rd src scA) + 1 + stepsOn (rd scB tgt)
    viaBA  = stepsOn (rd src scB) + 1 + stepsOn (rd scA tgt)

computePMF :: N -> K -> Pos -> Pos -> Pos -> Pos -> V.Vector Double
computePMF n k scA scB src tgt =
    let !q = qRing k
        !defs = primitivesToDefects q n k [EdgeAdd scA scB]
        gf z = firstPassageGF q n k defs src tgt z
        !sp = shortestPath n k src tgt scA scB
        !pmf = invertPGF 1500 gf
    in V.imap (\i x -> if i + 1 < sp then 0 else max 0 x) (V.tail pmf)

peakCount :: N -> K -> Pos -> Pos -> Pos -> Pos -> Int
peakCount n k scA scB src tgt = length $ mdPeaks (modality (computePMF n k scA scB src tgt))

tests :: TestTree
tests = testGroup "Regression"
    [ classificationTests
    , determinismTests
    , eigenvalueTests
    , matrixEntryTests
    , marris2023Tests
    , marris2025Tests
    , eigenvalueFormulaTests
    ]

classificationTests :: TestTree
classificationTests = testGroup "Classification"
    -- Two arrival populations are resolved only where the two routes carry
    -- separated timescales, which needs the circumference to be long against
    -- the local reach. These configurations were chosen against a rule under
    -- which the shortcut carried the whole holding probability and so was
    -- several times heavier than a ring bond; a genuine edge carries exactly
    -- what a ring bond carries, and at these sizes the direct route no longer
    -- separates. The sizes are therefore taken where it does.
    [ testGroup "Bimodal Snapshots"
        [ testCase "N=96 K=8 u=0 v=48" $
            assertBool "should have 2 peaks" $ peakCount 96 4 0 48 1 47 >= 2
        , testCase "N=96 K=6 u=0 v=48" $
            assertBool "should have 2 peaks" $ peakCount 96 3 0 48 1 47 >= 2
        ]
    , testGroup "Unimodal Snapshots"
        [ testCase "N=18 u=0 v=3" $
            assertBool "should have 1 peak" $ peakCount 18 2 0 3 4 16 <= 1
        , testCase "N=24 u=0 v=4" $
            assertBool "should have 1 peak" $ peakCount 24 2 0 4 5 20 <= 1
        ]
    ]

determinismTests :: TestTree
determinismTests = testGroup "Determinism"
    [ testCase "repeated QR gives same result" $ do
        let w = lazyTransition 18 2 0 9
            e1 = sort $ exactEigenvalues w
            e2 = sort $ exactEigenvalues w
        assertBool "eigenvalues match" $
            all (\(a,b) -> abs (a - b) < 1e-10) (zip e1 e2)

    , testCase "repeated full pipeline identical" $ do
        let pmf1 = computePMF 20 2 0 10 1 9
            pmf2 = computePMF 20 2 0 10 1 9
        assertBool "PMFs match" $
            V.all (\(a,b) -> abs (a - b) < 1e-12) (V.zip pmf1 pmf2)
    ]

eigenvalueTests :: TestTree
eigenvalueTests = testGroup "Eigenvalue"
    [ testCase "ring lam0 = 1.0 exactly" $
        ringEigenvalue 12 2 0 @?= 1.0

    , testCase "N=12 lam6 = 0.5 (N/2 mode, lazy)" $ do
        let expected = 0.5 + 0.5 * (cos (pi) + cos (2*pi)) / 2
        assertBool "lam6" $ abs (ringEigenvalue 12 2 6 - expected) < 1e-10

    , testCase "N=18 K=2: 18 eigenvalues" $
        length (exactEigenvalues $ lazyTransition 18 2 0 9) @?= 18

    , testCase "N=18 K=2: all in [0,1]" $
        assertBool "bounded" $ all (\e -> abs e <= 1 + 1e-6) $
            exactEigenvalues $ lazyTransition 18 2 0 9

    , testCase "transition N=12 K=2 shortcut(0,6): has lam=1" $
        assertBool "has 1" $ any (\e -> abs (e - 1) < 1e-6) $
            exactEigenvalues $ lazyTransition 12 2 0 6
    ]

matrixEntryTests :: TestTree
matrixEntryTests = testGroup "Matrix Entries"
    [ testCase "ring diagonal = 0.5 for K=2" $
        assertBool "diag" $ abs (matrixGet (ringMatrix (qRing 2) 12 2) 0 0 - 0.5) < 1e-10

    , testCase "ring neighbour = 0.125 for K=2" $
        assertBool "neigh" $ abs (matrixGet (ringMatrix (qRing 2) 12 2) 0 1 - 0.125) < 1e-10

    -- An added edge raises the degree of its endpoints and leaves the holding
    -- probability where it was: the mobility is shared over one more
    -- destination, and the walker is no less likely to stay than before.
    , testCase "shortcut endpoints keep their holding probability" $ do
        let w = lazyTransition 18 2 0 9
            q = qRing 2
        assertBool "W[0,0] = 1-q" $ abs (matrixGet w 0 0 - (1 - q)) < 1e-10

    , testCase "non-shortcut diagonal = 1-q" $ do
        let w = lazyTransition 18 2 0 9
            q = qRing 2
        assertBool "W[1,1]=1-q" $ abs (matrixGet w 1 1 - (1 - q)) < 1e-10
    ]

marris2023Tests :: TestTree
marris2023Tests = testGroup "Marris 2023"
    [ testGroup "Known Bimodal"
        [ testCase "N=96 d=48 K=8" $
            assertBool "bimodal" $ peakCount 96 4 0 48 1 47 >= 2
        , testCase "N=200 d=100 K=8" $
            assertBool "bimodal" $ peakCount 200 4 0 100 1 99 >= 2
        ]
    , testGroup "Known Unimodal"
        [ testCase "N=18 d=3" $ assertBool "unimodal" $ peakCount 18 2 0 3 4 16 <= 1
        , testCase "N=30 d=3" $ assertBool "unimodal" $ peakCount 30 2 0 3 4 28 <= 1
        ]
    ]

marris2025Tests :: TestTree
marris2025Tests = testGroup "Marris 2025"
    [ testCase "q(K) formula" $ do
        let q2 = qRing 2; q3 = qRing 3
        assertBool "q(2) = 0.5" $ abs (q2 - 0.5) < 1e-10
        assertBool "q(3) = 2/3" $ abs (q3 - 2/3) < 1e-10

    , testCase "higher K -> smaller |lam2| (faster mixing)" $ do
        let lam2k2 = ringEigQ (qRing 2) 20 2 1
            lam2k3 = ringEigQ (qRing 3) 20 3 1
        assertBool "|lam2(K=3)| < |lam2(K=2)|" $ abs lam2k3 < abs lam2k2

    , testCase "K=2 has 4 neighbours" $
        assertBool "4 nonzero off-diagonal" $
            length (filter (> 1e-10) [matrixGet (ringMatrix (qRing 2) 12 2) 0 j | j <- [1..11]]) == 4

    , testCase "K=3 has 6 neighbours" $
        assertBool "6 nonzero off-diagonal" $
            length (filter (> 1e-10) [matrixGet (ringMatrix (qRing 3) 12 3) 0 j | j <- [1..11]]) == 6

    , testCase "dense QR spectrum matches closed eigenvalue formula (N=24,K=2)" $ do
        let n = 24; k = 2; q = qRing k
            qrEigs = sort $ exactEigenvalues (ringMatrix q n k)
            formulaEigs = sort [ringEigQ q n k j | j <- [0..n-1]]
        assertBool "dense vs formula" $
            all (\(a,b) -> abs (a - b) < 1e-4) (zip qrEigs formulaEigs)
    ]

eigenvalueFormulaTests :: TestTree
eigenvalueFormulaTests = testGroup "Eigenvalue Formula"
    [ testCase "lam0 = 1 always" $
        assertBool "lam0" $ abs (ringEigQ (qRing 2) 20 2 0 - 1.0) < 1e-12

    , testCase "QR agreement N=18 K=2" $ do
        let n = 18; k = 2; q = qRing k
            w = ringMatrix (qRing k) n k
            qrEigs = sort $ exactEigenvalues w
            formulaEigs = sort [ringEigQ q n k j | j <- [0..n-1]]
        assertBool "close" $ all (\(a,b) -> abs (a - b) < 1e-4) (zip qrEigs formulaEigs)
    ]
