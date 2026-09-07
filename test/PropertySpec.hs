{-# LANGUAGE ScopedTypeVariables #-}
-- |
-- Module      : PropertySpec
-- Description : Universally quantified invariants for ring, matrix, GF, and
--               simulation.

module PropertySpec (tests) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC

import Data.Complex (Complex(..), realPart)
import qualified Data.Vector.Unboxed as V
import Data.List (sort)

import Types
import Ring
import Defect
import Passage
import Distribution
import Simulate

ringEigenvalue :: N -> K -> Int -> R
ringEigenvalue n k = ringEigQ (qRing k) n k

lazyTransition :: N -> K -> Pos -> Pos -> Matrix
lazyTransition n k u v = buildTransitionMatrix (qRing k) n k [EdgeAdd u v]

newtype RingSize = RingSize Int deriving (Show, Eq)
instance Arbitrary RingSize where
    arbitrary = RingSize <$> chooseInt (4, 50)

data MatrixCfg = MatrixCfg Int Int Int Int deriving Show
instance Arbitrary MatrixCfg where
    arbitrary = do
        n <- chooseInt (6, 25); k <- chooseInt (2, min 3 ((n-1) `div` 2))
        u <- chooseInt (0, n-1); v <- chooseInt (0, n-1) `suchThat` (/= u)
        return $ MatrixCfg n k u v

data FullCfg = FullCfg Int Int Int Int Int Int deriving Show
instance Arbitrary FullCfg where
    arbitrary = do
        n <- chooseInt (8, 30); k <- chooseInt (2, min 3 ((n-1) `div` 2))
        u <- chooseInt (0, n-1); v <- chooseInt (0, n-1) `suchThat` (/= u)
        src <- chooseInt (0, n-1); tgt <- chooseInt (0, n-1) `suchThat` (/= src)
        return $ FullCfg n k u v src tgt

validShortcut :: Int -> Int -> Int -> Int -> Bool
validShortcut n k u v = u /= v && ringDist n u v > k

tests :: TestTree
tests = testGroup "Property"
    [ ringProps
    , matrixProps
    , conservationProps
    , symmetryProps
    , gfProps
    , simulationProps
    ]

ringProps :: TestTree
ringProps = testGroup "Ring"
    [ QC.testProperty "eigenvalues in [0,1] for K=2" $ \(RingSize n) ->
        all (\j -> let e = ringEigenvalue n 2 j in e >= -1e-10 && e <= 1 + 1e-10)
            [0..n-1]

    , QC.testProperty "lam0 = 1 for all N" $ \(RingSize n) ->
        abs (ringEigenvalue n 2 0 - 1.0) < 1e-12

    , QC.testProperty "symmetry: lam_ell = lam_{N-ell}" $ \(RingSize n) ->
        all (\j -> abs (ringEigenvalue n 2 j - ringEigenvalue n 2 (n - j)) < 1e-10)
            [1 .. n `div` 2]

    , QC.testProperty "N eigenvalues returned" $ \(RingSize n) ->
        length [ringEigenvalue n 2 j | j <- [0..n-1]] == n

    , QC.testProperty "eigenvector has N components" $ \(RingSize n) ->
        length (ringEigenvector n 0) == n

    , QC.testProperty "lam0 = 1 for any q" $ forAll (choose (0.01, 1.0)) $ \q ->
        forAll (chooseInt (6, 30)) $ \n ->
            abs (ringEigQ q n 2 0 - 1.0) < 1e-12

    , QC.testProperty "lazy and non-lazy eigenvalues related by q" $ \(RingSize n) ->
        let k = 2; q = qRing k
        in all (\j ->
            let lazyE = ringEigQ q n k j
                nonLazyE = ringEigNonLazy n k j
                expected = (1 - q) + q * nonLazyE
            in abs (lazyE - expected) < 1e-10
            ) [0..n-1]
    ]

matrixProps :: TestTree
matrixProps = testGroup "Matrix"
    [ QC.testProperty "transition matrix stochastic" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            all (\s -> abs (s-1) < 1e-8) $ map sum $ matrixToLists $ lazyTransition n k u v

    , QC.testProperty "eigenvalues bounded by 1" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            all (\e -> abs e <= 1 + 1e-6) $ exactEigenvalues $ lazyTransition n k u v

    , QC.testProperty "eigenvalue count = N" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            length (exactEigenvalues $ lazyTransition n k u v) == n

    , QC.testProperty "perturbation rows sum to 0" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            let w0 = ringMatrix (qRing k) n k
                w  = lazyTransition n k u v
                vMat = matrixToLists w `zip` matrixToLists w0
            in all (\(wr, w0r) ->
                abs (sum wr - sum w0r) < 1e-10) vMat
    ]

conservationProps :: TestTree
conservationProps = testGroup "Conservation"
    [ QC.testProperty "trace = sum eigenvalues" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            let w = lazyTransition n k u v
                tr = sum [matrixGet w i i | i <- [0..n-1]]
                eigSum = sum $ exactEigenvalues w
            in abs (tr - eigSum) < 1e-6

    , testCase "ring columns sum to 1 (doubly stochastic)" $ do
        let m = ringMatrix (qRing 2) 12 2
        let colSums = [sum [matrixGet m i j | i <- [0..11]] | j <- [0..11]]
        assertBool "cols=1" $ all (\s -> abs (s-1) < 1e-10) colSums
    ]

symmetryProps :: TestTree
symmetryProps = testGroup "Symmetry"
    [ QC.testProperty "(u,v) <-> (v,u) same eigenvalues" $ \(MatrixCfg n k u v) ->
        validShortcut n k u v ==>
            let e1 = sort $ exactEigenvalues $ lazyTransition n k u v
                e2 = sort $ exactEigenvalues $ lazyTransition n k v u
            in all (\(a,b) -> abs (a - b) < 1e-8) (zip e1 e2)

    , QC.testProperty "ring eigenvalue symmetry" $ \(RingSize n) ->
        all (\j -> abs (ringEigenvalue n 2 j - ringEigenvalue n 2 (n - j)) < 1e-10)
            [1 .. n `div` 2]
    ]

gfProps :: TestTree
gfProps = testGroup "GF"
    [ QC.testProperty "pure ring PMF non-negative" $ \(FullCfg n k u v src tgt) ->
        validShortcut n k u v && src /= tgt ==>
            let q = qRing k
                gf z = pureRingFPGF q n k src tgt z
                pmf = invertPGF 200 gf
            in V.all (>= -1e-8) pmf

    , QC.testProperty "pure ring PMF sums close to 1" $ \(FullCfg n k u v src tgt) ->
        validShortcut n k u v && src /= tgt ==>
            let q = qRing k
                gf z = pureRingFPGF q n k src tgt z
                pmf = invertPGF 500 gf
                total = V.sum pmf
            in total > 0.85 && total < 1.01

    , QC.testProperty "P(T=0) = 0" $ \(FullCfg n k u v src tgt) ->
        validShortcut n k u v && src /= tgt ==>
            let q = qRing k
                gf z = pureRingFPGF q n k src tgt z
                pmf = invertPGF 200 gf
            in abs (pmf `V.unsafeIndex` 0) < 1e-6

    , QC.testProperty "MFPT from GF matches exact" $ \(RingSize n) ->
        n >= 6 ==>
            let k = 2; src = 0; tgt = n `div` 2
                exact = exactMFPT n k src tgt
                gf z = pureRingFPGF 1.0 n k src tgt z
                pmf = invertPGF 2000 gf
                gfMFPT = V.sum $ V.imap (\i p -> fromIntegral (i+1) * p) (V.tail pmf)
            in abs (gfMFPT - exact) < 0.5 + 0.02 * exact

    , QC.testProperty "network PMF non-negative" $
        forAll (chooseInt (8, 20)) $ \n ->
            let k = 2; q = qRing k
                u = 0; v = n `div` 2; src = 1; tgt = v - 1
                defs = primitivesToDefects q n k [EdgeAdd u v]
                gf z = firstPassageGF q n k defs src tgt z
                pmf = invertPGF 200 gf
            in V.all (>= -1e-6) pmf
    ]

simulationProps :: TestTree
simulationProps = testGroup "Simulation"
    [ QC.testProperty "walker accounting: absorbed + survived = total" $
        forAll (chooseInt (100, 500)) $ \w ->
            let sr = simulate $ defaultSimConfig
                    { simN = 12, simK = 2, simQ = 0.75
                    , simSrc = 0, simTgt = 6, simDefects = []
                    , simWalkers = w, simSeed = 42, simMaxT = 500 }
            in srAbsorbed sr + srSurvived sr == w

    , QC.testProperty "deterministic: same config same result" $
        forAll (chooseInt (50, 200)) $ \w ->
            let cfg = defaultSimConfig
                    { simN = 10, simK = 2, simQ = 0.75
                    , simSrc = 0, simTgt = 5, simDefects = []
                    , simWalkers = w, simSeed = 42, simMaxT = 300 }
                r1 = simulate cfg
                r2 = simulate cfg
            in srAbsorbed r1 == srAbsorbed r2

    , QC.testProperty "encounter: two walkers eventually meet" $
        forAll (chooseInt (6, 12)) $ \n ->
            let sr = simulate $ defaultSimConfig
                    { simN = n, simK = 1, simQ = 0.75
                    , simSrc = 0, simSrcB = n `div` 2
                    , simDefects = []
                    , simWalkers = 50, simSeed = 42, simMaxT = 5000
                    , simMode = Simulate.Encounter }
            in srAbsorbed sr > 0
    ]
