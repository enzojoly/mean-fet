{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : DecompositionSpec
-- Description : The per-contact-site decomposition of the encounter law. The
--               weights and the weighted conditional times must sum to the
--               splitting weight and to the mean respectively, the per-site
--               distributions must sum to the total distribution at every
--               step, and the two limits must agree with one another. The
--               later groups drive the production cell builder rather than the
--               ring shortcut, because the stationary-law factors are unity on
--               a lattice of equal degrees and can only be checked where they
--               are not.

module DecompositionSpec (tests) where

import qualified Data.Vector.Unboxed as V
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Types (R, Matrix, Primitive(..))
import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist)
import Cells (Dim(..), BC(..), cellMatrixWith)
import Encounter
    ( encounterMeanTwoQ
    , encounterPMFTwoQAcc
    , splittingProbsTwoQ
    , encounterSiteMeansTwoQ
    , encounterDecompositionTwoQ
    , encounterSitePMFTwoQAcc
    )

tests :: TestTree
tests = testGroup "Per-site decomposition"
    [ accountingIdentity
    , weightNormalisation
    , marginalRegression
    , limitAgreement
    , mirrorSymmetry
    , decoratedLattice
    , twoDimensions
    , relabelling
    ]

close :: Double -> Double -> Double -> Bool
close tol a b = abs (a - b) <= tol * max 1 (max (abs a) (abs b))

type Pieces = ([R], Matrix, [R], Matrix, [Double])

-- | The undecorated ring, by the same route the earlier groups use.
pieces :: Double -> Double -> Int -> Int -> Pieces
pieces qA qB n k =
    let !wA = buildTransitionMatrix qA n k []
        !wB = buildTransitionMatrix qB n k []
        (!eA, !vA) = exactEigensystem wA
        (!eB, !vB) = exactEigensystem wB
        !piV = stationaryDist wA
    in (eA, vA, eB, vB, piV)

-- | Any cell, by the builder the executable uses. A decorated or
-- two-dimensional cell reaches the decomposition only through this path.
cellPieces :: Dim -> BC -> Int -> Int -> Double -> Double -> [Primitive]
           -> Pieces
cellPieces dim bc size k qA qB prims =
    let !wA = cellMatrixWith dim bc size k qA prims
        !wB = cellMatrixWith dim bc size k qB prims
        (!eA, !vA) = exactEigensystem wA
        (!eB, !vB) = exactEigensystem wB
        !piV = stationaryDist wA
    in (eA, vA, eB, vB, piV)

-- | The assertions that pin the stationary-law factors together. The weights
-- summing to one is insensitive to a misplaced factor, being the
-- well-conditioned direction; the site contributions summing to the mean is
-- not. Neither alone is sufficient.
consistentWith :: String -> Pieces -> Double -> (Int, Int) -> TestTree
consistentWith label (eA, vA, eB, vB, piV) rho starts =
    testCase label $ do
        let (!m, !w, !sm) = encounterDecompositionTwoQ eA vA eB vB piV rho starts
            !m0 = encounterMeanTwoQ eA vA eB vB piV rho starts
        assertBool ("weights sum " ++ show (V.sum w))
            (close 1e-13 (V.sum w) 1.0)
        assertBool ("site means sum " ++ show (V.sum sm)
                    ++ " against own mean " ++ show m)
            (close 1e-12 (V.sum sm) m)
        assertBool ("own mean " ++ show m ++ " against evaluation route "
                    ++ show m0)
            (close 1e-7 m m0)
        assertBool "every weight is a probability"
            (V.all (\p -> p >= negate 1e-12 && p <= 1 + 1e-12) w)
        assertBool "every site contribution is non-negative"
            (V.all (>= negate 1e-9) sm)

accountingIdentity :: TestTree
accountingIdentity = testGroup "Site means sum to the mean"
    [ testCase "ring N=12 k=2 q=0.75/0.35 rho=1" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.75 0.35 12 2
            !m  = encounterMeanTwoQ eA vA eB vB piV 1.0 (1, 7)
            !sm = encounterSiteMeansTwoQ eA vA eB vB piV 1.0 (1, 7)
        assertBool ("sum " ++ show (V.sum sm) ++ " against mean " ++ show m)
            (close 1e-8 (V.sum sm) m)

    , testCase "ring N=10 k=1 q=0.9/0.2 rho=0.6" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.9 0.2 10 1
            !m  = encounterMeanTwoQ eA vA eB vB piV 0.6 (0, 5)
            !sm = encounterSiteMeansTwoQ eA vA eB vB piV 0.6 (0, 5)
        assertBool ("sum " ++ show (V.sum sm) ++ " against mean " ++ show m)
            (close 1e-8 (V.sum sm) m)

    , testCase "the bundled call reproduces each separate route" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.4 12 2
            (!m, !w, !sm) = encounterDecompositionTwoQ eA vA eB vB piV 1.0 (2, 8)
            !m0 = encounterMeanTwoQ eA vA eB vB piV 1.0 (2, 8)
            !w0 = splittingProbsTwoQ eA vA eB vB piV 1.0 (2, 8)
            !s0 = encounterSiteMeansTwoQ eA vA eB vB piV 1.0 (2, 8)
        assertBool ("mean " ++ show m ++ " against " ++ show m0)
            (close 1e-8 m m0)
        assertBool "weights"
            (and (zipWith (close 1e-8) (V.toList w) (V.toList w0)))
        assertBool "site means"
            (and (zipWith (close 1e-12) (V.toList sm) (V.toList s0)))

    , testCase "the pole expansion is self-consistent to machine precision" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.4 12 2
            (!m, !w, !sm) = encounterDecompositionTwoQ eA vA eB vB piV 1.0 (2, 8)
        assertBool ("weights sum " ++ show (V.sum w))
            (close 1e-14 (V.sum w) 1.0)
        assertBool ("site means sum " ++ show (V.sum sm) ++ " against " ++ show m)
            (close 1e-13 (V.sum sm) m)
    ]

weightNormalisation :: TestTree
weightNormalisation = testGroup "Weights sum to the splitting weight"
    [ testCase "ring N=12 k=2 rho=1 sums to one" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.75 0.35 12 2
            !w = splittingProbsTwoQ eA vA eB vB piV 1.0 (1, 7)
        assertBool ("sum " ++ show (V.sum w)) (close 1e-8 (V.sum w) 1.0)

    , testCase "every weight is a probability" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.75 0.35 12 2
            !w = splittingProbsTwoQ eA vA eB vB piV 1.0 (1, 7)
        assertBool "in the unit interval"
            (V.all (\p -> p >= 0 && p <= 1) w)
    ]

marginalRegression :: TestTree
marginalRegression = testGroup "Site distributions sum to the total"
    [ testCase "ring N=8 k=1 q=0.8/0.5 rho=1, at every step" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.5 8 1
            !tmax = 400
            !tot = encounterPMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (0, 4) tmax
            !per = encounterSitePMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (0, 4) tmax
            !summed = foldr (V.zipWith (+)) (V.replicate (tmax + 1) 0) per
            !worst = maximum
                [ abs (V.unsafeIndex summed t - V.unsafeIndex tot t)
                | t <- [0 .. tmax] ]
        assertBool ("worst departure " ++ show worst) (worst < 1e-10)

    , testCase "one distribution is returned per contact site" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.5 8 1
            !per = encounterSitePMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (0, 4) 200
        assertBool ("returned " ++ show (length per)) (length per == 8)

    , testCase "decorated ring, at every step" $ do
        let (eA, vA, eB, vB, piV) =
                cellPieces D1 Periodic 12 1 0.75 0.4 [EdgeAdd 0 6]
            !tmax = 400
            !tot = encounterPMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (1, 7) tmax
            !per = encounterSitePMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (1, 7) tmax
            !summed = foldr (V.zipWith (+)) (V.replicate (tmax + 1) 0) per
            !worst = maximum
                [ abs (V.unsafeIndex summed t - V.unsafeIndex tot t)
                | t <- [0 .. tmax] ]
        assertBool ("worst departure " ++ show worst) (worst < 1e-10)
    ]

limitAgreement :: TestTree
limitAgreement = testGroup "The two limits agree per site"
    [ testCase "summed site distribution matches the weight" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.5 8 1
            !tmax = 4000
            !per = encounterSitePMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (0, 4) tmax
            !w = splittingProbsTwoQ eA vA eB vB piV 1.0 (0, 4)
            !worst = maximum
                (zipWith (\s p -> abs (V.sum s - p)) per (V.toList w))
        assertBool ("worst departure " ++ show worst) (worst < 1e-6)

    , testCase "first moment of the site distribution matches the site mean" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.8 0.5 8 1
            !tmax = 4000
            !per = encounterSitePMFTwoQAcc 14.0 eA vA eB vB piV 1.0 (0, 4) tmax
            !sm = encounterSiteMeansTwoQ eA vA eB vB piV 1.0 (0, 4)
            firstMoment s = V.sum (V.imap (\i p -> fromIntegral i * p) s)
            !worst = maximum
                (zipWith (\s m -> abs (firstMoment s - m)) per (V.toList sm))
        assertBool ("worst departure " ++ show worst) (worst < 1e-3)
    ]

mirrorSymmetry :: TestTree
mirrorSymmetry = testGroup "A symmetric release gives mirror-equal sites"
    [ testCase "ring N=12 k=1 released antipodally" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.6 0.6 12 1
            !w = splittingProbsTwoQ eA vA eB vB piV 1.0 (0, 6)
            pairAt c = abs ( V.unsafeIndex w c
                           - V.unsafeIndex w ((12 - c) `mod` 12) )
            !worst = maximum [ pairAt c | c <- [1 .. 5] ]
        assertBool ("worst departure " ++ show worst) (worst < 1e-9)

    , testCase "site means share the same mirror symmetry" $ do
        let (eA, vA, eB, vB, piV) = pieces 0.6 0.6 12 1
            !sm = encounterSiteMeansTwoQ eA vA eB vB piV 1.0 (0, 6)
            pairAt c = abs ( V.unsafeIndex sm c
                           - V.unsafeIndex sm ((12 - c) `mod` 12) )
            !worst = maximum [ pairAt c | c <- [1 .. 5] ]
        assertBool ("worst departure " ++ show worst) (worst < 1e-9)
    ]

-- | A shortcut raises the degree of its endpoints, so the stationary law is no
-- longer uniform and the factors of Equation (B') stop being one. They enter
-- the pole expansion in the residue contraction and in the bordered source
-- vector, neither of which the undecorated groups above exercise at all.
decoratedLattice :: TestTree
decoratedLattice = testGroup "Decorated ring, stationary law not uniform"
    [ consistentWith "ring N=12 k=1 with one added edge, rho=1"
        (cellPieces D1 Periodic 12 1 0.75 0.4 [EdgeAdd 0 6]) 1.0 (1, 7)

    , consistentWith "ring N=12 k=1 with one added edge, rho=0.55"
        (cellPieces D1 Periodic 12 1 0.75 0.4 [EdgeAdd 0 6]) 0.55 (1, 7)

    , consistentWith "ring N=14 k=2 with two added edges"
        (cellPieces D1 Periodic 14 2 0.8 0.3 [EdgeAdd 0 7, EdgeAdd 3 10])
        1.0 (1, 8)

    , testCase "the stationary law is genuinely non-uniform" $ do
        let (_, _, _, _, piV) =
                cellPieces D1 Periodic 12 1 0.75 0.4 [EdgeAdd 0 6]
            !spread = maximum piV - minimum piV
        assertBool ("spread " ++ show spread) (spread > 1e-3)
    ]

-- | The square lattice, where the eigensystem is a Kronecker composition and
-- the spectrum carries heavy degeneracy. The linear size is kept small so the
-- contact set stays at twenty-five sites.
twoDimensions :: TestTree
twoDimensions = testGroup "Two dimensions"
    [ consistentWith "torus L=5 k=1 q=0.75/0.4 rho=1"
        (cellPieces D2 Periodic 5 1 0.75 0.4 []) 1.0 (6, 18)

    , consistentWith "torus L=5 k=1 q=0.75/0.4 rho=0.6"
        (cellPieces D2 Periodic 5 1 0.75 0.4 []) 0.6 (6, 18)

    , consistentWith "reflecting box L=5 k=1 q=0.75/0.4 rho=1"
        (cellPieces D2 Reflecting 5 1 0.75 0.4 []) 1.0 (6, 18)

    , consistentWith "reflecting box L=5 k=2 q=0.8/0.25 rho=1"
        (cellPieces D2 Reflecting 5 2 0.8 0.25 []) 1.0 (6, 18)

    , consistentWith "reflecting box L=4 k=1, asymmetric release"
        (cellPieces D2 Reflecting 4 1 0.75 0.1 []) 1.0 (5, 10)
    ]

-- | Relabelling the walkers exchanges both the mobilities and the starting
-- sites. The mean is already asserted invariant elsewhere; the decomposition
-- must be invariant component by component, which holds on decorated and
-- two-dimensional cells where mirror symmetry does not apply.
relabelling :: TestTree
relabelling = testGroup "Relabelling leaves the decomposition fixed"
    [ testCase "decorated ring N=12 k=1" $ do
        let (eA, vA, eB, vB, piV) =
                cellPieces D1 Periodic 12 1 0.75 0.4 [EdgeAdd 0 6]
            (fA, uA, fB, uB, piW) =
                cellPieces D1 Periodic 12 1 0.4 0.75 [EdgeAdd 0 6]
            (_, !w1, !s1) = encounterDecompositionTwoQ eA vA eB vB piV 1.0 (1, 7)
            (_, !w2, !s2) = encounterDecompositionTwoQ fA uA fB uB piW 1.0 (7, 1)
        assertBool "weights"
            (and (zipWith (close 1e-11) (V.toList w1) (V.toList w2)))
        assertBool "site means"
            (and (zipWith (close 1e-11) (V.toList s1) (V.toList s2)))

    , testCase "reflecting box L=5 k=1" $ do
        let (eA, vA, eB, vB, piV) = cellPieces D2 Reflecting 5 1 0.75 0.4 []
            (fA, uA, fB, uB, piW) = cellPieces D2 Reflecting 5 1 0.4 0.75 []
            (_, !w1, !s1) = encounterDecompositionTwoQ eA vA eB vB piV 1.0 (6, 18)
            (_, !w2, !s2) = encounterDecompositionTwoQ fA uA fB uB piW 1.0 (18, 6)
        assertBool "weights"
            (and (zipWith (close 1e-11) (V.toList w1) (V.toList w2)))
        assertBool "site means"
            (and (zipWith (close 1e-11) (V.toList s1) (V.toList s2)))
    ]
