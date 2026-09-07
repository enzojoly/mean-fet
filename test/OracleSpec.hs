{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : OracleSpec
-- Description : An independent solution of the two-walker problem, used to
--               judge the spectral engines rather than to agree with them.
--
--               The pair state is carried explicitly on the product of the two
--               single-walker state spaces, the diagonal is damped by the
--               absorption probability, and the mean is the solution of the
--               linear system it satisfies. Nothing here shares a line with the
--               machinery under test: no generating function is formed, no
--               spectrum is taken, no contour is traversed. The cost is the
--               cube of the pair dimension, which confines the method to small
--               lattices -- and that is precisely where it is wanted, since a
--               construction that is correct at one size and wrong at another
--               is not a construction but a coincidence.
--
--               This closes the standing gap in the two-dimensional results,
--               which until now rested on a single implementation with no
--               external reference of any kind.

module OracleSpec (tests) where

import qualified Numeric.LinearAlgebra as LA
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist)
import Encounter (encounterMeanRingTwoQ, encounterMeanTwoQ)
import Lattice (encounterMeanTorus, encounterMeanReflect)
import Types (Primitive(..))

tests :: TestTree
tests = testGroup "Independent pair-chain oracle"
    [ oracleAgainstItself
    , ringAgainstOracle
    , decoratedAgainstOracle
    , torusAgainstOracle
    , reflectAgainstOracle
    , absorptionAgainstOracle
    ]

-- | A single-walker step law as a sparse list of destinations and weights.
type StepLaw = Int -> [(Int, Double)]

-- | Mean time until the two walkers first occupy a common site and are
-- absorbed there, obtained by solving the pair chain directly.
--
-- Writing T for the mean from each joint state, one step costs unit time and
-- lands on a joint state that is retained with the full weight unless it lies
-- on the diagonal, where the surviving fraction is one minus the absorption
-- probability. Hence T = 1 + M T with M the damped joint kernel, and the mean
-- is the solution of (I - M) T = 1. Absorption is tested after the step, never
-- before, so a pair released from a common site is not absorbed at time zero.
pairMean :: Int -> StepLaw -> StepLaw -> Double -> (Int, Int) -> Double
pairMean !ns stepA stepB !rho (!a0, !b0) =
    let !n2 = ns * ns
        idx !a !b = a * ns + b

        !entries =
            [ ((idx a b, idx a' b'), pa * pb * damp a' b')
            | a <- [0 .. ns - 1]
            , b <- [0 .. ns - 1]
            , (a', pa) <- stepA a
            , (b', pb) <- stepB b
            ]
        damp !a' !b' = if a' == b' then 1 - rho else 1

        !zero = LA.konst 0 (n2, n2) :: LA.Matrix Double
        !kern = LA.accum zero (+) entries
        !sys  = LA.ident n2 - kern
        !rhs  = LA.konst 1 n2 :: LA.Vector Double
        !sol  = sys LA.<\> rhs
    in sol `LA.atIndex` idx a0 b0

ringStep :: Int -> Int -> Double -> StepLaw
ringStep !n !k !q = \s ->
    (s, 1 - q)
    : concat [ [ ((s + m) `mod` n, w), ((s - m + n * m) `mod` n, w) ]
             | m <- [1 .. k] ]
  where
    !w = q / (2 * fromIntegral k)

torusStep :: Int -> Double -> StepLaw
torusStep !l !q = \i ->
    let !x = i `mod` l
        !y = i `div` l
        nb = [ ((x + 1) `mod` l, y), ((x - 1 + l) `mod` l, y)
             , (x, (y + 1) `mod` l), (x, (y - 1 + l) `mod` l) ]
    in (i, 1 - q) : [ (u + l * v, q / 4) | (u, v) <- nb ]

-- | Reflection by refusal: a step that would leave the box is not taken, and
-- its weight is returned to the site of origin. This is the convention the
-- cosine eigenbasis diagonalises.
reflectStep :: Int -> Double -> StepLaw
reflectStep !l !q = \i ->
    let !x = i `mod` l
        !y = i `div` l
        cand = [ (x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1) ]
        inside (u, v) = u >= 0 && u < l && v >= 0 && v < l
        moves = [ (u + l * v, q / 4) | (u, v) <- cand, inside (u, v) ]
        !blocked = fromIntegral (length (filter (not . inside) cand)) * q / 4
    in (i, 1 - q + blocked) : moves

within :: Double -> Double -> Double -> Bool
within tol want got
    | want == 0 = abs got <= tol
    | otherwise = abs (got - want) / abs want <= tol

report :: String -> Double -> Double -> String
report lbl want got =
    lbl ++ ": oracle " ++ show want ++ ", engine " ++ show got
    ++ ", relative " ++ show (abs (got - want) / abs want)

-- | The oracle is checked against a value obtained without it before being
-- used to judge anything else. A ring whose mean is known in closed form
-- serves: if the pair construction, the damping convention or the timing of
-- absorption were wrong, this would not agree.
oracleAgainstItself :: TestTree
oracleAgainstItself = testGroup "The oracle itself"
    [ testCase "ring N=16 k=1 q=0.75/0.6 rho=1 matches the closed form" $ do
        let !o = pairMean 16 (ringStep 16 1 0.75) (ringStep 16 1 0.6) 1.0 (1, 9)
            !c = encounterMeanRingTwoQ 0.75 0.6 16 1 1.0 1 9
        assertBool (report "ring closed" o c) (within 1e-10 o c)
    , testCase "ring N=12 k=1 q=0.9/0.3 rho=0.5 matches the closed form" $ do
        let !o = pairMean 12 (ringStep 12 1 0.9) (ringStep 12 1 0.3) 0.5 (0, 6)
            !c = encounterMeanRingTwoQ 0.9 0.3 12 1 0.5 0 6
        assertBool (report "ring closed rho" o c) (within 1e-10 o c)
    ]

-- | A lattice whose degrees are unequal. Every cell checked above is regular,
-- so its stationary law is uniform and every correction for it is unity; a
-- decorated lattice is the only place those corrections can be wrong and be
-- seen to be wrong.
decoratedAgainstOracle :: TestTree
decoratedAgainstOracle = testGroup "Unequal degrees"
    [ testCase "ring N=10 with one added edge" $
        check 10 1 0.7 0.7 1.0 (0, 5) (1, 6)
    , testCase "ring N=12 with one added edge, unequal mobilities" $
        check 12 1 0.8 0.45 1.0 (0, 6) (2, 9)
    , testCase "ring N=10 with one added edge, imperfect absorption" $
        check 10 1 0.7 0.7 0.4 (0, 5) (1, 6)
    , testCase "ring N=12 k=2 with one added edge" $
        check 12 2 0.75 0.6 1.0 (0, 6) (1, 7)
    ]
  where
    check n k q1 q2 rho (u, v) (a, b) = do
        let !wA = decorated n k q1 u v
            !wB = decorated n k q2 u v
            !o = pairMean n (lawOf wA) (lawOf wB) rho (a, b)
            !e = encounterMeanDecorated n k q1 q2 rho u v a b
        assertBool (report ("N=" ++ show n) o e) (within 1e-7 o e)

-- | Row-stochastic transition law of a ring carrying one added edge, built here
-- rather than taken from the library, so that the oracle judges the library
-- rather than agreeing with it.
decorated :: Int -> Int -> Double -> Int -> Int -> LA.Matrix Double
decorated n k q u v = LA.accum (LA.konst 0 (n, n)) (+) entries
  where
    nbr c = [ (c + m) `mod` n | m <- [1 .. k] ]
         ++ [ (c - m + n * k) `mod` n | m <- [1 .. k] ]
    !bigK = fromIntegral (2 * k) :: Double
    deg c = if c == u || c == v then bigK + 1 else bigK
    entries = concat
        [ ((c, c), 1 - q)
          : [ ((c, j), q / deg c) | j <- nbr c ]
          ++ [ ((u, v), q / (bigK + 1)) | c == u ]
          ++ [ ((v, u), q / (bigK + 1)) | c == v ]
        | c <- [0 .. n - 1] ]

lawOf :: LA.Matrix Double -> StepLaw
lawOf w = \s -> [ (t, w `LA.atIndex` (s, t))
                | t <- [0 .. LA.cols w - 1]
                , w `LA.atIndex` (s, t) /= 0 ]

encounterMeanDecorated :: Int -> Int -> Double -> Double -> Double
                       -> Int -> Int -> Int -> Int -> Double
encounterMeanDecorated n k q1 q2 rho u v a b =
    let !wA = buildTransitionMatrix q1 n k [EdgeAdd u v]
        !wB = buildTransitionMatrix q2 n k [EdgeAdd u v]
        (!eA, !vA) = exactEigensystem wA
        (!eB, !vB) = exactEigensystem wB
        !piV = stationaryDist wA
    in encounterMeanTwoQ eA vA eB vB piV rho (a, b)

ringAgainstOracle :: TestTree
ringAgainstOracle = testGroup "Ring engine"
    [ testCase "N=20 k=2 homogeneous" $ do
        let !o = pairMean 20 (ringStep 20 2 0.7) (ringStep 20 2 0.7) 1.0 (2, 12)
            !e = encounterMeanRingTwoQ 0.7 0.7 20 2 1.0 2 12
        assertBool (report "ring k=2" o e) (within 1e-10 o e)
    , testCase "N=18 k=1 strongly asymmetric mobilities" $ do
        let !o = pairMean 18 (ringStep 18 1 0.95) (ringStep 18 1 0.15) 1.0 (0, 9)
            !e = encounterMeanRingTwoQ 0.95 0.15 18 1 1.0 0 9
        assertBool (report "ring asymmetric" o e) (within 1e-10 o e)
    ]

torusAgainstOracle :: TestTree
torusAgainstOracle = testGroup "Torus engine"
    [ testCase "L=4 heterogeneous with imperfect absorption" $ do
        let !o = pairMean 16 (torusStep 4 0.7) (torusStep 4 0.5) 0.6 (0, 2 + 4 * 1)
            !e = encounterMeanTorus 4 1 0.7 0.5 0.6 (0, 0) (2, 1)
        assertBool (report "torus L=4" o e) (within 1e-9 o e)
    , testCase "L=5 homogeneous" $ do
        let !o = pairMean 25 (torusStep 5 0.75) (torusStep 5 0.75) 1.0 (0, 2 + 5 * 2)
            !e = encounterMeanTorus 5 1 0.75 0.75 1.0 (0, 0) (2, 2)
        assertBool (report "torus L=5" o e) (within 1e-9 o e)
    , testCase "L=6 heterogeneous" $ do
        let !o = pairMean 36 (torusStep 6 0.8) (torusStep 6 0.35) 1.0 (0, 3 + 6 * 3)
            !e = encounterMeanTorus 6 1 0.8 0.35 1.0 (0, 0) (3, 3)
        assertBool (report "torus L=6" o e) (within 1e-9 o e)
    ]

-- | The two-dimensional reflecting results have until now had no external
-- reference of any kind, every cell having been produced by the same engine
-- that is being asked to confirm them.
reflectAgainstOracle :: TestTree
reflectAgainstOracle = testGroup "Reflecting box engine"
    [ testCase "L=4 homogeneous" $ do
        let !o = pairMean 16 (reflectStep 4 0.7) (reflectStep 4 0.7) 1.0 (0, 2 + 4 * 1)
            !e = encounterMeanReflect 4 1 0.7 0.7 1.0 (0, 0) (2, 1)
        assertBool (report "reflect L=4" o e) (within 1e-9 o e)
    , testCase "L=5 heterogeneous" $ do
        let !o = pairMean 25 (reflectStep 5 0.8) (reflectStep 5 0.4) 1.0 (0, 2 + 5 * 2)
            !e = encounterMeanReflect 5 1 0.8 0.4 1.0 (0, 0) (2, 2)
        assertBool (report "reflect L=5" o e) (within 1e-9 o e)
    , testCase "L=6 homogeneous, walkers at opposite corners" $ do
        let !o = pairMean 36 (reflectStep 6 0.75) (reflectStep 6 0.75) 1.0 (0, 5 + 6 * 5)
            !e = encounterMeanReflect 6 1 0.75 0.75 1.0 (0, 0) (5, 5)
        assertBool (report "reflect L=6" o e) (within 1e-9 o e)
    ]

-- | Imperfect absorption on both geometries. Each failed capture returns the
-- pair to circulation, so the cost of a lowered absorption probability is a
-- whole return excursion; an engine that mishandled the geometric series over
-- returns would agree at unit absorption and diverge below it.
absorptionAgainstOracle :: TestTree
absorptionAgainstOracle = testGroup "Imperfect absorption"
    [ testCase "torus L=4, rho = 0.25" $ do
        let !o = pairMean 16 (torusStep 4 0.75) (torusStep 4 0.75) 0.25 (0, 2 + 4 * 2)
            !e = encounterMeanTorus 4 1 0.75 0.75 0.25 (0, 0) (2, 2)
        assertBool (report "torus rho" o e) (within 1e-9 o e)
    , testCase "reflecting L=4, rho = 0.4" $ do
        let !o = pairMean 16 (reflectStep 4 0.75) (reflectStep 4 0.75) 0.4 (0, 2 + 4 * 2)
            !e = encounterMeanReflect 4 1 0.75 0.75 0.4 (0, 0) (2, 2)
        assertBool (report "reflect rho" o e) (within 1e-9 o e)
    , testCase "ring N=16, rho = 0.3" $ do
        let !o = pairMean 16 (ringStep 16 1 0.6) (ringStep 16 1 0.6) 0.3 (1, 9)
            !e = encounterMeanRingTwoQ 0.6 0.6 16 1 0.3 1 9
        assertBool (report "ring rho" o e) (within 1e-10 o e)
    ]
