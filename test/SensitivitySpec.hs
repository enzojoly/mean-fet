{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : SensitivitySpec
-- Description : The mobility sensitivity of the mean first-encounter time and
--               its decomposition by contact site, at perfect absorption. Six
--               groups. Mobility is a pure time change for one walker and is
--               not one for a pair, which is the zero and the signal of the
--               between-site contribution. The frozen partner reduces the
--               encounter mean to the single-walker mean exactly. The bare ring
--               carries a stationary point of its own, but only where a mode of
--               the relative spectrum approaches unity: it is a minimum, it is
--               confined to the last few per cent of the mobility range, it
--               requires a bipartite lattice, and it requires the fixed
--               mobility to exceed one half. The midpoint decomposition is an
--               algebraic identity at every step, and the between-site term
--               depends on the partition of the contact set while the total
--               does not.

module SensitivitySpec (tests) where

import Data.List (zipWith4)
import qualified Data.Vector.Unboxed as V
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Types (Primitive(..))
import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist,
               primitivesToDefects)
import Ring (ringEigQ)
import Passage (firstAbsorptionGF, meanFromGF, exactMFPTq)
import Encounter (encounterDecompositionTwoQ, encounterMeanRingTwoQ)

tests :: TestTree
tests = testGroup "Mobility sensitivity"
    [ timeChangeOneWalker
    , timeChangeFailsOnPair
    , frozenPartner
    , ringStationaryPoint
    , midpointExactness
    , partitionDependence
    ]

relClose :: Double -> Double -> Double -> Bool
relClose tol a b = abs (a - b) <= tol * max 1 (max (abs a) (abs b))

report :: String -> Double -> Double -> String
report lbl want got =
    lbl ++ ": want " ++ show want ++ " got " ++ show got
        ++ " absolute " ++ show (abs (want - got))

mobilities :: [Double]
mobilities = [0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 0.95]

nonLazyEig :: Int -> Int -> Int -> Double
nonLazyEig n k l = ringEigQ 1.0 n k l

ringWeight :: Int -> Int -> Int -> Double
ringWeight n d0 l =
    1 - cos (2 * pi * fromIntegral (l * d0) / fromIntegral n)

ringScaledMean :: Double -> Int -> Int -> Int -> Double
ringScaledMean q n k d0 = sum
    [ let !a = 1 - nonLazyEig n k l
      in ringWeight n d0 l / (a * (2 - q * a))
    | l <- [1 .. n - 1] ]

ringMeanDeriv :: Double -> Double -> Int -> Int -> Int -> Double
ringMeanDeriv q1 q2 n k d0 = sum
    [ let !a  = 1 - nonLazyEig n k l
          !l1 = ringEigQ q1 n k l
          !l2 = ringEigQ q2 n k l
          !gp = 1 - l1 * l2
      in negate (ringWeight n d0 l * l1 * a) / (gp * gp)
    | l <- [1 .. n - 1] ]

ringMeanDiff :: Double -> Double -> Int -> Int -> Int -> Double -> Double
ringMeanDiff q1 q2 n k d0 h =
    ( encounterMeanRingTwoQ q1 (q2 + h) n k 1.0 0 d0
    - encounterMeanRingTwoQ q1 (q2 - h) n k 1.0 0 d0 ) / (2 * h)

leastRelativeGap :: Double -> Double -> Int -> Int -> Double
leastRelativeGap q1 q2 n k = minimum
    [ abs (1 - ringEigQ q1 n k l * ringEigQ q2 n k l) | l <- [1 .. n - 1] ]

leastNonLazyEig :: Int -> Int -> Double
leastNonLazyEig n k = minimum [ nonLazyEig n k l | l <- [1 .. n - 1] ]

bisect :: (Double -> Double) -> Double -> Double -> Int -> Double
bisect f !lo !hi !i
    | i <= 0    = 0.5 * (lo + hi)
    | otherwise =
        let !mid = 0.5 * (lo + hi)
        in if f lo * f mid <= 0
           then bisect f lo mid (i - 1)
           else bisect f mid hi (i - 1)

singleWalkerMean :: Double -> Int -> Int -> [Primitive] -> Int -> Int -> Double
singleWalkerMean q n k prims n0 nt =
    meanFromGF (firstAbsorptionGF q n k (primitivesToDefects q n k prims)
                                  n0 nt 1.0)

type Decomp = (Double, V.Vector Double, V.Vector Double)

decompAt :: Int -> Int -> [Primitive] -> Double -> Double -> (Int, Int)
         -> Decomp
decompAt n k prims q1 q2 starts =
    let !wA = buildTransitionMatrix q1 n k prims
        !wB = buildTransitionMatrix q2 n k prims
        (!eA, !vA) = exactEigensystem wA
        (!eB, !vB) = exactEigensystem wB
        !piV = stationaryDist wA
    in encounterDecompositionTwoQ eA vA eB vB piV 1.0 starts

blockSums :: [[Int]] -> V.Vector Double -> [Double]
blockSums blocks v = [ sum [ v V.! j | j <- b ] | b <- blocks ]

guardedQuot :: Double -> Double -> Double
guardedQuot m p = if abs p < 1e-300 then 0 else m / p

midpointAB :: [[Int]] -> Decomp -> Decomp -> Double -> (Double, Double)
midpointAB blocks (_, wLo, mLo) (_, wHi, mHi) h =
    let pLo = blockSums blocks wLo
        pHi = blockSums blocks wHi
        eLo = zipWith guardedQuot (blockSums blocks mLo) pLo
        eHi = zipWith guardedQuot (blockSums blocks mHi) pHi
        aTerm pl ph el eh = ((ph - pl) / (2 * h)) * (0.5 * (eh + el))
        bTerm pl ph el eh = (0.5 * (ph + pl)) * ((eh - el) / (2 * h))
    in ( sum (zipWith4 aTerm pLo pHi eLo eHi)
       , sum (zipWith4 bTerm pLo pHi eLo eHi) )

fineBlocks :: Int -> [[Int]]
fineBlocks n = [ [j] | j <- [0 .. n - 1] ]

trivialBlock :: Int -> [[Int]]
trivialBlock n = [ [0 .. n - 1] ]

halfBlocks :: Int -> [[Int]]
halfBlocks n = [ [0 .. n `div` 2 - 1], [n `div` 2 .. n - 1] ]

timeChangeOneWalker :: TestTree
timeChangeOneWalker = testGroup "Mobility is a pure time change for one walker"
    [ testCase "frozen partner on the ring: q E is constant to 1e-12" $ do
        let n = 48
            k = 1
            vals = [ q * encounterMeanRingTwoQ q 0.0 n k 1.0 0 17
                   | q <- mobilities ]
            ref = head vals
        mapM_ (\(q, v) -> assertBool (report ("q=" ++ show q) ref v)
                                     (relClose 1e-12 ref v))
              (zip mobilities vals)

    , testCase "bare ring by the defect route: q T is constant to 1e-5" $
        invariant 24 1 [] 0 7

    , testCase "shortcut ring: q T is constant to 1e-5" $
        invariant 24 1 [EdgeAdd 0 12] 0 7

    , testCase "deleted edge: q T is constant to 1e-5" $
        invariant 24 1 [EdgeDel 5 6] 0 7

    , testCase "permeable barrier: q T is constant to 1e-5" $
        invariant 24 1 [Barrier 5 6 0.3] 0 7

    , testCase "range two: q T is constant to 1e-5" $
        invariant 24 2 [] 0 7
    ]
  where
    invariant n k prims n0 nt = do
        let vals = [ q * singleWalkerMean q n k prims n0 nt | q <- mobilities ]
            ref  = head vals
        mapM_ (\(q, v) -> assertBool (report ("q=" ++ show q) ref v)
                                     (relClose 1e-5 ref v))
              (zip mobilities vals)

timeChangeFailsOnPair :: TestTree
timeChangeFailsOnPair = testGroup "Mobility is not a time change for a pair"
    [ testCase "closed form matches the scaled spectral sum, N=48 k=1" $
        matchesScaled 48 1 24

    , testCase "closed form matches the scaled spectral sum, N=24 k=2" $
        matchesScaled 24 2 11

    , testCase "q E(q,q) is strictly increasing, N=48 k=1 d0=24" $
        increasing 48 1 24

    , testCase "q E(q,q) is strictly increasing, N=24 k=2 d0=11" $
        increasing 24 2 11

    , testCase "the departure from constancy exceeds five per cent" $ do
        let lo = 0.10 * encounterMeanRingTwoQ 0.10 0.10 48 1 1.0 0 24
            hi = 0.95 * encounterMeanRingTwoQ 0.95 0.95 48 1 1.0 0 24
        assertBool ("ratio " ++ show (hi / lo)) (hi / lo > 1.05)
    ]
  where
    matchesScaled n k d0 =
        mapM_ (\q ->
            let !got  = q * encounterMeanRingTwoQ q q n k 1.0 0 d0
                !want = ringScaledMean q n k d0
            in assertBool (report ("q=" ++ show q) want got)
                          (relClose 1e-11 want got))
        mobilities

    increasing n k d0 = do
        let vals = [ q * encounterMeanRingTwoQ q q n k 1.0 0 d0
                   | q <- mobilities ]
        mapM_ (\(p, (a, b)) ->
                assertBool ("step " ++ show p ++ ": " ++ show a
                            ++ " then " ++ show b)
                           (b > a))
              (zip [(0 :: Int) ..] (zip vals (tail vals)))

frozenPartner :: TestTree
frozenPartner = testGroup "A frozen partner reduces to single-walker passage"
    [ testCase "N=48 k=1: encounter mean equals the exact MFPT" $
        agrees 48 1 0 17

    , testCase "N=48 k=2: encounter mean equals the exact MFPT" $
        agrees 48 2 0 17

    , testCase "N=24 k=1, antipodal: encounter mean equals the exact MFPT" $
        agrees 24 1 0 12

    , testCase "weights concentrate on the partner's site as qB falls" $ do
        let n = 12
            (_, !wSlow, _) = decompAt n 1 [] 0.8 0.005 (0, 6)
            (_, !wFast, _) = decompAt n 1 [] 0.8 0.500 (0, 6)
            !slow = wSlow V.! 6
            !fast = wFast V.! 6
        assertBool ("weight at the partner's site, qB=0.005: " ++ show slow)
                   (slow > 0.5)
        assertBool ("qB=0.005 gives " ++ show slow
                    ++ " and qB=0.5 gives " ++ show fast)
                   (slow > fast)
        assertBool ("weights sum " ++ show (V.sum wSlow))
                   (relClose 1e-10 1.0 (V.sum wSlow))
    ]
  where
    agrees n k n0 nt =
        mapM_ (\q ->
            let !got  = encounterMeanRingTwoQ q 0.0 n k 1.0 n0 nt
                !want = exactMFPTq q n k n0 nt
            in assertBool (report ("q=" ++ show q) want got)
                          (relClose 1e-10 want got))
        mobilities

ringStationaryPoint :: TestTree
ringStationaryPoint = testGroup "The bare ring: parity, not geometry"
    [ testCase "analytic derivative matches a central difference, d0=24" $
        matchesDiff 24

    , testCase "analytic derivative matches a central difference, d0=13" $
        matchesDiff 13

    , testCase "k=1 d0=24: the derivative is negative below qB=0.95" $
        mapM_ (\q ->
            let !d = ringMeanDeriv 0.9 q 48 1 24
            in assertBool ("qB=" ++ show q ++ " derivative " ++ show d) (d < 0))
        [0.05, 0.25, 0.50, 0.75, 0.90, 0.95]

    , testCase "k=1 d0=24: the derivative is positive at qB=1" $ do
        let !d = ringMeanDeriv 0.9 1.0 48 1 24
        assertBool ("derivative " ++ show d) (d > 0)

    , testCase "k=1 d0=24: the stationary point lies above qB=0.95" $ do
        let f x = ringMeanDeriv 0.9 x 48 1 24
            !qs = bisect f 0.95 1.0 50
        assertBool ("stationary point at " ++ show qs
                    ++ " with derivative " ++ show (f qs))
                   (qs > 0.95 && qs < 1.0 && abs (f qs) < 1e-3)

    , testCase "k=1 d0=24: the stationary point is a minimum" $ do
        let f x = ringMeanDeriv 0.9 x 48 1 24
            !qs = bisect f 0.95 1.0 50
        assertBool ("below " ++ show (f (qs - 1e-3))
                    ++ " above " ++ show (f (qs + 1e-3)))
                   (f (qs - 1e-3) < 0 && f (qs + 1e-3) > 0)

    , testCase "k=1 d0=13: the stationary point lies lower than at d0=24" $ do
        let f d0 x = ringMeanDeriv 0.9 x 48 1 d0
            !q13 = bisect (f 13) 0.05 1.0 50
            !q24 = bisect (f 24) 0.95 1.0 50
        assertBool ("d0=13 at " ++ show q13 ++ ", d0=24 at " ++ show q24)
                   (q13 < q24)

    , testCase "the spectral gap collapses only where parity permits it" $ do
        assertBool ("k=1 least eigenvalue " ++ show (leastNonLazyEig 48 1))
                   (relClose 1e-12 (negate 1) (leastNonLazyEig 48 1))
        assertBool ("k=2 least eigenvalue " ++ show (leastNonLazyEig 48 2))
                   (leastNonLazyEig 48 2 > negate 0.99)
        assertBool ("k=1 least gap at unit mobility "
                    ++ show (leastRelativeGap 1.0 1.0 48 1))
                   (leastRelativeGap 1.0 1.0 48 1 < 1e-12)
        assertBool ("k=2 least gap at unit mobility "
                    ++ show (leastRelativeGap 1.0 1.0 48 2))
                   (leastRelativeGap 1.0 1.0 48 2 > 1e-3)

    , testCase "range two admits no stationary point at any mobility" $
        mapM_ (\q ->
            let !d = ringMeanDeriv 0.9 q 48 2 24
            in assertBool ("qB=" ++ show q ++ " derivative " ++ show d) (d < 0))
        [0.05, 0.25, 0.50, 0.75, 0.95, 1.00]

    , testCase "a fixed mobility of one half admits no stationary point" $
        mapM_ (\q ->
            let !d = ringMeanDeriv 0.5 q 48 1 24
            in assertBool ("qB=" ++ show q ++ " derivative " ++ show d) (d <= 0))
        [0.05, 0.25, 0.50, 0.75, 0.95, 1.00]

    , testCase "a fixed mobility below one half admits no stationary point" $
        mapM_ (\q ->
            let !d = ringMeanDeriv 0.4 q 48 1 24
            in assertBool ("qB=" ++ show q ++ " derivative " ++ show d) (d < 0))
        [0.05, 0.25, 0.50, 0.75, 0.95, 1.00]
    ]
  where
    matchesDiff d0 =
        mapM_ (\q ->
            let !want = ringMeanDeriv 0.9 q 48 1 d0
                !got  = ringMeanDiff 0.9 q 48 1 d0 1.0e-5
            in assertBool (report ("qB=" ++ show q) want got)
                          (relClose 1.0e-5 want got))
        [0.05, 0.25, 0.50, 0.75, 0.95]

midpointExactness :: TestTree
midpointExactness = testGroup "The midpoint decomposition is exact at every step"
    [ testCase "bare ring N=12 k=1, every step" $ exact 12 1 [] (1, 7)
    , testCase "shortcut ring N=12 k=1, every step" $
        exact 12 1 [EdgeAdd 0 6] (1, 7)
    , testCase "range two N=12 k=2, every step" $ exact 12 2 [] (1, 7)
    ]
  where
    steps = [1.0e-2, 1.0e-3, 1.0e-4]
    exact n k prims starts =
        mapM_ (\h -> do
            let lo = decompAt n k prims 0.75 (0.5 - h) starts
                hi = decompAt n k prims 0.75 (0.5 + h) starts
                (!a, !b) = midpointAB (fineBlocks n) lo hi h
                (mLo, _, _) = lo
                (mHi, _, _) = hi
                !slope = (mHi - mLo) / (2 * h)
                !scale = max 1 (abs a + abs b)
            assertBool ("h=" ++ show h ++ ": A+B " ++ show (a + b)
                        ++ " against slope " ++ show slope
                        ++ " residual " ++ show (abs (a + b - slope)))
                       (abs (a + b - slope) <= 1e-9 * scale))
        steps

partitionDependence :: TestTree
partitionDependence = testGroup "The total is partition-free and A is not"
    [ testCase "the trivial partition gives A = 0 exactly" $ do
        let (!a, !b) = atBlocks (trivialBlock n)
            !scale = max 1 (abs b)
        assertBool ("A " ++ show a) (abs a <= 1e-9 * scale)

    , testCase "every partition reproduces the same total" $ do
        let (!a0, !b0) = atBlocks (trivialBlock n)
            (!a1, !b1) = atBlocks (halfBlocks n)
            (!a2, !b2) = atBlocks (fineBlocks n)
        assertBool (report "half against trivial" (a0 + b0) (a1 + b1))
                   (relClose 1e-9 (a0 + b0) (a1 + b1))
        assertBool (report "fine against trivial" (a0 + b0) (a2 + b2))
                   (relClose 1e-9 (a0 + b0) (a2 + b2))

    , testCase "the between-site term differs between partitions" $ do
        let (!aFine, !bFine) = atBlocks (fineBlocks n)
            (!aHalf, _)      = atBlocks (halfBlocks n)
            !scale = max 1 (abs aFine + abs bFine)
        assertBool ("fine A " ++ show aFine ++ " half A " ++ show aHalf)
                   (abs (aFine - aHalf) > 1e-6 * scale)

    , testCase "the fine partition carries a non-vanishing between-site term" $ do
        let (!aFine, !bFine) = atBlocks (fineBlocks n)
            !scale = max 1 (abs aFine + abs bFine)
        assertBool ("A " ++ show aFine ++ " B " ++ show bFine)
                   (abs aFine > 1e-6 * scale)
    ]
  where
    n = 12
    h = 1.0e-3
    atBlocks blocks =
        let lo = decompAt n 1 [EdgeAdd 0 6] 0.75 (0.5 - h) (1, 7)
            hi = decompAt n 1 [EdgeAdd 0 6] 0.75 (0.5 + h) (1, 7)
        in midpointAB blocks lo hi h
