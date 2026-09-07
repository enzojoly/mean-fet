{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Encounter
-- Description : Two-walker first-encounter statistics on a shared environment.
--               One staged builder serves the generating function, the mean and
--               the splitting probabilities; the resolvent tensor is split into
--               real GEMMs; and the translation-invariant ring carries a closed
--               two-mobility relative-walk engine. The per contact site
--               decomposition is taken from an expansion about the stationary
--               pole rather than by evaluation near it, because the renewal
--               system is genuinely singular at the pole and its components are
--               fixed by the first and second order terms rather than by the
--               limit.

module Encounter
    ( encounterGF
    , encounterPMF
    , encounterGFTwoQ
    , encounterPMFTwoQ
    , encounterPMFTwoQAcc
    , encounterMeanTwoQ
    , splittingProbs
    , splittingProbsTwoQ
    , encounterSiteMeansTwoQ
    , encounterDecompositionTwoQ
    , encounterSitePMFTwoQAcc
    , encounterGFRingTwoQ
    , encounterPMFRingTwoQ
    , encounterMeanRingTwoQ
    , degenerateRingCell
    , degeneratePeriodicCell
    , countPeaks
    ) where

import Data.Complex (Complex(..))
import qualified Data.Vector.Unboxed as V
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra ((><), (#>))

import Types (R, Matrix)
import Ring (ringEigQ, stateCount)
import Distribution (invertPGF, invertPGFAcc, invertPGFVec, countPeaks)

data Stage = Stage
    { stN     :: !Int
    , stEigsA :: !(V.Vector Double)
    , stEigsB :: !(V.Vector Double)
    , stWA    :: !Matrix
    , stWB    :: !Matrix
    , stInitA :: !Matrix
    , stInitB :: !Matrix
    , stPi    :: !(V.Vector Double)
    , stRho   :: !Double
    , stSA    :: !Int
    , stSB    :: !Int
    , stOnes  :: !(LA.Vector Double)
    }

mkStage :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
        -> (Int, Int) -> Stage
mkStage eigsA evecsA eigsB evecsB piVec rho (startA, startB) =
    let !n      = length eigsA
        !eigsAV = V.fromList eigsA
        !eigsBV = V.fromList eigsB

        !wMatA = ((n * n) >< n)
            [ evecsA `LA.atIndex` (k, j) * evecsA `LA.atIndex` (l, j)
            | k <- [0 .. n - 1], l <- [0 .. n - 1], j <- [0 .. n - 1]
            ]

        !wMatB = ((n * n) >< n)
            [ evecsB `LA.atIndex` (k, j) * evecsB `LA.atIndex` (l, j)
            | k <- [0 .. n - 1], l <- [0 .. n - 1], j <- [0 .. n - 1]
            ]

        !initA = (n >< n)
            [ evecsA `LA.atIndex` (startA, j) * evecsA `LA.atIndex` (c, j)
            | c <- [0 .. n - 1], j <- [0 .. n - 1]
            ]

        !initB = (n >< n)
            [ evecsB `LA.atIndex` (startB, j) * evecsB `LA.atIndex` (c, j)
            | c <- [0 .. n - 1], j <- [0 .. n - 1]
            ]

        !onesN = LA.konst 1 n :: LA.Vector Double

    in Stage n eigsAV eigsBV wMatA wMatB initA initB
             (V.fromList piVec) rho startA startB onesN

crossOf :: Stage -> Matrix -> (LA.Vector Double, LA.Vector Double)
crossOf st !s =
    ( ((stWA st LA.<> s) * stWB st) #> stOnes st
    , ((stInitA st LA.<> s) * stInitB st) #> stOnes st
    )

-- | Renewal system assembled from contracted pair propagators, in the complex
-- arithmetic the contour evaluation requires.
buildRenewal :: Int -> V.Vector Double -> Double
             -> (Int -> Complex Double)
             -> (Int -> Int -> Complex Double)
             -> (Int -> Complex Double)
             -> Int -> Int
             -> (LA.Matrix (Complex Double), LA.Vector (Complex Double))
buildRenewal !n !piV !rho selfAt crossAt initAt !sA !sB =
    -- The eigensystem is that of the symmetrised operator, so the propagator in
    -- the original basis carries a factor sqrt(pi_j / pi_i) from i to j. On the
    -- diagonal contact set both walkers travel between the same pair of sites,
    -- so the two square roots combine to a whole power; the source vector runs
    -- from two distinct sites and carries their geometric mean instead. Both
    -- factors are unity on a lattice of equal degrees, which is why a decorated
    -- lattice is the only place they can be checked.
    let piRatio k l =
            ((piV `V.unsafeIndex` l) / (piV `V.unsafeIndex` k)) :+ 0

        !fM = (n >< n)
            [ if l == k
              then (((1 - rho) / rho) :+ 0) / selfAt k + 1
              else piRatio k l * crossAt k l / selfAt l
            | l <- [0 .. n - 1], k <- [0 .. n - 1]
            ]

        !piCorrA = sqrt (piV `V.unsafeIndex` sA)
        !piCorrB = sqrt (piV `V.unsafeIndex` sB)
        piC c = ((piV `V.unsafeIndex` c) / (piCorrA * piCorrB)) :+ 0

        !bV = LA.fromList [ piC c * initAt c / selfAt c | c <- [0 .. n - 1] ]

    in (fM, bV)

-- | Solution vector of the renewal system at a point on the contour. Its
-- components are the individual splitting generating functions, one per
-- contact site; their sum is the encounter generating function. Nothing sums
-- them here, so a caller wanting the per-site law pays no additional solve.
solveAtZ :: Stage -> Complex Double -> LA.Vector (Complex Double)
solveAtZ st (!zx :+ !zy) =
    let !n  = stN st
        !eA = stEigsA st
        !eB = stEigsB st

        sigPart f = (n >< n)
            [ let !g = (eA `V.unsafeIndex` l) * (eB `V.unsafeIndex` m)
                  !a = 1 - zx * g
                  !b = zy * g
                  !den = a * a + b * b
              in f a b den
            | l <- [0 .. n - 1], m <- [0 .. n - 1] ]

        !sr = sigPart (\a _ den -> a / den)
        !si = sigPart (\_ b den -> b / den)

        (!crossR, !initR) = crossOf st sr
        (!crossI, !initI) = crossOf st si

        !crossC = LA.toComplex (crossR, crossI)
        !initC  = LA.toComplex (initR, initI)

        selfAt k    = crossC `LA.atIndex` (k * n + k)
        crossAt k l = crossC `LA.atIndex` (k * n + l)
        initAt c    = initC  `LA.atIndex` c

        (!fMat, !bVec) = buildRenewal n (stPi st) (stRho st)
                             selfAt crossAt initAt (stSA st) (stSB st)

    in fMat LA.<\> bVec

evalAtZ :: Stage -> Complex Double -> Complex Double
evalAtZ st z = LA.sumElements (solveAtZ st z)

-- | Slope of the encounter generating function at z = 1 - delta, together with
-- the per-contact-site solution there. The renewal system is differentiated
-- analytically -- F x' = b' - F' x, reusing the factorisation of F -- so
-- nothing here is a finite difference; delta fixes only the point at which the
-- exact derivative is evaluated.
--
-- The limit cannot be taken at z = 1 itself. Every entry of the renewal matrix
-- tends to the same value there, so F degenerates to a rank-one matrix and the
-- solve is singular; its condition number grows as 1/delta on approach. The
-- summed quantities returned here lie in the well-conditioned direction of
-- that degeneracy and survive it. The individual components do not, and are
-- taken from expandAtPole instead.
slopeVecAndValues :: Stage -> Double -> (LA.Vector Double, LA.Vector Double)
slopeVecAndValues st !delta =
    let !n  = stN st
        !eA = stEigsA st
        !eB = stEigsB st
        !piV = stPi st
        !rho = stRho st
        !x = 1 - delta

        !sM = (n >< n)
            [ let !g = (eA `V.unsafeIndex` l) * (eB `V.unsafeIndex` m)
              in 1 / (1 - x * g)
            | l <- [0 .. n - 1], m <- [0 .. n - 1] ]

        !sM' = (n >< n)
            [ let !g = (eA `V.unsafeIndex` l) * (eB `V.unsafeIndex` m)
                  !d = 1 - x * g
              in g / (d * d)
            | l <- [0 .. n - 1], m <- [0 .. n - 1] ]

        (!cross,  !initV)  = crossOf st sM
        (!cross', !initV') = crossOf st sM'

        at :: LA.Vector Double -> Int -> Double
        at v i = v `LA.atIndex` i

        sAt  k = at cross  (k * n + k)
        sAt' k = at cross' (k * n + k)
        cAt  k l = at cross  (k * n + l)
        cAt' k l = at cross' (k * n + l)

        piRatio k l = (piV `V.unsafeIndex` l) / (piV `V.unsafeIndex` k)

        !fM = (n >< n)
            [ if l == k
              then (1 - rho) / (rho * sAt k) + 1
              else piRatio k l * cAt k l / sAt l
            | l <- [0 .. n - 1], k <- [0 .. n - 1]
            ] :: Matrix

        !fM' = (n >< n)
            [ if l == k
              then negate (1 - rho) / rho * sAt' k / (sAt k * sAt k)
              else piRatio k l * (cAt' k l * sAt l - cAt k l * sAt' l)
                   / (sAt l * sAt l)
            | l <- [0 .. n - 1], k <- [0 .. n - 1]
            ] :: Matrix

        !piCorrA = sqrt (piV `V.unsafeIndex` stSA st)
        !piCorrB = sqrt (piV `V.unsafeIndex` stSB st)
        piC c = (piV `V.unsafeIndex` c) / (piCorrA * piCorrB)

        !bV = LA.fromList
            [ piC c * at initV c / sAt c | c <- [0 .. n - 1] ]

        !bV' = LA.fromList
            [ piC c * (at initV' c * sAt c - at initV c * sAt' c)
              / (sAt c * sAt c)
            | c <- [0 .. n - 1] ]

        !lu  = LA.luPacked fM
        !xV  = LA.flatten (LA.luSolve lu (LA.asColumn bV))
        !rhs = bV' - (fM' #> xV)
        !xV' = LA.flatten (LA.luSolve lu (LA.asColumn rhs))

    in (xV', xV)

slopeAndValues :: Stage -> Double -> (Double, LA.Vector Double)
slopeAndValues st !delta =
    let (!s, !v) = slopeVecAndValues st delta
    in (LA.sumElements s, v)

-- | Evaluation point for the derivative, balancing truncation against the
-- conditioning of the renewal solve. Not a tuning parameter: the optimum is
-- flat to within a factor of a few, and the residual is of order 1e-9
-- relative across the configurations exercised by the reference assertions.
meanDelta :: Double
meanDelta = 1.0e-7

-- | Evaluation point for the splitting probabilities. These are values rather
-- than derivatives, so truncation dominates and a looser point is preferable
-- for conditioning.
splitDelta :: Double
splitDelta = 1.0e-6

-- | Expansion of the renewal system about the stationary pole.
--
-- The system is genuinely singular at z = 1, not merely awkward there. The
-- stationary term of the resolvent is the only singular one, and dividing it
-- out of numerator and denominator alike leaves the renewal matrix tending to
-- the matrix of ones and the source vector to the vector of ones. That system
-- does not determine the per-site splitting probabilities at all: every vector
-- summing to one solves it. Writing w = 1 - z and expanding,
--
--     F(w) = 1.1^T + w F1 + w^2 F2 + ...,   b(w) = 1 + w b1 + w^2 b2 + ...,
--     x(w) = x0 + w x1 + ...,
--
-- the order zero relation fixes only the sum of x0. Orders one and two give
-- two bordered systems sharing a single matrix,
--
--     [ F1  1 ] [ x0     ]   [ b1         ]      [ F1  1 ] [ x1     ]   [ b2 - F2 x0 ]
--     [ 1^T 0 ] [ sigma1 ] = [ 1          ] ,    [ 1^T 0 ] [ sigma2 ] = [ sigma1     ]
--
-- whose solutions are the splitting probabilities x0 and the per-site
-- contributions to the mean, -x1. Nothing is evaluated near the pole and
-- nothing is extrapolated, so the components carry the conditioning of the
-- bordered matrix alone, which is of order 1e6 and independent of any step.
-- Measured against a mirror-symmetric release the components are symmetric to
-- 1e-13, where evaluation at 1e-7 followed by extrapolation gives 1e-3.
--
-- The Taylor coefficients of each entry follow from series division,
-- S0 = P0/Q0, S1 = (P1 - S0 Q1)/Q0, S2 = (P2 - S0 Q2 - S1 Q1)/Q0, applied to
-- the ratio of residue-shifted numerator and denominator. Only the value and
-- first derivative of the analytic remainder are needed, so three contractions
-- serve where the evaluation route needed fourteen.
expandAtPole :: Stage -> (Double, V.Vector Double, V.Vector Double)
expandAtPole st
    | not stationaryPresent =
        error "expandAtPole: the cell is not measure-preserving, so there is \
              \no stationary pole to expand about."
    | doubled =
        error "expandAtPole: the relative spectrum touches unity away from the \
              \stationary mode. The renewal system carries a repeated \
              \singularity and the mean does not exist."
    | otherwise = (mean, weights, siteMeans)
  where
    !n   = stN st
    !eA  = stEigsA st
    !eB  = stEigsB st
    !piV = stPi st
    !rho = stRho st
    !ia  = V.maxIndex eA
    !ib  = V.maxIndex eB
    !gam = (1 - rho) / rho

    !stationaryPresent =
        abs (eA `V.unsafeIndex` ia - 1) < 1e-9
        && abs (eB `V.unsafeIndex` ib - 1) < 1e-9

    !doubled =
        or [ abs (1 - (eA `V.unsafeIndex` l) * (eB `V.unsafeIndex` m)) < 1e-10
           | l <- [0 .. n - 1], m <- [0 .. n - 1]
           , not (l == ia && m == ib) ]

    tensor :: (Double -> Double) -> Matrix
    tensor f = (n >< n)
        [ if l == ia && m == ib
          then 0
          else f ((eA `V.unsafeIndex` l) * (eB `V.unsafeIndex` m))
        | l <- [0 .. n - 1], m <- [0 .. n - 1] ]

    !residue = (n >< n)
        [ if l == ia && m == ib then 1 else 0
        | l <- [0 .. n - 1], m <- [0 .. n - 1] ] :: Matrix

    (!rC,  !rI)  = crossOf st residue
    (!t0C, !t0I) = crossOf st (tensor (\g -> 1 / (1 - g)))
    (!t1C, !t1I) = crossOf st
        (tensor (\g -> negate g / ((1 - g) * (1 - g))))

    at :: LA.Vector Double -> Int -> Double
    at v i = v `LA.atIndex` i

    rAt  k l = at rC  (k * n + l)
    t0At k l = at t0C (k * n + l)
    t1At k l = at t1C (k * n + l)

    quot3 :: (Double, Double, Double) -> (Double, Double, Double)
          -> (Double, Double, Double)
    quot3 (!p0, !p1, !p2) (!q0, !q1, !q2) =
        let !s0 = p0 / q0
            !s1 = (p1 - s0 * q1) / q0
            !s2 = (p2 - s0 * q2 - s1 * q1) / q0
        in (s0, s1, s2)

    denomAt l = (rAt l l, t0At l l, t1At l l)

    piRatio k l = (piV `V.unsafeIndex` l) / (piV `V.unsafeIndex` k)

    entry :: Int -> Int -> (Double, Double)
    entry l k
        | l == k =
            let (_, !s1, !s2) = quot3 (0, gam, 0) (denomAt l)
            in (s1, s2)
        | otherwise =
            let (_, !s1, !s2) =
                    quot3 (rAt k l, t0At k l, t1At k l) (denomAt l)
                !r = piRatio k l
            in (r * s1, r * s2)

    !fM1 = (n >< n)
        [ fst (entry l k) | l <- [0 .. n - 1], k <- [0 .. n - 1] ] :: Matrix

    !fM2 = (n >< n)
        [ snd (entry l k) | l <- [0 .. n - 1], k <- [0 .. n - 1] ] :: Matrix

    !piCorrA = sqrt (piV `V.unsafeIndex` stSA st)
    !piCorrB = sqrt (piV `V.unsafeIndex` stSB st)
    piC c = (piV `V.unsafeIndex` c) / (piCorrA * piCorrB)

    srcAt c =
        let (_, !s1, !s2) = quot3 (at rI c, at t0I c, at t1I c) (denomAt c)
        in (piC c * s1, piC c * s2)

    !bV1 = LA.fromList [ fst (srcAt c) | c <- [0 .. n - 1] ]
    !bV2 = LA.fromList [ snd (srcAt c) | c <- [0 .. n - 1] ]

    !bordered = ((n + 1) >< (n + 1))
        [ if r < n && c < n then fM1 `LA.atIndex` (r, c)
          else if r == n && c == n then 0
          else 1
        | r <- [0 .. n], c <- [0 .. n] ] :: Matrix

    !lu = LA.luPacked bordered

    solveBordered :: LA.Vector Double -> Double -> (LA.Vector Double, Double)
    solveBordered rhsTop !lastEntry =
        let !rhs = LA.vjoin [rhsTop, LA.fromList [lastEntry]]
            !sol = LA.flatten (LA.luSolve lu (LA.asColumn rhs))
        in (LA.subVector 0 n sol, sol `LA.atIndex` n)

    (!x0, !sigma1) = solveBordered bV1 1
    (!x1, _)       = solveBordered (bV2 - (fM2 #> x0)) sigma1

    !weights   = V.generate n $ \c -> x0 `LA.atIndex` c
    !siteMeans = V.generate n $ \c -> negate (x1 `LA.atIndex` c)
    !mean      = negate (LA.sumElements x1)

encounterGFTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
                -> (Int, Int) -> Complex Double -> Complex Double
encounterGFTwoQ eigsA evecsA eigsB evecsB piVec rho starts =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
    in evalAtZ st

encounterPMFTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
                 -> (Int, Int) -> Int -> V.Vector Double
encounterPMFTwoQ eigsA evecsA eigsB evecsB piVec rho starts tmax =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
    in invertPGF tmax (evalAtZ st)

encounterPMFTwoQAcc :: Double -> [R] -> Matrix -> [R] -> Matrix -> [Double]
                    -> Double -> (Int, Int) -> Int -> V.Vector Double
encounterPMFTwoQAcc acc eigsA evecsA eigsB evecsB piVec rho starts tmax =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
    in invertPGFAcc acc tmax (evalAtZ st)

encounterMeanTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
                  -> (Int, Int) -> Double
encounterMeanTwoQ eigsA evecsA eigsB evecsB piVec rho starts =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
        (!mD, _) = slopeAndValues st meanDelta
        (!mH, _) = slopeAndValues st (meanDelta / 2)
    in 2 * mH - mD

encounterGF :: [R] -> Matrix -> [Double] -> Double
            -> (Int, Int) -> Complex Double -> Complex Double
encounterGF eigs evecs piVec rho starts z =
    encounterGFTwoQ eigs evecs eigs evecs piVec rho starts z

encounterPMF :: [R] -> Matrix -> [Double] -> Double
             -> (Int, Int) -> Int -> V.Vector Double
encounterPMF eigs evecs piVec rho starts tmax =
    encounterPMFTwoQ eigs evecs eigs evecs piVec rho starts tmax

splittingProbsTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
                   -> (Int, Int) -> V.Vector Double
splittingProbsTwoQ eigsA evecsA eigsB evecsB piVec rho starts =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
        (_, !sD) = slopeAndValues st splitDelta
        (_, !sH) = slopeAndValues st (splitDelta / 2)
        !n = stN st
    in V.generate n $ \c ->
        max 0 (2 * (sH `LA.atIndex` c) - (sD `LA.atIndex` c))

splittingProbs :: [R] -> Matrix -> [Double] -> Double
              -> (Int, Int) -> V.Vector Double
splittingProbs eigs evecs piVec rho starts =
    splittingProbsTwoQ eigs evecs eigs evecs piVec rho starts

-- | Contribution of each contact site to the mean, that is the splitting
-- probability there multiplied by the mean encounter time conditional on
-- encounter occurring there. Taken from the second order term of the pole
-- expansion, so the components are exact rather than extrapolated and nothing
-- is clipped: a component returned negative would be a statement about the
-- conditioning of the bordered solve and must not be hidden.
encounterSiteMeansTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double] -> Double
                       -> (Int, Int) -> V.Vector Double
encounterSiteMeansTwoQ eigsA evecsA eigsB evecsB piVec rho starts =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
        (_, _, !sm) = expandAtPole st
    in sm

-- | The mean, the splitting probabilities and the per-site contributions to
-- the mean, from one expansion about the stationary pole. The three are
-- mutually consistent by construction: the weights sum to one and the site
-- contributions sum to the mean, both to machine precision rather than to the
-- accuracy of an evaluation point.
encounterDecompositionTwoQ :: [R] -> Matrix -> [R] -> Matrix -> [Double]
                           -> Double -> (Int, Int)
                           -> (Double, V.Vector Double, V.Vector Double)
encounterDecompositionTwoQ eigsA evecsA eigsB evecsB piVec rho starts =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
    in expandAtPole st

-- | The encounter law resolved by contact site, in the time domain. One
-- distribution per site, summing to the encounter distribution at every step.
-- The solve at each contour node is that of the summed route, so this costs
-- one transform per site and no additional solves.
encounterSitePMFTwoQAcc :: Double -> [R] -> Matrix -> [R] -> Matrix -> [Double]
                        -> Double -> (Int, Int) -> Int -> [V.Vector Double]
encounterSitePMFTwoQAcc acc eigsA evecsA eigsB evecsB piVec rho starts tmax =
    let !st = mkStage eigsA evecsA eigsB evecsB piVec rho starts
    in invertPGFVec acc tmax (stN st) (solveAtZ st)

foldSep :: Int -> Int -> Int
foldSep bigN d =
    let !d' = ((d `mod` bigN) + bigN) `mod` bigN
    in min d' (bigN - d')

encounterGFRingTwoQ :: Double -> Double -> Int -> Int -> Double -> Int
                    -> Complex Double -> Complex Double
encounterGFRingTwoQ qA qB bigN k rho d0 =
    let !nd   = fromIntegral bigN :: Double
        !mus  = V.generate bigN $ \l ->
                    ringEigQ qA bigN k l * ringEigQ qB bigN k l
        !cosD = V.generate bigN $ \l ->
                    cos (2 * pi * fromIntegral (l * d0) / nd)
        !rhoC = rho :+ 0
        !oneMinusRho = (1 - rho) :+ 0
        propD z = V.sum (V.zipWith
            (\mu c -> (c :+ 0) / (1 - z * (mu :+ 0))) mus cosD)
            / (nd :+ 0)
        prop0 z = V.sum (V.map
            (\mu -> 1 / (1 - z * (mu :+ 0))) mus)
            / (nd :+ 0)
    in \z -> rhoC * propD z / (oneMinusRho + rhoC * prop0 z)

encounterPMFRingTwoQ :: Double -> Double -> Int -> Int -> Double
                     -> Int -> Int -> Int -> V.Vector Double
encounterPMFRingTwoQ qA qB bigN k rho sA sB tmax =
    invertPGF tmax (encounterGFRingTwoQ qA qB bigN k rho (foldSep bigN (sB - sA)))

-- | True when the relative spectrum touches unity away from the stationary
-- mode, at which point the closed form diverges and encounter from an odd
-- separation is forbidden outright. This is the non-lazy nearest-neighbour
-- corner on an even periodic cell: 1 - mu = e*(S - P*e) vanishes at e = 2 only
-- when both mobilities are one, k = 1 and the side length is even, giving the
-- mode l = N/2. In two dimensions the axis eigenvalues must reach minus one
-- together, which is the same condition on the side length, so the predicate
-- serves the torus unchanged. A reflecting cell is exempt: the wall gives its
-- boundary sites a holding probability at unit mobility, which breaks the
-- parity alternation and keeps minus one out of the spectrum. The generic
-- route needs no predicate at all, expandAtPole inspecting the spectrum
-- directly.
degeneratePeriodicCell :: Double -> Double -> Int -> Int -> Bool
degeneratePeriodicCell qA qB size k =
    qA >= 1 && qB >= 1 && k == 1 && even size

degenerateRingCell :: Double -> Double -> Int -> Int -> Bool
degenerateRingCell = degeneratePeriodicCell

-- | Mean encounter time on the homogeneous ring, in closed form. The
-- stationary mode is removed analytically, leaving a finite sum over the
-- decaying spectrum plus the Kac correction for imperfect absorption, in
-- which the return cost is the number of states.
encounterMeanRingTwoQ :: Double -> Double -> Int -> Int -> Double
                      -> Int -> Int -> Double
encounterMeanRingTwoQ qA qB bigN k rho sA sB
    | degenerateRingCell qA qB bigN k =
        error "encounterMeanRingTwoQ: degenerate cell (qA = qB = 1, k = 1, \
              \N even) places a unit eigenvalue off the stationary mode; the \
              \mean does not exist. Use q < 1."
    | otherwise =
        let !nd  = fromIntegral bigN :: Double
            !d0  = foldSep bigN (sB - sA)
            !bulk = sum
                [ let !mu = ringEigQ qA bigN k l * ringEigQ qB bigN k l
                      !c  = cos (2 * pi * fromIntegral (l * d0) / nd)
                  in (1 - c) / (1 - mu)
                | l <- [1 .. bigN - 1] ]
            !kac = fromIntegral (stateCount bigN) * (1 - rho) / rho
        in bulk + kac
