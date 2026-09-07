{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Passage
-- Description : First passage, first return, and partial absorption generating
--               functions on the defected ring. All homogeneous propagator
--               lookups go through staged ring tables built once per
--               configuration, defects are consumed as a boxed vector, and the
--               modified determinant ratios are evaluated through a single LU
--               factorisation by the matrix determinant lemma.

module Passage
    ( firstPassageGF
    , firstAbsorptionGF
    , firstReturnGF
    , mkEigensystemGF
    , pureRingFPGF
    , exactMFPT
    , exactMFPTq
    , meanFromGF
    , gfMeanDelta
    ) where

import Data.Complex (Complex(..))
import Data.List (foldl')
import qualified Data.Vector as BV
import qualified Data.Vector.Unboxed as V
import qualified Numeric.LinearAlgebra as LA

import Types (R, Pos, Matrix, matrixGet, matrixSize)
import Ring (RingTables, mkRingTables, ringDenomAt, ringPropAt, ringQn)
import Defect (DefectEntry(..), buildHWith)

mkEigensystemGF :: [R] -> Matrix -> Pos -> Pos -> R -> R
                -> (Complex Double -> Complex Double)
mkEigensystemGF eigs evecs src tgt rho piRatio = gf
  where
    !n       = length eigs
    !eigV    = V.fromList eigs
    !cSrcTgt = V.generate n (\k -> matrixGet evecs src k * matrixGet evecs tgt k)
    !cTgtTgt = V.generate n (\k -> matrixGet evecs tgt k ^ (2 :: Int))
    !piR     = piRatio :+ 0

    spectralSum :: V.Vector Double -> Complex Double -> Complex Double
    spectralSum coeffs z = V.ifoldl' step 0 coeffs
      where
        step !acc !k !c =
            let !lam = eigV `V.unsafeIndex` k
            in acc + (c :+ 0) / (1 - z * (lam :+ 0))
    {-# INLINE spectralSum #-}

    gf :: Complex Double -> Complex Double
    gf z
        | rho >= 1.0 - 1e-15 =
            piR * spectralSum cSrcTgt z / spectralSum cTgtTgt z
        | otherwise =
            let !qST  = spectralSum cSrcTgt z
                !qTT  = spectralSum cTgtTgt z
                !rhoC = rho :+ 0
            in piR * rhoC * qST / (((1 - rho) :+ 0) + rhoC * qTT)

firstPassageGF :: Double -> Int -> Int -> [DefectEntry] -> Int -> Int
               -> Complex Double -> Complex Double
firstPassageGF q bigN k defects n0 n =
    let !tables = mkRingTables q bigN k
        !dvec   = BV.fromList defects
    in \z -> let (!sN0N, !sNN) = sPairAt tables dvec n0 n z
             in sN0N / sNN

firstAbsorptionGF :: Double -> Int -> Int -> [DefectEntry] -> Int -> Int
                  -> Double -> Complex Double -> Complex Double
firstAbsorptionGF q bigN k defects n0 n rho
    | rho >= 1.0 - 1e-15 = firstPassageGF q bigN k defects n0 n
    | otherwise =
        let !tables = mkRingTables q bigN k
            !dvec   = BV.fromList defects
            !rhoC   = rho :+ 0
            !oneMinusRho = (1 - rho) :+ 0
        in \z -> let (!sN0N, !sNN) = sPairAt tables dvec n0 n z
                 in rhoC * sN0N / (oneMinusRho + rhoC * sNN)

firstReturnGF :: Double -> Int -> Int -> [DefectEntry] -> Int
              -> Complex Double -> Complex Double
firstReturnGF q bigN k defects n =
    let !tables = mkRingTables q bigN k
        !dvec   = BV.fromList defects
    in \z -> let (_, !sNN) = sPairAt tables dvec n n z
             in 1 - recip sNN

pureRingFPGF :: Double -> Int -> Int -> Int -> Int
             -> Complex Double -> Complex Double
pureRingFPGF q bigN k n0 n z =
    ringQn q bigN k n0 n z / ringQn q bigN k n n z

sPairAt :: RingTables -> BV.Vector DefectEntry -> Int -> Int
        -> Complex Double -> (Complex Double, Complex Double)
sPairAt !tables !dvec !n0 !n !z =
    let !denom = ringDenomAt tables z
        qAt i j = ringPropAt tables denom i j
        !qN0N = qAt n0 n
        !qNN  = qAt n n
    in case BV.length dvec of
        0 -> (qN0N, qNN)
        1 ->
            let !d = dvec `BV.unsafeIndex` 0
                !eta = defXvu d :+ 0
                !u = defU d
                !v = defV d
                !h11 = eta * (qAt u u + qAt v v - qAt u v - qAt v u) - recip z
                !dqN  = qAt u n - qAt v n
                !brN0 = qAt n0 u - qAt n0 v
                !brNN = qAt n u - qAt n v
                !hmod1 = h11 - dqN * (eta * brN0)
                !hmod2 = h11 - dqN * (eta * brNN)
            in (qN0N - 1 + hmod1 / h11, qNN - 1 + hmod2 / h11)
        m ->
            let !hMat = buildHWith dvec qAt z
                !lu   = LA.luPacked hMat
                !vN   = LA.fromList
                    [ let !dj = dvec `BV.unsafeIndex` j
                      in qAt (defU dj) n - qAt (defV dj) n
                    | j <- [0 .. m - 1] ]
                bracketAt !x !i =
                    let !di = dvec `BV.unsafeIndex` i
                        !euv = defXuv di :+ 0
                        !evu = defXvu di :+ 0
                    in evu * qAt x (defU di) - euv * qAt x (defV di)
                !u1 = LA.fromList [ bracketAt n0 i | i <- [0 .. m - 1] ]
                !u2 = LA.fromList [ bracketAt n  i | i <- [0 .. m - 1] ]
                !sols = LA.luSolve lu (LA.fromColumns [u1, u2])
                !x1 = head (LA.toColumns sols)
                !x2 = LA.toColumns sols !! 1
                !r1 = 1 - LA.sumElements (vN * x1)
                !r2 = 1 - LA.sumElements (vN * x2)
            in (qN0N - 1 + r1, qNN - 1 + r2)

exactMFPT :: Int -> Int -> Int -> Int -> Double
exactMFPT !bigN !k !n0 !n = kd * foldl' (+) 0 (map term [1 .. bigN - 1])
  where
    !kd = fromIntegral k :: Double
    !nd = fromIntegral bigN :: Double
    !delta = n - n0
    term ell =
      let !l   = fromIntegral ell
          !num = (1 - cos (2 * pi * l * fromIntegral delta / nd))
               * sin (l * pi / nd)
          !den = kd * sin (l * pi / nd)
               - sin (kd * l * pi / nd) * cos ((kd + 1) * l * pi / nd)
      in if abs den < 1e-30 then 0 else num / den

exactMFPTq :: Double -> Int -> Int -> Int -> Int -> Double
exactMFPTq q bigN k n0 n
    | q <= 0    = error "exactMFPTq: q must be > 0"
    | otherwise = exactMFPT bigN k n0 n / q

-- | Evaluation point for the generating-function derivative on the real axis.
-- The same balance as elsewhere: the one-sided difference carries an error
-- linear in the step, removed to second order by extrapolation, against the
-- cancellation incurred by differencing values that both approach unity.
gfMeanDelta :: Double
gfMeanDelta = 1.0e-6

-- | Mean arrival time read directly from a first-passage or first-absorption
-- generating function, without inverting it. On a finite recurrent lattice the
-- generating function equals one at the boundary of the unit disc, so the mean
-- is the limit of [1 - F(1 - d)] / d as d vanishes. Differencing against the
-- known boundary value rather than against a second evaluation halves the
-- cancellation, and the leading error is linear in d and so removed by
-- extrapolating two steps.
--
-- This is the whole cost of a mean: two evaluations of the generating
-- function, against the thousands an inversion requires. Where only the mean
-- is wanted the distribution need never be formed.
meanFromGF :: (Complex Double -> Complex Double) -> Double
meanFromGF gf =
    let atZ !d = let (!re, _) = decompose (gf ((1 - d) :+ 0))
                 in (1 - re) / d
        !gD = atZ gfMeanDelta
        !gH = atZ (gfMeanDelta / 2)
    in 2 * gH - gD
  where
    decompose (a :+ b) = (a, b)
