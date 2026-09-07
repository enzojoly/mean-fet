{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Estimate
-- Description : Recovering an unbiased mean from a distribution observed only
--               up to a finite horizon. A sum truncated at t_max omits the
--               slowest arrivals, and although the omitted mass is small the
--               times it carries are large, so the deficit is a lever rather
--               than a rounding. The tail of a rational generating function is
--               asymptotically geometric, which makes the omitted contribution
--               a computable quantity rather than an unknown one, and the same
--               relation turns a mean taken by one route into a check on a mean
--               taken by another.

module Estimate
    ( truncatedSum
    , conditionalMean
    , censoredFraction
    , tailScale
    , correctMean
    , correctedFromPMF
    , selfConsistency
    , correctedEstimate
    , extrapolatedEstimate
    , sampledEstimate
    ) where

import qualified Data.Vector.Unboxed as V

import Types (Estimate(..), Provenance(..))

-- | Raw first moment over the tabulated range. This is a biased estimate of
-- the mean, low by the whole of the omitted tail, and additionally deflated by
-- the missing mass; it is recorded here only so that the bias it carries can
-- be named and corrected rather than inherited silently.
truncatedSum :: V.Vector Double -> Double
truncatedSum pmf = V.sum (V.imap (\t p -> fromIntegral t * p) pmf)

-- | Mean conditioned on arrival within the horizon. Renormalising by the
-- recovered mass removes the deflation but not the omission: the surviving
-- trajectories are precisely the slow ones, so this remains biased low.
conditionalMean :: V.Vector Double -> Double
conditionalMean pmf =
    let !m = V.sum pmf
    in if m <= 0 then 0 else truncatedSum pmf / m

censoredFraction :: V.Vector Double -> Double
censoredFraction pmf = max 0 (1 - V.sum pmf)

-- | Decay scale of the geometric tail, fitted by least squares to the
-- logarithm of the density over the closing stretch of the horizon. The
-- coefficients of a rational generating function are a finite sum of geometric
-- modes, so far enough out the slowest surviving pole dominates and the
-- logarithm is linear in time; the reciprocal slope is the mean residual life
-- of a trajectory that has not yet arrived.
--
-- The window opens at three fifths of the horizon, far enough for the
-- subdominant modes to have decayed and close enough to retain samples above
-- the noise floor. Where too few usable samples remain the scale falls back to
-- the conditional mean, which is the correct order of magnitude and errs
-- towards under-correction.
tailScale :: V.Vector Double -> Double
tailScale pmf
    | n < 8          = fallback
    | length pts < 4 = fallback
    | slope >= 0     = fallback
    | otherwise      = negate (1 / slope)
  where
    !n = V.length pmf
    !fallback = max 1 (conditionalMean pmf)
    !peak = if n == 0 then 0 else V.maximum pmf
    !floorV = peak * 1e-13
    !lo = (3 * n) `div` 5
    !pts = [ (fromIntegral t, log v)
           | t <- [lo .. n - 1]
           , let !v = pmf `V.unsafeIndex` t
           , v > floorV
           ]
    !np = fromIntegral (length pts) :: Double
    !sx = sum (map fst pts)
    !sy = sum (map snd pts)
    !sxx = sum (map (\(x, _) -> x * x) pts)
    !sxy = sum (map (\(x, y) -> x * y) pts)
    !den = np * sxx - sx * sx
    !slope = if abs den < 1e-30 then 0 else (np * sxy - sx * sy) / den

-- | Unbiased mean from a censored observation. Splitting the expectation over
-- whether arrival preceded the horizon,
--
--   E[T] = (1 - eps) E[T | T <= t_max] + eps E[T | T > t_max],
--
-- and using the memorylessness of a geometric tail to put
-- E[T | T > t_max] at t_max + tau. Every term is available from the
-- observation itself.
correctMean :: Double -> Double -> Int -> Double -> Double
correctMean !condMean !eps !tmax !tau =
    (1 - eps) * condMean + eps * (fromIntegral tmax + tau)

-- | The correction applied to a tabulated distribution, taking the horizon to
-- be the last tabulated index and fitting the tail scale from the distribution
-- itself.
correctedFromPMF :: V.Vector Double -> Double
correctedFromPMF pmf
    | V.length pmf < 2 = 0
    | otherwise = correctMean (conditionalMean pmf)
                              (censoredFraction pmf)
                              (V.length pmf - 1)
                              (tailScale pmf)

-- | Agreement between a mean obtained from the generating function and the
-- same mean reconstructed from the recovered distribution. The two are
-- computed by entirely separate routes -- one a derivative at the boundary of
-- the unit disc, the other a sum over inverted coefficients -- so their
-- agreement is evidence about both, and their disagreement localises to
-- whichever route the reference values do not support.
--
-- Returns the reconstruction and the relative residual.
selfConsistency :: Double -> V.Vector Double -> (Double, Double)
selfConsistency !gfMean pmf =
    let !recon = correctedFromPMF pmf
        !resid = if abs gfMean < 1e-300
                 then abs (recon - gfMean)
                 else abs (gfMean - recon) / abs gfMean
    in (recon, resid)

-- | The tail-corrected mean of a tabulated distribution, carrying the censored
-- fraction and fitted tail scale that produced it. The error bound is the
-- residual uncertainty in the correction: the omitted mass multiplied by the
-- uncertainty in where it arrives, taken conservatively as the tail scale
-- itself, since the geometric approximation to the tail is exact only
-- asymptotically.
correctedEstimate :: V.Vector Double -> Estimate
correctedEstimate pmf =
    let !v   = correctedFromPMF pmf
        !eps = censoredFraction pmf
        !tau = tailScale pmf
        !err = eps * tau
    in Estimate v err (Truncated eps tau)

-- | A quantity evaluated short of a limit and extrapolated. The bound combines
-- the residual truncation, which falls as the square of the step once the
-- leading term is removed, with the conditioning of the underlying solve,
-- which grows as its reciprocal; the coefficients are those measured on the
-- renewal system this codebase solves.
extrapolatedEstimate :: Double -> Double -> Estimate
extrapolatedEstimate !delta !v =
    let !truncErr = 8.5e4 * delta * delta * abs v
        !condErr  = 2.8e-17 / delta
    in Estimate v (truncErr + condErr) (Extrapolated delta)

sampledEstimate :: Int -> Double -> Double -> Estimate
sampledEstimate !n !se !v = Estimate v se (Sampled n se)
