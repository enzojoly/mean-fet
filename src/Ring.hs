{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Ring
-- Description : Homogeneous ring eigenstructure and propagator generating
--               functions, with staged tables that fold the spectrum to half
--               size and cache every cosine so no transcendental is evaluated
--               inside a z-loop.

module Ring
    ( ringDist
    , isRingNeighbour
    , qRing
    , ringEigNonLazy
    , ringEigQ
    , ringEigenvector
    , ringEigenvectorAt
    , ringQn
    , ringMatrix
    , RingTables
    , mkRingTables
    , ringDenomAt
    , ringPropAt
    , stateCount
    ) where

import Data.Complex (Complex(..), cis)
import Numeric.LinearAlgebra ((><))
import qualified Data.Vector as BV
import qualified Data.Vector.Unboxed as V

import Types (N, K, Pos, C, R, Matrix)

-- | Number of states in the walk, which is what the Kac return cost counts.
-- On the ring this coincides with the linear size; in two dimensions it does
-- not, so the two must not be conflated.
stateCount :: Int -> Int
stateCount bigN = bigN
{-# INLINE stateCount #-}

ringDist :: Int -> Int -> Int -> Int
ringDist bigN a b = min d (bigN - d)
  where !d = abs (a - b)

isRingNeighbour :: N -> K -> Pos -> Pos -> Bool
isRingNeighbour n k i j =
    let d = ringDist n i j
    in d > 0 && d <= k

qRing :: K -> R
qRing k = 1 - 1 / fromIntegral k

ringEigNonLazy :: Int -> Int -> Int -> Double
ringEigNonLazy !bigN !k !ell
    | ell == 0  = 1.0
    | otherwise = s / fromIntegral k
  where
    !nd = fromIntegral bigN :: Double
    !ld = fromIntegral ell  :: Double
    s = sum [ cos (2 * pi * ld * fromIntegral m / nd)
            | m <- [1 .. k] ]
{-# INLINE ringEigNonLazy #-}

ringEigQ :: Double -> Int -> Int -> Int -> Double
ringEigQ !q !bigN !k !ell = (1 - q) + q * ringEigNonLazy bigN k ell
{-# INLINE ringEigQ #-}

ringEigenvector :: N -> Int -> [C]
ringEigenvector n j = [ringEigenvectorAt n j x | x <- [0..n-1]]

ringEigenvectorAt :: N -> Int -> Pos -> C
ringEigenvectorAt n j x =
    cis (2 * pi * fromIntegral j * fromIntegral x / fromIntegral n)
    / sqrt (fromIntegral n)

ringQn :: Double
       -> Int
       -> Int
       -> Int
       -> Int
       -> Complex Double
       -> Complex Double
ringQn !q !bigN !k !n0 !n !z = s / fromIntegral bigN
  where
    !delta = ringDist bigN n0 n
    !nd    = fromIntegral bigN :: Double
    s = sum [ let !lam = ringEigQ q bigN k ell
                  !cosD = cos (2 * pi * fromIntegral ell * fromIntegral delta / nd)
              in (cosD :+ 0) / (1 - z * (lam :+ 0))
            | ell <- [0 .. bigN - 1] ]

ringMatrix :: R -> N -> K -> Matrix
ringMatrix q n k = (n >< n) [w i j | i <- [0..n-1], j <- [0..n-1]]
  where
    selfLoop = 1 - q
    neighbourWeight = q / fromIntegral (2 * k)
    w i j
        | i == j          = selfLoop
        | isRingNeighbour n k i j = neighbourWeight
        | otherwise       = 0

data RingTables = RingTables
    { rtSize   :: !Int
    , rtTop    :: !Int
    , rtLambda :: !(V.Vector Double)
    , rtWeight :: !(V.Vector Double)
    , rtCosine :: !(BV.Vector (V.Vector Double))
    }

mkRingTables :: Double -> Int -> Int -> RingTables
mkRingTables !q !bigN !k = RingTables bigN top lams ws coss
  where
    !top  = bigN `div` 2
    !nd   = fromIntegral bigN :: Double
    !lams = V.generate (top + 1) (ringEigQ q bigN k)
    !ws   = V.generate (top + 1) $ \l ->
                if l == 0 || (even bigN && l == top) then 1 else 2
    !coss = BV.generate (top + 1) $ \d ->
                V.generate (top + 1) $ \l ->
                    cos (2 * pi * fromIntegral (l * d) / nd)

ringDenomAt :: RingTables -> Complex Double -> V.Vector (Complex Double)
ringDenomAt (RingTables _ top lams ws _) !z =
    V.generate (top + 1) $ \l ->
        let !lam = lams `V.unsafeIndex` l
            !w   = ws   `V.unsafeIndex` l
        in (w :+ 0) / (1 - z * (lam :+ 0))

ringPropAt :: RingTables -> V.Vector (Complex Double) -> Int -> Int -> Complex Double
ringPropAt (RingTables bigN _ _ _ coss) !denom !n0 !n =
    let !d   = ringDist bigN n0 n
        !row = coss `BV.unsafeIndex` d
        !s   = V.sum (V.zipWith (\c dv -> (c :+ 0) * dv) row denom)
    in s / fromIntegral bigN
