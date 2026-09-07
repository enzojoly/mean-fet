{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Distribution
-- Description : Abate-Whitt inversion of probability generating functions on a
--               damped contour, and modality of the recovered distribution.
--               The sampler exploits conjugate symmetry of real coefficient
--               series and evaluates in parallel chunks, and the radix-two
--               transform reads its roots from a precomputed twiddle table. A
--               vector valued variant returns one distribution per component
--               of a vector valued generating function from a single contour,
--               which is what the per contact site transmission law needs.
--               Modality is the sign-change count of the first difference,
--               gated by the inversion error budget rather than by a chosen
--               threshold.

module Distribution
    ( invertPGF
    , invertPGFAcc
    , invertPGFWith
    , invertPGFVec
    , invertPGFVecWith
    , invertPGFNaive
    , nextPow2
    , noiseFloor
    , Modality(..)
    , modality
    , modalityAcc
    , countPeaks
    ) where

import Control.Monad (forM_)
import Control.Monad.ST (runST)
import Control.Parallel.Strategies (using, parListChunk, rdeepseq)
import Data.Bits (shiftL, shiftR, (.&.), (.|.))
import Data.Complex (Complex(..), mkPolar, conjugate, realPart)
import qualified Data.Vector.Storable as SV
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as V

invertPGF :: Int -> (Complex Double -> Complex Double) -> V.Vector Double
invertPGF !tmax gf = invertPGFWith tmax (nextPow2 (2 * tmax + 1)) 14.0 gf

invertPGFAcc :: Double -> Int -> (Complex Double -> Complex Double) -> V.Vector Double
invertPGFAcc !acc !tmax gf = invertPGFWith tmax (nextPow2 (2 * tmax + 1)) acc gf

invertPGFWith :: Int -> Int -> Double
              -> (Complex Double -> Complex Double) -> V.Vector Double
invertPGFWith !tmax !bigM !acc gf = V.generate (tmax + 1) pt
  where
    !md = fromIntegral bigM :: Double
    !r  = 10 ** (negate acc / (2 * md))

    !half = bigM `div` 2

    lower :: SV.Vector (Complex Double)
    !lower = SV.fromList
        ( [ gf (mkPolar r (2 * pi * fromIntegral j / md)) | j <- [0 .. half] ]
          `using` parListChunk 64 rdeepseq )

    node :: Int -> Complex Double
    node !j
        | j <= half = lower `SV.unsafeIndex` j
        | otherwise = conjugate (lower `SV.unsafeIndex` (bigM - j))
    {-# INLINE node #-}

    !transformed = fftForward (SV.generate bigM node)

    pt :: Int -> Double
    pt !t =
        let !rNegT = r ** fromIntegral (negate t)
        in realPart (transformed `SV.unsafeIndex` t) * rNegT / md

-- | One distribution per component of a vector valued generating function,
-- from a single pass over the contour. The nodes, the radius and the conjugate
-- symmetry are those of invertPGFWith, so the sum over components reproduces
-- the scalar inversion of the summed generating function to inversion
-- precision. The solve behind a vector valued generating function is performed
-- once per node whether or not its result is summed, so the additional cost
-- here is one transform per component and no additional solves.
invertPGFVec :: Double -> Int -> Int
             -> (Complex Double -> SV.Vector (Complex Double))
             -> [V.Vector Double]
invertPGFVec !acc !tmax !ncomp gfv =
    invertPGFVecWith tmax (nextPow2 (2 * tmax + 1)) acc ncomp gfv

invertPGFVecWith :: Int -> Int -> Double -> Int
                 -> (Complex Double -> SV.Vector (Complex Double))
                 -> [V.Vector Double]
invertPGFVecWith !tmax !bigM !acc !ncomp gfv = map component byComp
  where
    !md   = fromIntegral bigM :: Double
    !r    = 10 ** (negate acc / (2 * md))
    !half = bigM `div` 2

    lower :: [SV.Vector (Complex Double)]
    !lower =
        ( [ gfv (mkPolar r (2 * pi * fromIntegral j / md)) | j <- [0 .. half] ]
          `using` parListChunk 16 rdeepseq )

    byComp :: [SV.Vector (Complex Double)]
    !byComp =
        [ SV.fromList [ v `SV.unsafeIndex` c | v <- lower ]
        | c <- [0 .. ncomp - 1] ]

    component :: SV.Vector (Complex Double) -> V.Vector Double
    component !lc = V.generate (tmax + 1) pt
      where
        node :: Int -> Complex Double
        node !j
            | j <= half = lc `SV.unsafeIndex` j
            | otherwise = conjugate (lc `SV.unsafeIndex` (bigM - j))
        {-# INLINE node #-}

        !transformed = fftForward (SV.generate bigM node)

        pt :: Int -> Double
        pt !t =
            let !rNegT = r ** fromIntegral (negate t)
            in realPart (transformed `SV.unsafeIndex` t) * rNegT / md

invertPGFNaive :: Int -> Int -> Double
               -> (Complex Double -> Complex Double) -> V.Vector Double
invertPGFNaive !tmax !bigM !acc gf = V.generate (tmax + 1) pt
  where
    !md = fromIntegral bigM :: Double
    !r  = 10 ** (negate acc / (2 * md))

    gfReals :: V.Vector Double
    !gfReals = V.generate bigM $ \j ->
      realPart (gf (mkPolar r (2 * pi * fromIntegral j / md)))

    gfImags :: V.Vector Double
    !gfImags = V.generate bigM $ \j ->
      imagPart (gf (mkPolar r (2 * pi * fromIntegral j / md)))
      where imagPart (_ :+ b) = b

    pt :: Int -> Double
    pt !t = (rNegT / md) * V.foldl' (+) 0 (V.generate bigM termJ)
      where
        !rNegT = r ** fromIntegral (negate t)
        !td    = fromIntegral t :: Double
        termJ :: Int -> Double
        termJ !j =
          let !phi = 2 * pi * fromIntegral j * td / md
              !gr  = gfReals `V.unsafeIndex` j
              !gi  = gfImags `V.unsafeIndex` j
          in gr * cos phi + gi * sin phi

-- | Smallest feature height the inversion can resolve, as a fraction of the
-- peak. The Abate-Whitt contour leaves an aliasing residue of 10^(-acc/2) and
-- the r^(-t) prefactor amplifies rounding to about 10^(acc/4) times machine
-- epsilon; a local extremum below their maximum is not resolvable in
-- principle, whatever threshold is chosen. Applied to a per component
-- distribution the peak must be that of the summed distribution, since a floor
-- taken component by component erases structure at weakly visited sites and
-- admits rounding at strongly visited ones.
noiseFloor :: Double -> Double -> Double
noiseFloor !acc !peak =
    let !aliasing = 10 ** (negate acc / 2)
        !roundoff = 10 ** (acc / 4) * 2.220446049250313e-16
    in max aliasing roundoff * peak

-- | Modality of a recovered distribution. The peak count is the number of
-- interior maxima surviving the resolution floor; the sign-change count is the
-- number of sign changes of the first difference, which is one for a unimodal
-- distribution and three for a bimodal one. The late-population weight is the
-- mass beyond the deepest interior valley and is the physically meaningful
-- summary: it is a probability, needs no normalisation, and its sensitivity to
-- misplacing the valley is second order, since the density there is minimal.
data Modality = Modality
    { mdPeaks       :: ![(Int, Double)]
    , mdValley      :: !(Maybe (Int, Double))
    , mdW2          :: !Double
    , mdSignChanges :: !Int
    , mdFloor       :: !Double
    } deriving (Show, Eq)

modality :: V.Vector Double -> Modality
modality = modalityAcc 14.0

modalityAcc :: Double -> V.Vector Double -> Modality
modalityAcc !acc pmf
    | V.length pmf < 3 = Modality [] Nothing 0 0 0
    | peak <= 0        = Modality [] Nothing 0 0 0
    | otherwise        = Modality peaks valley w2 changes flr
  where
    !np   = V.length pmf
    !peak = if V.null pmf then 0 else V.maximum pmf
    !flr  = noiseFloor acc peak

    -- Differences below the resolution floor are held at zero, so a flat apex
    -- does not read as an absence of extremum and a ripple in the tail does
    -- not read as the presence of one.
    diffAt :: Int -> Double
    diffAt !i =
        let !d = pmf `V.unsafeIndex` (i + 1) - pmf `V.unsafeIndex` i
        in if abs d <= flr then 0 else d
    {-# INLINE diffAt #-}

    -- Sign changes of the first difference over the whole distribution,
    -- including the rise out of zero. A unimodal law changes sign once and a
    -- bimodal law three times. A law peaking at the origin never rises and
    -- changes sign not at all, which is the one case the count does not
    -- distinguish; every encounter law released from a nonzero separation
    -- rises out of zero.
    !signs = filter (/= 0) [ signum (diffAt i) | i <- [0 .. np - 2] ]
    !changes = length (filter id (zipWith (/=) signs (drop 1 signs)))

    valueAt :: Int -> Double
    valueAt !i = pmf `V.unsafeIndex` i
    {-# INLINE valueAt #-}

    -- Local maxima above the floor, with the endpoints admitted: a strict
    -- two-sided test would silently discard a late mode sitting at the
    -- truncation horizon.
    !rawPeaks =
        [ (i, v)
        | i <- [0 .. np - 1]
        , let !v = valueAt i
        , v > flr
        , let !before = if i == 0      then negativeInfinity else valueAt (i - 1)
        , let !after  = if i == np - 1 then negativeInfinity else valueAt (i + 1)
        , v >= before && v >= after
        ]

    -- A plateau is one mode, not one per sample.
    !distinctPeaks = collapse rawPeaks
      where
        collapse ((i, v) : rest@((j, w) : _))
            | j == i + 1 && w == v = collapse ((i, v) : drop 1 rest)
        collapse (p : rest) = p : collapse rest
        collapse []         = []

    -- Topographic prominence: the descent from the mode to the higher of the
    -- two saddles bounding it. Where the walk meets higher ground the saddle
    -- is the lowest point reached on the way; where it runs off the end of the
    -- support the lowest point reached still bounds the mode, since the
    -- distribution has descended that far. A side with no terrain at all --
    -- the mode sitting on the first or last sample -- imposes no bound, which
    -- is what admits a late mode at the truncation horizon while still
    -- rejecting a flat distribution, whose every point descends nowhere.
    saddle :: (Int -> Int) -> Int -> Double -> Double
    saddle step !start !v
        | j0 < 0 || j0 >= np = negativeInfinity
        | otherwise          = go j0 v
      where
        !j0 = step start
        go !j !lo
            | j < 0 || j >= np = lo
            | x > v            = lo
            | otherwise        = go (step j) (min lo x)
          where !x = valueAt j

    prominence :: Int -> Double -> Double
    prominence !i !v =
        let !key = max (saddle (subtract 1) i v) (saddle (+ 1) i v)
        in if isInfinite key then v else v - key

    !peaks = [ p | p@(i, v) <- distinctPeaks, prominence i v > flr ]

    -- The deepest point between the two leading modes. Adjacent modes leave no
    -- interior to search, in which case there is no valley to report.
    !valley = case peaks of
        (t1, _) : (t2, _) : _
            | t2 > t1 + 1 ->
                let !cands = [ (valueAt i, i) | i <- [t1 + 1 .. t2 - 1] ]
                    !lo    = minimum cands
                in Just (snd lo, fst lo)
        _ -> Nothing

    -- Mass carried by the late population. Its sensitivity to misplacing the
    -- valley is second order, the density there being minimal by construction.
    !w2 = case valley of
        Nothing      -> 0
        Just (tv, _) -> V.sum (V.drop (tv + 1) pmf)

negativeInfinity :: Double
negativeInfinity = negate (1 / 0)

-- | Retained for call sites that want only the number of resolved modes. The
-- accuracy argument is the same resolution floor the rest of the module uses.
countPeaks :: Double -> V.Vector Double -> [(Int, Double)]
countPeaks !acc = mdPeaks . modalityAcc acc

fftForward :: SV.Vector (Complex Double) -> SV.Vector (Complex Double)
fftForward !input = runST $ do
    let !n = SV.length input
        !logN = intLog2 n
        omega :: SV.Vector (Complex Double)
        !omega = SV.generate (max 1 (n `div` 2)) $ \j ->
            mkPolar 1 (negate 2 * pi * fromIntegral j / fromIntegral n)
    buf <- SM.new n
    forM_ [0 .. n - 1] $ \i -> do
        let !j = bitReverse logN i
        SM.unsafeWrite buf j (input `SV.unsafeIndex` i)
    forM_ [1 .. logN] $ \s -> do
        let !halfBlock = shiftL 1 (s - 1)
            !block     = shiftL 1 s
            !stride    = shiftR n s
        forM_ [0, block .. n - 1] $ \base ->
            forM_ [0 .. halfBlock - 1] $ \j -> do
                let !w   = omega `SV.unsafeIndex` (j * stride)
                    !top = base + j
                    !bot = top + halfBlock
                !u <- SM.unsafeRead buf top
                !v <- SM.unsafeRead buf bot
                let !wv = w * v
                SM.unsafeWrite buf top (u + wv)
                SM.unsafeWrite buf bot (u - wv)
    SV.unsafeFreeze buf
{-# NOINLINE fftForward #-}

bitReverse :: Int -> Int -> Int
bitReverse !bits !x = go bits x 0
  where
    go :: Int -> Int -> Int -> Int
    go 0 _ !acc = acc
    go !b !v !acc = go (b - 1) (shiftR v 1) (shiftL acc 1 .|. (v .&. 1))
{-# INLINE bitReverse #-}

intLog2 :: Int -> Int
intLog2 !x = go 0 1
  where
    go !acc !v
        | v >= x    = acc
        | otherwise = go (acc + 1) (shiftL v 1)
{-# INLINE intLog2 #-}

nextPow2 :: Int -> Int
nextPow2 x = go 1 where go !n = if n >= x then n else go (2 * n)
