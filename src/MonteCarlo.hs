{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : MonteCarlo
-- Description : The Monte Carlo reproducibility protocol. Work is always split
--               into exactly protocolBatches logical batches whose generator
--               streams derive only from the seed, never from the machine, so
--               results are identical at any core count; the batch count is a
--               protocol constant and must never be changed to track the number
--               of capabilities. Substreams are keyed by seed, batch, and trial
--               and never by any sweep parameter, so common random numbers hold
--               across parameter sweeps by construction.

module MonteCarlo
    ( protocolBatches
    , batchGens
    , batchSizes
    , runBatches
    , mergeWelford
    ) where

import Control.DeepSeq (NFData)
import Control.Parallel.Strategies (parMap, rdeepseq)
import System.Random (StdGen, mkStdGen, splitGen)

protocolBatches :: Int
protocolBatches = 64

batchGens :: Int -> [StdGen]
batchGens seed = go protocolBatches (mkStdGen seed)
  where
    go :: Int -> StdGen -> [StdGen]
    go 0 _  = []
    go !n g = let (!h, !g') = splitGen g in h : go (n - 1) g'

batchSizes :: Int -> [Int]
batchSizes total =
    let !b   = protocolBatches
        !per = total `div` b
    in replicate (b - 1) per ++ [total - per * (b - 1)]

runBatches :: NFData r => Int -> Int -> (Int -> StdGen -> r) -> [r]
runBatches seed total runOne =
    parMap rdeepseq (\(g, m) -> runOne m g) (zip (batchGens seed) (batchSizes total))

mergeWelford :: (Int, Double, Double) -> (Int, Double, Double) -> (Int, Double, Double)
mergeWelford (!nA, !meanA, !m2A) (!nB, !meanB, !m2B) = (nAB, meanAB, m2AB)
  where
    !nAB   = nA + nB
    !delta = meanB - meanA
    !meanAB = if nAB > 0
              then meanA + delta * fromIntegral nB / fromIntegral nAB
              else 0
    !m2AB  = m2A + m2B + delta * delta * fromIntegral nA * fromIntegral nB
                                       / max 1 (fromIntegral nAB)
