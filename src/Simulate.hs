{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Simulate
-- Description : Monte Carlo first-passage and first-encounter trials on the
--               defected ring. Trial physics and Welford merging are unchanged;
--               generator streams come from the fixed MonteCarlo batch
--               protocol, so results depend only on the seed and never on the
--               core count. Trials that have not arrived by the horizon are
--               reported as such rather than discarded: the mean over arrivals
--               alone is a conditional quantity, low by the whole of the
--               omitted tail, and the omitted trajectories are by construction
--               the slowest ones. Both the conditional mean and the corrected
--               estimate are returned, together with the censoring diagnostics
--               that relate them.

module Simulate
    ( SimMode(..)
    , Termination(..)
    , SimConfig(..)
    , SimResult(..)
    , defaultSimConfig
    , terminationHorizon
    , parityForbidden
    , simulate
    ) where

import Control.DeepSeq (NFData(..))
import Control.Monad.ST (runST)
import Data.List (foldl1')
import System.Random (StdGen, uniformR, splitGen)
import qualified Data.Vector as BV
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as MV
import qualified Numeric.LinearAlgebra as LA

import Types (Primitive, Matrix)
import Defect (buildTransitionMatrix)
import Estimate (correctMean, tailScale)
import qualified MonteCarlo as MC

data SimMode = Passage | Encounter
    deriving (Show, Eq)

-- | How a trial is allowed to end. A horizon of zero or less defers to the
-- configuration's own horizon field. Under a horizon the run is cut at a fixed
-- time and the survivors are censored, which is the only admissible policy
-- where arrival is not certain -- an absorbing domain can lose a walker
-- permanently, and there the censored fraction is itself the observable.
-- Running to absorption removes the censoring at the cost of an unbounded
-- trial, and carries a safety cap so that a configuration in which arrival is
-- impossible terminates and is reported rather than hanging.
data Termination
    = Horizon !Int
    | RunToAbsorption !Int
    deriving (Show, Eq)

terminationHorizon :: Termination -> Int
terminationHorizon (Horizon t)         = t
terminationHorizon (RunToAbsorption c) = c

-- | The horizon a run will actually use. A non-positive horizon defers to the
-- legacy field, so a configuration built before the policy existed behaves
-- exactly as it did and nothing is silently capped.
effectiveHorizon :: SimConfig -> Int
effectiveHorizon cfg = case simTermination cfg of
    Horizon t
        | t <= 0    -> max 1 (simMaxT cfg)
        | otherwise -> max 1 t
    RunToAbsorption c -> max 1 c

-- | Whether the sublattice parity of a nearest-neighbour walk forbids arrival
-- outright. Two non-lazy walkers on a bipartite lattice alternate between
-- sublattices in step, so an odd separation is preserved for ever and the
-- encounter probability is exactly zero. Running such a configuration to
-- absorption would never terminate.
parityForbidden :: SimConfig -> Bool
parityForbidden cfg =
    simMode cfg == Encounter
    && simK cfg == 1
    && simQ cfg >= 1 - 1e-15
    && simQB cfg >= 1 - 1e-15
    && odd (abs (simSrc cfg - simSrcB cfg))

data SimConfig = SimConfig
    { simN       :: !Int
    , simK       :: !Int
    , simQ       :: !Double
    , simQB      :: !Double
    , simSrc     :: !Int
    , simTgt     :: !Int
    , simSrcB    :: !Int
    , simDefects :: ![Primitive]
    , simWalkers :: !Int
    , simSeed    :: !Int
    , simMaxT    :: !Int
    , simMode    :: !SimMode
    , simRho     :: !Double
    , simTermination :: !Termination
    } deriving (Show)

defaultSimConfig :: SimConfig
defaultSimConfig = SimConfig
    { simN = 96, simK = 4, simQ = 0.75, simQB = 0.75
    , simSrc = 0, simTgt = 48, simSrcB = 48
    , simDefects = []
    , simWalkers = 10000, simSeed = 42, simMaxT = 2000
    , simMode = Passage, simRho = 1.0
    , simTermination = Horizon 0
    }

-- | Outcome of a Monte Carlo run. The mean carried by srMeanFPT is the
-- corrected estimate of the unconditional mean, so that a consumer reading the
-- obvious field reads the quantity the exact engines compute. The conditional
-- mean over arrivals alone is retained beside it, together with the censored
-- fraction and the fitted tail scale that relate the two, so the correction is
-- auditable rather than implicit.
data SimResult = SimResult
    { srFPTs             :: !(V.Vector Int)
    , srHistogram        :: !(V.Vector Double)
    , srMeanFPT          :: !Double
    , srConditionalMean  :: !Double
    , srCensoredFraction :: !Double
    , srTailScale        :: !Double
    , srCorrection       :: !Double
    , srStdErr           :: !Double
    , srAbsorbed         :: !Int
    , srSurvived         :: !Int
    } deriving (Show)

type TransRow = V.Vector Double

extractRows :: Matrix -> BV.Vector TransRow
extractRows w = BV.fromList [V.fromList (map realToFrac (LA.toList r)) | r <- LA.toRows w]

stepFromRow :: TransRow -> StdGen -> (Int, StdGen)
stepFromRow !row !gen =
    let (!r, !gen') = uniformR (0.0 :: Double, 1.0) gen
    in (scanRow row r 0 0.0, gen')

scanRow :: TransRow -> Double -> Int -> Double -> Int
scanRow !row !r !j !cumul
    | j >= V.length row - 1 = j
    | otherwise =
        let !cumul' = cumul + row `V.unsafeIndex` j
        in if r <= cumul' then j else scanRow row r (j + 1) cumul'

trialFP :: BV.Vector TransRow -> SimConfig -> StdGen -> Int
trialFP !rows !cfg !gen0 = go (simSrc cfg) 0 gen0
  where
    !maxT = effectiveHorizon cfg
    !tgt  = simTgt cfg
    !rho  = simRho cfg
    go !pos !t !gen
        | t >= maxT = 0
        | t > 0 && pos == tgt =
            if rho >= 1.0 - 1e-15
            then t
            else let (!r, !gen') = uniformR (0.0 :: Double, 1.0) gen
                 in if r <= rho then t
                    else let (!pos', !gen'') = stepFromRow (rows `BV.unsafeIndex` pos) gen'
                         in go pos' (t + 1) gen''
        | otherwise =
            let (!pos', !gen') = stepFromRow (rows `BV.unsafeIndex` pos) gen
            in go pos' (t + 1) gen'

trialFE :: BV.Vector TransRow -> BV.Vector TransRow -> SimConfig -> StdGen -> Int
trialFE !rowsA !rowsB !cfg !gen0 =
    let (!genA, !genB) = splitGen gen0
    in go (simSrc cfg) (simSrcB cfg) 0 genA genB
  where
    !maxT = effectiveHorizon cfg
    !rho  = simRho cfg
    go !posA !posB !t !gA !gB
        | t >= maxT             = 0
        | t > 0 && posA == posB =
            if rho >= 1.0 - 1e-15
            then t
            else let (!r, !gA') = uniformR (0.0 :: Double, 1.0) gA
                 in if r <= rho then t
                    else let (!posA', !gA'') = stepFromRow (rowsA `BV.unsafeIndex` posA) gA'
                             (!posB', !gB')  = stepFromRow (rowsB `BV.unsafeIndex` posB) gB
                         in go posA' posB' (t + 1) gA'' gB'
        | otherwise =
            let (!posA', !gA') = stepFromRow (rowsA `BV.unsafeIndex` posA) gA
                (!posB', !gB') = stepFromRow (rowsB `BV.unsafeIndex` posB) gB
            in go posA' posB' (t + 1) gA' gB'

data BatchResult = BatchResult
    { brHistCounts :: !(V.Vector Int)
    , brAbsorbed   :: !Int
    , brMean       :: !Double
    , brM2         :: !Double
    }

instance NFData BatchResult where
    rnf (BatchResult h a m m2) = rnf h `seq` rnf a `seq` rnf m `seq` rnf m2

runBatch :: (StdGen -> Int) -> StdGen -> Int -> Int -> BatchResult
runBatch !trialFn !gen0 !batchSize !maxT = runST $ do
    hist <- MV.replicate maxT (0 :: Int)
    let go !i !g !absorbed !mean !m2
            | i >= batchSize = do
                frozen <- V.unsafeFreeze hist
                return $! BatchResult frozen absorbed mean m2
            | otherwise = do
                let (!g1, !g2) = splitGen g
                    !fpt = trialFn g1
                if fpt > 0 && fpt <= maxT
                    then do
                        MV.unsafeModify hist (+ 1) (fpt - 1)
                        let !absorbed' = absorbed + 1
                            !delta  = fromIntegral fpt - mean
                            !mean'  = mean + delta / fromIntegral absorbed'
                            !delta2 = fromIntegral fpt - mean'
                            !m2'    = m2 + delta * delta2
                        go (i + 1) g2 absorbed' mean' m2'
                    else
                        go (i + 1) g2 absorbed mean m2
    go 0 gen0 0 0.0 0.0

mergeResults :: [BatchResult] -> Int -> SimResult
mergeResults batches totalWalkers = SimResult
    { srFPTs             = V.empty
    , srHistogram        = finalHist
    , srMeanFPT          = correctedMean
    , srConditionalMean  = finalMean
    , srCensoredFraction = censored
    , srTailScale        = tau
    , srCorrection       = correctedMean - finalMean
    , srStdErr           = finalStdErr
    , srAbsorbed         = totalAbsorbed
    , srSurvived         = totalWalkers - totalAbsorbed
    }
  where
    !horizon  = V.length rawCounts
    !censored = if totalWalkers <= 0 then 0
                else fromIntegral (totalWalkers - totalAbsorbed)
                     / fromIntegral totalWalkers
    !tau = tailScale finalHist
    !correctedMean = if totalAbsorbed <= 0 then 0
                     else correctMean finalMean censored horizon tau
    !totalAbsorbed = sum (map brAbsorbed batches)
    !rawCounts     = foldl1' (V.zipWith (+)) (map brHistCounts batches)
    !finalHist     = V.map (\c -> fromIntegral c / fromIntegral totalWalkers) rawCounts
    (!_, !finalMean, !finalM2) = foldl1' MC.mergeWelford
        [(brAbsorbed b, brMean b, brM2 b) | b <- batches]
    !finalVariance = if totalAbsorbed > 1
                     then finalM2 / fromIntegral (totalAbsorbed - 1)
                     else 0
    !finalStdErr   = if totalAbsorbed > 0
                     then sqrt finalVariance / sqrt (fromIntegral totalAbsorbed)
                     else 0

simulate :: SimConfig -> SimResult
simulate cfg
    | parityForbidden cfg =
        error "Simulate.simulate: sublattice parity forbids encounter for this \
              \configuration (two non-lazy nearest-neighbour walkers at odd \
              \separation); arrival has probability zero and no horizon will \
              \produce one. Use q < 1 or an even separation."
    | otherwise =
    let !total   = simWalkers cfg
        !maxT    = effectiveHorizon cfg
        !n       = simN cfg
        !k       = simK cfg
        !qA      = simQ cfg
        !qB      = simQB cfg
        !wA      = buildTransitionMatrix qA n k (simDefects cfg)
        !wB      = if abs (qA - qB) < 1e-15
                   then wA
                   else buildTransitionMatrix qB n k (simDefects cfg)
        !rowsA   = extractRows wA
        !rowsB   = extractRows wB
        trialFn  = case simMode cfg of
            Passage   -> trialFP rowsA cfg
            Encounter -> trialFE rowsA rowsB cfg
        !batches = MC.runBatches (simSeed cfg) total
                     (\m g -> runBatch trialFn g m maxT)
    in mergeResults batches total
