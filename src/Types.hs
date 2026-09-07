{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module      : Types
-- Description : Shared aliases, defect primitives, and the domain taxonomy
--               (dimension by boundary condition) used for engine routing.

module Types
    ( N, K, Pos, C, R
    , Provenance(..)
    , Estimate(..)
    , closedEstimate
    , provenanceTag
    , Primitive(..)
    , reversible
    , primitiveDegreeShift
    , Boundary(..)
    , Domain(..)
    , translationInvariant
    , EncounterEngine(..)
    , selectEncounterEngine
    , Matrix
    , Vector
    , matrixGet
    , matrixSize
    , matrixToLists
    , matrixFromLists
    ) where

import Data.Complex (Complex(..))
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra ((!))

type N = Int

type K = Int

type Pos = Int

type C = Complex Double

type R = Double

type Matrix = LA.Matrix R

type Vector = LA.Vector R

data Primitive
    = EdgeAdd      !Int !Int
    | DirectedAdd  !Int !Int
    | EdgeDel      !Int !Int
    | WattsStrogatz !Int !Int
    | Reweight     !Int !Int !Double
    | Barrier      !Int !Int !Double
    | Asymmetric   !Int !Int !Double
    | Teleport     !Int !Double
    deriving (Show, Eq)

data Boundary = Periodic | Reflecting | Absorbing
    deriving (Show, Eq)

-- | Whether a primitive leaves the walk reversible.
--
-- The spectral construction the exact engines use requires it: a reversible
-- chain is similar to a symmetric matrix under conjugation by the square root
-- of its stationary law, and it is that similarity which supplies the real
-- spectrum and orthonormal basis. Reversibility is not needed for the renewal
-- structure, which holds for any Markov chain, nor for simulation; it is needed
-- for the cheap and well-conditioned evaluation of it.
--
-- A directed edge, a biased bond and a resetting rule each break detailed
-- balance and are therefore admissible only to the simulator.
reversible :: Primitive -> Bool
reversible p = case p of
    EdgeAdd{}       -> True
    EdgeDel{}       -> True
    WattsStrogatz{} -> True
    Barrier{}       -> True
    Reweight{}      -> True
    DirectedAdd{}   -> False
    Asymmetric{}    -> False
    Teleport{}      -> False

-- | How a primitive changes the degree of the sites it names, as pairs of site
-- and increment. Used to construct the stationary law and to assert it.
primitiveDegreeShift :: Int -> Primitive -> [(Int, Int)]
primitiveDegreeShift n p = case p of
    EdgeAdd a b       -> [(a, 1), (b, 1)]
    EdgeDel a b       -> [(a, -1), (b, -1)]
    WattsStrogatz a b -> [((a + 1) `mod` n, -1), (b, 1)]
    DirectedAdd u _   -> [(u, 1)]
    _                 -> []

data Domain = Domain
    { dmDim       :: !Int
    , dmBoundary  :: !Boundary
    , dmDefects   :: ![Primitive]
    , dmShortcuts :: ![(Int, Int)]
    } deriving (Show, Eq)

translationInvariant :: Domain -> Bool
translationInvariant d =
    dmBoundary d == Periodic
    && null (dmDefects d)
    && null (dmShortcuts d)

data EncounterEngine = ClosedRelativeRing | ClosedRelativeTorus | GenericPair
    deriving (Show, Eq)

selectEncounterEngine :: Domain -> EncounterEngine
selectEncounterEngine d
    | translationInvariant d && dmDim d == 1 = ClosedRelativeRing
    | translationInvariant d && dmDim d == 2 = ClosedRelativeTorus
    | otherwise                              = GenericPair

matrixGet :: Matrix -> Int -> Int -> R
matrixGet m i j = m ! i ! j

matrixSize :: Matrix -> Int
matrixSize = LA.rows

matrixToLists :: Matrix -> [[R]]
matrixToLists = LA.toLists

matrixFromLists :: [[R]] -> Matrix
matrixFromLists = LA.fromLists

-- | How a numerical result was obtained, and therefore what may be believed of
-- it. A bare Double cannot say whether it is exact, extrapolated from a point
-- short of a limit, summed over a truncated range, or averaged over a sample,
-- and a value that cannot say what it is cannot be checked by anything that
-- consumes it. Carrying the provenance makes an unqualified claim
-- unrepresentable: a censored mean cannot be serialised without the fraction
-- that was censored, and an extrapolated one cannot be reported without the
-- step it was extrapolated from.
data Provenance
    = Closed
      -- ^ Exact up to floating-point rounding.
    | Extrapolated !Double
      -- ^ Evaluated short of a limit and extrapolated; carries the step.
    | Truncated !Double !Double
      -- ^ Corrected for a censored tail; carries the censored fraction and the
      --   fitted tail scale.
    | Sampled !Int !Double
      -- ^ Monte Carlo; carries the sample size and the standard error.
    deriving (Show, Eq)

-- | A value together with an absolute error bound and the route that produced
-- it. The bound is derived from the provenance rather than asserted, so it
-- moves when the method does.
data Estimate = Estimate
    { estValue      :: !Double
    , estError      :: !Double
    , estProvenance :: !Provenance
    } deriving (Show, Eq)

closedEstimate :: Double -> Estimate
closedEstimate v = Estimate v (8 * 2.220446049250313e-16 * abs v) Closed

provenanceTag :: Provenance -> String
provenanceTag Closed             = "closed"
provenanceTag (Extrapolated _)   = "extrapolated"
provenanceTag (Truncated _ _)    = "tail-corrected"
provenanceTag (Sampled _ _)      = "sampled"
