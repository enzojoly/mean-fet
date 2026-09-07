{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Cells
-- Description : The six cells of the programme behind one interface. Dimension
--               and boundary condition are the two axes along which the work is
--               organised, so they are the two arguments a caller supplies; the
--               choice of engine, of basis, and of what the geometry admits as
--               an observable follows from them rather than from the caller.
--
--               Two things vary with the cell and neither is a matter of
--               presentation. A measure-preserving domain returns walkers to
--               circulation, so the encounter law is proper and its mean is the
--               observable; an absorbing domain may destroy a walker first, so
--               the law is defective and the observable is the pair of a
--               splitting weight and a conditional mean. And a
--               translation-invariant domain admits a relative-coordinate
--               reduction that the others do not, which is a difference of cost
--               rather than of principle.
--
--               Reporting the pair everywhere removes the case distinction: a
--               measure-preserving cell is the weight-one member of the same
--               family. Where the contact renewal is used the result carries
--               the decomposition of the mean by contact site as well, which
--               the closed and direct routes cannot supply.

module Cells
    ( Dim(..)
    , BC(..)
    , CellResult(..)
    , cellMatrix
    , cellMatrixWith
    , cellStates
    , cellIndex
    , closedSpectrumExists
    , encounterCell
    ) where

import qualified Data.Vector.Unboxed as V
import qualified Numeric.LinearAlgebra as LA

import Types (Matrix, R, Primitive(..))
import Defect (exactEigensystem, stationaryDist)
import Absorbing (encounterAbsorbing)
import Encounter
    ( encounterPMFTwoQ
    , encounterMeanRingTwoQ
    , encounterPMFRingTwoQ
    , encounterDecompositionTwoQ
    , degenerateRingCell
    , degeneratePeriodicCell
    )

data Dim = D1 | D2
    deriving (Show, Eq)

data BC = Periodic | Reflecting | Absorbing
    deriving (Show, Eq)

-- | What a cell yields. The splitting weight is the probability that the pair
-- meets at all, which is one wherever the domain preserves measure; the mean is
-- conditioned on that event, and therefore coincides with the unconditional
-- mean exactly when the weight is one.
--
-- The two per-site vectors are the splitting probability at each contact site
-- and that site's contribution to the mean, which is its weight multiplied by
-- the mean conditional on encounter occurring there. They sum to the splitting
-- weight and to the mean respectively. Only the contact renewal produces them,
-- so the closed relative-coordinate route and the direct absorbing route
-- return them empty rather than fabricating them.
data CellResult = CellResult
    { crSplittingWeight :: !Double
    , crMean            :: !Double
    , crPMF             :: !(V.Vector Double)
    , crProvenance      :: !String
    , crWeights         :: !(V.Vector Double)
    , crSiteMeans       :: !(V.Vector Double)
    } deriving (Show)

-- | Number of states, which is what a Kac return costs and what the pair space
-- is built on. In one dimension it coincides with the linear size; in two it
-- does not, and conflating them understates a mean by the difference.
cellStates :: Dim -> Int -> Int
cellStates D1 n = n
cellStates D2 l = l * l

cellIndex :: Dim -> Int -> (Int, Int) -> Int
cellIndex D1 _ (x, _) = x
cellIndex D2 l (x, y) = x + l * y

-- | Whether the cell has a closed-form spectrum at the given range.
--
-- The periodic and reflecting cells have one at every range: the Fourier and
-- cosine bases diagonalise the wrapping and mirror walls exactly, whatever the
-- reach of a step. The absorbing cell does not. The odd extension that
-- annihilates the eigenfunction at the absorbing site leaves it non-zero one
-- site further out, and at unit range the walk cannot reach that site while at
-- greater range it can. That cell is therefore diagonalised numerically beyond
-- unit range, which costs a factor in construction and nothing in accuracy.
closedSpectrumExists :: BC -> Int -> Bool
closedSpectrumExists Absorbing k = k == 1
closedSpectrumExists _         _ = True

-- | The undecorated cell. Delegated rather than written out, so that the
-- boundary rule cannot differ between the bare lattice and the modified one.
cellMatrix :: Dim -> BC -> Int -> Int -> Double -> Matrix
cellMatrix dim bc size k q = cellMatrixWith dim bc size k q []

-- | The cell's transition matrix carrying a list of local modifications.
--
-- Built from the attempted steps rather than by adding a perturbation to the
-- undecorated matrix. The perturbations are written against a periodic
-- neighbourhood and adding one to a bounded cell would modify the wrong bonds:
-- the first site of an absorbing interval has one neighbour inside the domain
-- where a ring would give it two. Working from the attempts states the walk
-- once and lets each boundary dispose of those that leave.
--
-- The boundary rules are three and they are not interchangeable. A periodic
-- wall wraps. An absorbing wall keeps what leaves, which is why those rows fall
-- short of one. A reflecting wall returns it --- and how it returns it matters.
--
-- The rule used here is mirror reflection: an attempt that would land at
-- $-m$ arrives instead at $m-1$, reflected about the half-integer line outside
-- the domain. This is the wall the cosine basis diagonalises exactly, at every
-- range, and it is the wall for which the reflecting interval on N sites
-- carries the even sub-spectrum of a ring of 2N.
--
-- The alternative, under which a refused attempt returns to the site it came
-- from, coincides with this at unit range and diverges beyond it. It is also
-- worse physics: it gives boundary sites a raised and site-dependent holding
-- probability -- at range two on ten sites, 0.625 against 0.250 in the bulk --
-- which is a spatial heterogeneity introduced by the numerics rather than by
-- the model.
cellMatrixWith :: Dim -> BC -> Int -> Int -> Double -> [Primitive] -> Matrix
cellMatrixWith dim bc !size !k !q prims =
    LA.accum (LA.konst 0 (nn, nn)) (+) entries
  where
    !nn = cellStates dim size
    !attempts = decorate dim size prims (baseAttempts dim size k)
    entries = concat
        [ ((c, c), 1 - q) : row c | c <- [0 .. nn - 1] ]
    row c =
        let ts = attempts !! c
            !deg = fromIntegral (length ts) :: Double
        in [ ((c, j), q / deg) | Just j <- map resolve ts ]

    -- Where an attempted step actually lands, or nothing if it is lost.
    resolve (x, y) = case bc of
        Periodic   -> Just (idx (wrap x) (wrap y))
        Reflecting -> Just (idx (mirror x) (mirror y))
        Absorbing  -> if inRange x && inRange y then Just (idx x y) else Nothing

    idx x y = case dim of
        D1 -> x
        D2 -> x + size * y
    inRange v = v >= 0 && v < size
    wrap v = ((v `mod` size) + size) `mod` size
    -- Reflection about the half-integer lines just outside the domain. Valid
    -- while the range is shorter than the domain, which every cell here
    -- satisfies; a longer step would need repeated folding.
    mirror v
        | v < 0     = negate 1 - v
        | v >= size = 2 * size - 1 - v
        | otherwise = v

-- | The steps each site attempts, before any boundary rule is applied and
-- before any modification. Targets outside the domain are kept as they are, so
-- that the boundary rule may be applied once and in one place rather than being
-- folded into the neighbour list of every cell separately.
--
-- A site is carried as a pair, with the second coordinate zero in one
-- dimension, so that a reflection may act on each axis independently.
baseAttempts :: Dim -> Int -> Int -> [[(Int, Int)]]
baseAttempts D1 n k =
    [ [ (c + m, 0) | m <- [1 .. k] ] ++ [ (c - m, 0) | m <- [1 .. k] ]
    | c <- [0 .. n - 1] ]
baseAttempts D2 l k =
    [ concat [ [ (x + m, y), (x - m, y), (x, y + m), (x, y - m) ]
             | m <- [1 .. k] ]
    | y <- [0 .. l - 1], x <- [0 .. l - 1] ]

-- | Apply the adjacency-changing primitives. A weight-changing primitive has no
-- adjacency to change and is rejected here rather than silently misapplied.
decorate :: Dim -> Int -> [Primitive] -> [[(Int, Int)]] -> [[(Int, Int)]]
decorate dim size ps adj0 = foldl one adj0 ps
  where
    coord i = case dim of
        D1 -> (i, 0)
        D2 -> (i `mod` size, i `div` size)
    one adj p = case p of
        EdgeAdd a b       -> add b a (add a b adj)
        DirectedAdd u v   -> add u v adj
        EdgeDel a b       -> del b a (del a b adj)
        WattsStrogatz a b -> add b a (add a b (del aPlus a (del a aPlus adj)))
          where aPlus = (a + 1) `mod` length adj
        other -> error ("Cells.decorate: " ++ show other
                        ++ " modifies bond weights rather than adjacency and is \
                           \available on the periodic ring alone")
    add i j adj =
        [ if c == i then coord j : ns else ns | (c, ns) <- zip [0 ..] adj ]
    del i j adj =
        [ if c == i then filter (/= coord j) ns else ns
        | (c, ns) <- zip [0 ..] adj ]

-- | First-encounter statistics for a cell.
--
-- The route is chosen by the cell rather than by the caller. A
-- translation-invariant one-dimensional cell admits the relative-coordinate
-- closed form, which needs no limit and no linear solve. An absorbing cell has
-- no pole at the boundary of the unit disc, so its generating function is
-- evaluated there directly and both the weight and the conditional mean follow.
-- Everything else goes through the contact renewal, which additionally supplies
-- the decomposition of the mean by contact site.
--
-- The horizon governs the distribution only. Every mean returned here is
-- truncation-free, so a caller wanting the mean alone may pass a small horizon,
-- or zero to skip the inversion entirely.
encounterCell :: Dim -> BC -> Int -> Int -> Double -> Double -> Double
              -> ((Int, Int), (Int, Int)) -> [Primitive] -> Int -> CellResult
encounterCell dim bc !size !k !q1 !q2 !rho (startA, startB) prims !tmax =
    case (dim, bc) of
        -- The relative-coordinate reduction needs translation invariance, which
        -- a defect destroys, so the closed form serves the undecorated
        -- periodic ring alone and everything else goes through the renewal.
        (D1, Periodic)
            | not (null prims) -> generic
            | degenerateRingCell q1 q2 size k ->
                error "encounterCell: degenerate cell (q1 = q2 = 1, k = 1, \
                      \N even). The relative spectrum touches unity away from \
                      \the stationary mode and the mean does not exist."
            | otherwise ->
                let !sA = fst startA
                    !sB = fst startB
                    !m  = encounterMeanRingTwoQ q1 q2 size k rho sA sB
                    !p  = if tmax > 0
                          then encounterPMFRingTwoQ q1 q2 size k rho sA sB tmax
                          else V.empty
                in CellResult 1.0 m p "closed" V.empty V.empty

        -- The torus carries the same obstruction as the ring and by the same
        -- mechanism: at unit mobility and unit range parity is conserved, the
        -- axis eigenvalues reach minus one together when the side length is
        -- even, and two walkers released on opposite parity never meet. The
        -- renewal route detects this from the spectrum itself, but failing
        -- here names the configuration rather than the symptom.
        (D2, Periodic)
            | null prims && degeneratePeriodicCell q1 q2 size k ->
                error "encounterCell: degenerate torus (q1 = q2 = 1, k = 1, \
                      \L even). Parity is conserved, the relative spectrum \
                      \touches unity away from the stationary mode, and two \
                      \walkers released on opposite parity never meet."
            | otherwise -> generic

        (_, Absorbing) ->
            let !wA = cellMatrixWith dim bc size k q1 prims
                !wB = cellMatrixWith dim bc size k q2 prims
                !a  = cellIndex dim size startA
                !b  = cellIndex dim size startB
                (!phi, !cond) = encounterAbsorbing wA wB rho (a, b)
                !p = if tmax > 0 then pmfVia wA wB a b else V.empty
            in CellResult phi cond p "direct" V.empty V.empty

        _ -> generic
  where
    generic =
        let !wA = cellMatrixWith dim bc size k q1 prims
            !wB = cellMatrixWith dim bc size k q2 prims
            !a  = cellIndex dim size startA
            !b  = cellIndex dim size startB
            (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !piV = stationaryDist wA
            (!m, !w, !sm) = encounterDecompositionTwoQ eA vA eB vB piV rho (a, b)
            !p = if tmax > 0
                 then encounterPMFTwoQ eA vA eB vB piV rho (a, b) tmax
                 else V.empty
        in CellResult 1.0 m p "expanded" w sm

    pmfVia wA wB a b =
        let (!eA, !vA) = exactEigensystem wA
            (!eB, !vB) = exactEigensystem wB
            !piV = replicate (LA.rows wA) (1 / fromIntegral (LA.rows wA))
        in encounterPMFTwoQ eA vA eB vB piV rho (a, b) tmax
