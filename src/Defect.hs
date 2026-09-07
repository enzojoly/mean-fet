{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Defect
-- Description : Defect primitives, their perturbation matrices, the defect
--               technique H matrices, and dense eigensystem utilities. The
--               staged entry point buildHWith consumes defects as a boxed
--               vector and homogeneous propagators through a supplied lookup,
--               so nothing is rebuilt per contour node.

module Defect
    ( DefectEntry(..)
    , primitivesToDefects
    , buildHWith
    , scalarH
    , scalarHmod
    , buildTransitionMatrix
    , primitivePerturbation
    , exactEigenvalues
    , exactEigensystem
    , stationaryDist
    ) where

import Data.Complex (Complex(..))
import qualified Data.Map.Strict as Map
import qualified Data.Vector as BV
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra ((><))

import Types (R, N, K, Pos, Matrix, Primitive(..))
import Ring (ringQn, ringDist, ringMatrix, qRing)

data DefectEntry = DefectEntry
    { defU    :: !Int
    , defV    :: !Int
    , defXvu  :: !Double
    , defXuv  :: !Double
    } deriving Show

primitivesToDefects :: Double -> Int -> Int -> [Primitive] -> [DefectEntry]
primitivesToDefects q bigN k = concatMap (primitiveToDefects q bigN k)

primitiveToDefects :: Double -> Int -> Int -> Primitive -> [DefectEntry]
primitiveToDefects q bigN k prim = case prim of
    EdgeAdd a b      -> edgeAddDefects q bigN k a b
    DirectedAdd a b  -> directedAddDefects q bigN k a b
    EdgeDel a b      -> edgeDelDefects q bigN k a b
    WattsStrogatz a b -> wattsStrogatzDefects q bigN k a b
    Barrier a b p    -> barrierDefects q bigN k a b p
    Asymmetric a b d -> asymmetricDefects q bigN k a b d
    Reweight a b d   -> reweightDefects q bigN k a b d
    Teleport m r     -> teleportDefects q bigN k m r

-- | The same primitives again, as a list of bond modifications for the defect
-- technique. This representation and the matrix representation above describe
-- one object and must agree; they are held together by an assertion rather than
-- by a shared body, since one is a dense matrix and the other a sparse list.
--
-- A single DefectEntry carries a rank-one modification in the difference vector
-- of its two sites, which can move weight along one bond but cannot alter a
-- degree. Raising the degree of a site therefore takes one entry per bond it
-- owns, not one entry in total: the mobility spread over the enlarged
-- neighbourhood is a change to every bond at that site.
edgeAddDefects :: Double -> Int -> Int -> Int -> Int -> [DefectEntry]
edgeAddDefects q bigN k a b
    | ringDist bigN a b <= k = error $
        "EdgeAdd: " ++ show a ++ "," ++ show b ++ " already neighbours"
    | otherwise = ringA ++ ringB ++ [shortcut]
  where
    !bigK = fromIntegral (2 * k) :: Double
    -- Weight withdrawn from each existing bond, positive because a positive
    -- strength removes weight.
    !dRing = q / bigK - q / (bigK + 1)
    -- Weight supplied to the new bond, negative because it is an addition.
    !dSC   = negate (q / (bigK + 1))
    neighbours x = [ (x + d) `mod` bigN | d <- [negate k .. negate 1] ++ [1 .. k] ]
    ringA = [ DefectEntry a r dRing 0.0 | r <- neighbours a ]
    ringB = [ DefectEntry b s dRing 0.0 | s <- neighbours b ]
    shortcut = DefectEntry a b dSC dSC

directedAddDefects :: Double -> Int -> Int -> Int -> Int -> [DefectEntry]
directedAddDefects q bigN k u v = ringU ++ [shortcut]
  where
    !bigK = fromIntegral (2 * k) :: Double
    !dRing = q / bigK - q / (bigK + 1)
    !dSC   = negate (q / (bigK + 1))
    neighbours x = [ (x + d) `mod` bigN | d <- [negate k .. negate 1] ++ [1 .. k] ]
    ringU = [ DefectEntry u r dRing 0.0 | r <- neighbours u ]
    shortcut = DefectEntry u v 0.0 dSC

edgeDelDefects :: Double -> Int -> Int -> Int -> Int -> [DefectEntry]
edgeDelDefects q bigN k a b
    | ringDist bigN a b > k = error $
        "EdgeDel: " ++ show a ++ "," ++ show b ++ " are not neighbours"
    | k < 1 = []
    | otherwise = ringA ++ ringB ++ [cut]
  where
    !bigK = fromIntegral (2 * k) :: Double
    -- Weight supplied to each surviving bond, negative because it is an
    -- addition; the mobility is spread over one fewer destination.
    !dRing = negate (q / (bigK - 1) - q / bigK)
    !dCut  = q / bigK
    neighbours x = [ (x + d) `mod` bigN | d <- [negate k .. negate 1] ++ [1 .. k] ]
    ringA = [ DefectEntry a r dRing 0.0 | r <- neighbours a, r /= b ]
    ringB = [ DefectEntry b s dRing 0.0 | s <- neighbours b, s /= a ]
    cut = DefectEntry a b dCut dCut

-- | Rewiring, given directly rather than as a removal followed by an addition.
-- Composing those two would have the removal shrink a neighbourhood the
-- addition has already enlarged, since each is written against the undecorated
-- lattice. Here the source keeps its degree, having exchanged one neighbour for
-- another; the abandoned neighbour loses one; the new endpoint gains one.
wattsStrogatzDefects :: Double -> Int -> Int -> Int -> Int -> [DefectEntry]
wattsStrogatzDefects q bigN k u v
    | ringDist bigN u v <= k = error $
        "WattsStrogatz: " ++ show u ++ "," ++ show v ++ " already neighbours"
    | k < 1 = []
    | otherwise = [swapOut, swapIn] ++ ringAbandoned ++ ringNew
  where
    !uPlus = (u + 1) `mod` bigN
    !bigK = fromIntegral (2 * k) :: Double
    neighbours x = [ (x + d) `mod` bigN | d <- [negate k .. negate 1] ++ [1 .. k] ]
    -- The source exchanges one destination for another at unchanged weight.
    swapOut = DefectEntry u uPlus (q / bigK) 0.0
    swapIn  = DefectEntry u v (negate (q / bigK)) 0.0
    -- The abandoned neighbour spreads the same mobility over one fewer bond.
    ringAbandoned =
        [ DefectEntry uPlus r (negate (q / (bigK - 1) - q / bigK)) 0.0
        | r <- neighbours uPlus, r /= u ]
    -- The new endpoint spreads it over one more.
    ringNew =
        [ DefectEntry v s (q / bigK - q / (bigK + 1)) 0.0
        | s <- neighbours v ]

-- | A permeable barrier scales the weight of one existing bond by p, the
-- refused mass returning to the walker. One formula at every mobility: the bond
-- carries q/2k before and p times that after.
barrierDefects :: Double -> Int -> Int -> Int -> Int -> Double -> [DefectEntry]
barrierDefects q _bigN k a b p =
    let !w0 = q / fromIntegral (2 * k)
        !reduction = w0 * (1.0 - p)
    in [ DefectEntry a b reduction reduction ]

asymmetricDefects :: Double -> Int -> Int -> Int -> Int -> Double -> [DefectEntry]
asymmetricDefects _q _bigN _k a b d =
    [ DefectEntry a b (negate d) d ]

reweightDefects :: Double -> Int -> Int -> Int -> Int -> Double -> [DefectEntry]
reweightDefects _q _bigN _k a b d =
    [ DefectEntry a b (negate d) (negate d) ]

teleportDefects :: Double -> Int -> Int -> Int -> Double -> [DefectEntry]
teleportDefects _q bigN _k m r =
    [ DefectEntry i m (negate r) 0.0 | i <- [0 .. bigN - 1], i /= m ]

buildTransitionMatrix :: R -> N -> K -> [Primitive] -> Matrix
buildTransitionMatrix q n k prims = foldl LA.add w0 perturbations
  where
    !w0 = ringMatrix q n k
    perturbations = map (primitivePerturbation q n k) prims

-- | Perturbation matrix for one primitive.
--
-- Every primitive is a modification of the neighbour structure at a small
-- number of sites, and each is written here as a single expression valid at
-- every mobility. The previous formulation carried a separate branch at unit
-- mobility, because the rule it used tied the strength of a shortcut to the
-- probability of holding and therefore had nothing left to redirect once the
-- walker stopped holding. That is a different model from the one these names
-- describe, and it is discontinuous where the two meet.
--
-- The rule used throughout is the one the names mean: adding an edge raises the
-- degree of its endpoints and the mobility is shared over the enlarged
-- neighbourhood, so that a walker at a site of degree d steps to each neighbour
-- with probability q/d and holds with probability 1-q whatever d may be. The
-- holding probability is then a property of the walker and the degree a
-- property of the graph, and the two do not interfere.
--
-- One consequence is worth stating. A graph with unequal degrees has a
-- stationary distribution proportional to degree rather than uniform, so a
-- decorated lattice is reversible but not symmetric. That is correct, and the
-- eigensystem is obtained by conjugating with the square root of the stationary
-- law; it is only the undecorated cells whose stationary law is uniform.
primitivePerturbation :: R -> N -> K -> Primitive -> Matrix
primitivePerturbation q n k prim = (n >< n) [p i j | i <- [0..n-1], j <- [0..n-1]]
  where
    p i j = case prim of
        EdgeAdd a b       -> edgeAddPert q n k a b i j
        DirectedAdd u v   -> directedAddPert q n k u v i j
        EdgeDel a b       -> edgeDelPert q n k a b i j
        WattsStrogatz a b -> wattsStrogatzPert q n k a b i j
        Barrier a b bp    -> barrierPert a b
                                 (q / fromIntegral (2 * k) * (1 - bp)) i j
        Asymmetric a b d  -> asymPert q n k a b d i j
        Reweight a b d    -> reweightPert a b d i j
        Teleport m r      -> teleportPert q n k m r i j

-- | Ring neighbours of a site, on the lattice the caller actually has.
ringNbrs :: N -> K -> Int -> [Int]
ringNbrs n k c =
    [ (c + m) `mod` n | m <- [1 .. k] ] ++ [ (c - m + n * k) `mod` n | m <- [1 .. k] ]

isNbr :: N -> K -> Int -> Int -> Bool
isNbr n k c x = x `elem` ringNbrs n k c

-- | Every matrix here is row-stochastic: the entry at (i, j) is the probability
-- of stepping from i to j, and each row sums to one. A perturbation therefore
-- rewrites the rows of the sites whose neighbourhood it alters and leaves every
-- other row untouched, and each row it rewrites must sum to zero.
--
-- A two-way shortcut. Both endpoints gain one neighbour, so each spreads its
-- mobility over 2k+1 destinations rather than 2k; the holding probability is
-- untouched. Only the rows of the endpoints change, which is why the
-- perturbation is not symmetric: the endpoints have gained degree and their
-- neighbours have not.
edgeAddPert :: R -> N -> K -> Int -> Int -> Int -> Int -> R
edgeAddPert q n k a b i j
    | a == b        = 0
    | isNbr n k a b = 0   -- already adjacent; adding the edge is a no-op
    | i == a        = contribution b
    | i == b        = contribution a
    | otherwise     = 0
  where
    !bigK = fromIntegral (2 * k) :: R
    !shrink = q / (bigK + 1) - q / bigK
    contribution partner
        | j == partner  = q / (bigK + 1)
        | isNbr n k i j = shrink
        | otherwise     = 0

-- | A one-way shortcut. Only the source gains an out-edge, so only its row
-- changes and the destination keeps its degree.
--
-- This breaks detailed balance: the chain is not reversible, has no
-- symmetrising conjugation, and is therefore outside the spectral construction
-- the exact engines use. It is retained for simulation, where nothing requires
-- reversibility, and callers on the exact path must reject it.
directedAddPert :: R -> N -> K -> Int -> Int -> Int -> Int -> R
directedAddPert q n k u v i j
    | u == v          = 0
    | i /= u          = 0
    | j == v          = q / (bigK + 1)
    | isNbr n k u j   = q / (bigK + 1) - q / bigK
    | otherwise       = 0
  where
    !bigK = fromIntegral (2 * k) :: R

-- | Removal of an existing ring edge. Both endpoints lose a neighbour and
-- spread the same mobility over 2k-1 destinations.
edgeDelPert :: R -> N -> K -> Int -> Int -> Int -> Int -> R
edgeDelPert q n k a b i j
    | a == b              = 0
    | not (isNbr n k a b) = 0   -- no such edge; removal is a no-op
    | k < 1               = 0
    | i == a              = contribution b
    | i == b              = contribution a
    | otherwise           = 0
  where
    !bigK = fromIntegral (2 * k) :: R
    !grow = q / (bigK - 1) - q / bigK
    contribution partner
        | j == partner  = negate (q / bigK)
        | isNbr n k i j = grow
        | otherwise     = 0

-- | Rewiring in the sense of Watts and Strogatz: the edge from a to its
-- clockwise neighbour is cut and reattached to a distant site b.
--
-- This is not the sum of a removal and an addition computed independently.
-- Each of those is written against the undecorated lattice, and composing them
-- additively would have the removal shrink a neighbourhood the addition has
-- already enlarged. Rewiring is therefore given directly: the degree of a is
-- unchanged, since it exchanges one neighbour for another; the abandoned
-- neighbour loses one; and the new endpoint gains one.
wattsStrogatzPert :: R -> N -> K -> Int -> Int -> Int -> Int -> R
wattsStrogatzPert q n k a b i j
    | a == b                  = 0
    | not (isNbr n k a aPlus) = 0
    | isNbr n k a b           = 0
    | k < 1                   = 0
    | i == a                  = swapAtSource
    | i == aPlus              = lossAtAbandoned
    | i == b                  = gainAtNew
    | otherwise               = 0
  where
    !aPlus = (a + 1) `mod` n
    !bigK = fromIntegral (2 * k) :: R
    swapAtSource
        | j == aPlus        = negate (q / bigK)
        | j == b            = q / bigK
        | otherwise         = 0
    lossAtAbandoned
        | j == a            = negate (q / bigK)
        | isNbr n k aPlus j = q / (bigK - 1) - q / bigK
        | otherwise         = 0
    gainAtNew
        | j == a            = q / (bigK + 1)
        | isNbr n k b j     = q / (bigK + 1) - q / bigK
        | otherwise         = 0

barrierPert :: Int -> Int -> R -> Int -> Int -> R
barrierPert a b red i j
    | (i, j) == (a, b) || (i, j) == (b, a) = negate red
    | i == j && (i == a || i == b)           = red
    | otherwise                              = 0

-- | A directional bias on one bond, given as a fraction of the weight that bond
-- already carries: the step from a to b is favoured by d and the step back
-- disfavoured by the same fraction. The compensating weight is taken from and
-- given to the remaining bonds of each endpoint rather than to the diagonal,
-- so that the holding probability is untouched and the construction remains
-- admissible where the walker never holds.
--
-- The result is not reversible. A single biased bond on a ring leaves a cycle
-- with net circulation, which no stationary law can balance.
asymPert :: R -> N -> K -> Int -> Int -> R -> Int -> Int -> R
asymPert q n k a b d i j
    | a == b              = 0
    | not (isNbr n k a b) = 0
    | k < 1               = 0
    | i == a              = shift b (negate 1)
    | i == b              = shift a 1
    | otherwise           = 0
  where
    !bigK = fromIntegral (2 * k) :: R
    !amount = d * q / bigK
    shift partner sgn
        | j == partner      = negate sgn * amount
        | isNbr n k i j     = sgn * amount / (bigK - 1)
        | otherwise         = 0

-- | A bare reduction of one bond by an absolute amount, the refused weight
-- returning to the walker. Positive strength removes, as it does for a barrier,
-- so that the two primitives share a sign convention; the amount must not
-- exceed the weight the bond carries.
reweightPert :: Int -> Int -> R -> Int -> Int -> R
reweightPert a b d i j
    | (i, j) == (a, b) || (i, j) == (b, a) = negate d
    | i == j && (i == a || i == b)         = d
    | otherwise                            = 0

-- | Resetting to a fixed site. From anywhere but the target, the walker
-- abandons its step with probability r and reappears at m; otherwise it walks
-- as usual. The whole column is therefore scaled, not merely its diagonal,
-- which is what keeps the column summing to one.
--
-- Resetting has no reversible representation: it carries probability to one
-- site from everywhere and returns none, so detailed balance fails for every
-- stationary law. It is admissible to the simulator and not to the spectral
-- construction.
teleportPert :: R -> N -> K -> Int -> R -> Int -> Int -> R
teleportPert q n k m r i j
    | i == m    = 0
    | j == m    = r * (1 - w0 i m)
    | otherwise = negate r * w0 i j
  where
    w0 x y
        | x == y        = 1 - q
        | isNbr n k x y = q / fromIntegral (2 * k)
        | otherwise     = 0

type QCache = Map.Map (Int, Int) (Complex Double)

cachedQ :: QCache -> Int -> Int -> Complex Double
cachedQ cache src obs = case Map.lookup (src, obs) cache of
    Just v  -> v
    Nothing -> error $ "QCache miss: " ++ show (src, obs)

scalarH :: DefectEntry -> QCache -> Complex Double -> Complex Double
scalarH !d !cache !z =
    let !eta = defXvu d :+ 0
        !u = defU d
        !v = defV d
        !qUU = cachedQ cache u u
        !qVV = cachedQ cache v v
        !qUV = cachedQ cache u v
        !qVU = cachedQ cache v u
    in eta * (qUU + qVV - qUV - qVU) - recip z
{-# INLINE scalarH #-}

scalarHmod :: DefectEntry -> QCache -> Int -> Int
           -> Complex Double -> Complex Double -> Complex Double
scalarHmod !d !cache !n !n0 _z !h11base =
    let !eta = defXvu d :+ 0
        !u = defU d
        !v = defV d
        !dqN  = cachedQ cache u n  - cachedQ cache v n
        !dqN0 = cachedQ cache n0 u - cachedQ cache n0 v
    in h11base - dqN * (eta * dqN0)
{-# INLINE scalarHmod #-}

buildHWith :: BV.Vector DefectEntry -> (Int -> Int -> Complex Double)
           -> Complex Double -> LA.Matrix (Complex Double)
buildHWith defv qAt z = (m LA.>< m) [entry i j | i <- [0..m-1], j <- [0..m-1]]
  where
    !m = BV.length defv

    entry :: Int -> Int -> Complex Double
    entry !i !j =
      let !di = defv `BV.unsafeIndex` i
          !dj = defv `BV.unsafeIndex` j
          !ui = defU di; !vi = defV di
          !uj = defU dj; !vj = defV dj
          !euv = defXuv di :+ 0
          !evu = defXvu di :+ 0
          !dqU = qAt uj ui - qAt vj ui
          !dqV = qAt uj vi - qAt vj vi
          !kronZ = if i == j then recip z else 0
      in evu * dqU - euv * dqV - kronZ

stationaryDist :: Matrix -> [R]
stationaryDist m = LA.toList (LA.scale (1.0 / LA.sumElements piRaw) piRaw)
  where
    !n     = LA.rows m
    !wmI   = LA.tr m - LA.ident n
    !rList = LA.toRows wmI
    !aug   = LA.fromRows (take (n - 1) rList
                          ++ [LA.fromList (replicate n 1.0)])
    !rhs   = LA.fromList (replicate (n - 1) 0 ++ [1.0])
    !piRaw = aug LA.<\> rhs

matIsSymmetric :: Matrix -> Bool
matIsSymmetric m =
    maximum (map abs (LA.toList (LA.flatten (m - LA.tr m)))) < 1e-10

exactEigenvalues :: Matrix -> [R]
exactEigenvalues m = fst (exactEigensystem m)

exactEigensystem :: Matrix -> ([R], Matrix)
exactEigensystem m
    | matIsSymmetric m =
        let (eigVals, vecs) = LA.eigSH (LA.trustSym m)
        in (LA.toList eigVals, vecs)
    | otherwise = symmetrisedEigensystem m

symmetrisedEigensystem :: Matrix -> ([R], Matrix)
symmetrisedEigensystem m = (LA.toList eigVals, vecs)
  where
    piVec      = stationaryDist m
    sqrtPi     = LA.fromList (map sqrt piVec)
    invSqrtPi  = LA.fromList (map (\x -> 1.0 / sqrt x) piVec)
    !dMat      = LA.diag sqrtPi
    !dInvMat   = LA.diag invSqrtPi
    !raw       = dMat LA.<> m LA.<> dInvMat
    !sym       = LA.scale 0.5 (raw + LA.tr raw)
    (eigVals, vecs) = LA.eigSH (LA.trustSym sym)

