{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Diagnostics
-- Description : Interactive diagnostics for the mobility sensitivity of the
--               mean first-encounter time at perfect absorption. Every site
--               index in this module is one-based, matching the write-up rather
--               than the engine: sites run from one to N, and the single
--               conversion to the engine's zero-based indexing happens in
--               toEngine, applied once to the primitives and once to the
--               starting pair. Nothing else subtracts or adds one, and site
--               labels reported by siteDetail are one-based throughout.
--               Provides the exact decomposition at a mobility, the midpoint
--               between-site and within-site contributions over any partition
--               of the contact set, a step scan for locating the usable
--               differencing step, a mobility sweep carrying the least
--               splitting weight alongside the contributions, a site-resolved
--               table of the individual terms forming the aggregate balance,
--               and a survey reducing each of many configurations to a single
--               line so that a claim about one cell can be checked against a
--               family.

module Diagnostics
    ( Cell(..)
    , Turn(..)
    , ring
    , grid
    , decomposition
    , contributions
    , sites
    , whole
    , halves
    , leastEig
    , netAt
    , turningPoints
    , stepScan
    , sweep
    , blockSweep
    , siteDetail
    , survey
    ) where

import qualified Data.Vector.Unboxed as V
import Text.Printf (printf)

import Types (Primitive(..))
import Defect (buildTransitionMatrix, exactEigensystem, stationaryDist)
import Encounter (encounterDecompositionTwoQ)

data Cell = Cell
    { cellN      :: !Int
    , cellK      :: !Int
    , cellPrims  :: ![Primitive]
    , cellQ1     :: !Double
    , cellStarts :: !(Int, Int)
    }

data Turn = Turn
    { turnQ      :: !Double
    , turnRising :: !Bool
    , turnA      :: !Double
    , turnB      :: !Double
    }

ring :: Int -> Int -> Double -> (Int, Int) -> Cell
ring n k q1 starts = Cell n k [] q1 starts

grid :: Double -> Double -> Int -> [Double]
grid lo hi n =
    [ lo + (hi - lo) * fromIntegral i / fromIntegral n | i <- [0 .. n] ]

check :: Int -> String -> Int -> Int
check n what s
    | s >= 1 && s <= n = s - 1
    | otherwise = error (what ++ ": site " ++ show s
                    ++ " outside 1 to " ++ show n ++ "; this module is one-based")

toEngine :: Int -> Primitive -> Primitive
toEngine n p = case p of
    EdgeAdd a b       -> EdgeAdd (c a) (c b)
    DirectedAdd a b   -> DirectedAdd (c a) (c b)
    EdgeDel a b       -> EdgeDel (c a) (c b)
    WattsStrogatz a b -> WattsStrogatz (c a) (c b)
    Barrier a b x     -> Barrier (c a) (c b) x
    Asymmetric a b x  -> Asymmetric (c a) (c b) x
    Reweight a b x    -> Reweight (c a) (c b) x
    Teleport m x      -> Teleport (c m) x
  where
    c = check n "primitive"

type Decomp = (Double, V.Vector Double, V.Vector Double)

primitivesFor :: Cell -> [Primitive]
primitivesFor c = map (toEngine (cellN c)) (cellPrims c)

startsFor :: Cell -> (Int, Int)
startsFor c =
    let (a, b) = cellStarts c
        n = cellN c
    in (check n "start" a, check n "start" b)

decomposition :: Cell -> Double -> Decomp
decomposition c q2 =
    let !ps = primitivesFor c
        !wA = buildTransitionMatrix (cellQ1 c) (cellN c) (cellK c) ps
        !wB = buildTransitionMatrix q2 (cellN c) (cellK c) ps
        (!eA, !vA) = exactEigensystem wA
        (!eB, !vB) = exactEigensystem wB
        !piV = stationaryDist wA
    in encounterDecompositionTwoQ eA vA eB vB piV 1.0 (startsFor c)

leastEig :: Cell -> Double
leastEig c = minimum (fst (exactEigensystem
    (buildTransitionMatrix 1.0 (cellN c) (cellK c) (primitivesFor c))))

blockTotal :: [Int] -> V.Vector Double -> Double
blockTotal b v = sum [ v V.! (j - 1) | j <- b ]

blockMean :: [Int] -> V.Vector Double -> V.Vector Double -> Double
blockMean b w m =
    let !p = blockTotal b w
    in if p <= 0 then 0 else blockTotal b m / p

contributions :: Cell -> [[Int]] -> Double -> Double -> (Double, Double)
contributions c bs q2 h =
    let (_, !wLo, !mLo) = decomposition c (q2 - h)
        (_, !wHi, !mHi) = decomposition c (q2 + h)
        pieceA b = ((blockTotal b wHi - blockTotal b wLo) / (2 * h))
                 * 0.5 * (blockMean b wHi mHi + blockMean b wLo mLo)
        pieceB b = 0.5 * (blockTotal b wHi + blockTotal b wLo)
                 * ((blockMean b wHi mHi - blockMean b wLo mLo) / (2 * h))
    in (sum (map pieceA bs), sum (map pieceB bs))

sites :: Cell -> [[Int]]
sites c = [ [j] | j <- [1 .. cellN c] ]

whole :: Cell -> [[Int]]
whole c = [ [1 .. cellN c] ]

halves :: Cell -> [[Int]]
halves c =
    let !n = cellN c
    in [ [1 .. n `div` 2], [n `div` 2 + 1 .. n] ]

netAt :: Cell -> Double -> Double -> Double
netAt c h q2 =
    let (!a, !b) = contributions c (sites c) q2 h in a + b

refineTurn :: Cell -> Double -> (Double, Double) -> (Double, Double) -> Int
           -> Double
refineTurn c h (lo, vLo) (hi, vHi) i
    | i <= 0 || hi - lo <= 0 =
        lo + (hi - lo) * abs vLo / max 1e-300 (abs vLo + abs vHi)
    | otherwise =
        let !mid  = 0.5 * (lo + hi)
            !vMid = netAt c h mid
        in if vLo * vMid <= 0
           then refineTurn c h (lo, vLo) (mid, vMid) (i - 1)
           else refineTurn c h (mid, vMid) (hi, vHi) (i - 1)

turningPoints :: Cell -> Double -> [Double] -> [Turn]
turningPoints c h qs =
    let vals  = map (netAt c h) qs
        spans = zip (zip qs vals) (zip (drop 1 qs) (drop 1 vals))
    in [ mk lo hi | (lo, hi) <- spans, snd lo * snd hi < 0 ]
  where
    mk lo hi =
        let !q = refineTurn c h lo hi 14
            (!a, !b) = contributions c (sites c) q h
        in Turn q (snd lo < 0) a b

stepScan :: Cell -> Double -> [Double] -> IO ()
stepScan c q2 hs = do
    putStrLn (printf "%10s %16s %16s %16s" "h" "A" "B" "A+B" :: String)
    mapM_ row hs
  where
    row :: Double -> IO ()
    row h =
        let (!a, !b) = contributions c (sites c) q2 h
        in putStrLn (printf "%10.1e %16.9f %16.9f %16.9f"
                            h a b (a + b) :: String)

sweep :: Cell -> Double -> [Double] -> IO ()
sweep c h qs = do
    putStrLn (printf "%8s %12s %12s %12s %12s %12s"
                     "qB" "E" "minPhi" "A" "B" "A+B" :: String)
    mapM_ row qs
  where
    row :: Double -> IO ()
    row q2 =
        let (!e, !w, _) = decomposition c q2
            (!a, !b)    = contributions c (sites c) q2 h
        in putStrLn (printf "%8.4f %12.6f %12.3e %12.6f %12.6f %12.6f"
                            q2 e (V.minimum w) a b (a + b) :: String)

blockSweep :: Cell -> [[Int]] -> Double -> [Double] -> IO ()
blockSweep c bs h qs = do
    putStrLn (printf "%8s %14s %14s %14s" "qB" "A" "B" "A+B" :: String)
    mapM_ row qs
  where
    row :: Double -> IO ()
    row q2 =
        let (!a, !b) = contributions c bs q2 h
        in putStrLn (printf "%8.4f %14.6f %14.6f %14.6f"
                            q2 a b (a + b) :: String)

siteDetail :: Cell -> Double -> Double -> IO ()
siteDetail c q2 h = do
    let (!e, !w, _)     = decomposition c q2
        (_, !wLo, !mLo) = decomposition c (q2 - h)
        (_, !wHi, !mHi) = decomposition c (q2 + h)
    putStrLn (printf "qB %.4f   E %.6f   sites are one-based" q2 e :: String)
    putStrLn (printf "%6s %12s %14s %14s %14s"
                     "site" "Phi" "Ej-E" "a_j" "b_j" :: String)
    mapM_ (row e w wLo mLo wHi mHi) [1 .. cellN c]
  where
    row :: Double -> V.Vector Double -> V.Vector Double -> V.Vector Double
        -> V.Vector Double -> V.Vector Double -> Int -> IO ()
    row e w wLo mLo wHi mHi j =
        let !pLo  = wLo V.! (j - 1)
            !pHi  = wHi V.! (j - 1)
            !eLo  = blockMean [j] wLo mLo
            !eHi  = blockMean [j] wHi mHi
            !cond = 0.5 * (eHi + eLo)
            !aj   = ((pHi - pLo) / (2 * h)) * cond
            !bj   = 0.5 * (pHi + pLo) * ((eHi - eLo) / (2 * h))
        in putStrLn (printf "%6d %12.6f %14.6f %14.6f %14.6f"
                            j (w V.! (j - 1)) (cond - e) aj bj :: String)

survey :: Double -> [Double] -> [(String, Cell)] -> IO ()
survey h qs cells = do
    putStrLn (printf "%-22s %4s %3s %6s %8s %10s %10s %10s %7s %5s %9s %10s %10s"
                     "cell" "N" "k" "q1" "lam-" "E(lo)" "E(hi)" "minPhi"
                     "cancel" "turns" "q*" "A*" "B*" :: String)
    mapM_ row cells
  where
    row :: (String, Cell) -> IO ()
    row (label, c) =
        let !lam  = leastEig c
            means = [ let (e, _, _) = decomposition c q in e | q <- qs ]
            phis  = [ let (_, w, _) = decomposition c q in V.minimum w
                    | q <- qs ]
            !qLo  = case qs of { (x : _) -> x ; [] -> error "survey: empty grid" }
            (!aLo, !bLo) = contributions c (sites c) qLo h
            !canc = (abs aLo + abs bLo) / max 1e-300 (abs (aLo + bLo))
            turns = turningPoints c h qs
            !nT   = length turns
            lead  = printf "%-22s %4d %3d %6.2f %8.4f %10.4f %10.4f %10.2e %7.2f %5d"
                           label (cellN c) (cellK c) (cellQ1 c) lam
                           (head means) (last means) (minimum phis)
                           canc nT :: String
            rest  = case turns of
                []      -> printf " %9s %10s %10s" "-" "-" "-" :: String
                (t : _) -> printf " %9.5f %10.5f %10.5f  %s"
                                  (turnQ t) (turnA t) (turnB t)
                                  (if turnRising t then "min" else "MAX")
                           :: String
        in putStrLn (lead ++ rest)
