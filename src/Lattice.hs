{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Lattice
-- Description : Square-lattice first-encounter machinery. The periodic torus
--               carries closed relative-walk forms; the reflecting box uses the
--               abort convention, whose operator is a Kronecker sum of one
--               dimensional bands, so its eigensystem is one small
--               diagonalisation tensored up; the iterative pair solver applies
--               separable stencils with sparse shortcut corrections instead of
--               dense products; and every simulator runs on the fixed
--               MonteCarlo batch protocol, with shortcut-faithful row-sampling
--               variants.

module Lattice
    ( latIndex
    , latEigQ
    , encounterGFTorus
    , encounterPMFTorus
    , encounterMeanTorus
    , simulateTorus
    , simulateTorusSC
    , reflectMatrix
    , reflectBand1D
    , torusMatrix
    , reflectEigensystem
    , encounterPMFReflect
    , encounterPMFReflectIter
    , encounterPMFTorusIter
    , pairEncounterPMFIter
    , addShortcuts
    , encounterPMFReflectIterSC
    , encounterPMFTorusIterSC
    , encounterMeanReflect
    , simulateReflect
    , simulateReflectSC
    ) where

import Control.DeepSeq (NFData(..))
import Control.Monad (forM_)
import Control.Monad.ST (runST)
import Data.Complex (Complex(..))
import Data.List (foldl1')
import System.Random (StdGen, uniformR, splitGen)
import qualified Data.Vector as BV
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as MV
import Numeric.LinearAlgebra ((><))
import qualified Numeric.LinearAlgebra as LA

import Types (R, Matrix)
import Ring (ringEigNonLazy)
import Distribution (invertPGF)
import Encounter (encounterPMFTwoQ, encounterMeanTwoQ)
import qualified MonteCarlo as MC

latIndex :: Int -> (Int, Int) -> Int
latIndex l (x, y) = x + l * y

latEigQ :: Double -> Int -> Int -> (Int, Int) -> Double
latEigQ !q !l !k (!kx, !ky) =
    (1 - q) + q * 0.5 * (ringEigNonLazy l k kx + ringEigNonLazy l k ky)
{-# INLINE latEigQ #-}

torusRelProp :: Int -> Int -> Double -> Double -> (Int, Int)
             -> Complex Double -> Complex Double
torusRelProp !l !k !q1 !q2 (!dx, !dy) !z = go 0 0 0
  where
    !nn = fromIntegral (l * l) :: Double
    !ld = fromIntegral l :: Double
    go !acc !kx !ky
        | kx >= l   = acc / (nn :+ 0)
        | ky >= l   = go acc (kx + 1) 0
        | otherwise =
            let !lam  = latEigQ q1 l k (kx, ky) * latEigQ q2 l k (kx, ky)
                !cosD = cos (2 * pi * fromIntegral (kx * dx + ky * dy) / ld)
            in go (acc + (cosD :+ 0) / (1 - z * (lam :+ 0))) kx (ky + 1)

encounterGFTorus :: Int -> Int -> Double -> Double -> Double -> (Int, Int)
                 -> Complex Double -> Complex Double
encounterGFTorus l k q1 q2 rho (dx, dy) z =
    let !rhoC  = rho :+ 0
        !propD = torusRelProp l k q1 q2 (dx, dy) z
        !prop0 = torusRelProp l k q1 q2 (0, 0) z
    in rhoC * propD / (((1 - rho) :+ 0) + rhoC * prop0)

encounterPMFTorus :: Int -> Int -> Double -> Double -> Double
                  -> (Int, Int) -> (Int, Int) -> Int -> V.Vector Double
encounterPMFTorus l k q1 q2 rho (ax, ay) (bx, by) tmax =
    let !dx = (bx - ax) `mod` l
        !dy = (by - ay) `mod` l
    in invertPGF tmax (encounterGFTorus l k q1 q2 rho (dx, dy))

encounterMeanTorus :: Int -> Int -> Double -> Double -> Double
                   -> (Int, Int) -> (Int, Int) -> Double
encounterMeanTorus l k q1 q2 rho (ax, ay) (bx, by) =
    let !dx = (bx - ax) `mod` l
        !dy = (by - ay) `mod` l
        !nn = fromIntegral (l * l) :: Double
        !ld = fromIntegral l :: Double
        !s  = sum
            [ (1 - cosD) / (1 - lam)
            | kx <- [0 .. l - 1], ky <- [0 .. l - 1]
            , not (kx == 0 && ky == 0)
            , let !lam  = latEigQ q1 l k (kx, ky) * latEigQ q2 l k (kx, ky)
            , let !cosD = cos (2 * pi * fromIntegral (kx * dx + ky * dy) / ld)
            ]
    in s + nn * (1 - rho) / rho

reflectMatrix :: Int -> Int -> Double -> Matrix
reflectMatrix !bigL !k !q =
    let !n  = bigL * bigL
        !qf = q / (4 * fromIntegral k)
        coord i = (i `mod` bigL, i `div` bigL)
        outC v = max 0 (k - v) + max 0 (k - (bigL - 1 - v))
        hops (x, y) (x', y') =
            length [ () | m <- [1 .. k]
                   , (x + m == x' && y == y')
                     || (x - m == x' && y == y')
                     || (y + m == y' && x == x')
                     || (y - m == y' && x == x') ]
        entry i j =
            let (!x, !y) = coord i
                !h = hops (coord i) (coord j)
            in if i == j
               then (1 - q) + qf * fromIntegral (outC x + outC y)
               else qf * fromIntegral h
    in (n >< n) [ entry i j | i <- [0 .. n - 1], j <- [0 .. n - 1] ]

reflectBand1D :: Int -> Int -> Matrix
reflectBand1D !bigL !k =
    let outC v = max 0 (k - v) + max 0 (k - (bigL - 1 - v))
        entry x x'
            | x == x'                        = fromIntegral (outC x)
            | abs (x - x') <= k              = 1
            | otherwise                      = 0
    in (bigL >< bigL) [ entry x x' | x <- [0 .. bigL - 1], x' <- [0 .. bigL - 1] ]

reflectEigensystem :: Int -> Int -> Double -> ([R], Matrix)
reflectEigensystem !bigL !k !q =
    let !band = reflectBand1D bigL k
        (!bvals, !bvecs) = LA.eigSH (LA.trustSym band)
        !qf   = q / (4 * fromIntegral k)
        !vals = [ (1 - q) + qf * ((bvals `LA.atIndex` jA) + (bvals `LA.atIndex` jB))
                | jA <- [0 .. bigL - 1], jB <- [0 .. bigL - 1] ]
        !vecs = LA.kronecker bvecs bvecs
    in (vals, vecs)

uniformPi :: Int -> [Double]
uniformPi n = replicate n (1 / fromIntegral n)

encounterPMFReflect :: Int -> Int -> Double -> Double -> Double
                    -> (Int, Int) -> (Int, Int) -> Int -> V.Vector Double
encounterPMFReflect l k q1 q2 rho a b tmax =
    let (!eA, !vA) = reflectEigensystem l k q1
        (!eB, !vB) = reflectEigensystem l k q2
        !piV = uniformPi (l * l)
    in encounterPMFTwoQ eA vA eB vB piV rho (latIndex l a, latIndex l b) tmax

encounterMeanReflect :: Int -> Int -> Double -> Double -> Double
                     -> (Int, Int) -> (Int, Int) -> Double
encounterMeanReflect l k q1 q2 rho a b =
    let (!eA, !vA) = reflectEigensystem l k q1
        (!eB, !vB) = reflectEigensystem l k q2
        !piV = uniformPi (l * l)
    in encounterMeanTwoQ eA vA eB vB piV rho (latIndex l a, latIndex l b)

torusMatrix :: Int -> Int -> Double -> Matrix
torusMatrix !bigL !k !q =
    let !n  = bigL * bigL
        !qf = q / (4 * fromIntegral k)
        coord i = (i `mod` bigL, i `div` bigL)
        wrap d  = ((d `mod` bigL) + bigL) `mod` bigL
        hops (x, y) (x', y') =
            length [ () | m <- [1 .. k]
                   , (wrap (x + m) == x' && y == y')
                     || (wrap (x - m) == x' && y == y')
                     || (wrap (y + m) == y' && x == x')
                     || (wrap (y - m) == y' && x == x') ]
        entry i j =
            let !h = hops (coord i) (coord j)
            in if i == j then (1 - q) + qf * fromIntegral h
                         else qf * fromIntegral h
    in (n >< n) [ entry i j | i <- [0 .. n - 1], j <- [0 .. n - 1] ]

pairEncounterPMFIter :: Matrix -> Matrix -> Double -> Int -> Int -> Int
                     -> V.Vector Double
pairEncounterPMFIter !pa !pb !rho !sA !sB !tmax =
    let !n   = LA.rows pa
        !paT = LA.tr pa
        !u0  = (n >< n) [ if i == sA && j == sB then 1 else 0
                        | i <- [0 .. n - 1], j <- [0 .. n - 1] ]
        step !u =
            let !u1  = (paT LA.<> u) LA.<> pb
                !dg  = LA.takeDiag u1
                !met = LA.sumElements dg
                !u2  = u1 - LA.scale rho (LA.diag dg)
            in (u2, rho * met)
        go !t !u !acc
            | t > tmax  = reverse acc
            | otherwise = let (!u', !p) = step u in go (t + 1) u' (p : acc)
    in V.fromList (0 : go 1 u0 [])

addShortcuts :: Double -> Int -> [(Int, Int)] -> Matrix -> Matrix
addShortcuts !q !k edges !m0 =
    let !qf = q / (4 * fromIntegral k)
        !n  = LA.rows m0
        valid (a, b) = a /= b && a >= 0 && b >= 0 && a < n && b < n
        es  = filter valid edges
        deltas = concatMap (\(a, b) ->
                    [ ((a, b),  qf), ((a, a), -qf)
                    , ((b, a),  qf), ((b, b), -qf) ]) es
    in if null es then m0 else LA.accum m0 (+) deltas

data Geom = TorusG | ReflectG
    deriving Eq

scDeltas :: Double -> Int -> Int -> [(Int, Int)] -> [((Int, Int), Double)]
scDeltas !q !k !n edges =
    let !qf = q / (4 * fromIntegral k)
        valid (a, b) = a /= b && a >= 0 && b >= 0 && a < n && b < n
    in concatMap (\(a, b) ->
            [ ((a, b),  qf), ((a, a), -qf)
            , ((b, a),  qf), ((b, b), -qf) ]) (filter valid edges)

neighbourTable :: Geom -> Int -> Int -> (BV.Vector (V.Vector Int), V.Vector Int)
neighbourTable geom !l !k =
    let !n = l * l
        wrap d = ((d `mod` l) + l) `mod` l
        outC v = max 0 (k - v) + max 0 (k - (l - 1 - v))
        nbrsAt i =
            let (!x, !y) = (i `mod` l, i `div` l)
                torus = concat
                    [ [ wrap (x + m) + l * y
                      , wrap (x - m) + l * y
                      , x + l * wrap (y + m)
                      , x + l * wrap (y - m) ]
                    | m <- [1 .. k] ]
                refl = concat
                    [ [ x + m + l * y | x + m <= l - 1 ]
                      ++ [ x - m + l * y | x - m >= 0 ]
                      ++ [ x + l * (y + m) | y + m <= l - 1 ]
                      ++ [ x + l * (y - m) | y - m >= 0 ]
                    | m <- [1 .. k] ]
            in case geom of
                TorusG   -> V.fromList torus
                ReflectG -> V.fromList refl
        ocAt i =
            let (!x, !y) = (i `mod` l, i `div` l)
            in case geom of
                TorusG   -> 0
                ReflectG -> outC x + outC y
    in ( BV.generate n nbrsAt, V.generate n ocAt )

applySide :: Bool -> Int -> Double -> Double
          -> BV.Vector (V.Vector Int) -> V.Vector Int
          -> [((Int, Int), Double)]
          -> V.Vector Double -> V.Vector Double
applySide !isA !n !q !qf !nbrs !oc !deltas !u = V.create $ do
    out <- MV.new (n * n)
    let gather !ns !stride !off =
            let go !acc !i
                    | i >= V.length ns = acc
                    | otherwise =
                        go (acc + u `V.unsafeIndex`
                                ((ns `V.unsafeIndex` i) * stride + off)) (i + 1)
            in go 0 0
    if isA
        then forM_ [0 .. n - 1] $ \a -> do
            let !ns  = nbrs `BV.unsafeIndex` a
                !oca = fromIntegral (oc `V.unsafeIndex` a)
            forM_ [0 .. n - 1] $ \b -> do
                let !self = u `V.unsafeIndex` (a * n + b)
                    !s    = gather ns n b
                MV.unsafeWrite out (a * n + b)
                    ((1 - q) * self + qf * (s + oca * self))
        else forM_ [0 .. n - 1] $ \b -> do
            let !ns  = nbrs `BV.unsafeIndex` b
                !ocb = fromIntegral (oc `V.unsafeIndex` b)
            forM_ [0 .. n - 1] $ \a -> do
                let !self = u `V.unsafeIndex` (a * n + b)
                    !s    = gatherB ns a
                MV.unsafeWrite out (a * n + b)
                    ((1 - q) * self + qf * (s + ocb * self))
    forM_ deltas $ \((!r, !c), !v) ->
        if isA
            then forM_ [0 .. n - 1] $ \b ->
                MV.unsafeModify out (+ v * (u `V.unsafeIndex` (c * n + b))) (r * n + b)
            else forM_ [0 .. n - 1] $ \a ->
                MV.unsafeModify out (+ v * (u `V.unsafeIndex` (a * n + r))) (a * n + c)
    return out
  where
    gatherB !ns !a =
        let go !acc !i
                | i >= V.length ns = acc
                | otherwise =
                    go (acc + u `V.unsafeIndex`
                            (a * n + (ns `V.unsafeIndex` i))) (i + 1)
        in go 0 0

stencilPairPMF :: Geom -> Int -> Int -> Double -> Double
               -> [((Int, Int), Double)] -> [((Int, Int), Double)]
               -> Double -> Int -> Int -> Int -> V.Vector Double
stencilPairPMF geom !l !k !qA !qB !dA !dB !rho !sA !sB !tmax =
    let !n = l * l
        (!nbrs, !oc) = neighbourTable geom l k
        !qfA = qA / (4 * fromIntegral k)
        !qfB = qB / (4 * fromIntegral k)
        !u0 = V.generate (n * n) (\i -> if i == sA * n + sB then 1 else 0)
        step !u =
            let !v  = applySide True  n qA qfA nbrs oc dA u
                !w  = applySide False n qB qfB nbrs oc dB v
                !met = diagSum n w
                !u2 = scaleDiag n rho w
            in (u2, rho * met)
        go !t !u !acc
            | t > tmax  = reverse acc
            | otherwise = let (!u', !p) = step u in go (t + 1) u' (p : acc)
    in V.fromList (0 : go 1 u0 [])

diagSum :: Int -> V.Vector Double -> Double
diagSum !n !u =
    let go !acc !c
            | c >= n    = acc
            | otherwise = go (acc + u `V.unsafeIndex` (c * n + c)) (c + 1)
    in go 0 0

scaleDiag :: Int -> Double -> V.Vector Double -> V.Vector Double
scaleDiag !n !rho !u = V.modify
    (\mv -> forM_ [0 .. n - 1] $ \c ->
        MV.unsafeModify mv (* (1 - rho)) (c * n + c)) u

encounterPMFReflectIter :: Int -> Int -> Double -> Double -> Double
                        -> (Int, Int) -> (Int, Int) -> Int -> V.Vector Double
encounterPMFReflectIter !bigL !k !q1 !q2 !rho a b !tmax =
    stencilPairPMF ReflectG bigL k q1 q2 [] [] rho
        (latIndex bigL a) (latIndex bigL b) tmax

encounterPMFTorusIter :: Int -> Int -> Double -> Double -> Double
                      -> (Int, Int) -> (Int, Int) -> Int -> V.Vector Double
encounterPMFTorusIter !bigL !k !q1 !q2 !rho a b !tmax =
    stencilPairPMF TorusG bigL k q1 q2 [] [] rho
        (latIndex bigL a) (latIndex bigL b) tmax

encounterPMFReflectIterSC :: Int -> Int -> Double -> Double -> Double
                          -> (Int, Int) -> (Int, Int) -> [(Int, Int)] -> Int
                          -> V.Vector Double
encounterPMFReflectIterSC !bigL !k !q1 !q2 !rho a b scs !tmax =
    let !n = bigL * bigL
    in stencilPairPMF ReflectG bigL k q1 q2
        (scDeltas q1 k n scs) (scDeltas q2 k n scs) rho
        (latIndex bigL a) (latIndex bigL b) tmax

encounterPMFTorusIterSC :: Int -> Int -> Double -> Double -> Double
                        -> (Int, Int) -> (Int, Int) -> [(Int, Int)] -> Int
                        -> V.Vector Double
encounterPMFTorusIterSC !bigL !k !q1 !q2 !rho a b scs !tmax =
    let !n = bigL * bigL
    in stencilPairPMF TorusG bigL k q1 q2
        (scDeltas q1 k n scs) (scDeltas q2 k n scs) rho
        (latIndex bigL a) (latIndex bigL b) tmax

type GridStep = Int -> Int -> Double -> (Int, Int) -> StdGen -> ((Int, Int), StdGen)

torusStep :: GridStep
torusStep !l !k !q (!x, !y) !gen =
    let (!r, !gen') = uniformR (0.0 :: Double, 1.0) gen
    in if r < 1 - q
       then ((x, y), gen')
       else let !nb = 4 * k
                !s  = (r - (1 - q)) / q
                !j0 = floor (s * fromIntegral nb) :: Int
                !j  = if j0 >= nb then nb - 1 else j0
                !m  = j `div` 4 + 1
                !p  = case j `mod` 4 of
                        0 -> ((x + m) `mod` l, y)
                        1 -> ((x - m) `mod` l, y)
                        2 -> (x, (y + m) `mod` l)
                        _ -> (x, (y - m) `mod` l)
            in (p, gen')

reflectStep :: GridStep
reflectStep !l !k !q (!x, !y) !gen =
    let (!r, !gen') = uniformR (0.0 :: Double, 1.0) gen
    in if r < 1 - q
       then ((x, y), gen')
       else let !nb = 4 * k
                !s  = (r - (1 - q)) / q
                !j0 = floor (s * fromIntegral nb) :: Int
                !j  = if j0 >= nb then nb - 1 else j0
                !m  = j `div` 4 + 1
                !p  = case j `mod` 4 of
                        0 -> let !x' = x + m in if x' <= l - 1 then (x', y) else (x, y)
                        1 -> let !x' = x - m in if x' >= 0     then (x', y) else (x, y)
                        2 -> let !y' = y + m in if y' <= l - 1 then (x, y') else (x, y)
                        _ -> let !y' = y - m in if y' >= 0     then (x, y') else (x, y)
            in (p, gen')

trialGrid :: GridStep -> Int -> Int -> Double -> Double -> Double -> Int
          -> (Int, Int) -> (Int, Int) -> StdGen -> Int
trialGrid step !l !k !q1 !q2 !rho !maxT !pa0 !pb0 !gen0 =
    let (!gA, !gB) = splitGen gen0
    in go pa0 pb0 0 gA gB
  where
    go !pa !pb !t !gA !gB
        | t >= maxT         = 0
        | t > 0 && pa == pb =
            if rho >= 1.0 - 1e-15
            then t
            else let (!r, !gA') = uniformR (0.0 :: Double, 1.0) gA
                 in if r <= rho then t
                    else let (!pa', !gA'') = step l k q1 pa gA'
                             (!pb', !gB')  = step l k q2 pb gB
                         in go pa' pb' (t + 1) gA'' gB'
        | otherwise =
            let (!pa', !gA') = step l k q1 pa gA
                (!pb', !gB') = step l k q2 pb gB
            in go pa' pb' (t + 1) gA' gB'

data BatchT = BatchT
    { btHist :: !(V.Vector Int)
    , btCnt  :: !Int
    , btMean :: !Double
    , btM2   :: !Double
    }

instance NFData BatchT where
    rnf (BatchT h c m m2) = rnf h `seq` rnf c `seq` rnf m `seq` rnf m2

runBatchWith :: (StdGen -> Int) -> Int -> Int -> StdGen -> BatchT
runBatchWith trial !maxT !batch !gen0 = runST $ do
    hist <- MV.replicate maxT (0 :: Int)
    let go !i !g !cnt !mean !m2
            | i >= batch = do
                frozen <- V.unsafeFreeze hist
                return $! BatchT frozen cnt mean m2
            | otherwise = do
                let (!g1, !g2) = splitGen g
                    !fpt = trial g1
                if fpt > 0 && fpt <= maxT
                    then do
                        MV.unsafeModify hist (+ 1) (fpt - 1)
                        let !cnt'  = cnt + 1
                            !d     = fromIntegral fpt - mean
                            !mean' = mean + d / fromIntegral cnt'
                            !d2    = fromIntegral fpt - mean'
                            !m2'   = m2 + d * d2
                        go (i + 1) g2 cnt' mean' m2'
                    else go (i + 1) g2 cnt mean m2
    go 0 gen0 0 0.0 0.0

collectBatches :: Int -> [BatchT] -> (Double, Double, V.Vector Double)
collectBatches !walkers batches =
    let !rawCounts = foldl1' (V.zipWith (+)) (map btHist batches)
        !hist = V.map (\c -> fromIntegral c / fromIntegral walkers) rawCounts
        (!cnt, !mean, !m2) = foldl1' MC.mergeWelford
            [(btCnt b, btMean b, btM2 b) | b <- batches]
        !var = if cnt > 1 then m2 / fromIntegral (cnt - 1) else 0
        !se  = if cnt > 0 then sqrt var / sqrt (fromIntegral cnt) else 0
    in (mean, se, hist)

simulateGrid :: GridStep -> Int -> Int -> Double -> Double -> Double
             -> (Int, Int) -> (Int, Int) -> Int -> Int -> Int
             -> (Double, Double, V.Vector Double)
simulateGrid step !l !k !q1 !q2 !rho !pa0 !pb0 !walkers !seed !maxT =
    let !batches = MC.runBatches seed walkers
            (\m g -> runBatchWith (trialGrid step l k q1 q2 rho maxT pa0 pb0)
                                  maxT m g)
    in collectBatches walkers batches

simulateTorus :: Int -> Int -> Double -> Double -> Double
              -> (Int, Int) -> (Int, Int) -> Int -> Int -> Int
              -> (Double, Double, V.Vector Double)
simulateTorus = simulateGrid torusStep

simulateReflect :: Int -> Int -> Double -> Double -> Double
                -> (Int, Int) -> (Int, Int) -> Int -> Int -> Int
                -> (Double, Double, V.Vector Double)
simulateReflect = simulateGrid reflectStep

type Rows = BV.Vector (V.Vector Double)

rowsOf :: Matrix -> Rows
rowsOf w = BV.fromList [ V.fromList (LA.toList r) | r <- LA.toRows w ]

stepRow :: Rows -> Int -> StdGen -> (Int, StdGen)
stepRow !rows !pos !gen =
    let (!r, !gen') = uniformR (0.0 :: Double, 1.0) gen
        !row = rows `BV.unsafeIndex` pos
        scan !j !cumul
            | j >= V.length row - 1 = j
            | otherwise =
                let !cumul' = cumul + row `V.unsafeIndex` j
                in if r <= cumul' then j else scan (j + 1) cumul'
    in (scan 0 0.0, gen')

trialRows :: Rows -> Rows -> Double -> Int -> Int -> Int -> StdGen -> Int
trialRows !rowsA !rowsB !rho !maxT !ia0 !ib0 !gen0 =
    let (!gA, !gB) = splitGen gen0
    in go ia0 ib0 0 gA gB
  where
    go !ia !ib !t !gA !gB
        | t >= maxT         = 0
        | t > 0 && ia == ib =
            if rho >= 1.0 - 1e-15
            then t
            else let (!r, !gA') = uniformR (0.0 :: Double, 1.0) gA
                 in if r <= rho then t
                    else let (!ia', !gA'') = stepRow rowsA ia gA'
                             (!ib', !gB')  = stepRow rowsB ib gB
                         in go ia' ib' (t + 1) gA'' gB'
        | otherwise =
            let (!ia', !gA') = stepRow rowsA ia gA
                (!ib', !gB') = stepRow rowsB ib gB
            in go ia' ib' (t + 1) gA' gB'

simulateRows :: Rows -> Rows -> Double -> Int -> Int -> Int -> Int -> Int
             -> (Double, Double, V.Vector Double)
simulateRows !rowsA !rowsB !rho !ia !ib !walkers !seed !maxT =
    let !batches = MC.runBatches seed walkers
            (\m g -> runBatchWith (trialRows rowsA rowsB rho maxT ia ib)
                                  maxT m g)
    in collectBatches walkers batches

simulateTorusSC :: Int -> Int -> Double -> Double -> Double
                -> (Int, Int) -> (Int, Int) -> [(Int, Int)] -> Int -> Int -> Int
                -> (Double, Double, V.Vector Double)
simulateTorusSC !l !k !q1 !q2 !rho a b scs !walkers !seed !maxT =
    let !mA = addShortcuts q1 k scs (torusMatrix l k q1)
        !mB = addShortcuts q2 k scs (torusMatrix l k q2)
    in simulateRows (rowsOf mA) (rowsOf mB) rho
        (latIndex l a) (latIndex l b) walkers seed maxT

simulateReflectSC :: Int -> Int -> Double -> Double -> Double
                  -> (Int, Int) -> (Int, Int) -> [(Int, Int)] -> Int -> Int -> Int
                  -> (Double, Double, V.Vector Double)
simulateReflectSC !l !k !q1 !q2 !rho a b scs !walkers !seed !maxT =
    let !mA = addShortcuts q1 k scs (reflectMatrix l k q1)
        !mB = addShortcuts q2 k scs (reflectMatrix l k q2)
    in simulateRows (rowsOf mA) (rowsOf mB) rho
        (latIndex l a) (latIndex l b) walkers seed maxT
