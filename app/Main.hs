{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Main
-- Description : CLI entry point for encounters. First-passage and
--               first-encounter, exact and simulated, with two-mobility
--               (qA, qB) support, validated configuration, a wired accuracy
--               knob, closed relative-walk routing on the undecorated ring,
--               and the per contact site decomposition of the encounter law.

module Main (main) where

import System.Environment (getArgs)
import System.IO (hPutStrLn, stderr)
import qualified Data.Vector.Unboxed as V

import qualified Types as T
import Types (Primitive(..), Domain(..), EncounterEngine(..),
              selectEncounterEngine, reversible)
import Cells (Dim(..), BC(..), CellResult(..), encounterCell,
              closedSpectrumExists, cellMatrixWith, cellIndex)
import Defect (primitivesToDefects, buildTransitionMatrix, exactEigensystem,
               stationaryDist)
import Passage (firstAbsorptionGF, pureRingFPGF, exactMFPTq, meanFromGF)
import Distribution (invertPGFAcc, Modality(..), modality)
import Estimate (correctedFromPMF, selfConsistency, censoredFraction, tailScale, truncatedSum)
import Encounter (encounterPMFTwoQAcc, encounterDecompositionTwoQ,
                  encounterSitePMFTwoQAcc,
                  encounterPMFRingTwoQ, encounterMeanRingTwoQ)
import Simulate (SimMode(..), SimConfig(..), SimResult(..), defaultSimConfig,
                 simulate)
import Serialise
import Cli (readFlag, splitOn, checkRange, checkUnit, checkMaybe, orDie,
            writeOut)

data Engine = ExactOnly | SimOnly | Both
    deriving (Show, Eq)

-- | Which products a run is asked for. The mean is available directly from the
-- generating function at the cost of two evaluations, whereas the distribution
-- requires a contour quadrature over thousands of nodes; a sweep that plots
-- only the mean therefore has no reason to form the distribution at all, and
-- asking for it is the difference between seconds and quarter-hours per cell.
--
-- The per contact site distributions are a third level. They cost no
-- additional solves, the renewal system being solved at every node whether or
-- not its solution is summed, but they are one array per lattice site and so
-- multiply the output by the size of the contact set. The splitting
-- probabilities and the per-site contributions to the mean are not behind this
-- flag: they come from the pole expansion, need no contour at all, and are
-- written on every exact run.
data Want = WantMean | WantAll | WantSitePMF
    deriving (Show, Eq)

parseDim :: String -> Dim
parseDim "1" = D1
parseDim "2" = D2
parseDim x   = error ("Unknown --dim (expected 1|2): " ++ x)

parseBC :: String -> BC
parseBC "periodic"   = Periodic
parseBC "reflecting" = Reflecting
parseBC "absorbing"  = Absorbing
parseBC x = error
    ("Unknown --bc (expected periodic|reflecting|absorbing): " ++ x)

-- | A site, given as a scalar in one dimension and as X:Y in two.
parseSite :: String -> String -> (Int, Int)
parseSite flag s = case break (== ':') s of
    (x, ':' : y) -> (readFlag flag x, readFlag flag y)
    (x, _)       -> (readFlag flag x, 0)

cellTag :: Dim -> BC -> String
cellTag d b = dimTag d ++ bcTag b
  where
    dimTag D1 = "1"
    dimTag D2 = "2"
    bcTag Periodic   = "P"
    bcTag Reflecting = "R"
    bcTag Absorbing  = "A"

data Config = Config
    { cfgN       :: !Int
    , cfgK       :: !Int
    , cfgQ       :: !Double
    , cfgQB      :: !(Maybe Double)
    , cfgRho     :: !Double
    , cfgSrc     :: !Int
    , cfgTgt     :: !(Maybe Int)
    , cfgSrcB    :: !(Maybe Int)
    , cfgTmax    :: !Int
    , cfgAcc     :: !Double
    , cfgMode    :: !SimMode
    , cfgDefects :: ![Primitive]
    , cfgSim     :: !Int
    , cfgEngine  :: !Engine
    , cfgWant    :: !Want
    , cfgDim     :: !Dim
    , cfgBC      :: !BC
    , cfgSrcY    :: !Int
    , cfgSrcBY   :: !Int
    , cfgSeed    :: !Int
    , cfgOut     :: !(Maybe FilePath)
    }

defaultConfig :: Config
defaultConfig = Config
    { cfgN = 96, cfgK = 4, cfgQ = 0.75, cfgQB = Nothing
    , cfgRho = 1.0, cfgSrc = 0, cfgTgt = Nothing, cfgSrcB = Nothing
    , cfgTmax = 2000, cfgAcc = 14.0
    , cfgMode = Passage, cfgDefects = [], cfgSim = 0
    , cfgEngine = ExactOnly, cfgWant = WantAll
    , cfgDim = D1, cfgBC = Periodic, cfgSrcY = 0, cfgSrcBY = 0
    , cfgSeed = 42
    , cfgOut = Nothing
    }

parseWant :: String -> Want
parseWant "mean"    = WantMean
parseWant "all"     = WantAll
parseWant "sitepmf" = WantSitePMF
parseWant x = error ("Unknown --want (expected mean|all|sitepmf): " ++ x)

parseArgs :: [String] -> Config -> Config
parseArgs [] c = c
parseArgs ("--passage"   : rest) c = parseArgs rest c { cfgMode = Passage }
parseArgs ("--encounter" : rest) c = parseArgs rest c { cfgMode = Encounter }
parseArgs ("--n"     : v : rest) c = parseArgs rest c { cfgN = readFlag "--n" v }
parseArgs ("--k"     : v : rest) c = parseArgs rest c { cfgK = readFlag "--k" v }
parseArgs ("--q"     : v : rest) c = parseArgs rest c { cfgQ = readFlag "--q" v }
parseArgs ("--qB"    : v : rest) c = parseArgs rest c { cfgQB = Just (readFlag "--qB" v) }
parseArgs ("--rho"   : v : rest) c = parseArgs rest c { cfgRho = readFlag "--rho" v }
parseArgs ("--src"   : v : rest) c =
    let (x, y) = parseSite "--src" v
    in parseArgs rest c { cfgSrc = x, cfgSrcY = y }
parseArgs ("--tgt"   : v : rest) c = parseArgs rest c { cfgTgt = Just (readFlag "--tgt" v) }
parseArgs ("--srcB"  : v : rest) c =
    let (x, y) = parseSite "--srcB" v
    in parseArgs rest c { cfgSrcB = Just x, cfgSrcBY = y }
parseArgs ("--want"   : v : rest) c = parseArgs rest c { cfgWant = parseWant v }
parseArgs ("--dim"    : v : rest) c = parseArgs rest c { cfgDim = parseDim v }
parseArgs ("--bc"     : v : rest) c = parseArgs rest c { cfgBC = parseBC v }
parseArgs ("--L"      : v : rest) c = parseArgs rest c { cfgN = readFlag "--L" v }
parseArgs ("--no-pmf" : rest)     c = parseArgs rest c { cfgWant = WantMean }
parseArgs ("--tmax"  : v : rest) c = parseArgs rest c { cfgTmax = readFlag "--tmax" v }
parseArgs ("--acc"   : v : rest) c = parseArgs rest c { cfgAcc = readFlag "--acc" v }
parseArgs ("--sim"   : v : rest) c = parseArgs rest c { cfgSim = readFlag "--sim" v }
parseArgs ("--engine" : v : rest) c = parseArgs rest c { cfgEngine = parseEngine v }
parseArgs ("--seed"  : v : rest) c = parseArgs rest c { cfgSeed = readFlag "--seed" v }
parseArgs ("--out"   : v : rest) c = parseArgs rest c { cfgOut = Just v }
parseArgs ("--sc"       : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseEdgeAdd v] }
parseArgs ("--directed" : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseDirectedAdd v] }
parseArgs ("--ws"       : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseWS v] }
parseArgs ("--del"      : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseEdgeDel v] }
parseArgs ("--barrier"  : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseBarrier v] }
parseArgs ("--asym"     : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseAsym v] }
parseArgs ("--reset"    : v : rest) c = parseArgs rest c { cfgDefects = cfgDefects c ++ [parseTeleport v] }
parseArgs (x : _) _ = error $ "Unknown argument: " ++ x

parseEngine :: String -> Engine
parseEngine "exact" = ExactOnly
parseEngine "sim"   = SimOnly
parseEngine "both"  = Both
parseEngine x       = error $ "Unknown engine (want exact|sim|both): " ++ x

parseEdgeAdd :: String -> Primitive
parseEdgeAdd s = case break (== ':') s of
    (a, ':' : b) -> EdgeAdd (readFlag "--sc" a) (readFlag "--sc" b)
    _            -> error $ "Expected U:V, got: " ++ s

parseDirectedAdd :: String -> Primitive
parseDirectedAdd s = case break (== ':') s of
    (a, ':' : b) -> DirectedAdd (readFlag "--directed" a) (readFlag "--directed" b)
    _            -> error $ "Expected U:V, got: " ++ s

parseWS :: String -> Primitive
parseWS s = case break (== ':') s of
    (a, ':' : b) -> WattsStrogatz (readFlag "--ws" a) (readFlag "--ws" b)
    _            -> error $ "Expected U:V, got: " ++ s

parseEdgeDel :: String -> Primitive
parseEdgeDel s = case break (== ':') s of
    (a, ':' : b) -> EdgeDel (readFlag "--del" a) (readFlag "--del" b)
    _            -> error $ "Expected U:V, got: " ++ s

parseBarrier :: String -> Primitive
parseBarrier s = case splitOn ':' s of
    [a, b, p] -> Barrier (readFlag "--barrier" a) (readFlag "--barrier" b)
                         (readFlag "--barrier" p)
    _         -> error $ "Expected U:V:P, got: " ++ s

parseAsym :: String -> Primitive
parseAsym s = case splitOn ':' s of
    [a, b, d] -> Asymmetric (readFlag "--asym" a) (readFlag "--asym" b)
                            (readFlag "--asym" d)
    _         -> error $ "Expected U:V:D, got: " ++ s

parseTeleport :: String -> Primitive
parseTeleport s = case break (== ':') s of
    (m, ':' : r) -> Teleport (readFlag "--reset" m) (readFlag "--reset" r)
    _            -> error $ "Expected M:R, got: " ++ s

validateConfig :: Config -> Either String Config
validateConfig c = do
    _ <- checkRange "--n" 3 1000000 (cfgN c)
    _ <- checkRange "--k" 1 ((cfgN c - 1) `div` 2) (cfgK c)
    _ <- checkUnit "--q" (cfgQ c)
    _ <- checkMaybe (checkUnit "--qB") (cfgQB c)
    _ <- checkUnit "--rho" (cfgRho c)
    _ <- checkRange "--src" 0 (cfgN c - 1) (cfgSrc c)
    _ <- checkMaybe (checkRange "--tgt" 0 (cfgN c - 1)) (cfgTgt c)
    _ <- checkMaybe (checkRange "--srcB" 0 (cfgN c - 1)) (cfgSrcB c)
    _ <- checkRange "--tmax" 1 100000000 (cfgTmax c)
    _ <- checkRange "--acc" 1.0 30.0 (cfgAcc c)
    _ <- checkRange "--sim" 0 1000000000 (cfgSim c)
    mapM_ checkPrimitive (cfgDefects c)
    mapM_ checkReversible (cfgDefects c)
    _ <- checkRange "--src y" 0 (cfgN c - 1) (cfgSrcY c)
    _ <- checkRange "--srcB y" 0 (cfgN c - 1) (cfgSrcBY c)
    return c
  where
    -- The spectral construction requires a reversible chain: only then is the
    -- operator similar to a symmetric one under conjugation by the square root
    -- of the stationary law, which is what supplies the real spectrum and
    -- orthonormal basis it uses. A directed edge, a biased bond and a resetting
    -- rule each break detailed balance and belong to the simulator alone.
    checkReversible p
        | reversible p = Right ()
        | cfgEngine c == SimOnly = Right ()
        | otherwise = Left (show p ++ " is not reversible and has no spectral \
                            \representation; run it with --engine sim")
    inRange x = x >= 0 && x < cfgN c
    site name x
        | inRange x = Right x
        | otherwise = Left (name ++ " site " ++ show x ++ " outside [0, "
                            ++ show (cfgN c - 1) ++ "]")
    checkPrimitive p = case p of
        EdgeAdd a b       -> site "--sc" a >> site "--sc" b >> Right ()
        DirectedAdd a b   -> site "--directed" a >> site "--directed" b >> Right ()
        EdgeDel a b       -> site "--del" a >> site "--del" b >> Right ()
        WattsStrogatz a b -> site "--ws" a >> site "--ws" b >> Right ()
        Barrier a b _     -> site "--barrier" a >> site "--barrier" b >> Right ()
        Asymmetric a b _  -> site "--asym" a >> site "--asym" b >> Right ()
        Reweight a b _    -> site "reweight" a >> site "reweight" b >> Right ()
        Teleport m _      -> site "--reset" m >> Right ()

main :: IO ()
main = do
    args <- getArgs
    cfg <- orDie (validateConfig (parseArgs args defaultConfig))
    case cfgMode cfg of
        Passage   -> runPassage cfg
        Encounter -> runEncounter cfg

runPassage :: Config -> IO ()
runPassage cfg = do
    let !n    = cfgN cfg
        !k    = cfgK cfg
        !q    = cfgQ cfg
        !rho  = cfgRho cfg
        !src  = cfgSrc cfg
        !tgt  = maybe (n `div` 2) id (cfgTgt cfg)
        !tmax = cfgTmax cfg
        !acc  = cfgAcc cfg
        !defs = primitivesToDefects q n k (cfgDefects cfg)

    hPutStrLn stderr $ "Passage: N=" ++ show n ++ " K=" ++ show (2*k)
        ++ " q=" ++ show q ++ " rho=" ++ show rho
        ++ " src=" ++ show src ++ " tgt=" ++ show tgt

    let !wantPMF = cfgWant cfg /= WantMean
        !gf = firstAbsorptionGF q n k defs src tgt rho
        !pmf = if wantPMF then invertPGFAcc acc tmax gf else V.empty
        !ringGF = pureRingFPGF q n k src tgt
        !ringPMF = if wantPMF then invertPGFAcc acc tmax ringGF else V.empty
        !mfptNet = if wantPMF then correctedFromPMF pmf else meanFromGF gf
        !mfptResid = if wantPMF then snd (selfConsistency mfptNet pmf) else 0
        !mfptRing = exactMFPTq q n k src tgt
        !peaks = mdPeaks (modality pmf)

    simResult <- if cfgSim cfg > 0
        then do
            let sc = defaultSimConfig
                    { simN = n, simK = k, simQ = q, simQB = q
                    , simSrc = src, simTgt = tgt
                    , simDefects = cfgDefects cfg
                    , simWalkers = cfgSim cfg, simSeed = cfgSeed cfg
                    , simMaxT = tmax, simMode = Simulate.Passage
                    , simRho = rho
                    }
            let !sr = simulate sc
            hPutStrLn stderr $ "Sim MFPT (tail-corrected): " ++ show (srMeanFPT sr)
                ++ "  conditional: " ++ show (srConditionalMean sr)
                ++ "  censored: " ++ show (srCensoredFraction sr)
            return (Just sr)
        else return Nothing

    let bundle = PassageBundle
            { pbN = n, pbK = 2 * k, pbQ = q
            , pbSrc = src, pbTgt = tgt, pbRho = rho
            , pbDefects = show (cfgDefects cfg)
            , pbPMF = pmf, pbPureRingPMF = ringPMF
            , pbMFPTNetwork = mfptNet, pbMFPTRingExact = mfptRing
            , pbMFPTTruncated = if wantPMF then truncatedSum pmf else 0
            , pbCensoredFraction = if wantPMF then censoredFraction pmf else 0
            , pbTailScale = if wantPMF then tailScale pmf else 0
            , pbSelfConsistency = mfptResid
            , pbMeanProvenance =
                if wantPMF then "tail-corrected" else "extrapolated"
            , pbSplittingWeight = 1.0
            , pbPeaks = peaks
            , pbSimMFPT = fmap srMeanFPT simResult
            , pbSimCondMean = fmap srConditionalMean simResult
            , pbSimCensored = fmap srCensoredFraction simResult
            , pbSimTailScale = fmap srTailScale simResult
            , pbSimStdErr = fmap srStdErr simResult
            , pbSimHistogram = fmap srHistogram simResult
            }

    let json = renderJSONPretty 0 (exportPassage bundle)
    case cfgOut cfg of
        Nothing   -> putStrLn json
        Just path -> writeOut path json

runEncounter :: Config -> IO ()
runEncounter cfg = do
    let !n    = cfgN cfg
        !k    = cfgK cfg
        !q    = cfgQ cfg
        !qB   = maybe q id (cfgQB cfg)
        !rho  = cfgRho cfg
        !srcA = cfgSrc cfg
        !srcB = maybe (n `div` 2) id (cfgSrcB cfg)
        !tmax = cfgTmax cfg
        !acc  = cfgAcc cfg
        !domain = Domain 1 T.Periodic (cfgDefects cfg) []
        !engine = selectEncounterEngine domain
        !routed = engine == ClosedRelativeRing && srcA /= srcB

    hPutStrLn stderr $ "Encounter: N=" ++ show n ++ " K=" ++ show (2*k)
        ++ " q=" ++ show q ++ " qB=" ++ show qB ++ " rho=" ++ show rho
        ++ " srcA=" ++ show srcA ++ " srcB=" ++ show srcB
        ++ " cell=" ++ cellTag (cfgDim cfg) (cfgBC cfg)
        ++ " engine=" ++ show (cfgEngine cfg)
        ++ (if cfgDim cfg == D1 && cfgBC cfg == Periodic && routed
            then " route=closed-relative-ring" else " route=generic-pair")
        ++ (if closedSpectrumExists (cfgBC cfg) k
            then "" else " spectrum=numerical")

    let !runExact   = cfgEngine cfg /= SimOnly
        !runSim     = cfgEngine cfg /= ExactOnly || cfgSim cfg > 0
        !nWalkers   = if runSim && cfgSim cfg <= 0 then 1000000 else cfgSim cfg

    let !wantPMF     = cfgWant cfg /= WantMean
        !wantSitePMF = cfgWant cfg == WantSitePMF

    -- Cells other than the undecorated periodic ring are dispatched by
    -- dimension and boundary condition together, which is the pair the
    -- programme is organised by; the ring keeps its own route because the
    -- relative-coordinate reduction is available there and nowhere else.
    let !isRing = cfgDim cfg == D1 && cfgBC cfg == Periodic
        !startA = (srcA, cfgSrcY cfg)
        !startB = (srcB, cfgSrcBY cfg)
        !horizon = if wantPMF then tmax else 0

    -- The reported quantities, the decomposition by contact site, and the
    -- route that produced the mean. The provenance is carried out of the
    -- branch that computed it rather than reconstructed from the
    -- configuration afterwards, so a value cannot be labelled by a route it
    -- did not take.
    let (!pmf, !mfpt, !weight, !peaks, !splits, !siteMeans, !provenance)
            | runExact && isRing && routed =
                let !p  = if wantPMF
                          then encounterPMFRingTwoQ q qB n k rho srcA srcB tmax
                          else V.empty
                    !m  = encounterMeanRingTwoQ q qB n k rho srcA srcB
                    !pk = if wantPMF then mdPeaks (modality p) else []
                    !wMatA = buildTransitionMatrix q  n k []
                    !wMatB = buildTransitionMatrix qB n k []
                    (!eigsA, !evecsA) = exactEigensystem wMatA
                    (!eigsB, !evecsB) = exactEigensystem wMatB
                    !piVec = stationaryDist wMatA
                    (_, !sp, !sm) = encounterDecompositionTwoQ
                        eigsA evecsA eigsB evecsB piVec rho (srcA, srcB)
                in (p, m, 1.0, pk, sp, sm, "closed")
            | runExact && isRing =
                let !wMatA = buildTransitionMatrix q  n k (cfgDefects cfg)
                    !wMatB = buildTransitionMatrix qB n k (cfgDefects cfg)
                    (!eigsA, !evecsA) = exactEigensystem wMatA
                    (!eigsB, !evecsB) = exactEigensystem wMatB
                    !piVec = stationaryDist wMatA
                    !p  = if wantPMF
                          then encounterPMFTwoQAcc acc eigsA evecsA eigsB evecsB piVec rho (srcA, srcB) tmax
                          else V.empty
                    !pk = if wantPMF then mdPeaks (modality p) else []
                    (!m, !sp, !sm) = encounterDecompositionTwoQ
                        eigsA evecsA eigsB evecsB piVec rho (srcA, srcB)
                in (p, m, 1.0, pk, sp, sm, "expanded")
            | runExact =
                let !r = encounterCell (cfgDim cfg) (cfgBC cfg) n k q qB rho
                                       (startA, startB) (cfgDefects cfg) horizon
                    !p = crPMF r
                    !pk = if wantPMF then mdPeaks (modality p) else []
                in ( p, crMean r, crSplittingWeight r, pk
                   , crWeights r, crSiteMeans r, crProvenance r )
            | otherwise = (V.empty, 0, 1.0, [], V.empty, V.empty, "none")

    -- The encounter law resolved by contact site, in time. The renewal system
    -- is solved at every contour node whether or not its solution is summed,
    -- so this costs one transform per site and no additional solves; it is
    -- behind a flag because it is one array per site rather than because it is
    -- expensive to obtain. An absorbing cell is excluded: its killed kernel has
    -- no stationary law for the correction factors to use.
    let !sitePMFs
            | runExact && wantSitePMF && cfgBC cfg /= Absorbing =
                let !wMatA = cellMatrixWith (cfgDim cfg) (cfgBC cfg) n k q  (cfgDefects cfg)
                    !wMatB = cellMatrixWith (cfgDim cfg) (cfgBC cfg) n k qB (cfgDefects cfg)
                    (!eigsA, !evecsA) = exactEigensystem wMatA
                    (!eigsB, !evecsB) = exactEigensystem wMatB
                    !piVec = stationaryDist wMatA
                    !a = cellIndex (cfgDim cfg) n startA
                    !b = cellIndex (cfgDim cfg) n startB
                in encounterSitePMFTwoQAcc acc eigsA evecsA eigsB evecsB
                       piVec rho (a, b) tmax
            | otherwise = []

    if runExact
        then do
            hPutStrLn stderr "Exact path built, encounter statistics ready."
            hPutStrLn stderr $ "Peaks: " ++ show peaks
            hPutStrLn stderr $ "Encounter MFPT (" ++ provenance ++ "): " ++ show mfpt
            hPutStrLn stderr $ "Contact sites: " ++ show (V.length splits)
                ++ "  weights sum " ++ show (V.sum splits)
                ++ "  site means sum " ++ show (V.sum siteMeans)
        else hPutStrLn stderr "Engine=sim: skipping exact eigensystem."

    simResult <- if runSim
        then do
            let sc = defaultSimConfig
                    { simN = n, simK = k, simQ = q, simQB = qB
                    , simSrc = srcA, simSrcB = srcB
                    , simDefects = cfgDefects cfg
                    , simWalkers = nWalkers, simSeed = cfgSeed cfg
                    , simMaxT = tmax, simMode = Simulate.Encounter
                    , simRho = rho
                    }
            let !sr = simulate sc
            hPutStrLn stderr $ "Sim MFPT (tail-corrected): " ++ show (srMeanFPT sr)
                ++ "  conditional: " ++ show (srConditionalMean sr)
                ++ "  censored: " ++ show (srCensoredFraction sr)
            return (Just sr)
        else return Nothing

    let bundle = EncounterBundle
            -- Degree, not range: a site of a one-dimensional lattice at range
            -- k has 2k neighbours and one of a square lattice has 4k. The
            -- reported quantity is the degree in both, so that a figure
            -- comparing the two dimensions compares like with like.
            { ebN = n
            , ebK = (case cfgDim cfg of D1 -> 2; D2 -> 4) * k
            , ebQ = q, ebQB = qB
            , ebSrcA = srcA, ebSrcB = srcB, ebRho = rho
            , ebDefects = show (cfgDefects cfg)
            , ebPMF = pmf, ebMFPT = mfpt
            , ebCensoredFraction =
                if wantPMF && weight > 0.999 then censoredFraction pmf else 0
            , ebTailScale = if wantPMF then tailScale pmf else 0
            -- The censoring correction assumes a proper law: it treats the
            -- mass beyond the horizon as a tail that will arrive. On an
            -- absorbing domain the missing mass was never censored, it was
            -- lost at the boundary, and those pairs never meet at all. The
            -- residual is therefore withheld there rather than reported as a
            -- large number, and the corresponding check is that the tabulated
            -- mass equals the splitting weight.
            , ebSelfConsistency =
                if wantPMF && weight > 0.999
                then snd (selfConsistency mfpt pmf)
                else 0
            , ebMeanProvenance = provenance
            , ebSplittingWeight = weight
            , ebPeaks = peaks
            , ebSplitting = splits
            , ebSiteMeans = siteMeans
            , ebSitePMF = sitePMFs
            , ebSimMFPT = fmap srMeanFPT simResult
            , ebSimCondMean = fmap srConditionalMean simResult
            , ebSimCensored = fmap srCensoredFraction simResult
            , ebSimTailScale = fmap srTailScale simResult
            , ebSimStdErr = fmap srStdErr simResult
            , ebSimHistogram = fmap srHistogram simResult
            }

    let json = renderJSONPretty 0 (exportEncounter bundle)
    case cfgOut cfg of
        Nothing   -> putStrLn json
        Just path -> writeOut path json
