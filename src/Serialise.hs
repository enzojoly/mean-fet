{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Serialise
-- Description : JSON rendering and result bundles. Output bytes are identical
--               to the previous String renderer; assembly now goes through
--               ByteString Builder so concatenation is linear. The encounter
--               bundle carries the per contact site decomposition alongside the
--               summed quantities, the site distributions being optional
--               because they are one array per lattice site.

module Serialise
    ( JValue(..)
    , renderJSON
    , renderJSONPretty
    , PassageBundle(..)
    , EncounterBundle(..)
    , exportPassage
    , exportEncounter
    ) where

import Data.List (intersperse)
import Numeric (showFFloat)
import qualified Data.ByteString.Builder as B
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Vector.Unboxed as V

data JValue
    = JNum  !Double
    | JInt  !Int
    | JStr  !String
    | JBool !Bool
    | JNull
    | JArr  [JValue]
    | JObj  [(String, JValue)]

renderJSON :: JValue -> String
renderJSON = BL.unpack . B.toLazyByteString . buildJSON

renderJSONPretty :: Int -> JValue -> String
renderJSONPretty ind = BL.unpack . B.toLazyByteString . buildPretty ind

numText :: Double -> String
numText d
    | isNaN d      = "null"
    | isInfinite d = "null"
    | d == fromIntegral (round d :: Integer) && abs d < 1e15
                   = show (round d :: Integer)
    | otherwise    = showFFloat Nothing d ""

buildJSON :: JValue -> B.Builder
buildJSON (JNum d)   = B.stringUtf8 (numText d)
buildJSON (JInt i)   = B.stringUtf8 (show i)
buildJSON (JStr s)   = B.charUtf8 '"' <> B.stringUtf8 (escapeJSON s) <> B.charUtf8 '"'
buildJSON (JBool b)  = B.stringUtf8 (if b then "true" else "false")
buildJSON JNull      = B.stringUtf8 "null"
buildJSON (JArr xs)  =
    B.charUtf8 '[' <> mconcat (intersperse (B.charUtf8 ',') (map buildJSON xs))
                   <> B.charUtf8 ']'
buildJSON (JObj kvs) =
    B.charUtf8 '{' <> mconcat (intersperse (B.charUtf8 ',') (map kv kvs))
                   <> B.charUtf8 '}'
  where
    kv (k, v) = B.charUtf8 '"' <> B.stringUtf8 (escapeJSON k)
             <> B.stringUtf8 "\":" <> buildJSON v

buildPretty :: Int -> JValue -> B.Builder
buildPretty = go
  where
    go _ (JNum d)  = B.stringUtf8 (numText d)
    go _ (JInt i)  = B.stringUtf8 (show i)
    go _ (JStr s)  = buildJSON (JStr s)
    go _ (JBool b) = buildJSON (JBool b)
    go _ JNull     = B.stringUtf8 "null"
    go _ (JArr []) = B.stringUtf8 "[]"
    go indent (JArr xs) =
        B.stringUtf8 "[\n"
        <> mconcat (intersperse (B.stringUtf8 ",\n")
                    (map (\x -> pad (indent + 2) <> go (indent + 2) x) xs))
        <> B.charUtf8 '\n' <> pad indent <> B.charUtf8 ']'
    go _ (JObj []) = B.stringUtf8 "{}"
    go indent (JObj kvs) =
        B.stringUtf8 "{\n"
        <> mconcat (intersperse (B.stringUtf8 ",\n") (map (kvp (indent + 2)) kvs))
        <> B.charUtf8 '\n' <> pad indent <> B.charUtf8 '}'

    kvp ind (k, v) =
        pad ind <> B.charUtf8 '"' <> B.stringUtf8 (escapeJSON k)
        <> B.stringUtf8 "\": " <> go ind v

    pad n = B.stringUtf8 (replicate n ' ')

escapeJSON :: String -> String
escapeJSON = concatMap esc
  where
    esc '"'  = "\\\""
    esc '\\' = "\\\\"
    esc '\n' = "\\n"
    esc '\r' = "\\r"
    esc '\t' = "\\t"
    esc c
        | c < ' '   = "\\u" ++ padHex (showHex' (fromEnum c))
        | otherwise  = [c]
    padHex s = replicate (4 - length s) '0' ++ s
    showHex' n
        | n < 16    = [hexDigit n]
        | otherwise = showHex' (n `div` 16) ++ [hexDigit (n `mod` 16)]
    hexDigit n
        | n < 10    = toEnum (n + fromEnum '0')
        | otherwise = toEnum (n - 10 + fromEnum 'a')

data PassageBundle = PassageBundle
    { pbN             :: !Int
    , pbK             :: !Int
    , pbQ             :: !Double
    , pbSrc           :: !Int
    , pbTgt           :: !Int
    , pbRho           :: !Double
    , pbDefects       :: !String
    , pbPMF           :: !(V.Vector Double)
    , pbPureRingPMF   :: !(V.Vector Double)
    , pbMFPTNetwork   :: !Double
    , pbMFPTRingExact :: !Double
    , pbMFPTTruncated :: !Double
    , pbCensoredFraction :: !Double
    , pbTailScale     :: !Double
    , pbSelfConsistency :: !Double
    , pbMeanProvenance :: !String
    , pbSplittingWeight :: !Double
    , pbPeaks         :: ![(Int, Double)]
    , pbSimMFPT       :: !(Maybe Double)
    , pbSimCondMean   :: !(Maybe Double)
    , pbSimCensored   :: !(Maybe Double)
    , pbSimTailScale  :: !(Maybe Double)
    , pbSimStdErr     :: !(Maybe Double)
    , pbSimHistogram  :: !(Maybe (V.Vector Double))
    }

-- | The encounter result. Beside the summed quantities the bundle carries the
-- decomposition by contact site: the splitting probability at each site, that
-- site's contribution to the mean, and optionally the site's own distribution
-- in time.
--
-- The contribution is the weight multiplied by the conditional mean, not the
-- conditional mean itself, because the quotient is undefined where the weight
-- vanishes and the product is what sums to the mean. The site distributions
-- are one array per lattice site and are written only when asked for.
data EncounterBundle = EncounterBundle
    { ebN             :: !Int
    , ebK             :: !Int
    , ebQ             :: !Double
    , ebQB            :: !Double
    , ebSrcA          :: !Int
    , ebSrcB          :: !Int
    , ebRho           :: !Double
    , ebDefects       :: !String
    , ebPMF           :: !(V.Vector Double)
    , ebMFPT          :: !Double
    , ebCensoredFraction :: !Double
    , ebTailScale     :: !Double
    , ebSelfConsistency :: !Double
    , ebMeanProvenance :: !String
    , ebSplittingWeight :: !Double
    , ebPeaks         :: ![(Int, Double)]
    , ebSplitting     :: !(V.Vector Double)
    , ebSiteMeans     :: !(V.Vector Double)
    , ebSitePMF       :: ![V.Vector Double]
    , ebSimMFPT       :: !(Maybe Double)
    , ebSimCondMean   :: !(Maybe Double)
    , ebSimCensored   :: !(Maybe Double)
    , ebSimTailScale  :: !(Maybe Double)
    , ebSimStdErr     :: !(Maybe Double)
    , ebSimHistogram  :: !(Maybe (V.Vector Double))
    }

exportPassage :: PassageBundle -> JValue
exportPassage pb = JObj
    [ ("mode", JStr "passage")
    , ("config", JObj
        [ ("N",   JInt (pbN pb))
        , ("K",   JInt (pbK pb))
        , ("q",   JNum (pbQ pb))
        , ("src", JInt (pbSrc pb))
        , ("tgt", JInt (pbTgt pb))
        , ("rho", JNum (pbRho pb))
        , ("defects", JStr (pbDefects pb))
        ])
    , ("pmf", JArr [JNum p | p <- V.toList (pbPMF pb)])
    , ("pure_ring_pmf", JArr [JNum p | p <- V.toList (pbPureRingPMF pb)])
    , ("mfpt", JObj
        [ ("network", JNum (pbMFPTNetwork pb))
        , ("ring_exact", JNum (pbMFPTRingExact pb))
        , ("network_truncated_sum", JNum (pbMFPTTruncated pb))
        ])
    , ("splitting_weight", JNum (pbSplittingWeight pb))
    , ("diagnostics", JObj
        [ ("censored_fraction", JNum (pbCensoredFraction pb))
        , ("tail_scale", JNum (pbTailScale pb))
        , ("self_consistency_residual", JNum (pbSelfConsistency pb))
        , ("mean_provenance", JStr (pbMeanProvenance pb))
        ])
    , ("peaks", JArr [JObj [("t", JInt t), ("height", JNum h)] | (t, h) <- pbPeaks pb])
    , ("simulation", JObj
        [ ("mfpt", maybe JNull JNum (pbSimMFPT pb))
        , ("conditional_mean", maybe JNull JNum (pbSimCondMean pb))
        , ("censored_fraction", maybe JNull JNum (pbSimCensored pb))
        , ("tail_scale", maybe JNull JNum (pbSimTailScale pb))
        , ("stderr", maybe JNull JNum (pbSimStdErr pb))
        , ("histogram", case pbSimHistogram pb of
            Nothing -> JNull
            Just h  -> JArr [JNum v | v <- V.toList h])
        ])
    ]

exportEncounter :: EncounterBundle -> JValue
exportEncounter eb = JObj
    [ ("mode", JStr "encounter")
    , ("config", JObj
        [ ("N",   JInt (ebN eb))
        , ("K",   JInt (ebK eb))
        , ("q",   JNum (ebQ eb))
        , ("qB",  JNum (ebQB eb))
        , ("srcA", JInt (ebSrcA eb))
        , ("srcB", JInt (ebSrcB eb))
        , ("rho", JNum (ebRho eb))
        , ("defects", JStr (ebDefects eb))
        ])
    , ("encounter_pmf", JArr [JNum p | p <- V.toList (ebPMF eb)])
    , ("encounter_mfpt", JNum (ebMFPT eb))
    , ("splitting_weight", JNum (ebSplittingWeight eb))
    , ("diagnostics", JObj
        [ ("censored_fraction", JNum (ebCensoredFraction eb))
        , ("tail_scale", JNum (ebTailScale eb))
        , ("self_consistency_residual", JNum (ebSelfConsistency eb))
        , ("mean_provenance", JStr (ebMeanProvenance eb))
        ])
    , ("peaks", JArr [JObj [("t", JInt t), ("height", JNum h)] | (t, h) <- ebPeaks eb])
    , ("splitting_probabilities", JArr [JNum p | p <- V.toList (ebSplitting eb)])
    , ("site_mean_contributions", JArr [JNum p | p <- V.toList (ebSiteMeans eb)])
    , ("site_pmf", JArr [ JArr [JNum v | v <- V.toList s] | s <- ebSitePMF eb ])
    , ("simulation", JObj
        [ ("mfpt", maybe JNull JNum (ebSimMFPT eb))
        , ("conditional_mean", maybe JNull JNum (ebSimCondMean eb))
        , ("censored_fraction", maybe JNull JNum (ebSimCensored eb))
        , ("tail_scale", maybe JNull JNum (ebSimTailScale eb))
        , ("stderr", maybe JNull JNum (ebSimStdErr eb))
        , ("histogram", case ebSimHistogram eb of
            Nothing -> JNull
            Just h  -> JArr [JNum v | v <- V.toList h])
        ])
    ]
