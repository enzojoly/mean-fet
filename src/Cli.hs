-- |
-- Module      : Cli
-- Description : Shared command line plumbing for both executables: flag value
--               parsing with named errors, range validation returning Either,
--               and output path resolution into the json directory.

module Cli
    ( readFlag
    , splitOn
    , parsePair
    , checkRange
    , checkUnit
    , checkMaybe
    , orDie
    , resolveOut
    , writeOut
    ) where

import System.Directory (createDirectoryIfMissing)
import System.Exit (exitFailure)
import System.FilePath (takeDirectory)
import System.IO (hPutStrLn, stderr)
import Text.Read (readMaybe)

readFlag :: Read a => String -> String -> a
readFlag flag v = case readMaybe v of
    Just x  -> x
    Nothing -> error ("Bad value for " ++ flag ++ ": " ++ v)

splitOn :: Char -> String -> [String]
splitOn _ [] = [""]
splitOn c s  = case break (== c) s of
    (w, [])       -> [w]
    (w, _ : rest) -> w : splitOn c rest

parsePair :: String -> (Int, Int)
parsePair s = case break (== ':') s of
    (a, ':' : b) -> (readFlag "coordinate pair" a, readFlag "coordinate pair" b)
    _            -> error ("Expected X:Y, got: " ++ s)

checkRange :: (Ord a, Show a) => String -> a -> a -> a -> Either String a
checkRange name lo hi x
    | x < lo || x > hi =
        Left (name ++ " = " ++ show x ++ " outside [" ++ show lo ++ ", " ++ show hi ++ "]")
    | otherwise = Right x

checkUnit :: String -> Double -> Either String Double
checkUnit name x
    | x <= 0 || x > 1 = Left (name ++ " = " ++ show x ++ " outside (0, 1]")
    | otherwise       = Right x

checkMaybe :: (a -> Either String a) -> Maybe a -> Either String (Maybe a)
checkMaybe _ Nothing  = Right Nothing
checkMaybe f (Just x) = fmap Just (f x)

orDie :: Either String a -> IO a
orDie (Right x) = return x
orDie (Left e)  = do
    hPutStrLn stderr ("Config error: " ++ e)
    exitFailure

resolveOut :: FilePath -> FilePath
resolveOut path
    | '/' `elem` path = path
    | otherwise       = "json/" ++ path

writeOut :: FilePath -> String -> IO ()
writeOut path s = do
    let full = resolveOut path
        dir  = takeDirectory full
    createDirectoryIfMissing True dir
    writeFile full s
    hPutStrLn stderr ("Wrote " ++ full)
