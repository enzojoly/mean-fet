{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : DefectSpec
-- Description : The defect correction judged against the object it claims to
--               compute. A localised modification of the transition matrix is
--               a finite-rank perturbation, so the perturbed resolvent follows
--               from the unperturbed one by a correction whose size is set by
--               the number of modified bonds rather than by the lattice. The
--               claim is checked by inverting the perturbed matrix outright at
--               a size where doing so is cheap.
--
--               The multi-bond path had no assertion of any kind before this,
--               every stored result having carried a single shortcut.

module DefectSpec (tests) where

import Data.Complex (Complex(..), magnitude, realPart)
import qualified Data.Vector as BV
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra ((><))
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Defect (DefectEntry(..), buildHWith)

tests :: TestTree
tests = testGroup "Defect technique"
    [ singleBond
    , twoBonds
    , signDiagnostic
    ]

ringMatrix :: Int -> Int -> Double -> LA.Matrix Double
ringMatrix !n !k !q = LA.accum (LA.konst 0 (n, n)) (+) entries
  where
    !w = q / (2 * fromIntegral k)
    entries =
        [ ((s, s), 1 - q) | s <- [0 .. n - 1] ]
        ++ concat
        [ [ (((s + m) `mod` n, s), w), (((s - m + n * m) `mod` n, s), w) ]
        | s <- [0 .. n - 1], m <- [1 .. k] ]

-- | A bond modification of strength eta on the pair (u, v). Positive strength
-- withdraws weight from the bond and negative strength supplies it, so an
-- added shortcut carries a negative value. The convention is not incidental:
-- reversing it produces a matrix that is not a probability kernel.
bondPerturbation :: Int -> [(Int, Int, Double)] -> LA.Matrix Double
bondPerturbation !n defs = LA.accum (LA.konst 0 (n, n)) (+) entries
  where
    entries = concat
        [ [ ((u, u), eta), ((v, v), eta)
          , ((u, v), negate eta), ((v, u), negate eta) ]
        | (u, v, eta) <- defs ]

resolventAt :: LA.Matrix Double -> Complex Double -> LA.Matrix (Complex Double)
resolventAt w z =
    let !n = LA.rows w
    in LA.inv (LA.ident n - LA.konst z (n, n) * LA.complex w)

within :: Double -> Double -> Double -> Bool
within tol want got
    | abs want < 1e-14 = abs got <= tol
    | otherwise        = abs (got - want) / abs want <= tol

singleBond :: TestTree
singleBond = testGroup "One bond"
    [ testCase "correction reproduces the perturbed resolvent" $ do
        let !z = 0.7 :+ 0
            !eta = negate 0.25
            !w0 = ringMatrix 11 1 0.6
            !wd = w0 + bondPerturbation 11 [(1, 6, eta)]
            !qt = resolventAt w0 z
            !st = resolventAt wd z
            qAt i j = qt `LA.atIndex` (i, j)
            !e = eta :+ 0
            !h = e * (qAt 1 1 + qAt 6 6 - qAt 1 6 - qAt 6 1) - recip z
            !got = realPart (qAt 0 5 - (qAt 1 5 - qAt 6 5) * (e * (qAt 0 1 - qAt 0 6)) / h)
            !want = realPart (st `LA.atIndex` (5, 0))
        assertBool ("want " ++ show want ++ " got " ++ show got)
            (within 1e-10 want got)
    ]

-- | Two bonds at once. The correction is no longer a scalar fraction but a
-- bilinear form in the inverse defect matrix, and the bonds interact through
-- its off-diagonal entries; treating them independently would agree with
-- neither this nor the resolvent.
twoBonds :: TestTree
twoBonds = testGroup "Two bonds"
    [ testCase "bilinear correction reproduces the perturbed resolvent" $ do
        let !z = 0.7 :+ 0
            !bonds = [(1, 6, negate 0.25), (3, 9, negate 0.40)]
            !w0 = ringMatrix 11 1 0.6
            !wd = w0 + bondPerturbation 11 bonds
            !qt = resolventAt w0 z
            !st = resolventAt wd z
            qAt i j = qt `LA.atIndex` (i, j)
            !entries = BV.fromList [ DefectEntry u v e e | (u, v, e) <- bonds ]
            !hMat = buildHWith entries qAt z
            !aVec = LA.fromList
                [ qAt u 5 - qAt v 5 | (u, v, _) <- bonds ]
            !cVec = LA.fromList
                [ (e :+ 0) * qAt 0 u - (e :+ 0) * qAt 0 v | (u, v, e) <- bonds ]
            !x = LA.flatten (LA.luSolve (LA.luPacked hMat) (LA.asColumn cVec))
            !got = realPart (qAt 0 5 - LA.sumElements (aVec * x))
            !want = realPart (st `LA.atIndex` (5, 0))
        assertBool ("want " ++ show want ++ " got " ++ show got)
            (within 1e-10 want got)

    , testCase "buildHWith agrees with the defect matrix written out" $ do
        let !z = 0.7 :+ 0
            !bonds = [(1, 6, negate 0.25), (3, 9, negate 0.40)]
            !w0 = ringMatrix 11 1 0.6
            !qt = resolventAt w0 z
            qAt i j = qt `LA.atIndex` (i, j)
            !entries = BV.fromList [ DefectEntry u v e e | (u, v, e) <- bonds ]
            !viaCode = buildHWith entries qAt z
            !m = length bonds
            !direct = (m >< m)
                [ let (ui, vi, ei) = bonds !! i
                      (uj, vj, _)  = bonds !! j
                      !e = ei :+ 0
                  in e * (qAt uj ui - qAt vj ui)
                     - e * (qAt uj vi - qAt vj vi)
                     - (if i == j then recip z else 0)
                | i <- [0 .. m - 1], j <- [0 .. m - 1] ]
            !dev = LA.maxElement (LA.cmap magnitude (viaCode - direct))
        assertBool ("max deviation " ++ show dev) (dev < 1e-12)
    ]

-- | A bond strength of the wrong sign does not merely give a wrong number: it
-- places a real zero of the defect matrix inside the unit interval, and every
-- inversion downstream then carries a spurious pole. The zero is cheap to look
-- for, and its absence is a standing guarantee about the configuration.
signDiagnostic :: TestTree
signDiagnostic = testGroup "Bond strength sign"
    [ testCase "an added shortcut leaves no zero on the unit interval" $
        assertBool "expected no sign change" (crossings (negate 0.25) == 0)
    , testCase "the reversed sign plants a zero on the unit interval" $
        assertBool "expected a sign change" (crossings 0.25 > 0)
    ]
  where
    crossings eta =
        let !w0 = ringMatrix 9 1 0.6
            vals =
                [ let !z = t :+ 0
                      !qt = resolventAt w0 z
                      qAt i j = qt `LA.atIndex` (i, j)
                  in realPart ((eta :+ 0)
                        * (qAt 1 1 + qAt 6 6 - qAt 1 6 - qAt 6 1) - recip z)
                | t <- [0.05, 0.10 .. 0.95] ]
        in length (filter id (zipWith (\a b -> a * b < 0) vals (drop 1 vals)))
