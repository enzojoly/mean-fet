{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : SlopeSpec
-- Description : The slope of the mean first-encounter time at zero mobility of
--               walker B, checked two ways: the derivative of the pair-space
--               resolvent taken directly, and the closed form of Appendix E
--               built from walker A's hitting times and its occupation before
--               reaching walker B's release. The pair space is solved densely,
--               so the lattices are small; the closed form is size-independent.

module SlopeSpec (tests) where

import qualified Numeric.LinearAlgebra as LA
import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit (testCase, assertBool)

import Cells (Dim(..), BC(..), cellMatrix, cellStates)

tests :: TestTree
tests = testGroup "Slope of the mean at zero mobility of walker B"
    [ testCase "interval, mid-domain release, K = 4" $
        check D1 Reflecting 24 2 0.75 1 11
    , testCase "ring, mid-domain release, K = 4" $
        check D1 Periodic 24 2 0.75 1 11
    , testCase "ring, adjacent release, K = 2" $
        check D1 Periodic 24 1 0.75 1 2
    , testCase "interval, near-wall release, K = 4" $
        check D1 Reflecting 24 2 0.75 1 2
    , testCase "interval, far-wall release, K = 4" $
        check D1 Reflecting 24 2 0.75 1 22
    , testCase "2D box, interior release, K = 4" $
        check D2 Reflecting 6 1 0.75 (0 + 6 * 1) (3 + 6 * 3)
    ]

check :: Dim -> BC -> Int -> Int -> Double -> Int -> Int -> IO ()
check dim bc size k q1 a b =
    let !direct = slopeResolvent dim bc size k q1 a b
        !closed = slopeClosedForm dim bc size k q1 a b
        !scale  = max 1 (abs direct)
    in assertBool (report direct closed) (abs (direct - closed) / scale < 1e-9)

report :: Double -> Double -> String
report direct closed =
    "resolvent " ++ show direct ++ ", closed form " ++ show closed

slopeResolvent :: Dim -> BC -> Int -> Int -> Double -> Int -> Int -> Double
slopeResolvent dim bc size k q1 a b =
    let !ns    = cellStates dim size
        !wA    = cellMatrix dim bc size k q1
        !jB    = cellMatrix dim bc size k 1.0
        !eye   = LA.ident ns
        !off   = [ x * ns + y | x <- [0 .. ns - 1], y <- [0 .. ns - 1], x /= y ]
        !sel   = LA.Pos (LA.idxs off)
        !q0    = LA.kronecker wA eye LA.?? (sel, sel)
        !dq    = LA.kronecker wA (jB - eye) LA.?? (sel, sel)
        !sys   = LA.ident (length off) - q0
        !t0    = sys LA.<\> LA.konst 1 (length off)
        !dt    = sys LA.<\> (dq LA.#> t0)
        !place = length (takeWhile (/= a * ns + b) off)
    in dt LA.! place

slopeClosedForm :: Dim -> BC -> Int -> Int -> Double -> Int -> Int -> Double
slopeClosedForm dim bc size k q1 a b =
    let !ns   = cellStates dim size
        !wA   = cellMatrix dim bc size k q1
        !jB   = cellMatrix dim bc size k 1.0
        !hb   = hitting wA ns b
        !occ  = LA.accum (occupation wA ns b a) (+) [(a, -1)]
        hop j = (jB LA.! b) LA.! j
        term j =
            let !hj    = hitting wA ns j
                !inner = sum [ (occ LA.! x) * ((hj LA.! x) - (hb LA.! x))
                             | x <- [0 .. ns - 1], x /= b ]
            in hop j * (inner + hj LA.! b)
    in sum [ term j | j <- [0 .. ns - 1], j /= b, hop j > 0 ]

hitting :: LA.Matrix Double -> Int -> Int -> LA.Vector Double
hitting w ns j =
    let !keep = [ x | x <- [0 .. ns - 1], x /= j ]
        !sel  = LA.Pos (LA.idxs keep)
        !sub  = w LA.?? (sel, sel)
        !sol  = (LA.ident (ns - 1) - sub) LA.<\> LA.konst 1 (ns - 1)
    in LA.assoc ns 0 (zip keep (LA.toList sol))

occupation :: LA.Matrix Double -> Int -> Int -> Int -> LA.Vector Double
occupation w ns b a =
    let !keep = [ x | x <- [0 .. ns - 1], x /= b ]
        !sel  = LA.Pos (LA.idxs keep)
        !sub  = w LA.?? (sel, sel)
        !fund = LA.inv (LA.ident (ns - 1) - sub)
        !row  = length (takeWhile (/= a) keep)
    in LA.assoc ns 0 (zip keep (LA.toList (fund LA.! row)))
