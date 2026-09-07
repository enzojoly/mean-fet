{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : Absorbing
-- Description : Domains whose boundary destroys the walker. A step that would
--               leave the region is taken and the walker is lost, so the
--               transition matrix conserves nothing and the walk terminates
--               with probability one whether or not the two walkers ever meet.
--
--               Two consequences follow, and both change what may be reported.
--               There is no stationary state and hence no pole at the boundary
--               of the unit disc, so the generating function may be evaluated
--               there directly and every delicacy of the recurrent case falls
--               away. But the encounter law is defective: its total mass is the
--               probability that the pair meets before either is lost, and the
--               unconditional mean is infinite whenever that probability falls
--               short of one. The honest pair of numbers is that probability
--               and the mean conditioned upon it, and this module returns them
--               together so that neither can be quoted without the other.

module Absorbing
    ( killedRing1D
    , killedBox2D
    , absorbingSpectrum1D
    , encounterAbsorbing
    ) where

import qualified Numeric.LinearAlgebra as LA

import Types (Matrix, R)

-- | Transition matrix of a lazy k-range walk on the sites 1..N with absorbing
-- sites at 0 and N+1. Weight directed at a site outside the interval is not
-- returned to the walker, as it would be at a reflecting wall, but simply
-- leaves: the column sums fall short of one by exactly the probability of
-- being lost on that step.
killedRing1D :: Int -> Int -> Double -> Matrix
killedRing1D !n !k !q = LA.accum (LA.konst 0 (n, n)) (+) entries
  where
    !w = q / (2 * fromIntegral k)
    entries =
        [ ((s, s), 1 - q) | s <- [0 .. n - 1] ]
        ++ [ ((t, s), w)
           | s <- [0 .. n - 1], m <- [1 .. k], t <- [s + m, s - m]
           , t >= 0, t < n ]

-- | The two-dimensional counterpart on an L by L region with every edge
-- absorbing. A step leaving the region on either axis is lost.
killedBox2D :: Int -> Int -> Double -> Matrix
killedBox2D !l !k !q = LA.accum (LA.konst 0 (n, n)) (+) entries
  where
    !n = l * l
    !w = q / (4 * fromIntegral k)
    idx !x !y = x + l * y
    inside !x !y = x >= 0 && x < l && y >= 0 && y < l
    entries =
        [ ((idx x y, idx x y), 1 - q) | x <- [0 .. l - 1], y <- [0 .. l - 1] ]
        ++ [ ((idx x' y', idx x y), w)
           | x <- [0 .. l - 1], y <- [0 .. l - 1], m <- [1 .. k]
           , (x', y') <- [ (x + m, y), (x - m, y), (x, y + m), (x, y - m) ]
           , inside x' y' ]

-- | Closed-form spectrum of the killed one-dimensional walk, valid only at
-- unit range.
--
-- The sine basis diagonalises the killed operator because the odd extension
-- that annihilates the eigenfunction at the absorbing site also annihilates it
-- one step beyond, and at unit range the walk cannot reach further. At range
-- two and above it can: the image at the second site outside the interval does
-- not vanish, nothing cancels it, and the basis fails. This is not a gap in
-- the derivation but a statement about the operator, so the caller is given
-- the closed form where it exists and directed to diagonalise numerically
-- where it does not.
absorbingSpectrum1D :: Int -> Int -> Double -> Maybe [R]
absorbingSpectrum1D !n !k !q
    | k /= 1    = Nothing
    | otherwise = Just
        [ (1 - q) + q * cos (pi * fromIntegral l / fromIntegral (n + 1))
        | l <- [1 .. n] ]

-- | First-encounter statistics on an absorbing domain, as the pair of numbers
-- the geometry admits: the probability that the walkers meet before either is
-- lost, and the mean meeting time conditioned on their doing so.
--
-- The generating function is evaluated at the boundary of the unit disc
-- directly. Nothing there is singular, since a substochastic kernel has every
-- eigenvalue inside the disc and the renewal matrix stays well conditioned;
-- the derivative is taken by a short extrapolation only because the function
-- is available as values rather than in closed form.
--
-- Returns the splitting weight and the conditional mean. A weight at or below
-- zero means the pair cannot meet at all, and the conditional mean is then
-- undefined rather than large.
encounterAbsorbing :: Matrix -> Matrix -> Double -> (Int, Int)
                   -> (Double, Double)
encounterAbsorbing !wA !wB !rho (!sA, !sB) =
    let (!eA, !vA) = symEigen wA
        (!eB, !vB) = symEigen wB
        !n = LA.rows wA

        atZ :: Double -> Double
        atZ !z = totalAt n eA vA eB vB rho sA sB z

        !phi = atZ 1.0
        !d   = 1.0e-6
        slope !dd = (phi - atZ (1 - dd)) / dd
        !deriv = 2 * slope (d / 2) - slope d
    in if phi <= 0 then (0, 0) else (phi, deriv / phi)

symEigen :: Matrix -> ([R], Matrix)
symEigen w =
    let !sym = LA.trustSym (LA.scale 0.5 (w + LA.tr w))
        (!vals, !vecs) = LA.eigSH sym
    in (LA.toList vals, vecs)

-- | Total encounter generating function at a real argument, assembled from the
-- two single-walker spectra. The kernel couples one mode of each walker and
-- carries the whole of the dependence on the pair; the renewal over the
-- contact sites is the same system solved everywhere else, in real arithmetic
-- because the argument is real.
totalAt :: Int -> [R] -> Matrix -> [R] -> Matrix -> Double
        -> Int -> Int -> Double -> Double
totalAt !n eA vA eB vB !rho !sA !sB !z =
    let !ea = LA.fromList eA
        !eb = LA.fromList eB

        sigma !j !m =
            1 / (1 - z * (ea `LA.atIndex` j) * (eb `LA.atIndex` m))

        vAat i j = vA `LA.atIndex` (i, j)
        vBat i j = vB `LA.atIndex` (i, j)

        cross !p !r = sum
            [ vAat p j * vAat r j * vBat p m * vBat r m * sigma j m
            | j <- [0 .. n - 1], m <- [0 .. n - 1] ]

        initTo !c = sum
            [ vAat sA j * vAat c j * vBat sB m * vBat c m * sigma j m
            | j <- [0 .. n - 1], m <- [0 .. n - 1] ]

        !diag = [ cross c c | c <- [0 .. n - 1] ]
        diagAt c = diag !! c

        !fM = (LA.><) n n
            [ if l == k
              then (1 - rho) / (rho * diagAt k) + 1
              else cross k l / diagAt l
            | l <- [0 .. n - 1], k <- [0 .. n - 1] ]

        !bV = LA.fromList [ initTo c / diagAt c | c <- [0 .. n - 1] ]
        !x  = fM LA.<\> bV
    in LA.sumElements x
