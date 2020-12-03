module Data.CG.Minus.Bezier where

import Data.CG.Minus.Types

-- | Three-point bezier curve interpolation.  The index /mu/ is
--   in the range zero to one.
bezier3 :: Num n => Pt n -> Pt n -> Pt n -> n -> Pt n
bezier3 (Pt x1 y1) (Pt x2 y2) (Pt x3 y3) mu =
    let a = mu*mu
        b = 1 - mu
        c = b*b
        x = x1*c + 2*x2*b*mu + x3*a
        y = y1*c + 2*y2*b*mu + y3*a
    in Pt x y

-- | Four-point bezier curve interpolation.  The index /mu/ is
--   in the range zero to one.
bezier4 :: Num n => Pt n -> Pt n -> Pt n -> Pt n -> n -> Pt n
bezier4 (Pt x1 y1) (Pt x2 y2) (Pt x3 y3) (Pt x4 y4) mu =
    let a = 1 - mu
        b = a*a*a
        c = mu*mu*mu
        x = b*x1 + 3*mu*a*a*x2 + 3*mu*mu*a*x3 + c*x4
        y = b*y1 + 3*mu*a*a*y2 + 3*mu*mu*a*y3 + c*y4
    in Pt x y

-- | Generic variant, given a scaling function, ie. 'pt_scale' or 'p3_scale'.
bezier4' :: (Num a, Num t) => (a -> t -> t) -> t -> t -> t -> t -> a -> t
bezier4' f p0 p1 p2 p3 t =
    let u = 1 - t
        tt = t * t
        ttt = tt * t
        uu = u * u
        uuu = uu * u
    in (uuu `f` p0) + ((3 * uu * t) `f` p1) + ((3 * u * tt) `f` p2) + (ttt `f` p3)

unwrap_mu :: Fractional a => a -> a -> a -> a
unwrap_mu x0 x3 mu = let d = x3 - x0 in (mu - x0) / d

next_mu :: Fractional a => a -> a -> a
next_mu x0 x3 = x0 + ((x3 - x0) / 2)

-- | Requires that x is monotonic (increasing), so that a simple
-- binary search will work.  Initial /mu/ allows fast forwarding to
-- last result for left to right lookup (unused in variants here).
bezier4_y_mt' :: (Ord a, Fractional a) => a -> Pt a -> Pt a -> Pt a -> Pt a -> a -> a -> (a,a)
bezier4_y_mt' dx (Pt x0 y0) (Pt x1 y1) (Pt x2 y2) (Pt x3 y3) i_mu x =
    let f = bezier4 (Pt x0 y0) (Pt x1 y1) (Pt x2 y2) (Pt x3 y3)
        recur l r mu =
            let (Pt x' y') = f (unwrap_mu x0 x3 mu)
            in if abs (x - x') <= dx
               then (mu,y')
               else let (l',r') = if x' < x then (mu,r) else (l,mu)
                        mu' = next_mu l' r'
                    in recur l' r' mu'
    in recur x0 x3 i_mu

-- | Variant with initial /mu/ at mid-point.
--
-- > import Sound.SC3.Plot
--
-- > let f = bezier4_y_mt 0.0001 (Pt 0 0) (Pt 0.1 0.3) (Pt 0.5 0.2) (Pt 1 0.5)
-- > plotTable1 (map (snd . f) [0.0,0.01 .. 1.0])
--
-- > let f = bezier4_y_mt 0.0001 (Pt 0 0) (Pt 0.4 1.3) (Pt 0.6 1.3) (Pt 1 0)
-- > plotTable1 (map (snd . f) [0.0,0.01 .. 1.0])
--
-- > let f = bezier4_y_mt 0.0001 (Pt 1 0) (Pt 1.4 1.3) (Pt 1.6 1.3) (Pt 2 0)
-- > plotTable1 (map (snd . f) [1.0,1.01 .. 2.0])
bezier4_y_mt :: (Ord a, Fractional a) => a -> Pt a -> Pt a -> Pt a -> Pt a -> a -> (a, a)
bezier4_y_mt dx p0 c1 c2 p3 x =
    let x0 = pt_x p0
        x3 = pt_x p3
        mu = x0 + ((x3 - x0) / 2)
    in bezier4_y_mt' dx p0 c1 c2 p3 mu x

-- | Variant that scans a list [p0,c1,c2,p3,c4,c5,p6...] to resolve /x/.
bezier4_seq_y' :: (Ord a, Fractional a) => a -> [Pt a] -> a -> Maybe (a,a)
bezier4_seq_y' dx pt_l x =
    case pt_l of
      p0:c1:c2:p3:pt_l' -> if x >= pt_x p0 && x <= pt_x p3
                           then Just (bezier4_y_mt dx p0 c1 c2 p3 x)
                           else bezier4_seq_y' dx (p3 : pt_l') x
      _ -> Nothing

-- | Variant giving only /y/ and zero for out of range /x/.
bezier4_seq_y0 :: (Ord a, Fractional a) => a -> [Pt a] -> a -> a
bezier4_seq_y0 dx pt = maybe 0 snd . bezier4_seq_y' dx pt

-- | Variant to generate a /wavetable/, ie. /n/ equally spaced /x/ across range of /pt/.
--
-- > let pt = [Pt 0 0,Pt 0.4 (-1.3),Pt 0.6 (-1.3),Pt 1 0,Pt 1.4 1.3,Pt 1.6 1.3,Pt 2 0]
-- > plotTable1 (bezier4_seq_wt 1e-4 pt 1024)
bezier4_seq_wt :: (Ord a, Integral i, Fractional a, Enum a) => a -> [Pt a] -> i -> [a]
bezier4_seq_wt dx pt n =
    let l = pt_x (head pt)
        r = pt_x (last pt)
        l' = l + ((r - l) / fromIntegral n)
        ix = [l,l' .. r]
    in map (bezier4_seq_y0 dx pt) ix
