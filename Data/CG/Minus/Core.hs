module Data.CG.Minus.Core where

import qualified Data.Fixed as F {- base -}
import Data.List {- base -}
import Data.Maybe {- base -}
import System.Random {- random -}
import Text.Printf {- base -}

import Data.CG.Minus.Types
import Data.CG.Minus.Geometry

pt_eq_by :: (t -> t -> Bool) -> Pt t -> Pt t -> Bool
pt_eq_by f (Pt x1 y1) (Pt x2 y2) = f x1 x2 && f y1 y2

pt_eq_approx :: (Ord n,Floating n) => Pt n -> Pt n -> Bool
pt_eq_approx = pt_eq_by (~=)

-- * R(eal) functions

-- | Epsilon.
epsilon :: Floating n => n
epsilon = 0.000001

-- | Is absolute difference less than 'epsilon'.
(~=) :: (Floating a, Ord a) => a -> a -> Bool
p ~= q = abs (p - q) < epsilon

-- | Degrees to radians.
--
-- > map r_to_radians [-180,-90,0,90,180] == [-pi,-pi/2,0,pi/2,pi]
r_to_radians :: R -> R
r_to_radians x = (x / 180) * pi

-- | Radians to degrees, inverse of 'r_to_radians'.
--
-- > map r_from_radians [-pi,-pi/2,0,pi/2,pi] == [-180,-90,0,90,180]
-- > r_from_radians (pi/8) == 22.5
r_from_radians :: R -> R
r_from_radians x = (x / pi) * 180

-- | 'R' modulo within range.
--
-- > map (r_constrain (3,5)) [2.75,5.25] == [4.75,3.25]
r_constrain :: (R,R) -> R -> R
r_constrain (l,r) =
    let down n i x = if x > i then down n i (x - n) else x
        up n i x = if x < i then up n i (x + n) else x
        both n i j x = up n i (down n j x)
    in both (r - l) l r

-- | Sum of squares.
mag_sq :: Num a => a -> a -> a
mag_sq x y = x * x + y * y

-- | 'sqrt' of 'mag_sq'.
mag :: Floating c => c -> c -> c
mag x = sqrt . mag_sq x

-- * Pt functions

-- | Tuple constructor.
mk_pt :: (a,a) -> Pt a
mk_pt (x,y) = Pt x y

-- | Tuple accessor.
pt_xy :: Pt t -> (t,t)
pt_xy (Pt x y) = (x,y)

-- | 'Pt' of (0,0).
--
-- > pt_origin == Pt 0 0
pt_origin :: Num a => Pt a
pt_origin = Pt 0 0

-- | 'Pt' at /(n,n)/.
--
-- > pt_from_scalar 1 == Pt 1 1
pt_from_scalar :: a -> Pt a
pt_from_scalar a = Pt a a

-- | Clip /x/ and /y/ to lie in /(0,n)/.
--
-- > pt_clipu 1 (Pt 0.5 1.5) == Pt 0.5 1
pt_clipu :: (Ord a,Num a) => a -> Pt a -> Pt a
pt_clipu u =
    let f n = if n < 0 then 0 else if n > u then u else n
    in pt_uop f

-- | Swap /x/ and /y/ coordinates at 'Pt'.
--
-- > pt_swap (Pt 1 2) == Pt 2 1
pt_swap :: Pt a -> Pt a
pt_swap (Pt x y) = Pt y x

-- | Negate /y/ element of 'Pt'.
--
-- > pt_negate_y (Pt 1 1) == Pt 1 (-1)
pt_negate_y :: (Num a) => Pt a -> Pt a
pt_negate_y (Pt x y) = Pt x (negate y)

-- | 'Pt' variant of 'r_to_radians'.
--
-- > pt_to_radians (Pt 90 270) == Pt (pi/2) (pi*(3/2))
pt_to_radians :: Pt R -> Pt R
pt_to_radians = pt_uop r_to_radians

-- | 'Pt' variant of 'r_from_radians'.
pt_from_radians :: Pt R -> Pt R
pt_from_radians = pt_uop r_from_radians

-- | Cartesian (x,y) to polar (distance,angle) coordinate form.
--
-- > map pt_to_polar [Pt 1 0,Pt 0 pi] == [Pt 1 0,Pt pi (pi/2)]
pt_to_polar :: Pt R -> Pt R
pt_to_polar = mk_pt . rectangular_to_polar . pt_xy

{- | Polar to cartesian, inverse of 'pt_to_polar'.

<http://mathworld.wolfram.com/PolarCoordinates.html>

rho (magnitude) is the radial distance from the origin,
theta (phase )is the counterclockwise angle from the x-axis

> pt_eq_by (~=) (pt_from_polar (Pt pi (pi/2))) (Pt 0 pi)
> pt_eq_by (~=) (pt_from_polar (Pt (sqrt 2) (r_to_radians 45))) (Pt 1 1)

-}
pt_from_polar :: Pt R -> Pt R
pt_from_polar = mk_pt . polar_to_rectangular . pt_xy

mk_pt_polar :: (R,R) -> Pt R
mk_pt_polar = mk_pt . polar_to_rectangular

pt_polar :: R -> R -> Pt R
pt_polar = curry mk_pt_polar

-- | Scalar 'Pt' '+'.
--
-- > pt_offset 1 pt_origin == Pt 1 1
pt_offset :: Num a => a -> Pt a -> Pt a
pt_offset = pt_uop . (+)

-- | Scalar 'Pt' '*'.
--
-- > pt_scale 2 (Pt 1 2) == Pt 2 4
pt_scale :: Num a => a -> Pt a -> Pt a
pt_scale = pt_uop . (*)

{-
infixl 7 |*, *|

(|*) :: Num a => a -> Pt a -> Pt a
(|*) = pt_scale

(*|) :: Num a => Pt a -> a -> Pt a
(*|) = flip pt_scale
-}

-- | Pointwise 'min'.
pt_min :: (Ord a) => Pt a -> Pt a -> Pt a
pt_min = pt_binop min

-- | Pointwise 'max'.
pt_max :: (Ord a) => Pt a -> Pt a -> Pt a
pt_max = pt_binop max

-- | Apply function to /x/ and /y/ fields of three 'Pt'.
pt_ternary_f :: (a->a->b->b->c->c->d) -> Pt a -> Pt b -> Pt c -> d
pt_ternary_f f (Pt x0 y0) (Pt x1 y1) (Pt x2 y2) = f x0 y0 x1 y1 x2 y2

-- | Given a /(minima,maxima)/ pair, expand so as to include /p/.
--
-- > pt_minmax (Pt 0 0,Pt 1 1) (Pt (-1) 2) == (Pt (-1) 0,Pt 1 2)
pt_minmax :: Ord a => (Pt a,Pt a) -> Pt a -> (Pt a,Pt a)
pt_minmax (p0,p1) p =
    let f x0 y0 x1 y1 x y =
            (Pt (min x x0) (min y y0)
            ,Pt (max x x1) (max y y1))
    in pt_ternary_f f p0 p1 p

-- | 'Pt' variant of 'constrain'.
pt_constrain :: (Pt R,Pt R) -> Pt R -> Pt R
pt_constrain (p0,p1) p =
    let f x0 y0 x1 y1 x y =
            let x' = r_constrain (x0,x1) x
                y' = r_constrain (y0,y1) y
            in Pt x' y'
    in pt_ternary_f f p0 p1 p

-- | Angle to origin.  By convention the angle is zero rightwards
-- along the X axis, and is measured counter-clockwise.
--
-- > map pt_angle_o [Pt 0 1,Pt 1 0,Pt 0 (-1),Pt (-1) 0] == [pi / 2,0,-pi / 2,pi]
-- > r_from_radians (pt_angle_o (Pt 2 1)) ~= 26.565051177077986
pt_angle_o :: Pt R -> R
pt_angle_o (Pt x y) = atan2 y x

-- | Angle from /p/ to /q/.
--
-- > pt_angle (Pt 0 (-1)) (Pt 0 1) == pi/2
-- > pt_angle (Pt 1 0) (Pt 0 1) == pi * 3/4
-- > pt_angle (Pt 0 1) (Pt 0 1) == 0
pt_angle :: Pt R -> Pt R -> R
pt_angle p q = pt_angle_o (q - p)

-- | Pointwise '+'.
--
-- pt_translate (Vc 1 1) (Pt 0 0) == pt 1 1
pt_translate :: Num a => Vc a -> Pt a -> Pt a
pt_translate (Vc dx dy) (Pt x y) = Pt (x + dx) (y + dy)

-- | 'pt_uop' 'fromIntegral'.
pt_from_i :: (Integral i,Num a) => Pt i -> Pt a
pt_from_i = pt_uop fromIntegral

-- | 'mag_sq' of /x/ /y/ of 'Pt'.
pt_mag_sq :: Num a => Pt a -> a
pt_mag_sq (Pt x y) = mag_sq x y

-- | 'mag' of /x/ /y/ of 'Pt'.
pt_mag :: Floating a => Pt a -> a
pt_mag (Pt x y) = mag x y

-- | Distance from 'Pt' /p/ to 'Pt' /q/.
--
-- > pt_distance (Pt 0 0) (Pt 0 1) == 1
-- > pt_distance (Pt 0 0) (Pt 1 1) == sqrt 2
pt_distance :: Floating a => Pt a -> Pt a -> a
pt_distance p1 p2 = pt_mag (p2 - p1)

-- | Are /x/ and /y/ of 'Pt' /p/ in range (0,1).
--
-- > map pt_is_normal [Pt 0 0,Pt 1 1,Pt 2 2] == [True,True,False]
pt_is_normal :: (Ord a,Num a) => Pt a -> Bool
pt_is_normal (Pt x y) = x >= 0 && x <= 1 && y >= 0 && y <= 1

-- | Rotate 'Pt' /n/ radians.
--
-- > pt_rotate pi (Pt 1 0) ~= Pt (-1) 0
pt_rotate :: Floating a => a -> Pt a -> Pt a
pt_rotate a (Pt x y) =
    let s = sin a
        c = cos a
    in Pt (x * c - y * s) (y * c + x * s)

pt_rotate_about :: Floating a => a -> Pt a -> Pt a -> Pt a
pt_rotate_about r p0 p1 = pt_rotate r (p1 - p0) + p0

pt_lift3 :: ((a,a) -> (b,b) -> (c,c) -> (d,d)) -> Pt a -> Pt b -> Pt c -> Pt d
pt_lift3 f (Pt x1 y1) (Pt x2 y2) (Pt x3 y3) =
    let (x,y) = f (x1,y1) (x2,y2) (x3,y3)
    in Pt x y

-- | Copy /x/ and /y/ from 'Pt' to 'Vc'.
pt_to_vc :: Pt a -> Vc a
pt_to_vc (Pt x y) = Vc x y

-- | The reflection of 'Pt' across a vertical line at /x/.
pt_reflect_x :: Num a => a -> Pt a -> Pt a
pt_reflect_x rx (Pt x y) = Pt (rx + (rx - x)) y

-- | The reflection of 'Pt' across a horizontal line at /y/.
pt_reflect_y :: Num a => a -> Pt a -> Pt a
pt_reflect_y ry (Pt x y) = Pt x (ry + (ry - y))

-- | The reflection of 'Pt' across the minimum distance to (/x/,/y/).
--
-- > pt_reflect_xy (1,1) (Pt (-1) 0) == Pt 3 2
pt_reflect_xy :: Num a => (a,a) -> Pt a -> Pt a
pt_reflect_xy (rx,ry) (Pt x y) = Pt (rx + (rx - x)) (ry + (ry - y))

-- * Vc functions

mk_vc :: (R,R) -> Vc R
mk_vc = uncurry Vc

mk_vc_polar :: (R,R) -> Vc R
mk_vc_polar = mk_vc . polar_to_rectangular

vc_polar :: R -> R -> Vc R
vc_polar = curry mk_vc_polar

vc_to_polar :: Vc R -> Vc R
vc_to_polar = mk_vc . rectangular_to_polar . vc_xy

vc_xy :: Vc t -> (t,t)
vc_xy (Vc x y) = (x,y)

-- | 'mag_sq' of 'Vc'.
vc_mag_sq :: Floating c => Vc c -> c
vc_mag_sq (Vc dx dy) = mag_sq dx dy

-- | 'mag' of 'Vc'.
vc_mag :: Floating c => Vc c -> c
vc_mag (Vc dx dy) = mag dx dy

-- | Multiply 'Vc' pointwise by scalar.
--
-- > vc_scale 2 (Vc 3 4) == Vc 6 8
vc_scale :: Num a => a -> Vc a -> Vc a
vc_scale n (Vc x y) = Vc (x * n) (y * n)

-- | 'Vc' dot product.
--
-- > vc_dot (Vc 1 2) (Vc 3 4) == 11
vc_dot :: Num a => Vc a -> Vc a -> a
vc_dot (Vc x y) (Vc x' y') = (x * x') + (y * y')

-- | Scale 'Vc' to have unit magnitude (to within tolerance).
--
-- > let x = (sqrt 2) / 2
-- > vc_mag_sq (Vc x x) ~= 1.0
-- > vc_unit (Vc x x) == Vc x x
-- > vc_unit (Vc 0 0) == Vc 0 0
-- > vc_unit (Vc 1 1) ~= Vc x x
vc_unit :: (Ord a, Floating a) => Vc a -> Vc a
vc_unit v =
    if abs (vc_mag_sq v - 1) < epsilon
    then v
    else if vc_mag_sq v == 0
         then v
         else let m = vc_mag v in Vc (vc_x v / m) (vc_y v / m)

-- | The angle between two vectors on a plane. The angle is from v1 to
-- v2, positive anticlockwise.  The result is in (-pi,pi)
vc_angle :: Vc R -> Vc R -> R
vc_angle (Vc x1 y1) (Vc x2 y2) =
    let t1 = atan2 y1 x1
        t2 = atan2 y2 x2
    in r_constrain (-pi,pi) (t2 - t1)

-- * Line functions

-- | Variant on 'Ln' which takes 'Pt' co-ordinates as duples.
--
-- > mk_ln ((0,0),(1,1)) == Ln (Pt 0 0) (Pt 1 1)
mk_ln :: ((a,a),(a,a)) -> Ln a
mk_ln ((x1,y1),(x2,y2)) = Ln (Pt x1 y1) (Pt x2 y2)

{-
ln :: (Pt a, Pt a) -> Ln a
ln = uncurry Ln
-}

ln' :: (a, a) -> (a, a) -> Ln a
ln' = curry mk_ln

-- | 'Vc' that 'pt_translate's start 'Pt' to end 'Pt' of 'Ln'.
--
-- > let l = Ln (Pt 0 0) (Pt 1 1)
-- > in ln_start l `pt_translate` ln_vc l == Pt 1 1
ln_vc :: Num a => Ln a -> Vc a
ln_vc (Ln p q) = let Pt x y = q - p in Vc x y

-- | 'ln_start' and 'ln_vc'.
ln_pt_vc :: Num a => Ln a -> (Pt a,Vc a)
ln_pt_vc l = (ln_start l,ln_vc l)

-- | 'Pt' UOp at 'Ln'.
ln_uop :: (Pt a -> Pt b) -> Ln a -> Ln b
ln_uop f (Ln l r) = Ln (f l) (f r)

-- | 'pt_translate' at 'Ln'.
ln_translate :: Num n => Vc n -> Ln n -> Ln n
ln_translate v (Ln a b) = Ln (pt_translate v a) (pt_translate v b)

-- | Variant with 'Pt' in place of 'Vc'.
ln_translate_pt :: Num n => Pt n -> Ln n -> Ln n
ln_translate_pt v (Ln a b) = Ln (a + v) (b + v)

ln_multiply :: Num n => Pt n -> Ln n -> Ln n
ln_multiply p (Ln a b) = Ln (a * p) (b * p)

-- | 'pt_scale' at 'Ln'.
ln_scale :: Num b => b -> Ln b -> Ln b
ln_scale m = ln_uop (pt_scale m)

-- | The angle, in /radians/, anti-clockwise from the /x/-axis.
--
-- > ln_angle (mk_ln ((0,0),(0,0))) == 0
-- > ln_angle (mk_ln ((0,0),(1,1))) == pi/4
-- > ln_angle (mk_ln ((0,0),(0,1))) == pi/2
-- > ln_angle (mk_ln ((0,0),(-1,1))) == pi * 3/4
ln_angle :: Ln R -> R
ln_angle ln =
    let Vc dx dy = ln_vc ln
    in if dx == 0 && dy == 0 then 0 else atan2 dy dx

-- | Start and end points of 'Ln'.
--
-- > ln_pt (Ln (Pt 1 0) (Pt 0 0)) == (Pt 1 0,Pt 0 0)
ln_pt :: Ln a -> (Pt a,Pt a)
ln_pt (Ln s e) = (s,e)

-- | Variant of 'ln_pt' giving co-ordinates as duples.
--
-- > ln_elem (Ln (Pt 1 0) (Pt 0 0)) == ((1,0),(0,0))
ln_elem :: Ln a -> ((a,a),(a,a))
ln_elem (Ln (Pt x1 y1) (Pt x2 y2)) = ((x1,y1),(x2,y2))

-- | Midpoint of a 'Ln'.
--
-- > ln_midpoint (Ln (Pt 0 0) (Pt 2 1)) == Pt 1 (1/2)
ln_midpoint :: Fractional a => Ln a -> Pt a
ln_midpoint (Ln (Pt x1 y1) (Pt x2 y2)) =
    let x = (x1 + x2) / 2
        y = (y1 + y2) / 2
    in Pt x y

-- | Variant on 'ln_midpoint'.
--
-- > cc_midpoint (Just (Pt 0 0),Nothing) == Pt 0 0
-- > cc_midpoint (Nothing,Just (Pt 2 1)) == Pt 2 1
-- > cc_midpoint (Just (Pt 0 0),Just (Pt 2 1)) == Pt 1 (1/2)
cc_midpoint :: (Maybe (Pt R), Maybe (Pt R)) -> Pt R
cc_midpoint cc =
    case cc of
      (Nothing,Nothing) -> Pt 0 0
      (Just p,Nothing) -> p
      (Nothing, Just q) -> q
      (Just p, Just q) -> ln_midpoint (Ln p q)

-- | Magnitude of 'Ln', ie. length of line.
--
-- > ln_magnitude (Ln (Pt 0 0) (Pt 1 1)) == sqrt 2
-- > pt_x (pt_to_polar (Pt 1 1)) == sqrt 2
ln_magnitude :: Ln R -> R
ln_magnitude = vc_mag . ln_vc

-- | Order 'Pt' at 'Ln' so that /p/ is to the left of /q/.  If /x/
-- fields are equal, sort on /y/.
--
-- > ln_sort (Ln (Pt 1 0) (Pt 0 0)) == Ln (Pt 0 0) (Pt 1 0)
-- > ln_sort (Ln (Pt 0 1) (Pt 0 0)) == Ln (Pt 0 0) (Pt 0 1)
ln_sort :: Ord a => Ln a -> Ln a
ln_sort ln =
    let Ln p q = ln
        Pt x1 y1 = p
        Pt x2 y2 = q
    in case compare x1 x2 of
         LT -> ln
         EQ -> if y1 <= y2 then ln else Ln q p
         GT -> Ln q p

-- | Adjust 'Ln' to have equal starting 'Pt' but magnitude 'R'.
--
-- > ln_adjust (sqrt 2) (Ln (Pt 0 0) (Pt 2 2)) == Ln (Pt 0 0) (Pt 1 1)
ln_adjust :: (Floating a, Ord a) => a -> Ln a -> Ln a
ln_adjust z ln =
    let Ln p _ = ln
        v = vc_scale z (vc_unit (ln_vc ln))
    in Ln p (pt_translate v p)

-- | Extend 'Ln' by 'R', ie. 'ln_adjust' with /n/ added to
-- 'ln_magnitude'.
--
-- > ln_extend (sqrt 2) (Ln (Pt 0 0) (Pt 1 1)) ~= Ln (Pt 0 0) (Pt 2 2)
ln_extend :: R -> Ln R -> Ln R
ln_extend n l = Ln (ln_start l) (pt_linear_extension n l)

-- | Variant definition of 'ln_extend'.
--
-- > ln_extend_ (sqrt 2) (Ln (Pt 0 0) (Pt 1 1)) == Ln (Pt 0 0) (Pt 2 2)
ln_extend_ :: R -> Ln R -> Ln R
ln_extend_ n l = ln_adjust (n + ln_magnitude l) l

-- | Calculate the point that extends a line by length 'n'.
--
-- > pt_linear_extension (sqrt 2) (Ln (Pt 1 1) (Pt 2 2)) ~= Pt 3 3
-- > pt_linear_extension 1 (Ln (Pt 1 1) (Pt 1 2)) ~= Pt 1 3
pt_linear_extension :: R -> Ln R -> Pt R
pt_linear_extension n (Ln p q) =
    let Pt mg ph = pt_to_polar (q - p)
    in pt_from_polar (Pt (mg + n) ph) + p

-- | Does 'Pt' /p/ lie on 'Ln' (inclusive).
--
-- > let {f = pt_on_line (Ln (Pt 0 0) (Pt 1 1))
-- >     ;r = [True,False,False,True]}
-- > in map f [Pt 0.5 0.5,Pt 2 2,Pt (-1) (-1),Pt 0 0] == r
pt_on_line :: Ln R -> Pt R -> Bool
pt_on_line l r =
    let (p,q) = ln_pt l
        Pt i j = pt_to_polar (q - p)
        Pt i' j' = pt_to_polar (r - p)
    in r == p || r == q || (j == j' && i' <= i)

-- | Vertical line.
ln_x_aligned :: a -> (a,a) -> Ln a
ln_x_aligned x (y0,y1) = Ln (Pt x y0) (Pt x y1)

-- | Horizontal line.
ln_y_aligned :: a -> (a,a) -> Ln a
ln_y_aligned y (x0,x1) = Ln (Pt x0 y) (Pt x1 y)

-- | Minimum Distance between a Point and a Line
-- <http://paulbourke.net/geometry/pointlineplane/>
--
-- > map (pt_ln_intersect_md_raw ((0,0),(10,10))) [(10,5),(15,0)] == [(7.5,7.5),(7.5,7.5)]
pt_ln_intersect_md_raw :: Fractional t => ((t,t),(t,t)) -> (t,t) -> (t,t,t)
pt_ln_intersect_md_raw ((x1,y1),(x2,y2)) (x3,y3) =
    let xd = x2 - x1
        yd = y2 - y1
        u = ((x3 - x1) * xd + (y3 - y1) * yd) / (xd * xd + yd * yd)
        x4 = x1 + u * (x2 - x1)
        y4 = y1 + u * (y2 - y1)
    in (u,x4,y4)

-- | Wrapper for 'pt_ln_intersect_md_raw'.
pt_ln_intersect_md :: Fractional a => Ln a -> Pt a -> (a,Pt a)
pt_ln_intersect_md ln pt =
    let (u,x,y) = pt_ln_intersect_md_raw (ln_elem ln) (pt_xy pt)
    in (u,Pt x y)

-- | 'pt_reflect_xy' about 'pt_ln_intersect_md'.
--
-- > map (pt_ln_reflect_md (Ln (Pt 0 0) (Pt 1 1))) [Pt 1 0,Pt 0.25 1] == [Pt 0 1,Pt 1 0.25]
pt_ln_reflect_md :: Fractional a => Ln a -> Pt a -> Pt a
pt_ln_reflect_md ln p =
    let (_,Pt x y) = pt_ln_intersect_md ln p
    in pt_reflect_xy (x,y) p

-- | Swap start and end points.
ln_reverse :: Ln a -> Ln a
ln_reverse (Ln p q) = Ln q p

-- * Intersection

-- | Intersection of two infinite lines, given as Pt and Vc.
ln_intersect_sg :: (Eq a,Fractional a) => (Pt a,Vc a) -> (Pt a,Vc a) -> Maybe (a,a)
ln_intersect_sg (Pt x1 y1,Vc dx1 dy1) (Pt x2 y2,Vc dx2 dy2) =
    let a = (dx2 * dy1) - (dx1 * dy2)
        t' = ((dx1 * (y2 - y1)) - (dy1 * (x2 - x1))) / a
        t = ((dx2 * (y1 - y2)) - (dy2 * (x1 - x2))) / (negate a)
    in if a == 0 then Nothing else Just (t,t')

ln_intersect :: (Eq t, Fractional t) => Ln t -> Ln t -> Maybe (t,t)
ln_intersect l1 l2 = ln_intersect_sg (ln_start l1,ln_vc l1) (ln_start l2,ln_vc l2)

-- | The 'Pt' at /z/ along 'Ln', 0 is the start of the line and 1 is the end.
ln_pt_along :: Num n => n -> Ln n -> Pt n
ln_pt_along z ln =
    let v = vc_scale z (ln_vc ln)
        Ln p _ = ln
    in pt_translate v p

-- | Do two 'Ln's intersect, and if so at which 'Pt'.
--
-- > ln_intersection (ln' (0,0) (5,5)) (ln' (5,0) (0,5)) == Just (Pt 2.5 2.5)
-- > ln_intersection (ln' (1,3) (9,3)) (ln' (0,1) (2,1)) == Nothing
-- > ln_intersection (ln' (1,5) (6,8)) (ln' (0.5,3) (6,4)) == Nothing
-- > ln_intersection (ln' (1,2) (3,6)) (ln' (2,4) (4,8)) == Nothing
-- > ln_intersection (ln' (2,3) (7,9)) (ln' (1,2) (5,7)) == Nothing
-- > ln_intersection (ln' (0,0) (1,1)) (ln' (0,0) (1,0)) == Just (Pt 0 0)
ln_intersection :: (Ord a,Fractional a) => Ln a -> Ln a -> Maybe (Pt a)
ln_intersection l0 l1 =
    case ln_intersect l0 l1 of
      Nothing -> Nothing
      Just (i,j) -> if i >= 0 && i <= 1 && j >= 0 && j <= 1
                    then Just (ln_pt_along i l0)
                    else Nothing

-- | Variant definition of 'ln_intersection', using algorithm at
-- <http://paulbourke.net/geometry/lineline2d/>.
--
-- > let comp_f x y = (ln_intersection x y,ln_intersection_ x y)
-- > comp_f (ln' (1,2) (3,6)) (ln' (2,4) (4,8)) == (Nothing,Nothing)
-- > comp_f (ln' (0,0) (1,1)) (ln' (0,0) (1,0)) == (Just (Pt 0 0),Just (Pt 0 0))
-- > comp_f (ln' (0,0) (1,2)) (ln' (0,1) (2,1)) == (Just (Pt 0.5 1),Just (Pt 0.5 1))
ln_intersection_ :: (Ord a,Fractional a) => Ln a -> Ln a -> Maybe (Pt a)
ln_intersection_ l0 l1 =
    let ((x1,y1),(x2,y2)) = ln_elem l0
        ((x3,y3),(x4,y4)) = ln_elem l1
        d = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
        ua' = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)
        ub' = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)
    in if d == 0
       then Nothing
       else if ua' == 0 && ub' == 0
            then Just (Pt x1 y1)
            else let ua = ua' / d
                     ub = ub' / d
                 in if in_range 0 1 ua && in_range 0 1 ub
                    then let x = x1 + ua * (x2 - x1)
                             y = y1 + ua * (y2 - y1)
                         in Just (Pt x y)
                    else Nothing

-- | Predicate variant of 'ln_intersection'.
--
-- > ln_intersect_p (ln' (1,1) (3,8)) (ln' (0.5,2) (4,7)) == True
-- > ln_intersect_p (ln' (3.5,9) (3.5,0.5)) (ln' (3,1) (9,1)) == True
ln_intersect_p :: (Ord a, Fractional a) => Ln a -> Ln a -> Bool
ln_intersect_p l = isJust . ln_intersection l

{- | Distances along a line, given as Pt and Vc, that it intersects with a circle.

> let o = Pt 0 0
> let c = (o,1)
> let z = sqrt 2 / 2
> map (\v -> ln_circle_intersection (o,v) c) [Vc 1 0,Vc 1 1] == [Just (1,-1),Just (z,-z)]

-}
ln_circle_intersection :: (Ord a,Floating a) => (Pt a,Vc a) -> (Pt a,a) -> Maybe (a,a)
ln_circle_intersection (Pt lx ly,Vc dx dy) (Pt cx cy,r) =
    let a = (dx * dx + dy * dy)
        b = 2 * ((lx - cx) * dx + (ly - cy) * dy)
        c = (lx - cx) * (lx - cx) + (ly - cy) * (ly - cy) - r * r
        z = b * b - 4 * a * c
    in if z < 0
       then Nothing
       else if a == 0
            then if c == 0 then Just (0,0) else Nothing
            else Just ((-b + sqrt z) / (2 * a)
                      ,(-b - sqrt z) / (2 * a))

-- | Variant with input as 'Ln' and output as list of actual intersection points.
ln_circle_intersection_set :: (Ord t, Floating t) => Ln t -> (Pt t,t) -> [Pt t]
ln_circle_intersection_set l c =
    case ln_circle_intersection (ln_pt_vc l) c of
      Nothing -> []
      Just (p,q) -> let f n = if n >= 0 && n <= 1 then Just n else Nothing
                    in map (flip ln_pt_along l) (catMaybes [f p,f q])

-- * Line slope

-- | Slope of 'Ln' or 'Nothing' if /vertical/.
--
-- > let l = map (ln' (0,0)) [(1,0),(2,1),(1,1),(0,1),(-1,1)]
-- > in map ln_slope l == [Just 0,Just (1/2),Just 1,Nothing,Just (-1)]
ln_slope :: (Fractional a,Eq a) => Ln a -> Maybe a
ln_slope l =
    let ((x1,y1),(x2,y2)) = ln_elem l
    in case x2 - x1 of
         0 -> Nothing
         dx -> Just ((y2 - y1) / dx)

-- | Are 'Ln's parallel, ie. have equal 'ln_slope'.  Note that the
-- direction of the 'Ln' is not relevant, ie. this is not equal to
-- 'ln_same_direction'.
--
-- > ln_parallel (ln' (0,0) (1,1)) (ln' (2,2) (1,1)) == True
-- > ln_parallel (ln' (0,0) (1,1)) (ln' (2,0) (1,1)) == False
-- > ln_parallel (ln' (1,2) (3,6)) (ln' (2,4) (4,8)) == True
-- > map ln_slope [ln' (2,2) (1,1),ln' (2,0) (1,1)] == [Just 1,Just (-1)]
ln_parallel :: (Ord a,Fractional a) => Ln a -> Ln a -> Bool
ln_parallel p q = ln_slope p == ln_slope q

-- | Are 'Ln's parallel, ie. have equal 'ln_angle'.
--
-- > ln_parallel_ (ln' (0,0) (1,1)) (ln' (2,2) (1,1)) == True
ln_parallel_ :: Ln R -> Ln R -> Bool
ln_parallel_ p q = ln_angle (ln_sort p) == ln_angle (ln_sort q)

-- | Are two vectors are in the same direction (to within a small
-- tolerance).
vc_same_direction :: (Ord a, Floating a) => Vc a -> Vc a -> Bool
vc_same_direction v w =
    let Vc dx1 dy1 = vc_unit v
        Vc dx2 dy2 = vc_unit w
    in abs (dx2 - dx1) < epsilon && abs (dy2 - dy1) < epsilon

-- | Do 'Ln's have same direction (within tolerance).
--
-- > ln_same_direction (ln' (0,0) (1,1)) (ln' (0,0) (2,2)) == True
-- > ln_same_direction (ln' (0,0) (1,1)) (ln' (2,2) (0,0)) == False
ln_same_direction :: (Ord a, Floating a) => Ln a -> Ln a -> Bool
ln_same_direction p q = ln_vc p `vc_same_direction` ln_vc q

-- | Are 'Ln's parallel, ie. does 'ln_vc' of each equal 'ln_same_direction'.
--
-- > ln_parallel__ (ln' (0,0) (1,1)) (ln' (2,2) (1,1)) == True
ln_parallel__ :: Ln R -> Ln R -> Bool
ln_parallel__ p q = ln_vc (ln_sort p) `vc_same_direction` ln_vc (ln_sort q)

-- | Is 'Ln' horizontal, ie. is 'ln_slope' zero.
--
-- > ln_horizontal (ln' (0,0) (1,0)) == True
-- > ln_horizontal (ln' (1,0) (0,0)) == True
ln_horizontal :: (Fractional a,Eq a) => Ln a -> Bool
ln_horizontal = (== Just 0) . ln_slope

-- | Is 'Ln' vertical, ie. is 'ln_slope' 'Nothing'.
--
-- > ln_vertical (ln' (0,0) (0,1)) == True
ln_vertical :: (Fractional a,Eq a) => Ln a -> Bool
ln_vertical = (== Nothing) . ln_slope

ln_minmax :: Ord a => Ln a -> (Pt a, Pt a)
ln_minmax (Ln (Pt x1 y1) (Pt x2 y2)) = (Pt (min x1 x2) (min y1 y2),Pt (max x1 x2) (max y1 y2))

ln_wn :: (Num n,Ord n) => Ln n -> Wn n
ln_wn = wn_from_extent . ln_minmax

-- | 'pt_rotate_about' at 'Ln'.
ln_rotate_about :: R -> Pt R -> Ln R -> Ln R
ln_rotate_about r c (Ln p q) = let f = pt_rotate_about r c in Ln (f p) (f q)

-- | Give translation and rotation from /p/ to /q/.
--  Magnitude is not considered.
--
-- > let l0 = Ln (Pt 0 0) (Pt 1 0)
-- > let l1 = Ln (Pt 1 1) (Pt 1 2)
-- > let (tr,(c,r)) = ln_align l0 l1
-- > ln_rotate_about r c (CG.ln_translate_pt tr l0) == l1
ln_align :: Ln R -> Ln R -> (Pt R, (Pt R, R))
ln_align p q =
    let (pt0,vc0) = ln_pt_vc p
        (pt1,vc1) = ln_pt_vc q
        (Vc _ ph0) = vc_to_polar vc0
        (Vc _ ph1) = vc_to_polar vc1
    in (pt1 - pt0,(pt1,ph1 -ph0))

-- * Ln sets

-- | 'pt_minmax' for set of 'Ln'.
lns_minmax :: Ord n => [Ln n] -> (Pt n,Pt n)
lns_minmax = ls_minmax . Ls . concatMap (\(Ln l r) -> [l,r])

lns_translate :: Num n => Vc n -> [Ln n] -> [Ln n]
lns_translate v = map (ln_translate v)

-- | Variant with 'Pt' not 'Vc'.
lns_translate_pt :: Num n => Pt n -> [Ln n] -> [Ln n]
lns_translate_pt p = map (ln_translate_pt p)

lns_multiply :: Num n => Pt n -> [Ln n] -> [Ln n]
lns_multiply p = map (ln_multiply p)

-- | Normalise to (0,m).
lns_normalise :: (Fractional n,Ord n) => n -> [Ln n] -> [Ln n]
lns_normalise m l =
    let w = wn_from_extent (lns_minmax l)
    in map (ln_scale m . ln_normalise_w w) l

-- * L(ine) s(egment) functions

-- | Alias for 'Ls'.
ls :: [Pt a] -> Ls a
ls = Ls

-- | Variant 'Ls' constructor from 'Pt' co-ordinates as duples.
mk_ls :: [(a,a)] -> Ls a
mk_ls = ls . map (uncurry Pt)

list_close :: [a] -> [a]
list_close l =
    case l of
      e:_ -> l ++ [e]
      _ -> error "list_close"

ls_close :: Ls t -> Ls t
ls_close = ls_elem_f list_close

ls_null :: Ls a -> Bool
ls_null = null . ls_elem

ls_map :: (Pt t -> Pt a) -> Ls t -> Ls a
ls_map f (Ls l) = Ls (map f l)

ls_multiply :: Num n => Pt n -> Ls n -> Ls n
ls_multiply p = ls_map (* p)

-- | Negate /y/ elements.
ls_negate_y :: (Num a) => Ls a -> Ls a
ls_negate_y = ls_map pt_negate_y

pts_minmax :: Ord a => [Pt a] -> (Pt a, Pt a)
pts_minmax s =
    case s of
      [] -> error "pts_minmax"
      p:ps -> foldl pt_minmax (p,p) ps

-- | Generate /minima/ and /maxima/ 'Point's from 'Ls'.
ls_minmax :: Ord a => Ls a -> (Pt a,Pt a)
ls_minmax = pts_minmax . ls_elem

-- | Separate 'Ls' at points where the 'Vc' from one element to the
-- next exceeds the indicated distance.
--
-- > map length (ls_separate (Vc 2 2) (map (uncurry Pt) [(0,0),(1,1),(3,3)])) == [2,1]
ls_separate :: (Ord a,Num a) => Vc a -> Ls a -> [Ls a]
ls_separate (Vc dx dy) (Ls l) =
    let f (Pt x0 y0) (Pt x1 y1) = abs (x1 - x0) < dx &&
                                  abs (y1 - y0) < dy
    in map Ls (segment_f f l)

-- | Delete 'Pt' from 'Ls' so that no two 'Pt' are within a tolerance
-- given by 'Vc'.
ls_tolerate :: (Ord a,Num a) => Vc a -> Ls a -> Ls a
ls_tolerate (Vc x y) (Ls l) =
    let too_close (Pt x0 y0) (Pt x1 y1) =
            let dx = abs (x1 - x0)
                dy = abs (y1 - y0)
            in dx < x && dy < y
    in Ls (delete_f too_close l)

-- | Variant of 'ls_tolerate' where 'Vc' is optional, and 'Nothing' gives 'id'.
ls_tolerate_maybe :: (Ord a,Num a) => Maybe (Vc a) -> Ls a -> Ls a
ls_tolerate_maybe i =
    case i of
      Nothing -> id
      Just i' -> ls_tolerate i'

-- | Test if point 'Pt' lies inside polygon 'Ls'.
--
-- > ls_pt_inside (mk_ls [(0,0),(1,0),(1,1),(0,1)]) (Pt 0.5 0.5) == True
ls_pt_inside :: Ls R -> Pt R -> Bool
ls_pt_inside (Ls s) (Pt x y) =
    case s of
      [] -> undefined
      l0:l -> let xs = pairs ((l0:l)++[l0])
                  f (Pt x1 y1,Pt x2 y2) =
                      and [y > min y1 y2
                          ,y <= max y1 y2
                          ,x <= max x1 x2
                          ,y1 /= y2
                          ,x1 == x2 ||
                           x <= (y-y1)*(x2-x1)/(y2-y1)+x1]
              in odd (length (filter id (map f xs)))

-- | Variant that counts points at vertices as inside.
--
-- > ls_pt_inside_or_vertex (mk_ls [(0,0),(1,0),(1,1),(0,1)]) (Pt 0 1) == True
ls_pt_inside_or_vertex :: Ls R -> Pt R -> Bool
ls_pt_inside_or_vertex l p = p `elem` ls_elem l || ls_pt_inside l p

-- | Check all 'Pt' at 'Ls' are 'pt_is_normal'.
ls_check_normalised :: (Ord a,Num a) => Ls a -> Bool
ls_check_normalised (Ls s) = all pt_is_normal s

-- | Line co-ordinates as /x/,/y/ list.
--
-- > ls_xy [Pt 0 0,Pt 1 1] == [0,0,1,1]
ls_xy :: Ls a -> [a]
ls_xy = concatMap (\(Pt x y) -> [x,y]) . ls_elem

-- | 'Ls' average.
ls_centroid :: Fractional t => Ls t -> Pt t
ls_centroid (Ls l) =
    let (x,y) = unzip (map pt_xy l)
        length' = fromIntegral . length
    in Pt (sum x / length' x) (sum y / length' y)

-- | 'ls_map' of 'pt_rotate_about'.
ls_rotate_about :: Floating t => t -> Pt t -> Ls t -> Ls t
ls_rotate_about r c = ls_map (pt_rotate_about r c)

-- | 'ls_rotate_about' of 'ls_centroid'.
ls_rotate_about_centroid :: Floating t => t -> Ls t -> Ls t
ls_rotate_about_centroid r l = ls_rotate_about r (ls_centroid l) l

ls_elem_f :: ([Pt t] -> [Pt u]) -> Ls t -> Ls u
ls_elem_f f = Ls . f . ls_elem

ls_unlift :: (Ls t -> Ls u) -> [Pt t] -> [Pt u]
ls_unlift f = ls_elem . f . Ls

-- | Midpoints of line segments.
ls_midpoints :: Ls R -> Ls R
ls_midpoints =
    let adj l = zip l (tail l)
    in ls_elem_f (map (\(p,q) -> ln_midpoint (Ln p q)) . adj)

-- | 'pt_x' of 'ls_minmax'.
ls_minmax_x :: Ord t => Ls t -> (t,t)
ls_minmax_x = bimap1 pt_x . ls_minmax

-- | 'pt_y' of 'ls_minmax'.
ls_minmax_y :: Ord t => Ls t -> (t,t)
ls_minmax_y = bimap1 pt_y . ls_minmax

ls_to_ln_set :: Ls t -> [Ln t]
ls_to_ln_set (Ls l) = zipWith Ln l (tail l)

-- | Points at which /ln/ intersects /ls/.
ls_intersections :: Ls R -> Ln R -> [Pt R]
ls_intersections p ln = mapMaybe (ln_intersection ln) (ls_to_ln_set p)

ls_intersect_p :: (Fractional a, Ord a) => Ls a -> Ls a -> Bool
ls_intersect_p p q = any id [ln_intersect_p r s | r <- ls_to_ln_set p, s <- ls_to_ln_set q]

-- | @y@ co-ordinates of intersection with vertical line given by @(y0,y1)@ and @x@.
ls_intersections_y :: (R,R) -> Ls R -> R -> [R]
ls_intersections_y y p x = map pt_y (ls_intersections p (ln_x_aligned x y))

-- | @x@ co-ordinates of intersection with horizontal line given by @(x0,x1)@ and @y@.
ls_intersections_x :: (R,R) -> Ls R -> R -> [R]
ls_intersections_x x p y = map pt_x (ls_intersections p (ln_y_aligned y x))

-- | Shift each 'Pt' set in a sequence, the amount and direction are given
-- by /get_mm/ and /mk_vc/.
pts_sequence :: Num a => ([Pt a] -> (a,a)) -> (a -> Vc a) -> [[Pt a]] -> [[Pt a]]
pts_sequence get_mm get_vc l =
    let f st e = let (x0,x1) = get_mm e
                 in (st + (x1 - x0)
                    ,map (pt_translate (get_vc (st - x0))) e)
    in case l of
         l0:l' -> l0 : snd (mapAccumL f (snd (get_mm l0)) l')
         [] -> l

pts_sequence_right :: (Num t,Ord t) => t -> [[Pt t]] -> [[Pt t]]
pts_sequence_right k = pts_sequence (bimap id (+ k) . bimap1 pt_x . pts_minmax) (\x -> Vc x 0)

pts_sequence_above :: (Num t,Ord t) => t -> [[Pt t]] -> [[Pt t]]
pts_sequence_above k = pts_sequence (bimap id (+ k) . bimap1 pt_y . pts_minmax) (\y -> Vc 0 y)

ls_sequence :: Num a => (Ls a -> (a,a)) -> (a -> Vc a) -> [Ls a] -> [Ls a]
ls_sequence get_mm get_vc l = map Ls (pts_sequence (get_mm . Ls) get_vc (map ls_elem l))

-- | Shift each 'Ls' so that it is /k/ to the right of the rightmost point of the preceding 'Ls'
ls_sequence_right :: (Num t,Ord t) => t -> [Ls t] -> [Ls t]
ls_sequence_right k = ls_sequence (bimap id (+ k) . ls_minmax_x) (\x -> Vc x 0)

-- | Vertical variant of 'ls_sequence_right'.
ls_sequence_above :: (Num t,Ord t) => t -> [Ls t] -> [Ls t]
ls_sequence_above k = ls_sequence (bimap id (+ k) . ls_minmax_y) (\y -> Vc 0 y)

-- * Window

-- | Variant 'Wn' constructor.
wn' :: (a,a) -> (a,a) -> Wn a
wn' (x,y) (i,j) = Wn (Pt x y) (Vc i j)

-- | Extract /(x,y)/ and /(dx,dy)/ pairs.
--
-- > wn_extract (Wn (Pt 0 0) (Vc 1 1)) == ((0,0),(1,1))
wn_extract :: Wn a -> ((a,a),(a,a))
wn_extract (Wn (Pt x y) (Vc dx dy)) = ((x,y),(dx,dy))

-- | Show function for window with fixed precision of 'n'.
--
-- > wn_show 1 (Wn (Pt 0 0) (Vc 1 1)) == "((0.0,0.0),(1.0,1.0))"
wn_show :: Int -> Wn R -> String
wn_show n (Wn (Pt x0 y0) (Vc dx dy)) =
    let fs = printf "((%%.%df,%%.%df),(%%.%df,%%.%df))" n n n n
    in printf fs x0 y0 dx dy

-- | Unit window, lower left at origin.
wn_unit :: Num n => Wn n
wn_unit = Wn pt_origin (Vc 1 1)

-- | Square bounding circle at 'Pt' with radius.
--
-- > wn_square (Pt 0 0) 1 == wn' (-1,-1) (2,2)
wn_square :: Num n => Pt n -> n -> Wn n
wn_square (Pt x y) n = Wn (Pt (x - n) (y - n)) (Vc (n * 2) (n * 2))

-- | Square window, center at origin.  /n/ is half the length of each
-- side of the square.
--
-- > wn_square_o 1 == wn' (-1,-1) (2,2)
wn_square_o :: Num n => n -> Wn n
wn_square_o = wn_square (Pt 0 0)

-- | Is 'Pt' within 'Wn' exclusive of edge.
--
-- > map (pt_in_window (wn' (0,0) (1,1))) [Pt 0.5 0.5,Pt 1 1] == [True,False]
pt_in_window :: (Ord a,Num a) => Wn a -> Pt a -> Bool
pt_in_window (Wn (Pt lx ly) (Vc dx dy)) (Pt x y) =
    let (ux,uy) = (lx+dx,ly+dy)
    in x > lx && x < ux && y > ly && y < uy

-- | 'Wn' from /(lower-left,upper-right)/ extent.
wn_from_extent :: Num a => (Pt a,Pt a) -> Wn a
wn_from_extent (Pt x0 y0,Pt x1 y1) = Wn (Pt x0 y0) (Vc (x1 - x0) (y1 - y0))

-- | 'Wn' containing 'Ls'.
--
-- > ls_window (mk_ls [(0,0),(1,1),(2,0)]) == wn' (0,0) (2,1)
ls_window :: (Num a,Ord a) => Ls a -> Wn a
ls_window = wn_from_extent . ls_minmax

pts_window :: (Num t,Ord t) => [Pt t] -> Wn t
pts_window = wn_from_extent . pts_minmax

-- | A 'Wn' that encompasses both input 'Wn's.
wn_join :: (Num a,Ord a) => Wn a -> Wn a -> Wn a
wn_join (Wn (Pt x0 y0) (Vc dx0 dy0)) (Wn (Pt x1 y1) (Vc dx1 dy1)) =
    let x = min x0 x1
        y = min y0 y1
        dx = max (x0+dx0) (x1+dx1) - x
        dy = max (y0+dy0) (y1+dy1) - y
    in Wn (Pt x y) (Vc dx dy)

wn_join_l :: (Num a,Ord a) => [Wn a] -> Wn a
wn_join_l = foldl1 wn_join

-- | Predictate to determine if two 'Wn's intersect.
wn_intersect :: (Num a,Ord a) => Wn a -> Wn a -> Bool
wn_intersect w0 w1  =
    let (Wn (Pt x0 y0) (Vc dx0 dy0)) = w0
        (Wn (Pt x1 y1) (Vc dx1 dy1)) = w1
    in not (x0 > x1+dx1 || x1 > x0+dx0 || y0 > y1+dy1 || y1 > y0+dy0)

-- | Are all points at 'Ls' within the 'Wn'.
ls_in_window :: Wn R -> Ls R -> Bool
ls_in_window w = all (pt_in_window w). ls_elem

-- | Are any points at 'Ls' within the window 'Wn'.
ls_enters_window :: Wn R -> Ls R -> Bool
ls_enters_window w = any (pt_in_window w) . ls_elem

-- | Are all points at 'Ls' outside the 'Wn'.
ls_not_in_window :: Wn R -> Ls R -> Bool
ls_not_in_window w = all (not . pt_in_window w) . ls_elem

-- | Break 'Ls' into segments that are entirely within the 'Wn'.
ls_segment_window :: Wn R -> Ls R -> [Ls R]
ls_segment_window w =
    let g [] = []
        g xs = let (i,xs') = span (pt_in_window w) xs
               in i : g (dropWhile (not . pt_in_window w) xs')
    in map Ls . filter (not . null) . g . ls_elem

-- | Normalisation function for 'Wn', ie. map 'Pt' to lie within (0,1).
wn_normalise_f :: (Ord n,Fractional n) => Wn n -> Pt n -> Pt n
wn_normalise_f (Wn (Pt x0 y0) (Vc dx dy)) (Pt x y) =
    let z = max dx dy
    in Pt ((x - x0) / z) ((y - y0) / z)

-- | Given 'Wn' normalise the 'Ls'.
ls_normalise_w :: (Ord n,Fractional n) => Wn n -> Ls n -> Ls n
ls_normalise_w w = ls_elem_f (map (wn_normalise_f w))

-- | Given 'Wn' normalise 'Ln'.
ln_normalise_w :: (Ord n,Fractional n) => Wn n -> Ln n -> Ln n
ln_normalise_w w (Ln p q) =
    let f = wn_normalise_f w
    in Ln (f p) (f q)

-- | Variant of 'ls_normalise_w', the window is determined by the
-- extent of the 'Ls'.
ls_normalise :: (Ord n,Fractional n) => Ls n -> Ls n
ls_normalise l = ls_normalise_w (wn_from_extent (ls_minmax l)) l

ls_concat :: [Ls a] -> Ls a
ls_concat = Ls . concat . map ls_elem

-- | Normalise a set of line segments using composite window.
ls_normalise_set :: (Ord n,Fractional n) => [Ls n] -> [Ls n]
ls_normalise_set l =
    let w = wn_from_extent (ls_minmax (ls_concat l))
    in map (ls_normalise_w w) l

-- | Shift lower left 'Pt' of 'Wn' by indicated 'Pt'.
pt_shift_w :: Num a => Pt a -> Wn a -> Wn a
pt_shift_w p (Wn dp ex) = Wn (p + dp) ex

-- | Negate /y/ field of lower left 'Pt' of 'Wn'.
wn_negate_y :: Num a => Wn a -> Wn a
wn_negate_y (Wn p v) = Wn (pt_negate_y p) v

-- | 'Wn' to 'Ls' (CCW), open.
wn_to_ls :: Num t => Wn t -> Ls t
wn_to_ls (Wn (Pt x y) (Vc dx dy)) =
    let (x',y') = (x + dx,y + dy)
    in Ls [Pt x y,Pt x y',Pt x' y',Pt x' y]

-- | Closed form.
wn_to_ls_closed :: Num t => Wn t -> Ls t
wn_to_ls_closed = ls_close .  wn_to_ls

-- | Vector giving width and height of each cell of an (r,c) grid at w.
--
-- > wn_grid_vc div (wn_square_o 200) (8,8) == Vc 50 50
wn_grid_vc :: (t -> t -> t) -> Wn t -> (t,t) -> Vc t
wn_grid_vc div_f (Wn _ (Vc dx dy)) (r,c) = Vc (dx `div_f` c) (dy `div_f` r)

-- | Grid dividing window into indicated number of rows and columns.
-- Row order, lowest row first.  Points at lower left.
--
-- > wn_grid_pt_ll div (wn_square_o 200) (8,8)
wn_grid_pt_ll :: (Enum t,Num t) => (t -> t -> t) -> Wn t -> (t,t) -> [[Pt t]]
wn_grid_pt_ll div_f w g =
    let Vc ix iy = wn_grid_vc div_f w g
        Wn (Pt x y) (Vc dx dy) = w
        f y' = map (\x' -> Pt x' y') [x,x + ix .. x + dx - ix]
    in map f [y,y + iy .. y + dy - iy]

-- | Variant with points at center.
--
-- > wn_grid_pt_c div (wn_square_o 200) (8,8)
wn_grid_pt_c :: (Enum t,Num t) => (t -> t -> t) -> Wn t -> (t,t) -> [[Pt t]]
wn_grid_pt_c div_f w g =
    let Vc ix iy = wn_grid_vc div_f w g
        v = Vc (ix `div_f` 2) (iy `div_f` 2)
    in map (map (\p -> pt_translate v p)) (wn_grid_pt_ll div_f w g)

-- | Variant with derived 'Wn' for each cell.
wn_grid_wn :: (Enum t,Num t) => (t -> t -> t) -> t -> Wn t -> (t,t) -> [[Wn t]]
wn_grid_wn div_f sc w g =
    let v = vc_scale sc (wn_grid_vc div_f w g)
        v' = vc_uop (\n -> negate (n `div_f` 2)) v
    in map (map (\p -> Wn (pt_translate v' p) v)) (wn_grid_pt_c div_f w g)

-- | Variant with derived closed 'Ls' for each cell.
--
-- > concat $ wn_grid_ls (/) 0.75 (wn_square_o 200) (8,8)
wn_grid_ls :: (Enum t,Num t) => (t -> t -> t) -> t -> Wn t -> (t,t) -> [[Ls t]]
wn_grid_ls div_f sc w = map (map wn_to_ls_closed) . wn_grid_wn div_f sc w

-- * Random

-- | Generate a random 'Pt' within 'Wn'.
--
-- > pt_random (wn_square_o 1.0)
pt_random :: (Num a,Random a) => Wn a -> IO (Pt a)
pt_random (Wn (Pt x y) (Vc dx dy)) = do
  x' <- randomRIO (x,x + dx)
  y' <- randomRIO (y,y + dy)
  return (Pt x' y')

-- | Generate a random 'P3' within cube.
--
-- > p3_random_u (-1.0,1.0)
p3_random :: Random a => (a,a) -> IO (P3 a)
p3_random (l,r) = do
  x <- randomRIO (l,r)
  y <- randomRIO (l,r)
  z <- randomRIO (l,r)
  return (P3 x y z)

-- | Generate a random 'Ln' within 'Wn'.
--
-- > ln_random (wn_square_o 1.0)
ln_random :: (Num a,Random a) => Wn a -> IO (Ln a)
ln_random w = do
  p <- pt_random w
  q <- pt_random w
  return (Ln p q)

-- * Matrix

-- | A translation matrix with independent x and y offsets.
mx_translation :: Num n => n -> n -> Matrix n
mx_translation = Matrix 1 0 0 1

-- | A scaling matrix with independent x and y scalars.
mx_scaling :: Num n => n -> n -> Matrix n
mx_scaling x y = Matrix x 0 0 y 0 0

-- | A rotation matrix through the indicated angle (in radians).
mx_rotation :: Floating n => n -> Matrix n
mx_rotation a =
    let c = cos a
        s = sin a
        t = negate s
    in Matrix c s t c 0 0

-- | The identity matrix.
mx_identity :: Num n => Matrix n
mx_identity = Matrix 1 0 0 1 0 0

mx_translate :: Num n => n -> n -> Matrix n -> Matrix n
mx_translate x y m = m * (mx_translation x y)

mx_scale :: Num n => n -> n -> Matrix n -> Matrix n
mx_scale x y m = m * (mx_scaling x y)

mx_rotate :: Floating n => n -> Matrix n -> Matrix n
mx_rotate r m = m * (mx_rotation r)

mx_scalar_multiply :: Num n => n -> Matrix n -> Matrix n
mx_scalar_multiply scalar = mx_uop (* scalar)

mx_adjoint :: Num n => Matrix n -> Matrix n
mx_adjoint (Matrix a b c d x y) =
    Matrix d (-b) (-c) a (c * y - d * x) (b * x - a * y)

mx_invert :: Fractional n => Matrix n -> Matrix n
mx_invert m =
    let Matrix xx yx xy yy _ _ = m
        d = xx * yy - yx * xy
    in mx_scalar_multiply (recip d) (mx_adjoint m)

mx_list :: Matrix n -> [n]
mx_list (Matrix a b c d e f) = [a,b,c,d,e,f]

-- | Apply a transformation matrix to a point.
pt_transform :: Num n => Matrix n -> Pt n -> Pt n
pt_transform (Matrix a1 a2 b1 b2 c1 c2) (Pt x y) =
    let x' = x * a1 + y * b1 + c1
        y' = x * a2 + y * b2 + c2
    in Pt x' y'

-- * Polygon

-- | {A,B,C,D,E} are vertices.
-- /A/ is the interior angle at vertex /A/.
-- /a/ is the distance from /E/ to /A/.
--
-- > let f p = map (pt_xy . pt_uop round) . polygon_unfold (Pt 0 0,0) p . map degrees_to_radians
-- > let rp = replicate
-- > let g n p q = f (rp n p) (rp n q)
-- > g 3 100 60 == [(0,0),(100,0),(50,87)]
-- > g 4 100 90 == [(0,0),(100,0),(100,100),(0,100)]
-- > g 5 100 108 == [(0,0),(100,0),(131,95),(50,154),(-31,95)]
-- > g 6 100 120 == [(0,0),(100,0),(150,87),(100,173),(0,173),(-50,87)]
-- > g 7 100 (900/7) == [(0,0),(100,0),(162,78),(140,176),(50,219),(-40,176),(-62,78)]
-- > g 8 100 135 == [(0,0),(100,0),(171,71),(171,171),(100,241),(0,241),(-71,171),(-71,71)]
-- > f [30,82.24,73.27,60,60] [90,90,101.5,101.5,157] == [(0,0),(30,0),(30,82),(-43,82),(-55,23)]
polygon_unfold :: (Pt R,R) -> [R] -> [R] -> [Pt R]
polygon_unfold i p q =
    let f (pt,ph0) (mg,ph) =
            let r = pt_translate (mk_vc_polar (mg,ph0)) pt
            in ((r,ph0 + (pi - ph)),r)
        rem_last = reverse . tail . reverse -- hmm...
    in fst i : rem_last (snd (mapAccumL f i (zip p q)))

-- | Inverse of 'polygon_unfold'.
--
-- > let f = polygon_param . map mk_pt
-- > let g (p,q) = (map round p,map (round . radians_to_degrees) q)
-- > let h = g . f
-- > h [(0,0),(100,0),(50,87)] == ([100,100,100],[60,60,60])
-- > h [(0,0),(100,0),(131,95),(50,154),(-31,95)] == ([100,100,100,100,100],[108,108,108,108,108])
-- > h [(0,0),(30,0),(30,82),(-43,82),(-55,23)] == ([30,82,73,60,60],[157,90,90,101,101])
polygon_param :: [Pt R] -> ([R],[R])
polygon_param p =
    let adj l = zip l (tail (cycle l))
        vc = map (vc_to_polar . ln_vc . uncurry Ln) (adj p)
        f z ph = let r = pi - (ph - z) in (ph,r `F.mod'` two_pi)
        ph0 = vc_y (last vc)
    in (map vc_x vc
       ,snd (mapAccumL f ph0 (map vc_y vc)))

-- * P3 functions

-- | Tuple constructor.
p3' :: (a,a,a) -> P3 a
p3' (x,y,z) = P3 x y z

-- | Tuple accessor.
p3_xyz :: P3 t -> (t,t,t)
p3_xyz (P3 x y z) = (x,y,z)

-- | 'P3' of (0,0,0).
--
-- > p3_origin == P3 0 0 0
p3_origin :: Num a => P3 a
p3_origin = P3 0 0 0

-- | 'P3' at /(n,n,n)/.
--
-- > p3_from_scalar 1 == P3 1 1 1
p3_from_scalar :: a -> P3 a
p3_from_scalar a = P3 a a a

-- | Scalar 'P3' '+'.
--
-- > p3_offset 1 p3_origin == P3 1 1 1
p3_offset :: Num a => a -> P3 a -> P3 a
p3_offset = p3_uop . (+)

-- | Scalar 'P3' '*'.
--
-- > p3_scale 2 (P3 1 2 3) == P3 2 4 6
p3_scale :: Num a => a -> P3 a -> P3 a
p3_scale = p3_uop . (*)

-- * Ord

-- | Given /left/ and /right/, is /x/ in range (inclusive).
--
-- > map (in_range 0 1) [-1,0,1,2] == [False,True,True,False]
in_range :: Ord a => a -> a -> a -> Bool
in_range l r x = l <= x && x <= r

-- * List

-- | Split list at element where predicate /f/ over adjacent elements
-- first holds.
--
-- > split_f (\p q -> q - p < 3) [1,2,4,7,11] == ([1,2,4],[7,11])
split_f :: (a -> a -> Bool) -> [a] -> ([a],[a])
split_f f =
    let go i [] = (reverse i,[])
        go i [p] = (reverse (p:i), [])
        go i (p:q:r) =
            if f p q
            then go (p:i) (q:r)
            else (reverse (p:i),q:r)
    in go []

-- | Variant on 'split_f' that segments input.
--
-- > segment_f (\p q -> abs (q - p) < 3) [1,3,7,9,15] == [[1,3],[7,9],[15]]
segment_f :: (a -> a -> Bool) -> [a] -> [[a]]
segment_f f xs =
    let (p,q) = split_f f xs
    in if null q
       then [p]
       else p : segment_f f q

-- | Delete elements of a list using a predicate over the
-- previous and current elements.
delete_f :: (a -> a -> Bool) -> [a] -> [a]
delete_f f =
    let go [] = []
        go [p] = [p]
        go (p:q:r) =
            if f p q
            then go (p:r)
            else p : go (q:r)
    in go

-- | All adjacent pairs of a list.
--
-- > pairs [1..5] == [(1,2),(2,3),(3,4),(4,5)]
pairs :: [x] -> [(x,x)]
pairs l =
    case l of
      x:y:z -> (x,y) : pairs (y:z)
      _ -> []

-- * Bifunctor

-- > bimap id succ (0,0) == (0,1)
-- > fmap succ (0,0) == (0,1)
bimap :: (a -> b) -> (c -> d) -> (a,c) -> (b,d)
bimap f g (p,q) = (f p,g q)

bimap1 :: (c -> d) -> (c,c) -> (d,d)
bimap1 f (p,q) = (f p,f q)

