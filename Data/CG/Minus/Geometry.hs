{- | Simple geometric functions.

Polygon variables:

n = degree, a = side length, r,ir = inradius, R,cr = circumradius, A,ar = area, s = sagitta

-}
module Data.CG.Minus.Geometry where

import Data.Complex {- base -}

-- * Math

two_pi :: Floating n => n
two_pi = 2 * pi

-- | Square, inverse of 'sqrt'.
--
-- > sqrt (sqr 3) == 3
sqr :: Num a => a -> a
sqr n = n * n

-- | Secant.
sec :: Floating a => a -> a
sec z = 1 / cos z

-- | Cotangent.
cot :: Floating a => a -> a
cot z = 1 / tan z

-- | Degrees to radians.
--
-- > degrees_to_radians 60 == (pi/3)
degrees_to_radians :: Floating n => n -> n
degrees_to_radians = (* pi) . (/ 180)

-- > radians_to_degrees (pi/3) == 60
radians_to_degrees :: Floating n => n -> n
radians_to_degrees = (* 180) . (/ pi)

polar_to_rectangular :: Floating t => (t,t) -> (t,t)
polar_to_rectangular (mg,ph) =
    let c = mkPolar mg ph
    in (realPart c,imagPart c)

rectangular_to_polar :: RealFloat a => (a, a) -> (a,a)
rectangular_to_polar (x,y) = polar (x :+ y)

-- * Triangle

-- > triangle_semiperimeter_sss 1 1 1 == 1.5
triangle_semiperimeter_sss :: Fractional a => a -> a -> a -> a
triangle_semiperimeter_sss a b c = (1/2) * (a + b + c)

-- > triangle_area_herons_sss 1 1 1
triangle_area_herons_sss :: Floating a => a -> a -> a -> a
triangle_area_herons_sss a b c =
    let s = triangle_semiperimeter_sss a b c
    in sqrt (s * (s - a) * (s - b) * (s - c))

-- > triangle_area_aas (radians 120) (radians 30) 1
triangle_area_aas :: Floating a => a -> a -> a -> a
triangle_area_aas alpha beta a =
    let gamma = pi - alpha - beta
        n = sqr a * sin beta * sin gamma
        d = 2 * sin alpha
    in n / d

-- > triangle_area_asa (radians 30) 1 (radians 30)
triangle_area_asa :: Floating a => a -> a -> a -> a
triangle_area_asa alpha c beta =
    let n = sqr c
        d = 2 * (cot alpha + cot beta)
    in n / d

-- > triangle_side_asa (radians 30) 1 (radians 30)
triangle_side_asa :: Floating a => a -> a -> a -> a
triangle_side_asa alpha c beta =
    let gamma = pi - alpha - beta
    in (sin alpha / sin gamma) * c

-- | <http://mathworld.wolfram.com/LawofCosines.html>
law_of_cosines :: Floating a => a -> a -> a -> a
law_of_cosines s a s' = sqrt (sqr s + sqr s' - 2 * s * s' * cos a)

-- > triangle_side_sas 0.5 (radians 120) 0.5
triangle_side_sas :: Floating a => a -> a -> a -> a
triangle_side_sas = law_of_cosines

triangle_centroid :: Fractional t => (t,t) -> (t,t) -> (t,t) -> (t,t)
triangle_centroid (x1,y1) (x2,y2) (x3,y3) =
    let x = (x1 + x2 + x3) / 3
        y = (y1 + y2 + y3) / 3
    in (x,y)

-- * Pentagon

-- > pentagon_circumradius_a 1
pentagon_circumradius_a :: Floating a => a -> a
pentagon_circumradius_a a = 0.1 * sqrt (50 + 10 * sqrt 5) * a

-- > pentagon_inradius_a 1
pentagon_inradius_a :: Floating a => a -> a
pentagon_inradius_a a = 0.1 * sqrt (25 + 10 * sqrt 5) * a

-- > pentagon_sagitta_a 1
pentagon_sagitta_a :: Floating a => a -> a
pentagon_sagitta_a a = 0.1 * sqrt (25 - 10 * sqrt 5) * a

-- > pentagon_area_a 1
pentagon_area_a :: Floating a => a -> a
pentagon_area_a a =  0.25 * sqrt (25 + 10 * sqrt 5) * a * a

-- * Hexagon

-- > hexagon_inradius_a 1 == 0.8660254037844386
hexagon_inradius_a :: Floating a => a -> a
hexagon_inradius_a a = 0.5 * sqrt 3 * a

-- | For hexagons, the side length (a) and the circumradius (R) are equal.
--
-- > hexagon_circumradius_a 1 == 1
hexagon_circumradius_a :: a -> a
hexagon_circumradius_a = id

-- > hexagon_sagitta_a 1 == 0.1339745962155614
hexagon_sagitta_a :: Floating a => a -> a
hexagon_sagitta_a a = 0.5 * (2 - sqrt 3) * a

-- > hexagon_area_a 1 == 2.598076211353316
hexagon_area_a :: Floating a => a -> a
hexagon_area_a a = 1.5 * sqrt 3 * sqr a

-- * Polygon

{- | <http://mathworld.wolfram.com/PolygonArea.html>

Note that the area of a convex polygon is defined to be positive if
the points are arranged in a counterclockwise order, and negative
if they are in clockwise order (Beyer 1987).

> let u = [(0,0),(1,0),(1,1),(0,1)]
> polygon_signed_area u  == 1
> polygon_signed_area (reverse u) == -1

-}
polygon_signed_area :: Fractional t => [(t,t)] -> t
polygon_signed_area p =
    let q = zip p (tail (cycle p))
        f ((x1,y1),(x2,y2)) = x1 * y2 - x2 * y1
    in sum (map f q) / 2

-- > regular_polygon_side_length 5 1
regular_polygon_side_length :: Floating a => a -> a -> a
regular_polygon_side_length n cr = 2 * cr * sin (pi / n)

-- > map (\n -> regular_polygon_inradius n 1) [3,4,5,6]
regular_polygon_inradius :: Floating a => a -> a -> a
regular_polygon_inradius n cr = cr * cos (pi / n)

-- > regular_polygon_circumradius 5 1
regular_polygon_circumradius :: Floating a => a -> a -> a
regular_polygon_circumradius n r = r * sec (pi / n)

-- > regular_polygon_area 6 1
regular_polygon_area :: Floating a => a -> a -> a
regular_polygon_area n cr = 0.5 * n * sqr cr * sin (two_pi / n)

-- > regular_polygon_sagitta 6 1
regular_polygon_sagitta :: Floating a => a -> a -> a
regular_polygon_sagitta n cr = 2 * cr * sqr (sin (pi / (2 * n)))
