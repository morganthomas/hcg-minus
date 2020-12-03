-- | Compass bearings.
module Data.CG.Minus.Bearing where

import Data.CG.Minus

-- | Enumeration of compass bearings (32 divisions).
data Bearing = N | NbE | NNE | NEbN | NE | NEbE | ENE | EbN
             | E | EbS | ESE | SEbE | SE | SEbS | SSE | SbE
             | S | SbW | SSW | SWbS | SW | SWbW | WSW | WbS
             | W | WbN | WNW | NWbW | NW | NWbN | NNW | NbW
               deriving (Eq,Enum,Bounded,Show)

-- | The 4, 8 and 16 element subsets of bearings, and all 32 bearings.
bearings_4,bearings_8,bearings_16,bearings_32 :: [Bearing]
bearings_4 = [N,E,S,W]
bearings_8 = [N,NE,E,SE,S,SW,W,NW]
bearings_16 = [N,NNE,NE,ENE,E,ESE,SE,SSE,S,SSW,SW,WSW,W,WNW,NW,NNW]
bearings_32 = [N .. NbW]

-- | The point on the unit circle of 'Bearing'.
--
-- > map bearing_pt_u [N,E,S,W]
-- > let n = sqrt 0.5 in pt_eq_by (~=) (bearing_pt_u NE) (Pt n n)
bearing_pt_u :: Bearing -> Pt R
bearing_pt_u b = pt_from_polar (Pt 1 (bearing_angle b))

-- | 'E' is zero, direction is counter-clockwise, in radians, in (-pi,pi).
--
-- > map (r_from_radians . bearing_angle) bearings_8 == [90,45,0,-45,-90,-135,180,135]
bearing_angle :: Bearing -> R
bearing_angle b =
    let ph = pi / 16
        n = 8 - fromIntegral (fromEnum b)
        a = n * ph
    in if a <= -pi then a + 2 * pi else a

bearing_n :: Integral i => i -> Pt R -> Pt R -> i
bearing_n dv p q =
    let dv' = fromIntegral dv
        a = negate (pt_angle p q) + pi
    in round ((a * ((dv' / 2) / pi)) - (dv' / 4)) `mod` dv

-- | Bearing from 'Pt' /p/ to /q/.
--
-- > let f (x,y) = bearing (Pt 0 0) (Pt x y)
-- > map f [(0,1),(1,1),(1,0),(1,-1)] == [N,NE,E,SE]
-- > map f [(0,-1),(-1,-1),(-1,0),(-1,1)] == [S,SW,W,NW]
-- > map f [(1/2,1),(1,1/2),(1,-1/2),(1/2,-1)] == [NNE,ENE,ESE,SSE]
-- > map f [(1/4,1),(1,1/4),(1,-1/4),(1/4,-1)] == [NbE,EbN,EbS,SbE]
-- > map f [(-1/2,-1),(-1,-1/2),(-1,1/2),(-1/2,1)] == [SSW,WSW,WNW,NNW]
-- > map f [(-1/4,-1),(-1,-1/4),(-1,1.4),(-1/4,1)] == [SbW,WbS,NWbN,NbW]
bearing :: Pt R -> Pt R -> Bearing
bearing p = toEnum . bearing_n 32 p

-- | Bearing to nearest eight point compass bearing
--
-- > let f (x,y) = bearing_8 (Pt 0 0) (Pt x y)
-- > map f [(1/4,1),(1,1/4),(1,-1/4),(1/4,-1)] == [N,E,E,S]
bearing_8 :: Pt R -> Pt R -> Bearing
bearing_8 p = toEnum . (* 4) . bearing_n 8 p

-- | Predicate that is 'True' if bearings are opposite.
--
-- > bearing_opposite (NW,SE) == True
-- > map bearing_opposite (zip [N,E,S,W] [S,W,N,E]) == [True,True,True,True]
bearing_opposite :: (Bearing,Bearing) -> Bool
bearing_opposite (p, q) =
    let n = (fromEnum p - fromEnum q) `mod` 32
    in n == 16
