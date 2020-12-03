-- | <http://afriggeri.github.io/RYB/> &
--   <http://web.siat.ac.cn/~baoquan/papers/InfoVis_Paint.pdf>
module Data.CG.Minus.Colour.RYB where

import Data.List {- base -}

import qualified Data.Fixed as F {- base -}

import qualified Data.Colour.SRGB as SRGB {- colour -}

type R = Double
type RYB = (R,R,R)
type RGB = (R,R,R)
type RGB_U8 = (Int,Int,Int)

t3_to_list :: (t,t,t) -> [t]
t3_to_list (p,q,r) = [p,q,r]

t3_from_list :: [t] -> (t, t, t)
t3_from_list l =
    case l of
      [p,q,r] -> (p,q,r)
      _ -> error "t3_from_list"

ryb_clr :: [(String,[R])]
ryb_clr =
    [("white",  [1,1,1])
    ,("red",    [1,0,0])
    ,("yellow", [1,1,0])
    ,("blue",   [0.163,0.373,0.6])
    ,("violet", [0.5,0,0.5])
    ,("green",  [0,0.66,0.2])
    ,("orange", [1,0.5,0])
    ,("black",  [0.2,0.094,0.0])]

ryb_clr_ix :: String -> Int -> R
ryb_clr_ix nm ix = maybe (error "ryb_clr_ix") (!! ix) (lookup nm ryb_clr)

unit_to_u8 :: R -> Int
unit_to_u8 = floor . (* 255.0)

t3_unit_to_u8 :: (R,R,R) -> (Int,Int,Int)
t3_unit_to_u8 (p,q,r) = (unit_to_u8 p,unit_to_u8 q,unit_to_u8 r)

rgb_to_u8 :: RGB -> RGB_U8
rgb_to_u8 = t3_unit_to_u8

-- > map (ryb_to_u8 . ryb_to_rgb) [(0,0,0),(1,1,1),(1,0,0),(0,1,0),(0,0,1)]
-- > ryb_to_rgb (0,1,1) == (0.0,0.66,0.2) -- green
-- > ryb_to_rgb (0,0.5,1) -- cyan = blue+green = (0,1,2)
ryb_to_rgb :: RYB -> RGB
ryb_to_rgb (r, y, b) =
    let f ix =
            ryb_clr_ix "white" ix  * (1-r) * (1 - b) * (1 - y) +
            ryb_clr_ix "red" ix    *    r  * (1 - b) * (1 - y) +
            ryb_clr_ix "blue" ix   * (1-r) *      b  * (1 - y) +
            ryb_clr_ix "violet" ix *    r  *      b  * (1 - y) +
            ryb_clr_ix "yellow" ix * (1-r) * (1 - b) *      y  +
            ryb_clr_ix "orange" ix *     r * (1 - b) *      y  +
            ryb_clr_ix "green" ix  * (1-r) *      b  *      y  +
            ryb_clr_ix "black" ix  *     r *      b  *      y
    in t3_from_list (map f [0..2])

-- > map (euclidian_distance [0,1,0]) [[1,1,0],[0,1,0],[1,1,1],[0.5,0.5,0.5]] == [1,0,2,0.75]
euclidian_distance :: Floating t => [t] -> [t] -> t
euclidian_distance p1 =
    let f x1 x2 = (x2 - x1) ** 2
    in sum . zipWith f p1

euclidian_distance_t3 :: Floating t => (t,t,t) -> (t,t,t) -> t
euclidian_distance_t3 (p1,p2,p3) (q1,q2,q3) =
    let f x1 x2 = (x2 - x1) ** 2
    in f p1 q1 + f p2 q2 + f p3 q3

euclidian_distance_t3_set :: Floating t => [(t,t,t)] -> (t,t,t) -> t
euclidian_distance_t3_set l x = sum (map (euclidian_distance_t3 x) l)

int_to_r :: Int -> R
int_to_r = fromIntegral

-- > map n_triples [8,27,64,125,216,343,512]
n_triples :: Int -> (R,R)
n_triples k =
    let fceil = int_to_r . ceiling
        b = fceil (int_to_r k ** (1/3))
        n = (b ** 3)
    in (b,n)

-- > gen_triples 8 == [(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]
gen_triples :: Int -> [RYB]
gen_triples k =
    let ffloor = int_to_r . floor
        (base,base_n) = n_triples k
        (%) = F.mod'
        f n = (ffloor (n / (base * base)) / (base - 1)
              ,ffloor ((n / base) % base) / (base - 1)
              ,ffloor (n % base) / (base - 1))
    in map f [0 .. base_n - 1]

-- > most_distant_set [(0,0,0)] (gen_triples 8)
most_distant_set :: (Ord t, Floating t) => [(t,t,t)] -> [(t,t,t)] -> ((t,t,t),[(t,t,t)])
most_distant_set x l =
    let d = map (euclidian_distance_t3_set x) l
        z = maximum d
    in case find ((== z) . fst) (zip d l) of
         Just (_,e) -> (e,delete e l)
         Nothing -> error "most_distant_set"

distance_step :: (Ord t, Floating t) => [(t,t,t)] -> [(t,t,t)] -> [(t,t,t)]
distance_step lhs rhs =
    case rhs of
      [] -> []
      _ -> let (e,rhs') = most_distant_set lhs rhs
           in e : distance_step (e:lhs) rhs'

-- > map (rgb_to_u8 . ryb_to_rgb) (distance_sort (gen_triples 8))
distance_sort :: (Ord t, Floating t) => [(t,t,t)] -> [(t,t,t)]
distance_sort l =
    case l of
      e:l' -> e : distance_step [e] l'
      _ -> error "distance_sort"

-- > ryb_colour_gen 8 == [(0,0,0),(1,1,1),(0,0,1),(1,1,0),(0,1,0),(1,0,1),(0,1,1),(1,0,0)]
ryb_colour_gen :: Int -> [RYB]
ryb_colour_gen = distance_sort . gen_triples

rgb_colour_gen :: Int -> [RGB]
rgb_colour_gen = map ryb_to_rgb . ryb_colour_gen

-- > rgb_u8_colour_gen 27
rgb_u8_colour_gen :: Int -> [RGB_U8]
rgb_u8_colour_gen = map rgb_to_u8 . rgb_colour_gen

colour_gen :: Int -> [SRGB.Colour R]
colour_gen = map (\(r,g,b) -> SRGB.sRGB r g b) . rgb_colour_gen
