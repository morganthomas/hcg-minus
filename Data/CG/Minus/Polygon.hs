module Data.CG.Minus.Polygon where

import Data.CG.Minus.Core
import Data.CG.Minus.Types

-- | Polygon.
type Polygon = [Pt R]

-- | Degree is the number of points.
--
-- > let tri = [Pt 0 0,Pt 1 0,Pt 0 1]
-- > polygon_degree tri == 3
polygon_degree :: Polygon -> Int
polygon_degree = length

-- | List of /k/ edges of /k/ polygon.
--
-- > polygon_edges tri == zipWith Ln tri (tail (cycle tri))
polygon_edges :: Polygon -> [Ln R]
polygon_edges =
    let adj_cyc l = zip l (tail (cycle l))
    in map (uncurry Ln) . adj_cyc

-- | Index /n/th point.
polygon_pt :: Polygon -> Int -> Pt R
polygon_pt p k = p !! k

-- | Index /n/th edge.
polygon_edge :: Polygon -> Int -> Ln R
polygon_edge p k = polygon_edges p !! k

polygon_align_ln :: Polygon -> Int -> (Int,Bool) -> Polygon
polygon_align_ln p k0 (k1,k1_rev) =
    let l0 = polygon_edge p k0
        l1 = (if k1_rev then ln_reverse else id) (polygon_edge p k1)
        (tr,(c,r)) = ln_align l0 l1
        f = pt_rotate_about r c . pt_translate (pt_to_vc tr)
    in map f p

polygon_reflect_ln_md :: Ln R -> Polygon -> Polygon
polygon_reflect_ln_md ln = map (pt_ln_reflect_md ln)

polygon_reflect_xy :: (R,R) -> Polygon -> Polygon
polygon_reflect_xy xy = map (pt_reflect_xy xy)

polygon_edge_reflect :: Polygon -> Int -> Polygon
polygon_edge_reflect p k = polygon_reflect_ln_md (polygon_edge p k) p

polygon_pt_reflect :: Polygon -> Int -> Polygon
polygon_pt_reflect p k = polygon_reflect_xy (pt_xy (polygon_pt p k)) p

polygon_reflect_x_right :: [Pt R] -> [Pt R]
polygon_reflect_x_right l =
    let m = maximum (map pt_x l)
    in map (pt_reflect_x m) l

polygon_reflect_y_up :: [Pt R] -> [Pt R]
polygon_reflect_y_up l =
    let m = maximum (map pt_y l)
    in map (pt_reflect_y m) l
