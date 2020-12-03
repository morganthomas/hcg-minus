-- | Very simple picture model.
module Data.CG.Minus.Picture where

import Data.List {- base -}
import Data.Maybe {- base -}

import Data.CG.Minus.Types {- hcg-minus -}
import qualified Data.CG.Minus as CG {- hcg-minus -}

type Line_Width = R

type Dash = ([R],R)

data Pen = Pen Line_Width Ca Dash
           deriving (Eq,Show)

data Mark = Line Pen (Ln R)
          | Polygon (Either Pen Ca) [Pt R]
          | Circle (Either Pen Ca) (Pt R,R)
          | Dot Ca (Pt R,R)
            deriving (Eq,Show)

type Picture = [Mark]

no_dash :: Dash
no_dash = ([],0)

line_seq :: Pen -> [Pt R] -> [Mark]
line_seq pen =
    let adj l = zip l (tail l)
    in map (Line pen . uncurry Ln) . adj

polygon_l :: Pen -> [Pt R] -> Mark
polygon_l pen = Polygon (Left pen)

polygon_f :: Ca -> [Pt R] -> Mark
polygon_f clr = Polygon (Right clr)

circle_l :: Pen -> (Pt R,R) -> Mark
circle_l pen = Circle (Left pen)

circle_f :: Ca -> (Pt R,R) -> Mark
circle_f clr = Circle (Right clr)

mark_wn :: Mark -> Wn R
mark_wn m =
    case m of
      Line _ ln -> CG.ln_wn ln
      Polygon _ p -> CG.pts_window p
      Circle _ (c,r) -> CG.wn_square c r
      Dot _ (c,r) -> CG.wn_square c r

mark_normal :: Mark -> Mark
mark_normal m =
    case m of
      Line p ln -> Line p (CG.ln_sort ln)
      Polygon _ _ -> m -- should ensure CCW
      Circle _ _ -> m
      Dot _ _ -> m

mark_pt_set :: Mark -> [Pt R]
mark_pt_set m =
    case m of
      Line _ (Ln p q) -> [p,q]
      Polygon _ p -> p
      Circle _ (p,_) -> [p]
      Dot _ (p,_) -> [p]

mark_ln :: Mark -> Maybe (Ln R)
mark_ln m =
    case m of
      Line _ l -> Just l
      _ -> Nothing

mark_circle :: Mark -> Maybe (Pt R,R)
mark_circle m =
    case m of
      Circle _ c -> Just c
      _ -> Nothing

picture_pt_set :: Picture -> [Pt R]
picture_pt_set = concatMap mark_pt_set

picture_ln_set :: Picture -> [Ln R]
picture_ln_set = mapMaybe mark_ln

picture_ln_intersections :: Picture -> [Pt R]
picture_ln_intersections p =
    let l = picture_ln_set p
    in catMaybes [CG.ln_intersection l0 l1 | l0 <- l, l1 <- l, l0 /= l1]

picture_ln_circle_intersections :: Picture -> [Pt R]
picture_ln_circle_intersections p =
    let l_set = picture_ln_set p
        c_set = mapMaybe mark_circle p
    in concat [CG.ln_circle_intersection_set l c | l <- l_set, c <- c_set]

picture_normalise :: Picture -> Picture
picture_normalise = nub . map mark_normal

picture_wn :: Picture -> Wn R
picture_wn = foldl1 CG.wn_join . map mark_wn

