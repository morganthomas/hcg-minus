-- | Arrows
module Data.CG.Minus.Arrow where

import Data.CG.Minus

-- | Given the arrow body 'Ln' and the arrow length and arrow angle
-- (in radians) 'R' calculate the 'Pt' of each arrow tip.
--
-- > arrow_coord (Ln (Pt 0 0) (Pt 1 1)) 0.1 (pi/9)
arrow_coord :: Ln R -> R -> R -> (Pt R,Pt R)
arrow_coord l n a =
    let ((x0,y0),(x1,y1)) = ln_elem l
        a' = atan2 (y1 - y0) (x1 - x0) + pi
        x2 = x1 + n * cos (a' - a)
        y2 = y1 + n * sin (a' - a)
        x3 = x1 + n * cos (a' + a)
        y3 = y1 + n * sin (a' + a)
    in (Pt x2 y2,Pt x3 y3)
