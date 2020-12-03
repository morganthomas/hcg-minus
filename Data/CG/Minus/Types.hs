-- | Types
module Data.CG.Minus.Types where

import qualified Data.Colour as C {- colour -}

-- | Two-dimensional point.
--
-- Pt are 'Num', pointwise, ie:
--
-- > Pt 1 2 + Pt 3 4 == Pt 4 6
-- > Pt 1 2 * Pt 3 4 == Pt 3 8
-- > negate (Pt 0 1) == Pt 0 (-1)
-- > abs (Pt (-1) 1) == Pt 1 1
-- > signum (Pt (-1/2) (1/2)) == Pt (-1) 1
data Pt a = Pt {pt_x :: a,pt_y :: a} deriving (Eq,Ord,Show)

-- | Unary operator at 'Pt', ie. basis for 'Num' instances.
pt_uop :: (a -> b) -> Pt a -> Pt b
pt_uop f (Pt x y) = Pt (f x) (f y)

-- | Binary operator at 'Pt', ie. basis for 'Num' instances.
pt_binop :: (a -> b -> c) -> Pt a -> Pt b -> Pt c
pt_binop f (Pt x1 y1) (Pt x2 y2) = Pt (x1 `f` x2) (y1 `f` y2)

instance Num a => Num (Pt a) where
    (+) = pt_binop (+)
    (-) = pt_binop (-)
    (*) = pt_binop (*)
    negate = pt_uop negate
    abs = pt_uop abs
    signum = pt_uop signum
    fromInteger n = let n' = fromInteger n in Pt n' n'

data P3 a = P3 {p3_x :: a,p3_y :: a,p3_z :: a} deriving (Eq,Ord,Show)

-- | Unary operator at 'P3', ie. basis for 'Num' instances.
p3_uop :: (a -> b) -> P3 a -> P3 b
p3_uop f (P3 x y z) = P3 (f x) (f y) (f z)

-- | Binary operator at 'P3', ie. basis for 'Num' instances.
p3_binop :: (a -> b -> c) -> P3 a -> P3 b -> P3 c
p3_binop f (P3 x1 y1 z1) (P3 x2 y2 z2) = P3 (x1 `f` x2) (y1 `f` y2) (z1 `f` z2)

instance Num a => Num (P3 a) where
    (+) = p3_binop (+)
    (-) = p3_binop (-)
    (*) = p3_binop (*)
    negate = p3_uop negate
    abs = p3_uop abs
    signum = p3_uop signum
    fromInteger n = let n' = fromInteger n in P3 n' n' n'

-- | Two-dimensional vector.  Vector are 'Num' in the same manner as 'Pt'.
data Vc a = Vc {vc_x :: a,vc_y :: a} deriving (Eq,Ord,Show)

-- | Unary operator at 'Vc', ie. basis for 'Num' instances.
vc_uop :: (a -> b) -> Vc a -> Vc b
vc_uop f (Vc x y) = Vc (f x) (f y)

-- | Binary operator at 'Vc', ie. basis for 'Num' instances.
vc_binop :: (a -> b -> c) -> Vc a -> Vc b -> Vc c
vc_binop f (Vc x1 y1) (Vc x2 y2) = Vc (x1 `f` x2) (y1 `f` y2)

instance Num a => Num (Vc a) where
    (+) = vc_binop (+)
    (-) = vc_binop (-)
    (*) = vc_binop (*)
    negate = vc_uop negate
    abs = vc_uop abs
    signum = vc_uop signum
    fromInteger n = let n' = fromInteger n in Vc n' n'

-- | Two-dimensional line.
--
-- > ln_start (Ln (Pt 0 0) (Pt 1 1)) == Pt 0 0
-- > ln_end (Ln (Pt 0 0) (Pt 1 1)) == Pt 1 1
data Ln a = Ln {ln_start :: Pt a,ln_end :: Pt a} deriving (Eq,Ord,Show)

-- | Line segments.
data Ls a = Ls {ls_elem :: [Pt a]} deriving (Eq,Show)

-- | Window, given by a /lower left/ 'Pt' and an /extent/ 'Vc'.
data Wn a = Wn {wn_ll :: Pt a,wn_ex :: Vc a} deriving (Eq,Show)

-- | Real number, synonym for 'Double'.
type R = Double

-- | Transformation matrix data type (Matrix a1 a2 b1 b2 c1 c2)
--  x' =  a1 * x + b1 * y + c1
--  y' =  a2 * x + b2 * y + c2
data Matrix n = Matrix n n n n n n deriving (Eq,Show)

-- | See 'pt_transform' for typed variant.
matrix_apply_raw :: Num t => (t,t,t,t,t,t) -> (t,t) -> (t,t)
matrix_apply_raw (a1,a2,b1,b2,c1,c2) (x,y) =
    let x' = x * a1 + y * b1 + c1
        y' = x * a2 + y * b2 + c2
    in (x',y')

-- | Enumeration of 'Matrix' indices.
data Matrix_Index = I0 | I1 | I2

mx_row :: Num n => Matrix n -> Matrix_Index -> (n,n,n)
mx_row (Matrix a b c d e f) i =
    case i of
      I0 -> (a,b,0)
      I1 -> (c,d,0)
      I2 -> (e,f,1)

mx_col :: Num n => Matrix n -> Matrix_Index -> (n,n,n)
mx_col (Matrix a b c d e f) i =
    case i of
      I0 -> (a,c,e)
      I1 -> (b,d,f)
      I2 -> (0,0,1)

mx_multiply :: Num n => Matrix n -> Matrix n -> Matrix n
mx_multiply a b =
    let f i j = let (r1,r2,r3) = mx_row a i
                    (c1,c2,c3) = mx_col b j
                in r1 * c1 + r2 * c2 + r3 * c3
    in Matrix (f I0 I0) (f I0 I1) (f I1 I0) (f I1 I1) (f I2 I0) (f I2 I1)

-- | Pointwise unary operator.
mx_uop :: (n -> n) -> Matrix n -> Matrix n
mx_uop g (Matrix a b c d e f) =
    Matrix (g a) (g b) (g c) (g d) (g e) (g f)

-- | Pointwise binary operator.
mx_binop :: (n -> n -> n) -> Matrix n -> Matrix n -> Matrix n
mx_binop g (Matrix a b c d e f) (Matrix a' b' c' d' e' f') =
    Matrix (g a a') (g b b') (g c c') (g d d') (g e e') (g f f')

instance Num n => Num (Matrix n) where
    (*) = mx_multiply
    (+) = mx_binop (+)
    (-) = mx_binop (-)
    abs = mx_uop abs
    signum = mx_uop signum
    fromInteger n = let n' = fromInteger n
                    in Matrix n' 0 0 n' 0 0

-- | Opaque colour.
type C = C.Colour R

-- | Colour with /alpha/ channel.
type Ca = C.AlphaColour R

