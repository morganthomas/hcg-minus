module Data.CG.Minus.Colour.VGA where

import Data.Bits {- base -}
import Data.Function {- base -}
import Data.List {- base -}

type T3 t = (t,t,t)

t3_map :: (t -> u) -> T3 t -> T3 u
t3_map f (p,q,r) = (f p,f q,f r)

type T6 t = (t,t,t,t,t,t)

t6_from_list :: [t] -> T6 t
t6_from_list l = case l of {[a,b,c,d,e,f] -> (a,b,c,d,e,f);_ -> error "t6_from_list"}

-- | Table in binary form (r,g,b,R,G,B), ie. {0,1}.
ega_color_table_t6 :: [T6 Int]
ega_color_table_t6 =
    let f i = map (fromEnum . testBit i) [5::Int,4 .. 0]
    in map (t6_from_list . f) [0::Int .. 63]

t6_one_bit_to_t3_two_bit :: T6 Int -> T3 Int
t6_one_bit_to_t3_two_bit (r,g,b,r',g',b') = (r + 2 * r',g + 2 * g',b + 2 * b')

-- | Table in ternary form (r,g,b), ie. {0,1,2}.
--
-- > length ega_color_table == 64
-- > map (ega_color_table !!) [31,32] == [(2,3,3),(1,0,0)]
-- > let f = t3_map (* 0x55) in map f ega_color_table
ega_color_table :: [T3 Int]
ega_color_table = map t6_one_bit_to_t3_two_bit ega_color_table_t6

-- | Indices for the default 16-colour table.
--
-- > map (ega_color_table !!) ega_default_16_indices == ega_default_16
ega_default_16_indices :: [Int]
ega_default_16_indices = [0 .. 5] ++ [20,7] ++ [56 .. 63]

-- | (index,name,two-bit-rgb)
cga_16_color_palette :: [(Int,String,T3 Int)]
cga_16_color_palette =
    [(0,"black",(0,0,0)),(8,"gray",(1,1,1))
    ,(1,"blue",(0,0,2)),(9,"light blue",(1,1,3))
    ,(2,"green",(0,2,0)),(10,"light green",(1,3,1))
    ,(3,"cyan",(0,2,2)),(11,"light cyan",(1,3,3))
    ,(4,"red",(2,0,0)),(12,"light red",(3,1,1))
    ,(5,"magenta",(2,0,2)),(13,"light magenta",(3,1,3))
    ,(6,"brown",(2,1,0)),(14,"yellow",(3,3,1))
    ,(7,"light gray",(2,2,2)),(15,"white (high intensity)",(3,3,3))]

-- > let f = t3_map (* 0x55) in map f cga_16_color_table
cga_16_color_table :: [T3 Int]
cga_16_color_table =
    let sq = map (\(i,_,c) -> (i,c)) cga_16_color_palette
    in map snd (sortBy (compare `on` fst) sq)

-- | Synonym for 'cga_16_color_table'.
ega_default_16 :: [T3 Int]
ega_default_16 = cga_16_color_table
