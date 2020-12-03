-- | Colour names, non-specific.
module Data.CG.Minus.Colour.Names where

import Data.Maybe {- base -}

import qualified Data.Colour as C {- colour -}

import Data.CG.Minus.Colour {- hcg-minus -}

hex_to_c :: Int -> C.Colour Double
hex_to_c = rgb24_to_c . hex_to_rgb24

colour_tbl :: [(String,C.Colour Double)]
colour_tbl =
    let f (_,nm,hex) = (nm,hex_to_c hex)
    in map f colour_names_table

-- > map colour_lookup_err ["Cerulean","Fern green","Candlelight yellow"]
colour_lookup_err :: String -> C.Colour Double
colour_lookup_err = fromMaybe (error "colour_lookup") . flip lookup colour_tbl

-- > map (\(_,_,n) -> hex_to_rgb24 n) colour_names_table
colour_names_table :: Num n => [(String,String,n)]
colour_names_table =
    [("red","Venetian red",0xC80815)
    ,("blue","Swedish azure blue",0x005B99)
    ,("orange","Safety orange",0xFF6600)
    ,("magenta","Dye magenta",0xCA1F7B)
    ,("magenta","Process magenta",0xFF0090)
    ,("yellow","Candlelight yellow",0xFCD116)
    ,("cyan","Cyan additive secondary",0x00FFFF)
    ,("cyan","Cyan subtractive primary",0x00B7EB)
    ,("green","Fern green",0x009246)
    ,("brown","Sepia brown",0x704214)
    ,("green","Verdigris",0x43B3AE)
    ,("cyan","Viridian",0x40826D)
    ,("cyan","Cerulean",0x007BA7)
    ]
