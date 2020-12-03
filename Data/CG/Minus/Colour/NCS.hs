-- | Natural Color System
module Data.CG.Minus.Colour.NCS where

import qualified Data.Colour as C {- colour -}

import Data.CG.Minus.Colour {- hcg-minus -}

ncs_table :: Num n => [(String,n,n,n,n)]
ncs_table =
    [("White",0xFFFFFF,255,255,255)
    ,("Black",0x000000,0,0,0)
    ,("Green",0x009F6B,0,159,107) -- 0x00A368
    ,("Red",0xC40233,196,2,51)
    ,("Yellow",0xFFD300,255,211,0)
    ,("Blue",0x0087BD,0,135,189) -- 0x0088BF
    ]

ncs_colour :: (Int,Int,Int) -> C.Colour Double
ncs_colour = rgb24_to_c

ncs_red,ncs_green,ncs_yellow,ncs_blue :: C.Colour Double
ncs_red = ncs_colour (196,2,51)
ncs_green = ncs_colour (0,159,107)
ncs_yellow = ncs_colour (255,211,0)
ncs_blue = ncs_colour (0,135,189)
