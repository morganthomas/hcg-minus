-- | Colour related functions
module Data.CG.Minus.Colour where

import Data.Bits {- base -}
import qualified Data.Colour as C {- colour -}
import qualified Data.Colour.SRGB as SRGB {- colour -}
import qualified Data.Colour.RGBSpace.HSL as HSL {- colour -}

-- | Grey 'Colour'.
mk_grey :: (Ord a,Floating a) => a -> C.Colour a
mk_grey x = SRGB.sRGB x x x

-- | Reduce 'Colour' to grey.  Constants are @0.3@, @0.59@ and @0.11@.
to_greyscale :: (Ord a,Floating a) => C.Colour a -> a
to_greyscale c =
    let (SRGB.RGB r g b) = SRGB.toSRGB c
    in r * 0.3 + g * 0.59 + b * 0.11

-- | 'mk_grey' '.' 'to_greyscale'.
to_greyscale_c :: (Ord a,Floating a) => C.Colour a -> C.Colour a
to_greyscale_c = mk_grey . to_greyscale

-- | Discard /alpha/ channel, if possible.
pureColour :: (Ord a, Fractional a) => C.AlphaColour a -> C.Colour a
pureColour c =
    let a = C.alphaChannel c
    in if a > 0
       then C.darken (recip a) (c `C.over` C.black)
       else error "hcg-: transparent has no pure colour"

-- * Tuples

-- > hex_to_rgb24 0xC80815 == (0xC8,0x08,0x15)
hex_to_rgb24 :: (Bits t, Num t) => t -> (t,t,t)
hex_to_rgb24 n = (shiftR (n .&. 0xFF0000) 16,shiftR (n .&. 0x00FF00) 8,n .&. 0x0000FF)

srgb_components :: SRGB.RGB t -> (t,t,t)
srgb_components c = (SRGB.channelRed c,SRGB.channelGreen c,SRGB.channelBlue c)

-- | 'C' to /(red,green,blue)/ tuple.
c_to_rgb :: (Ord t,Floating t) => C.Colour t -> (t,t,t)
c_to_rgb = srgb_components . SRGB.toSRGB

-- | Tuple to 'C', inverse of 'c_to_rgb'.
rgb_to_c :: (Ord t,Floating t) => (t,t,t) -> C.Colour t
rgb_to_c (r,g,b) = SRGB.sRGB r g b

-- > rgb24_to_c ((0xC8,0x08,0x15) :: (Int,Int,Int))
rgb24_to_c :: (Real r,Ord t,Floating t) => (r, r, r) -> C.Colour t
rgb24_to_c (r,g,b) = rgb_to_c (realToFrac r / 255,realToFrac g / 255,realToFrac b / 255)

c_to_hsl :: (Ord t,Floating t) => C.Colour t -> (t,t,t)
c_to_hsl c = HSL.hslView (SRGB.toSRGB c)

hsl_to_c :: (RealFrac t,Floating t) => (t,t,t) -> C.Colour t
hsl_to_c (h,s,l) = rgb_to_c (srgb_components (HSL.hsl h s l))

-- > rgb_to_hsl (1,0,0) == (0,1,0.5)
rgb_to_hsl :: (Ord t,Floating t) => (t,t,t) -> (t,t,t)
rgb_to_hsl = c_to_hsl . rgb_to_c

-- > hsl_to_rgb (0,1,0.5) == (1,0,0)
hsl_to_rgb :: (RealFrac t,Floating t) => (t,t,t) -> (t,t,t)
hsl_to_rgb = c_to_rgb . hsl_to_c

-- | Tuple to 'Ca', inverse of 'c_to_rgba'.
rgba_to_ca :: (Ord t,Floating t) => (t,t,t,t) -> C.AlphaColour t
rgba_to_ca (r,g,b,a) = rgb_to_c (r,g,b) `C.withOpacity` a

rgb_to_ca :: (Ord t,Floating t) => (t,t,t) -> C.AlphaColour t
rgb_to_ca (r,g,b) = rgb_to_c (r,g,b) `C.withOpacity` 1

-- | 'Ca' to /(red,green,blue,alpha)/ tuple
ca_to_rgba :: (Ord t,Floating t) => C.AlphaColour t -> (t,t,t,t)
ca_to_rgba x =
    let x' = SRGB.toSRGB (pureColour x)
    in (SRGB.channelRed x'
       ,SRGB.channelGreen x'
       ,SRGB.channelBlue x'
       ,C.alphaChannel x)

-- | Is the /alpha/ channel zero.
ca_is_transparent :: (Ord t, Num t) => C.AlphaColour t -> Bool
ca_is_transparent x = not (C.alphaChannel x > 0)
