-- | Planck radiation equation.
module Data.CG.Minus.Colour.Planck where

import qualified Data.Colour as C {- colour -}
import qualified Data.CG.Minus.Colour as CG {- hcg-minus -}

-- | Given wavelength (in microns) and temperature (in degrees Kelvin)
-- solve Planck's radiation equation.
--
-- > planck_rad_eq 0.7 2600 == 8.22656629154115e7
planck_rad_eq :: Floating a => a -> a -> a
planck_rad_eq l t =
    let k0 = 3.7403e10
        k1 = -5.0
        k2 = 2.7182818284590452354
        k3 = 14384.0
        pow = (**)
        n0 = k0 * pow l k1
        n1 = pow k2 (k3 / (l * t)) - 1.0
    in n0 / n1

-- | Return the color of a black body emitting light at a given
-- temperature.  The Planck radiation equation is solved directly for
-- the @R@, @G@, and @B@ wavelengths defined for the CIE 1931 Standard
-- Colorimetric Observer.  The colour temperature is specified in
-- degrees Kelvin.  Typical constraints for star temperatures are @>=@
-- 2600@K@ (/S Cephei, R Andromedae/) and @<=@ 28,000@K@ (/Spica/).
--
-- > let h (r,g,b) = let f = floor . (*) 255 in (f r,f g,f b)
-- > in map (h . k_to_rgb) [2600,28000] == [(255,95,22),(49,118,254)]
k_to_rgb :: (Floating t, Ord t) => t -> (t,t,t)
k_to_rgb k =
    let r = planck_rad_eq 0.7000 k
        g = planck_rad_eq 0.5461 k
        b = planck_rad_eq 0.4358 k
        s = 1.0 / max r (max g b)
    in (r * s,g * s,b * s)

-- | 'toC' '.' 'k_to_rgb'.
k_to_colour :: Double -> C.Colour Double
k_to_colour = CG.rgb_to_c . k_to_rgb
