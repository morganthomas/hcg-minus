-- | RGB to greyscale conversion functions.
module Data.CG.Minus.Colour.Grey where

-- | (R,G,B) triple.
type RGB a = (a,a,a)

-- | Grey luminance.
type GREY a = a

-- | Apply binary function at triple, leftmost first.
f2_t3_left :: (t -> t -> t) -> (t,t,t) -> t
f2_t3_left f (p,q,r) = f (f p q) r

-- | Triple variant of 'min'.
min3 :: Ord t => (t,t,t) -> t
min3 = f2_t3_left min

-- | Triple variant of 'max'.
max3 :: Ord t => (t,t,t) -> t
max3 = f2_t3_left max

-- | Desaturation.
rgb_to_gs_lightness :: (Fractional a, Ord a) => RGB a -> GREY a
rgb_to_gs_lightness c = (min3 c + max3 c) / 2

-- | Simple average.
rgb_to_gs_average :: Fractional a => RGB a -> GREY a
rgb_to_gs_average (r,g,b) = (r + g + b) / 3

-- | Luminosity coefficents, ie. (R,G,B) multipliers.
type COEF a = (a,a,a)

rgb_to_gs_luminosity :: Num a => COEF a -> RGB a -> GREY a
rgb_to_gs_luminosity (rm,gm,bm) (r,g,b) = r * rm + g * gm + b * bm

luminosity_coef_rec_709 :: Fractional a => COEF a
luminosity_coef_rec_709 = (0.2126,0.7152,0.0722)

luminosity_coef_rec_601 :: Fractional a => COEF a
luminosity_coef_rec_601 = (0.299,0.587,0.114)

luminosity_coef_smpte_240m :: Fractional a => COEF a
luminosity_coef_smpte_240m = (0.212,0.701,0.087)

rgb_to_gs_decompose_min :: Ord a => RGB a -> GREY a
rgb_to_gs_decompose_min = f2_t3_left min

rgb_to_gs_decompose_max :: Ord a => RGB a -> GREY a
rgb_to_gs_decompose_max = f2_t3_left max

rgb_to_gs_r :: RGB t -> GREY t
rgb_to_gs_r (r,_,_) = r

rgb_to_gs_g :: RGB t -> GREY t
rgb_to_gs_g (_,g,_) = g

rgb_to_gs_b :: RGB t -> GREY t
rgb_to_gs_b (_,_,b) = b

-- | Requires R and G and B components to be equal.
rgb_to_gs_eq :: Eq a => RGB a -> Maybe (GREY a)
rgb_to_gs_eq (r,g,b) = if r == g && r == b then Just r else Nothing

-- | 'error'ing variant.
rgb_to_gs_eq' :: Eq a => RGB a -> GREY a
rgb_to_gs_eq' = maybe (error "rgb_to_gs_eq: not equal") id . rgb_to_gs_eq
