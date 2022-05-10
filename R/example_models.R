# 
# 
# This file contains models 
# 

# Kubo's forest model (1996)
# 
# See also Génin et al. 2018 for model definition 
# 
# Kubo, Takuya, Yoh Iwasa, and Naoki Furumoto. 1996. “Forest Spatial Dynamics with Gap
# Expansion: Total Gap Area and Gap Size Distribution.” Journal of Theoretical Biology 
# 180 (3): 229–46.
# 
#'@export
forestgap <- function(parms = list(d = 0.125, 
                                   delta = 0.5, 
                                   alpha = 0.2)) { 
  camodel( 
    transition(from = "TREE", 
               to   = "EMPTY", 
               prob = ~ d + delta * q["EMPTY"] ), 
    transition(from = "EMPTY", 
               to   = "TREE", 
               prob = ~ alpha * p["TREE"]), 
    parms = parms, 
    all_states = c("EMPTY", "TREE")
  )
}

# Guichard's Mussel Bed model (2003)
# 
# NOTE: this is not exactly the original model as the effect of having disturbed 
# neighbors is proportional to the number of neighbors in that state, instead of being 
# the max as soon as one neighbors is disturbed. 
# 
# See also Génin et al. 2018 for model definition 
# 
# Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge, B.A. (2003). Mussel
# disturbance dynamics: signatures of oceanographic forcing from local interactions. The
# American Naturalist, 161, 889–904. doi: 10.1086/375300
# 
#'@export
musselbed <- function(parms = list(d = 0.1, 
                                   delta = 0.25, 
                                   alpha = 0.4)) { 
  camodel( 
    transition("EMPTY",   "MUSSEL",  ~ alpha * q["MUSSEL"]), 
    transition("DISTURB", "EMPTY",   ~ 1), 
    transition("MUSSEL",  "DISTURB", ~ d + delta * q["DISTURB"]), 
    parms = parms, 
    all_states = c("MUSSEL", "EMPTY", "DISTURB")
  )
}



# Kéfi's vegetation model 
# 
# Notation taken from Schneider and Kéfi 2016
# 
#'@export
aridvege <- function(parms = list(r = 0.01, 
                                  f = 0.9, 
                                  delta = 0.1, 
                                  c = 0.2, 
                                  d = 0.1, 
                                  m0 = 0.05, 
                                  g0 = 0.1,   # g0 and b in the bistable range
                                  b = 0.4)) { # 
  camodel(
    transition("DEGR", "VEGE", 
               ~ r + q["VEGE"] * f), 
    transition("EMPTY", "VEGE", 
               ~ ( delta * p["VEGE"] + ( 1 - delta ) * q["VEGE"]) * ( b - c * p["VEGE"] )), 
    transition("EMPTY", "DEGR", 
               ~ d), 
    transition("VEGE", "EMPTY", 
               ~ m0 + g0 * ( 1 - q["VEGE"] )), 
    parms = parms, 
    all_states = c("DEGR", "EMPTY", "VEGE")
  )
}

