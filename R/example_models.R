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
forestgap <- function(parms = list(d = 0.125, 
                                   delta = 0.5, 
                                   alpha = 0.2)) { 
  camodel( 
    transition(from = "TREE", 
               to   = "EMPTY", 
               prob = ~ d + delta * q["EMPTY"] ), 
    transition(from = "EMPTY", 
               to   = "TREE", 
               prob = ~ alpha), 
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


