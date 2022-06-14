# 
# 
# This file contains models 
# 


ca_library <- function(model, 
                       parms = NULL, 
                       neighbors = NULL, 
                       wrap = TRUE) { 
  mod <- NULL 
  
  if ( length(model) != 1 ) { 
    stop("Bad model name specified")
  }
  
# Kubo's forest model (1996)
# 
# See also Génin et al. 2018 for model definition 
# 
# Kubo, Takuya, Yoh Iwasa, and Naoki Furumoto. 1996. “Forest Spatial Dynamics with Gap
# Expansion: Total Gap Area and Gap Size Distribution.” Journal of Theoretical Biology 
# 180 (3): 229–46.
# 
  if ( model == "forestgap" || model == "forest-gap") { 
    if ( is.null(parms) ) { 
      parms <- list(d = 0.125, 
                    delta = 0.5, 
                    alpha = 0.2)
    }
    #TODO: what is the default number of neighbors of the forestgap model ?
    neighbors <- ifelse(is.null(neighbors), 4, neighbors)
    
    mod <- camodel(
      transition(from = "TREE", 
                 to   = "EMPTY", 
                 prob = ~ d + delta * q["EMPTY"] ), 
      transition(from = "EMPTY", 
                 to   = "TREE", 
                 prob = ~ alpha * p["TREE"]), 
      neighbors = neighbors, 
      wrap = wrap, 
      parms = parms, 
      all_states = c("EMPTY", "TREE")
    )
  
  }
  
# Guichard's Mussel Bed model (2003)
# 
# See also Génin et al. 2018 for model definition 
# 
# Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge, B.A. (2003). Mussel
# disturbance dynamics: signatures of oceanographic forcing from local interactions. The
# American Naturalist, 161, 889–904. doi: 10.1086/375300
# 
  if ( model == "musselbed" || model == "mussel-bed" ) { 
    if ( is.null(parms) ) { 
      parms <- list(d = 0.1, 
                    delta = 0.25, 
                    alpha = 0.4)
    }
    #TODO: what is the default number of neighbors of the musselbed model ?
    neighbors <- ifelse(is.null(neighbors), 4, neighbors)
    
    mod <- camodel( 
      transition("EMPTY",   "MUSSEL",  ~ alpha * q["MUSSEL"]), 
      transition("DISTURB", "EMPTY",   ~ 1), 
      transition("MUSSEL",  "DISTURB", ~ d + delta * sign(q["DISTURB"])), 
      neighbors = neighbors, 
      wrap = wrap, 
      parms = parms, 
      all_states = c("MUSSEL", "EMPTY", "DISTURB")
    )
  }
  
  
  
  # Conway's Game of life 
  # https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life
  if ( model == "gameoflife" || model == "game-of-life" ) { 
    
    neighbors <- ifelse(is.null(neighbors), 8, neighbors)
    
    mod <- camodel( 
      transition("LIVE", "DEAD", ~ q["LIVE"] < (2/nb) | q["LIVE"] > (3/nb)), 
      transition("DEAD", "LIVE", ~ q["LIVE"] == (3/nb)), 
      wrap = wrap, 
      parms = parms, 
      parms = list(nb = 8), 
      all_states = c("DEAD", "LIVE")
    )
  }
  
  
  
  
  if ( is.null(mod) ) { 
    all_models <- c("forestgap", "musselbed", "gameoflife")
    stop(paste0(model, " is an unknown model. Available models:", 
                paste(all_models, collapse = " ")))
  }
  
  return(mod)
}



# Genin's coral reef model (2022?) 
# 
#'@export
coralreef <- function(parms = list(r_a = 1.79, 
                                   l_a = 0.02, 
                                   i_a = 1.79e-3, 
                                   r_c = 1.31e-4, 
                                   d_0 = 4, 
                                   m_a = 0.011, 
                                   h_u = 1, 
                                   g = 0.18, 
                                   theta_b = 1.11, 
                                   theta_c = 1.05, 
                                   m_c     = 10^(-3.5)), 
                      neighbors = 4, 
                      wrap = TRUE) { 
  
  camodel(
    transition(from = "BARE", to = "CORAL", 
               ~ r_c * ( 1 + d_0 * q["CORAL"] )), 
    transition(from = "CORAL", to = "BARE", 
               ~ m_c), 
    transition(from = "BARE", to = "ALGAE", 
               ~ i_a + r_a * p["ALGAE"] + l_a * q["ALGAE"]), 
    transition(from = "ALGAE", to = "BARE", 
               ~ m_a + h_u * g * ( theta_b * q["BARE"] + theta_c * q["CORAL"])), 
    neighbors = neighbors, 
    wrap = wrap, 
    parms = parms, 
    all_states = c("BARE", "ALGAE", "CORAL")
  )
}

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
                                   alpha = 0.2), 
                      neighbors = 4, 
                      wrap = TRUE) { 
  message("forestgap is deprecated, use ca_library(model = 'forestgap')")
  
  camodel( 
    transition(from = "TREE", 
               to   = "EMPTY", 
               prob = ~ d + delta * q["EMPTY"] ), 
    transition(from = "EMPTY", 
               to   = "TREE", 
               prob = ~ alpha * p["TREE"]), 
    neighbors = neighbors, 
    wrap = wrap, 
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

