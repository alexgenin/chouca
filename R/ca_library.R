# 
# 
# This file contains models 
# 

#' @title Library of probabilistic cellular automata 
#'
#' @description Get one of the PCA model included in \code{chouca}
#' 
#' @param model The model to return, as a string. See Details for the full list of models
#'   included with \code{chouca}.
#' 
#' @param parms The model parameters to use, as a named list. If unset, the model default
#'   parameters will be used. 
#' 
#' @param neighbors The number of neighbors to use in the cellular automaton (4 for 4-way 
#'   or von-Neumann neghborhood, or 8 for a Moore neighborhood). If unset, the model 
#'   default neighborhood will be used. 
#' 
#' @param wrap Whether the 2D grid should wrap around at the edges. 
#' 
#' @examples 
#' 
#' # Import a model, create an initial landscape and run it for ten iterations
#' forestgap_model <- ca_library("forestgap")
#' im <- generate_initmat(forestgap_model, c(0.5, 0.5), nr = 64, nc = 100) 
#' run_camodel(forestgap_model, im, niter = 10) 
#' 
#'@export
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
  if ( model == "musselbed" || model == "mussel-bed" || model == "mussel bed") { 
    if ( is.null(parms) ) { 
      parms <- list(d = 0.1, 
                    delta = 0.25, 
                    r = 0.4)
    }
    # There is 8 neighrbors by default in MB model (caspr implementation: 
    # https://github.com/fdschneider/caspr/blob/master/R/musselbed.R)
    neighbors <- ifelse(is.null(neighbors), 8, neighbors)
    
    mod <- camodel( 
      transition("EMPTY",   "MUSSEL",  ~ r * q["MUSSEL"]), 
      transition("DISTURB", "EMPTY",   ~ 1), 
      transition("MUSSEL",  "DISTURB", ~ delta + d * (q["DISTURB"] > 0)), 
      neighbors = neighbors, 
      wrap = wrap, 
      parms = parms, 
      all_states = c("MUSSEL", "EMPTY", "DISTURB")
    )
  }
  
# Kéfi's Arid vegetation model (2007), extended by Schneider (2016) 
# 
# 
# Kéfi, Sonia, Max Rietkerk, Concepción L. Alados, Yolanda Pueyo, Vasilios P. 
# Papanastasis, Ahmed ElAich, and Peter C. de Ruiter. 2007. “Spatial Vegetation 
# Patterns and Imminent Desertification in Mediterranean Arid Ecosystems.” 
# Nature 449 (7159): 213–17. doi: 10.1038/nature06111.
# 
# Schneider, Florian D., and Sonia Kefi. 2016. “Spatially Heterogeneous Pressure Raises
# Risk of Catastrophic Shifts.” Theoretical Ecology 9 (2): 207-17. 
# doi: 10.1007/s12080-015-0289-1.
# 
# Notation follows what is in Schneider et al. (2016)
  if ( model == "aridvege" || model == "arid vege" || model == "arid-vege" ) { 
    
    if ( is.null(parms) ) { 
    parms <- list(r = 0.01, 
                  f = 0.9, 
                  delta = 0.1, 
                  c = 0.2, 
                  d = 0.1, 
                  m0 = 0.05, 
                  g0 = 0.1, # g0 and b in the bistable range
                  b = 0.4)
    }
    
    wrap <- ifelse(is.null(wrap), TRUE, wrap)
    neighbors <- ifelse(is.null(neighbors), 4, neighbors) 
    
    mod <- camodel(
      transition("DEGR", "VEGE", 
                 ~ r + q["VEGE"] * f), 
      transition("EMPTY", "VEGE", 
                 ~ ( delta * p["VEGE"] + ( 1 - delta ) * 
                   q["VEGE"]) * ( b - c * p["VEGE"] )), 
      transition("EMPTY", "DEGR", 
                 ~ d), 
      transition("VEGE", "EMPTY", 
                 ~ m0 + g0 * ( 1 - q["VEGE"] )), 
      parms = parms, 
      wrap = wrap, 
      neighbors = neighbors, 
      all_states = c("DEGR", "EMPTY", "VEGE")
    )
    
  }
  
  
  # Genin's coral reef model (2022?) 
  # 
  if ( model == "coralreef" || model == "coral-reef" || model == "coral reef") { 
    if ( is.null(parms) ) { 
      parms <- list(r_a = 1.792032617, 
                    l_a = 0.016357009, 
                    i_a = 0.001792033, 
                    r_c = 0.001308055, 
                    d_0 = 1, 
                    m_a = 0.011323493, 
                    h_u = 2.000000000, 
                    g = 0.178133479, 
                    theta_b = 1.111260092, 
                    theta_c = 1.051433267, 
                    m_c     = 0.000100000)
    }
    
    # Default neighbors = 4
    neighbors <- ifelse(is.null(neighbors), 4, neighbors)
    wrap <- ifelse(is.null(wrap), TRUE, wrap)
    
    mod <- camodel(
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
  
  
  # Conway's Game of life 
  # https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life
  if ( model == "gameoflife" || model == "game-of-life" ) { 
    
    neighbors <- ifelse(is.null(neighbors), 8, neighbors)
    if ( ! is.null(parms) ) { 
      warning("Game of life model takes no parameters, parms will be ignored")
    }
    
    mod <- camodel( 
      transition("LIVE", "DEAD", ~ q["LIVE"] < (2/8) | q["LIVE"] > (3/8)), 
      transition("DEAD", "LIVE", ~ q["LIVE"] == (3/8)), 
      wrap = wrap, 
      neighbors = neighbors, 
      all_states = c("DEAD", "LIVE")
    )
  }
  
  
  # Rock-Paper-Scissor model (not sure where it comes from)
  # More than two neighbors that beat it -> switch
  # https://www.youtube.com/watch?v=TvZI6Xc0J1Y
  if ( model == "rock-paper-scissor" || model == "rockpaperscissor" || 
       model == "rpc") { 
    
    neighbors <- ifelse(is.null(neighbors), 8, neighbors)
    
    if ( is.null(parms) ) { 
      parms <- list(prob = 1)
    }
    
    mod <- camodel(
      transition(from = "r", to = "p", ~ prob * ( q["p"] > (1/8)*2) ), 
      transition(from = "p", to = "c", ~ prob * ( q["c"] > (1/8)*2) ), 
      transition(from = "c", to = "r", ~ prob * ( q["r"] > (1/8)*2) ), 
      parms = parms, 
      wrap = wrap, 
      neighbors = neighbors
    )
    
  }
  
  
  if ( is.null(mod) ) { 
    all_models <- c("forestgap", "musselbed", "gameoflife", "rockpaperscissor")
    stop(paste0(model, " is an unknown model. Available models:", 
                paste(" - ", all_models, collapse = "\n")))
  }
  
  return(mod)
}

