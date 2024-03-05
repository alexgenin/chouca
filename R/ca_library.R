#
#
# This file contains models
#

#' @title Library of stochastic cellular automata
#'
#' @description Get one of the SCA models included in \code{chouca}
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
#' @param wrap Whether the 2D grid should wrap around at the edges. Default it to wrap
#'   around the edges of the landscape.
#'
#' @details
#'
#'   This function gives access to different stochastic cellular automata models. You
#'     can provide a named list of parameters, or set the number of neighbor or wrapping
#'     options, but default are chosen if left unspecified. This function provides
#'     the following models (the string represents the name of the model, as passed
#'     using the 'model' argument):
#'
#'   \enumerate{
#'     \item \code{"forestgap"} Kubo's forest gap model (1996), which describes how
#'       gaps form in a forest and expand through disturbances.
#'
#'     \item \code{"musselbed"} A model of mussel beds, in which disturbance by waves
#'       occurs at the edge of mussel patches (Guichard et al. 2003)
#'
#'     \item \code{"aridvege"} A model of arid vegetation, in which facilitation between
#'       neighboring plants occur, along with grazing. The original model is to be
#'       found in Kéfi et al. (2007), with extensions in Schneider et al. (2016)
#'
#'     \item \code{"aridvege-danet"} An extension of the previous model to two species
#'       with assymetric facilitation (Danet et al. 2021)
#'
#'     \item \code{"coralreef"} A model of coral reef with local feedbacks of corals and
#'       macroalgae (Génin, in prep)
#'
#'     \item \code{"gameoflife"} The famous Game of Life by Conway, a deterministic
#'       model.
#'
#'     \item \code{"rockpaperscissor"} A rock-paper-scissor model with three states, in
#'       which a cell changes state depending on its neighbors according the game rules
#'       (e.g. "rock beats scissors"). This deterministic model produces nice spirals.
#'   }
#' 
#' @returns 
#'   
#NOTE: Also update the same section in ca_library() when changing this
#' 
#' This function returns a \code{\link{list}} object with class \code{ca_model}, with 
#' the following named components. Please note that most are for internal use and may 
#' change with package updates. 
#' 
#' \describe{
#'   \item{\code{transitions}}{the list of transitions of the model, as returned 
#'       by \code{\link{transition}} }
#' 
#'   \item{\code{nstates}}{the number of states of the model}
#' 
#'   \item{\code{parms}}{the parameter values used for the model}
#' 
#'   \item{\code{beta_0},\code{beta_q}, \code{beta_pp}, \code{beta_pq}, \code{beta_qq}}{ 
#'     internal tables used to represent probabilities of transitions when running
#'     simulations, these tables are for internal use and probably not interesting for 
#'     end users, but more information is provided in the package source code}
#'   
#'   \item{\code{wrap}}{Whether the model uses a toric space that wraps around the edge}
#' 
#'   \item{\code{neighbors}}{The type of neighborhood (4 or 8)}
#'   
#'   \item{\code{epsilon}}{The \code{epsilon} values used in the model definition, below 
#'     which transition probabilities are assumed to be zero}
#'   
#'   \item{\code{xpoints}}{(for internal use only) The number of values used to 
#'      represent the proportion of neighbors of a cell in each state}
#'   
#'   \item{\code{max_error}, \code{max_rel_error}}{vector of numeric values containing 
#'     the maximum error and maximum relative error on each transition probability}
#' 
#'   \item{\code{fixed_neighborhood}}{flag equal to \code{TRUE} when cells have
#'     a fixed number of neighbors}
#' 
#' }
#' 
#' @references
#'
#'  Danet, Alain, Florian Dirk Schneider, Fabien Anthelme, and Sonia Kéfi. 2021.
#'  "Indirect Facilitation Drives Species Composition and Stability in Drylands."
#'  Theoretical Ecology 14 (2): 189–203. \doi{10.1007/s12080-020-00489-0}.
#'
#'  Genin, A., S. A. Navarrete, A. Garcia-Mayor, and E. A. Wieters. in
#'  press (2023). Emergent spatial patterns can indicate upcoming regime
#'  shifts in a realistic model of coral community. The American Naturalist.
#'
#'  Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge, B.A. (2003). Mussel
#'  disturbance dynamics: signatures of oceanographic forcing from local interactions.
#'  The American Naturalist, 161, 889–904. \doi{10.1086/375300}
#'
#'  Kefi, Sonia, Max Rietkerk, Concepción L. Alados, Yolanda Pueyo, Vasilios P.
#'  Papanastasis, Ahmed ElAich, and Peter C. de Ruiter. 2007. "Spatial Vegetation
#'  Patterns and Imminent Desertification in Mediterranean Arid Ecosystems."
#'  Nature 449 (7159): 213–17. \doi{10.1038/nature06111}.
#'
#'  Kubo, Takuya, Yoh Iwasa, and Naoki Furumoto. 1996. "Forest Spatial Dynamics with Gap
#'  Expansion: Total Gap Area and Gap Size Distribution." Journal of Theoretical Biology
#'  180 (3): 229–46.
#'
#'  Schneider, Florian D., and Sonia Kefi. 2016. "Spatially Heterogeneous Pressure Raises
#'  Risk of Catastrophic Shifts." Theoretical Ecology 9 (2): 207-17.
#'  \doi{10.1007/s12080-015-0289-1}.
#'
#' @examples
#'
#' # Import a model, create an initial landscape and run it for ten iterations
#' forestgap_model <- ca_library("forestgap")
#' im <- generate_initmat(forestgap_model, c(0.5, 0.5), nrow = 64, ncol = 100)
#' run_camodel(forestgap_model, im, times = seq(0,100))
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
# Kubo, Takuya, Yoh Iwasa, and Naoki Furumoto. 1996. "Forest Spatial Dynamics with Gap
# Expansion: Total Gap Area and Gap Size Distribution." Journal of Theoretical Biology
# 180 (3): 229–46.
#
  if ( model == "forestgap" || model == "forest-gap") {
    if ( is.null(parms) ) {
      parms <- list(d = 0.125,
                    delta = 0.05,
                    alpha = 0.2)
    }

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
      all_states = c("EMPTY", "TREE"),
      check_model = "quick"
    )

  }

# Guichard's Mussel Bed model (2003)
#
# See also Génin et al. 2018 for model definition
#
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
      all_states = c("MUSSEL", "EMPTY", "DISTURB"),
      check_model = "quick"
    )
  }

# Kéfi's Arid vegetation model (2007), extended by Schneider (2016)
#
#
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
                    b = 0.4,
                    pr = 1)
    }

    wrap <- ifelse(is.null(wrap), TRUE, wrap)
    neighbors <- ifelse(is.null(neighbors), 4, neighbors)

    mod <- camodel(
      transition("DEGR", "EMPTY",
                 ~ r + q["VEGE"] * f),
      transition("EMPTY", "VEGE",
                 ~ ( delta * p["VEGE"] + ( 1 - delta ) * q["VEGE"] ) *
                     ( b - c * p["VEGE"] ) ),
      transition("EMPTY", "DEGR",
                 ~ d),
      transition("VEGE", "EMPTY",
                 ~ m0 + g0 * ( 1 - pr * q["VEGE"] )),
      parms = parms,
      wrap = wrap,
      neighbors = neighbors,
      all_states = c("DEGR", "EMPTY", "VEGE"),
      check_model = "quick"
    )

  }

  if ( model == "aridvege-danet" ) {
    if ( is.null(parms) ) {
      parms <- list(delta = 0.1,
                    b = 0.8,
                    c = 0.2,
                    gamma = 0.1,
                    g = 0.1,
                    u = 5,
                    m = 0.1, # or .2 ???
                    d = 0.1,
                    r = 0.01,
                    f = 0.9)
    }

    mod <- camodel(
      transition(from = "0", to = "N",
        ~ ( delta * p["N"] + ( 1 - delta ) * q["N"] ) *
            ( b - c * ( p["P"] + p["N"] ) - gamma )
      ),
      transition(from = "0", to = "P",
        ~ ( delta * p["P"] + ( 1 - delta ) * q["P"] ) *
            ( b - c * ( p["P"] + p["N"] ) -
              g * ( 1 - ( 1 - exp( - u * q["N"] ) ) ) )
      ),
      transition(from = "N", to = "0", ~ m),
      transition(from = "P", to = "0", ~ m),
      transition(from = "0", to = "-", ~ d),
      transition(from = "-", to = "0", ~ r + f * ( q["N"] + q["P"] ) ),
      wrap = TRUE,
      parms = parms,
      neighbors = 4
    )
  }


  # Genin's coral reef model (2023?)
  #
  if ( model == "coralreef" || model == "coral-reef" || model == "coral reef") {
    if ( is.null(parms) ) {
      parms <- list(r_a = 1.813037, 
                    l_a = 0.029767236894167, 
                    alpha = 0.01, 
                    m_a = 0.0792644531506507, 
                    r_c = 0.00130805533933221, 
                    d_0 = 0, 
                    m_c = 1e-04, 
                    g = 0.203809923511696, 
                    theta_b = 0.965518113386164, 
                    theta_c = 1.04159074812267, 
                    h_u = 2, 
                    rate_bl = 0, 
                    m_b = 0.0953125028179528)
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
                 ~ r_a * ( alpha + ( 1 - alpha ) * p["ALGAE"] ) + l_a * q["ALGAE"]),
      transition(from = "ALGAE", to = "BARE",
                 ~ m_a + h_u * g * ( theta_b * q["BARE"] + theta_c * q["CORAL"])),
      neighbors = neighbors,
      wrap = wrap,
      parms = parms,
      all_states = c("BARE", "ALGAE", "CORAL"),
      check_model = "quick"
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
      neighbors = neighbors,
      check_model = "quick"
    )

  }


  if ( is.null(mod) ) {
    all_models <- c("forestgap", "musselbed",  "coralreef", "aridvege",
                    "aridvege-danet",
                    "gameoflife", "rockpaperscissor")
    stop(paste0(model, " is an unknown model. Available models:\n",
                paste(" - ", all_models, collapse = "\n")))
  }

  return(mod)
}

