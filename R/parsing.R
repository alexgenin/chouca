# 
# This file contains the function that compiles a declarative model definition into 
# something that is acceptable for the CA engines. 
# 

#' @title Definition of a stochastic cellular automaton 
#'
#' @description  High-level definition of a stochastic cellular automaton
#' 
#' @param ... a number of transition descriptions, as built by the
#'   \code{\link{transition}} function (see Details and Examples)
#' 
#' @param neighbors The number of neighbors to use in the cellular automaton (4 for 4-way 
#'   or von-Neumann neghborhood, or 8 for an 8-way or Moore neighborhood)
#' 
#' @param wrap If \code{TRUE}, then the 2D grid on which the model is run wraps around 
#'   at the edges (the top/leftmost cells will be considered neighbors of the 
#'   bottom/rightmost cells)
#' 
#' @param continuous If \code{TRUE}, the model definition should be treated as a 
#'   continuous stochastic cellular automaton
#' 
#' @param parms a named list of parameters, which should be all numeric, single values
#' 
#' @param all_states the complete list of states (a character vector). If unspecified,
#'   it will be guessed from the transition rules, but it is a good idea pass it here 
#'   to make sure the model definition is correct. 
#' 
#' @param verbose whether information should be printed when parsing the model
#'   definition. 
#' 
#' @param check_model A check of the model definition is done to make sure there 
#'   are no issues with it (e.g. probabilities outside the [1,0] interval, or an 
#'   unsupported model definition). A quick check that should catch most problems is
#'   performed if check_model is "quick", an extensive check that tests all 
#'   neighborhood configurations is done with "full", and no check is performed with 
#'   "none".
#' 
#' @param epsilon A small value under which coefficient values are considered to be 
#'   equal to zero
#' 
#' @param fixed_neighborhood When not using wrapping around the edges (\code{wrap = TRUE},
#'   the number of neighbors per cell is variable, which can slow down the simulation. 
#'   Set this option to \code{TRUE} to consider that the number of neighbors is always 
#'   four or eight, regardless of the position of the cell in the landscape, at the cost 
#'   of approximate dynamics on the edge of the landscape.
#' 
#' @details
#' 
#' This function allows defining a stochastic cellular automaton model by its set of
#' transition rules. These are defined by a set of calls to the \code{transition()}
#' function. Each of these calls defines the two states of the transition, and the
#' probability (as a one-sided formula involving constants and special values p, q, 
#' etc.). 
#' 
#' \code{transition()} calls takes three arguments: the state from which the transition 
#' is defined, the state to which the transition goes, and a transition probability, 
#' defined as a one-sided formula. This formula can include numerical constants, 
#' parameters defined in the named list \code{parms}, and any combination of p['a'] 
#' and q['b'], which respectively represent the proportion of cells in a landscape in 
#' state 'a', and the proportion of neighbors of a given cell in state 'b' ('a', and 
#' 'b' being, of course, any of the possible states defined in the model). See 
#' section Examples for examples of model implementations. 
#' 
#' It is important to remember when using this function that \code{chouca} only 
#' supports models where the probabilities follow the following functional form: 
#' 
#' P = \deqn{\beta_0 + \sum_{k=1}^S f(q_k) + \sum{k=1}^S \beta^p_k p_k + \beta^{pq}_1 p^2q^2 \dots \beta^{pq}_I p^5q^5 + \beta^{pp}_J p^2p^2 \dots \beta^{pp}_J p^5p^5 + + \beta^{qq}_K q^2q^2 \dots \beta^{qq}_K q^5q^5}
#' 
#' where p_{k} and q_{k} are the proportions of cells in state k in the landscape, and 
#' in the cell neighborhood, respectively, and the various \eqn{\beta} coefficient 
#' are estimated by \code{chouca} internally. When \code{check_model} is "quick" or 
#' "full", a check is performed to make sure the functional form above is able to 
#' accurately represent probabilities of transitions in the model, with "full" enabling 
#' more extensive tesing, and "none" removing it entirely.  
#' 
#' When using chouca, very small coefficients in the formula above may be rounded down 
#' to zero. This may be an issue if your transition probabilities are very low: in 
#' this case, consider reducing \code{epsilon} to a smaller value. 
#' 
#' To run a model once it is defined, the function \code{\link{run_camodel}} can be 
#' used, or \code{\link{run_meanfield}} for a mean-field approximation. An initial 
#' landscape for a simulation can be created using \code{\link{generate_initmat}}. 
#' 
#' @seealso run_camodel, generate_initmat, run_meanfield, update.ca_model
#' 
#' @examples 
#' 
#' # Redefine Kubo's 1996 forest gap model 
#' kubo <- camodel( 
#'   transition(from = "TREE", 
#'              to   = "EMPTY", 
#'              prob = ~ d + delta * q["EMPTY"] ), 
#'   transition(from = "EMPTY", 
#'              to   = "TREE", 
#'              prob = ~ alpha), 
#'   parms = list(d = 0.125, 
#'                delta = 0.5, 
#'                alpha = 0.2), 
#'   all_states = c("EMPTY", "TREE"), 
#'   neighbors = 4, 
#'   wrap = TRUE, 
#'   continuous = FALSE
#' )
#' 
#' # A fun plant model 
#' mod <- camodel(
#'   transition("plant", "empty", ~ death * ( 1 - (2*q["plant"]-1)^2) ), 
#'   transition("empty", "plant", ~ q["plant"]^2 ), 
#'   all_states = c("empty", "plant"), 
#'   wrap = TRUE, 
#'   neighbors = 4, 
#'   continuous = FALSE, 
#'   parms = list(death = 0.2496)
#' )
#' 
#' # Conway's Game of Life 
#' mod <- camodel( 
#'   transition("LIVE", "DEAD", ~ q["LIVE"] < (2/8) | q["LIVE"] > (3/8)), 
#'   transition("DEAD", "LIVE", ~ q["LIVE"] == (3/8)), 
#'   wrap = TRUE, 
#'   neighbors = 8, 
#'   continuous = FALSE, 
#'   all_states = c("DEAD", "LIVE")
#' )
#' 
#' # A spiral-generating rock-paper-scissor model
#' mod <- camodel(
#'   transition(from = "r", to = "p", ~ prob * ( q["p"] > (1/8)*2) ), 
#'   transition(from = "p", to = "c", ~ prob * ( q["c"] > (1/8)*2) ), 
#'   transition(from = "c", to = "r", ~ prob * ( q["r"] > (1/8)*2) ), 
#'   parms = list(prob = 1), 
#'   wrap = TRUE, 
#'   continuous = FALSE, 
#'   neighbors = 8
#' )
#' 
#'@export
camodel <- function(..., 
                    neighbors, 
                    wrap, 
                    continuous = FALSE, 
                    parms = list(), 
                    all_states = NULL, 
                    check_model = "quick", 
                    verbose = FALSE, 
                    epsilon = sqrt(.Machine[["double.eps"]]), 
                    fixed_neighborhood = FALSE) { 
  
  if ( ( ! identical(parms, list()) ) && 
        ( ! is.list(parms) || is.null(names(parms)) || any(names(parms) == "") ) ) { 
    stop("parms must be a list all elements named")
  }
  
  if ( any( c("p", "q") %in% names(parms) ) ) { 
    stop("no parameters must be named 'q' or 'p', these are reserved to refer to densities")
  }
  
  if ( is.logical(check_model) && isFALSE(check_model) ) { 
    check_model <- "none"
  }
  
  if ( ! check_model %in% c("none", "full", "quick") ) { 
    stop("'check_model' must be one of 'none', 'quick' or 'full'")
  }
  
  # Read all possible states 
  msg <- function(txt) if ( verbose ) message(txt, appendLF = FALSE)
  msg("Parsing CA model definition\n")

  # Read transition objects
  transitions <- list(...)
  
  # Check that they are all transitions 
  if ( ! all(sapply(transitions, inherits, "camodel_transition")) ) { 
    stop("One of the model transition was not processed by the transition() function. Please check your function arguments and the ?camodel help on how to define model state transitions.")
  }
  
  states <- plyr::ldply(transitions, function(o) as.data.frame(o[c("from", "to")]) )
  uniqstates <- unique(c(states[ ,"from"], states[ ,"to"]))
  if ( is.null(all_states) ) { 
    nstates <- length(uniqstates)
    states <- uniqstates
  } else { 
    nstates <- length(all_states) 
    states <- all_states
    if ( ! all(uniqstates %in% states) ) { 
      stop("Some cell states defined in transitions are not defined in 'all_states'")
    }
  }
  
  # Make sure transitions are unique
  trcodes <- plyr::laply(transitions, function(tr) { 
    paste(tr[["from"]], tr[["to"]], sep = " -> ")
  })
  if ( any(duplicated(trcodes)) ) { 
    msg <- paste("Duplicated transition(s) in model definition:\n", 
                 paste("  ", trcodes[duplicated(trcodes)], collapse = "\n"))
    stop(msg)
  }
  
  msg(sprintf("Using %s cell states\n", nstates))
  msg(sprintf("Using %s transitions\n", length(transitions)))
  
  # Parse transitions
  # In the worst case scenario, we need 25 points, so that the number of neighbors is 
  # always divisible by 2, 3, 4 or 8. (25 because zero included)
  xpoints <- 1 + ifelse(wrap, neighbors, ifelse(neighbors == 8, 120, 24))
  
  # If we want a fixed neighborhood, regardless of whether if wrap or not, then 
  # overwrite xpoints with the fixed value. 
  xpoints <- ifelse(fixed_neighborhood, 1+neighbors, xpoints)
  
  transitions_parsed <- lapply(transitions, parse_transition, states, parms, xpoints,
                               epsilon, neighbors, check_model)
  
  # Pack transitions into matrices
  beta_q  <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_q")
  beta_pp <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_pp")
  beta_pq <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_pq")
  beta_qq <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_qq")
  beta_0  <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_0")
  
  # Extract mean/max errors
  max_error <- plyr::laply(transitions_parsed, function(o) o[["max_error"]])
  max_rel_error <- plyr::laply(transitions_parsed, function(o) o[["max_rel_error"]])
  
  # Identify absorbing states 
  tr <- lapply(transitions, function(o) cbind(o[["from"]], o[["to"]]))
  tr <- do.call(rbind, tr)
  is_absorb_state <- tr[ ,2][ ! (tr[ ,2] %in% tr[ ,1]) ]
  is_fixed_state <- states[ ! states %in% tr ]
  absorb_states <- states[states %in% c(is_absorb_state, is_fixed_state)]
  

  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, levels = states), # convert explicitely to factor
                parms = parms, 
                beta_0 = beta_0, 
                beta_q = beta_q, 
                beta_pp = beta_pp, 
                beta_pq = beta_pq, 
                beta_qq = beta_qq, 
                absorbing_states = absorb_states, 
                wrap = wrap, 
                neighbors = neighbors, 
                continuous = continuous, 
                epsilon = epsilon, 
                xpoints = xpoints, 
                max_error = max_error, 
                max_rel_error = max_rel_error, 
                fixed_neighborhood = fixed_neighborhood)
  
  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

#'@describeIn camodel
#' 
#' @param from The state from which the transition is defined 
#' 
#' @param to The state to which the transition is defined 
#' 
#' @param prob a one-sided formula describing the probability of transition between the two 
#'  states (see Details section for more information).
#' 
#'@export
transition <- function(from, to, prob) { 
  if ( length(from) != 1 ) { 
    stop("from is not a single-length vector")
  }
  if ( length(to) != 1 ) { 
    stop("from is not a single-length vector")
  }
  if ( ! inherits(prob, "formula") ) { 
    stop("prob must be a one-sided formula ('~ <expr>')")
  }
  
  o <- list(from = from, 
            to   = to, 
            prob = prob)
  class(o) <- c("camodel_transition", "list")
  
  return(o)
}

# Helper function to pack transitions into model-level tables
pack_table_fromto <- function(tr, table) { 
  
  if ( nrow(tr[[table]]) == 0 ) { 
    data.frame(from = integer(0), to = integer(0), tr[[table]])
  } else { 
    data.frame(from = tr[["from"]], to = tr[["to"]], tr[[table]])
  }
}

# Update a ca_model with new arguments 
# first argument needs to be 'object' to respect S3 method naming
#' @title Update a cellular automaton 
#'
#' @description Update the definition of a stochastic cellular automaton 
#'   (SCA), using new parameters, type of wrapping, or any other parameters 
#'   entering in the definition of the model.
#'
#' @param object The SCA object (returned by \code{\link{camodel()}}) 
#' 
#' @param parms a named list of parameters, which should be all numeric, 
#'   single values
#' 
#' @param neighbors The number of neighbors to use in the cellular automaton 
#'   (4 for 4-way or von-Neumann neghborhood, or 8 for an 8-way or Moore 
#'   neighborhood)
#' 
#' @param wrap If \code{TRUE}, then the 2D grid on which the model is run wraps 
#'   around at the edges (the top/leftmost cells will be considered neighbors 
#'   of the bottom/rightmost cells)
#' 
#' @param continuous If \code{TRUE}, the model definition should be treated as 
#'   a continuous stochastic cellular automaton
#' 
#' @param fixed_neighborhood When not using wrapping around the edges (\code{wrap = TRUE},
#'   the number of neighbors per cell is variable, which can slow down the simulation. 
#'   Set this option to \code{TRUE} to consider that the number of neighbors is always 
#'   four or eight, regardless of the position of the cell in the landscape, at the cost 
#'   of approximate dynamics on the edge of the landscape.
#' 
#' @param check_model A check of the model definition is done to make sure there 
#'   are no issues with it (e.g. probabilities outside the [1,0] interval, or an 
#'   unsupported model definition). A quick check that should catch most problems is
#'   performed if check_model is "quick", an extensive check that tests all 
#'   neighborhood configurations is done with "full", and no check is performed with 
#'   "none".
#' 
#' @param verbose whether information should be printed when parsing the model
#'   definition. 
#' 
#' @param ... extra arguments are ignored 
#' 
#' @details 
#'  This function allows you to update some aspects of a pre-defined celullar 
#'    automaton, such as parameter values, the type of neighborhood, whether 
#'    to wrap around the edge of space, etc. It is handy when running multiple 
#'    simulations over a gradient of values a given parameter. 
#'
#' @seealso camodel, run_camodel
#'
#'@examples 
#' 
#' # Update the parameters of a model 
#' mussels <- ca_library("musselbed")
#' mussels[["parms"]] # old parameters 
#' mussels_new <- update(mussels, parms = list(d = 0.2, delta = 0.1, r = 0.8))
#' mussels_new[["parms"]] # updated parameters 
#' 
#' # Update the type of neighborhood, wrapping around the edges, and 
#' # the parameters
#' mussels_new <- update(mussels, 
#'                       parms = list(d = 0.2, delta = 0.1, r = 0.8), 
#'                       wrap = TRUE, 
#'                       neighbors = 8)
#' mussels_new 
#' 
#' # Run the model for different values of d, the death rate of mussels
#' ds <- seq(0, 0.25, length.out = 12)
#' initmat <- generate_initmat(mussels, c(0.5, 0.5, 0), nr = 64, nc = 64)
#' results <- lapply(ds, function(this_dvalue) { 
#'   musselmod <- update(mussels, parms = list(d = this_dvalue))
#'   run <- run_camodel(musselmod, initmat, times = seq(0, 128))
#'   data.frame(d = this_dvalue, 
#'              as.data.frame(tail(run[["output"]][["covers"]], 1)))
#' })
#' results <- do.call(rbind, results)
#' plot(results[ ,"d"], results[ ,"MUSSEL"], type = "b", 
#'      xlab = "d", ylab = "Mussel cover")
#' 
#'@export
update.ca_model <- function(object, 
                            parms = NULL, 
                            neighbors = NULL, 
                            wrap = NULL, 
                            continuous = NULL, 
                            fixed_neighborhood = NULL, 
                            check_model = "quick", 
                            verbose = FALSE, 
                            ...) { 
  
  if ( is.null(wrap) ) { 
    wrap <- object[["wrap"]]
  }
  if ( is.null(neighbors) ) { 
    neighbors <- object[["neighbors"]]
  }
  if ( is.null(parms) ) { 
    parms <- object[["parms"]]
  }
  if ( is.null(continuous) ) { 
    continuous <- object[["continuous"]]
  }
  if ( is.null(fixed_neighborhood) ) { 
    fixed_neighborhood <- object[["fixed_neighborhood"]]
  }
  
  if ( ! is.null(parms) ) { 
    parms_orig <- object[["parms"]]
    if ( is.null(names(parms)) || any( names(parms) == "" ) ) { 
      stop("Not all elements are named in the parameter list")   
    }
    if ( ! all( names(parms) %in% names(parms_orig) ) ) { 
      stop("New parameter names do not match original parameters")
    }
    parms_orig[names(parms)] <- parms
    parms <- parms_orig
  }
  
  # Extract model parameters, and do the call
  newcall <- c(object[["transitions"]], # always a list, so result of c() is a list
               neighbors = neighbors, 
               wrap = wrap, 
               fixed_neighborhood = fixed_neighborhood, 
               continuous = continuous, 
               parms = list(parms), 
               all_states = list(object[["states"]]), 
               check_model = check_model, 
               verbose = verbose, 
               epsilon = object[["epsilon"]])
  
  do.call(camodel, newcall)
}
