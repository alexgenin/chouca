# 
# This file contains the function that compiles a declarative model definition into 
# something that is acceptable for the CA engines. 
# 

#' @title Definition of a probabilistic cellular automaton 
#'
#' @description  High-level definition of a probabilistic cellular automaton
#' 
#' @param ... a number of transition descriptions, as built by the
#'   \code{\link{transition}} function (see Details)
#' 
#' @param neighbors The number of neighbors to use in the cellular automaton (4 for 4-way 
#'   or von-Neumann neghborhood, or 8 for a Moore neighborhood)
#' 
#' @param wrap Whether the 2D grid should wrap around at the edges
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
#' @param check_model A quick check of the model definition is done to make sure there 
#'   are no issues (e.g. probabilities outside the [1,0] interval, or an unsupported 
#'   model definition). Setting this argument to \code{FALSE} skips this checks, which is 
#'   a bit faster.
#' 
#' @param epsilon A small value under which coefficient values are considered to be 
#'   equal to zero
#' 
#' @param fixed_neighborhood When not using wrapping around the edges, the number of 
#'   neighbors per cell is variable, which can slow down the simulation. Set this 
#'   option to \code{TRUE} to consider that the number of neighbors is always four or 
#'   eight, regardless of the position of the cell in the landscape. 
#' 
#' @details
#' 
#' This function allows defining a stochastic cellular automaton model by its set of
#' transition rules. These are defined by a set of calls to the \code{transition()}
#' function. Each of these calls defines the two states of the transition, and the
#' probability (as a one-sided formula involving constants and special values p, q, 
#' etc.). 
#' 
#' \code{transition()} takes three arguments: the state from which the transition is 
#' defined, the state to which the transition goes, and a transition probability, defined
#' as a one-sided formula. This formula can include numerical constants, parameters 
#' defined in the named list \code{parms}, and any combination of p['a'] and q['b'], 
#' which respectively represent the proportion of cells in a landscape in state 'a', and 
#' the proportion of neighbors of a given cell in state 'b' ('a', and 'b' being, of 
#' course, any of the possible states defined in the model). See section Examples for 
#' examples of model implementations. 
#' 
#' It is important to remember that \code{chouca} only supports models where the 
#' transition probabilities are polynomials of p and q with maximum degree 2. In other
#' words, for a model with n states, the transition probabilities must be of the form: 
#' 
#' P = \eqn{a_0 + \sum_{k=1}^{n} b_k p_{k} + c_k p^2_{k} + d_k q_{k} + e_k q^2_{k} }
#' 
#' where p_{k} and q_{k} are the proportions of cells in state k in the landscape, and 
#' in the cell neighborhood, respectively. \eqn{a_0}, along with the coeffients 
#' \eqn{b_k}, \eqn{c_k}, \eqn{d_k} and \eqn{e_k} must be all constants. However, they 
#' can themselve depend on values defined in the named list \code{parms}. This 
#' allows defining the model in terms of literal constants, and compute internal
#' coefficients on the fly by picking values from the \code{parms} list. 
#' 
#' \code{camodel()} will run a few checks on your model definition to make sure 
#' transition probabilities do not go above one or below zero, and that the transition 
#' probabilities fit the formula above. These checks should catch most errors, but are 
#' not infaillible. They can be disabled using \code{check_model = TRUE}.
#' 
#' When using compiled code used to run the model, very small coefficients in the formula
#' above are rounded down to zero. This may be an issue if your transition 
#' probabilities are very low: in this case, consider reducing \code{epsilon} to a value
#' closer to zero.
#' 
#' To run a model once it is defined, the function \code{\link{run_camodel}} can be used.
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
#'   wrap = TRUE
#' )
#' 
#' # A fun plant model 
#' mod <- camodel(
#'   transition("plant", "empty", ~ death * ( 1 - (2*q["plant"]-1)^2) ), 
#'   transition("empty", "plant", ~ q["plant"]^2 ), 
#'   all_states = c("empty", "plant"), 
#'   wrap = TRUE, 
#'   neighbors = 4, 
#'   parms = list(death = 0.2496)
#' )
#' 
#' # Conway's Game of Life 
#' mod <- camodel( 
#'   transition("LIVE", "DEAD", ~ q["LIVE"] < (2/8) | q["LIVE"] > (3/8)), 
#'   transition("DEAD", "LIVE", ~ q["LIVE"] == (3/8)), 
#'   wrap = TRUE, 
#'   neighbors = 8, 
#'   all_states = c("DEAD", "LIVE")
#' )
#' 
#' # A spiral-generating rock-paper-scissor model (run with substeps = 1)
#' mod <- camodel(
#'   transition(from = "r", to = "p", ~ prob * ( q["p"] > (1/8)*2) ), 
#'   transition(from = "p", to = "c", ~ prob * ( q["c"] > (1/8)*2) ), 
#'   transition(from = "c", to = "r", ~ prob * ( q["r"] > (1/8)*2) ), 
#'   parms = list(prob = 1), 
#'   wrap = TRUE, 
#'   neighbors = 8
#' )
#' 
#'@export
camodel <- function(..., 
                    neighbors, # default to von-neumann neighborhood
                    wrap, 
                    parms = list(), 
                    all_states = NULL, 
                    check_model = TRUE, 
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
  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, levels = states), # convert explicitely to factor
                transitions_defs = list(...), 
                parms = parms, 
                beta_0 = beta_0, 
                beta_q = beta_q, 
                beta_pp = beta_pp, 
                beta_pq = beta_pq, 
                beta_qq = beta_qq, 
                wrap = wrap, 
                neighbors = neighbors, 
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
#'@export
update.ca_model <- function(object, 
                            parms = NULL, 
                            neighbors = NULL, 
                            wrap = NULL, 
                            check_model = TRUE, 
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
  
  # Extract model parameters, and do the call
  newcall <- c(object[["transitions_defs"]], # always a list, so result of c() is a list
               neighbors = neighbors, 
               wrap = wrap, 
               parms = list(parms), 
               all_states = list(object[["states"]]), 
               check_model = check_model, 
               verbose = verbose, 
               epsilon = object[["epsilon"]])
  
  do.call(camodel, newcall)
}
