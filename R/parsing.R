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
                    epsilon = sqrt(.Machine[["double.eps"]])) { 
  
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
  
  # TODO: check if some transitions are declared twice 
  
  msg(sprintf("Using %s cell states\n", nstates))
  msg(sprintf("Using %s transitions\n", length(transitions)))

  # Parse transitions
  # In the worst case scenario, we need 25 points, so that the number of neighbors is 
  # always divisible by 2, 3, 4 or 8. (25 because zero included)
  xpoints <- 1 + ifelse(wrap, neighbors, ifelse(neighbors == 8, 120, 24))
  transitions_parsed <- lapply(transitions, parse_transition, states, parms, xpoints,
                               epsilon)
  
  # Pack transitions into matrices
  qmat <- plyr::ldply(transitions_parsed, function(tr) { 
    data.frame(from = tr[["from"]], 
               to = tr[["to"]], 
               tr[["beta_q"]])
  })
  
  pmat <- plyr::ldply(transitions_parsed, function(tr) { 
    data.frame(from = tr[["from"]], 
               to = tr[["to"]], 
               tr[["beta_p"]])
  })
  
  alpha <- plyr::ldply(transitions_parsed, function(tr) { 
    data.frame(from = tr[["from"]], 
               to = tr[["to"]], 
               a0 = tr[["alpha0"]])
  })
  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, levels = states), # convert explicitely to factor
                transitions_defs = list(...), 
                parms            = parms, 
                alpha = alpha, 
                pmat = pmat, 
                qmat = qmat, 
                wrap = wrap, 
                neighbors = neighbors, 
                epsilon = epsilon, 
                xpoints = xpoints)
  
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

parse_transition <- function(tr, state_names, parms, xpoints, epsilon) { 
  
  if ( ! inherits(tr, "camodel_transition") ) { 
    m <- paste("The transition definition has not been defined using transition().", 
               "Please do using transition(from = ..., to = ..., prob = ~ ...)")
    stop(m)
  }
  ns <- length(state_names)
  
  # Make sure the environment for the formula is set to the empty environment. Otherwise
  # it trips up hash computing of the model, which is used to determine whether we have 
  # to recompile or not when using the 'compiled' engine.
  environment(tr[["prob"]]) <- emptyenv()
  
  pexpr <- as.expression( as.list(tr[["prob"]])[[2]] )
  zero <- rep(0, ns)
  names(zero) <- state_names
  
  # Constant probability component (when all the p and q are zero)
  warn <- character()
  prob_with <- function(p, q) { 
    p <- as.numeric( eval(pexpr, envir = c(parms, list(p = p, q = q))) ) 
    if ( p < 0 ) { 
      m <- paste0("Transition probabilities may go below zero. Problematic transition:\n", 
                 "  ", tr[["from"]], " -> ", tr[["to"]])
      warn <<- c(warn, m)
    }
    p
  }
  
  # Get intercept for this transition
  alpha0 <- prob_with(p = zero, q = zero)
  
  # Get coefficient/exponent table for q. We evaluate q at 'xpoints' points for vectors 
  # of q with zeros everywhere except one state. 
  beta_q <- plyr::ldply(state_names, function(s) { 
    q <- zero
    xs <- seq(0, 1, length.out = xpoints) 
    nqs <- seq_along(xs) - 1
    ys <- sapply(xs, function(x) { 
      q[s] <- x
      prob_with(p = zero, q = q) - alpha0
    })
    
    data.frame(state = s, qs = nqs, ys = ys)
  })
  
  beta_p <- plyr::ldply(state_names, function(s) { 
    p <- zero
    # In the worst case scenario, we need 9 points, in the case of 8 possible neighbors
    xs <- seq(0, 1, length.out = 25) 
    ys <- sapply(xs, function(x) { 
      p[s] <- x
      prob_with(p = p, q = zero) - alpha0
    })
    
    # Fit polynomial to this 
    poly <- fitpoly(xs, ys, epsilon)
    data.frame(state = s, coef = poly, expo = seq_along(poly))
  })
  
  # Print all the warnings accumulated so far
  if ( length(warn) > 0 ) { 
    # Some of them may be repeated
    warn <- unique(warn) 
    lapply(warn, warning, call. = FALSE)
  }
  
  list(from = tr[["from"]], 
       to   = tr[["to"]], 
       alpha0 = alpha0, 
       beta_q = beta_q, 
       beta_p = beta_p, 
       prob = pexpr)
  
}

# Update a ca_model with new arguments 
#'@export
update.ca_model <- function(object, parms, 
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

fitpoly <- function(x, y, epsilon) { 
  
  force(x)
  force(y)
  
  # If 
  if ( all(y < epsilon) ) { 
    return( 0 ) 
  }
  
  # We fit a polynomial
  polyc <- function(x, betas) { 
    sapply(x, function(i) { 
      sum(betas * i ^ seq_along(betas))
    })
  }

  ss <- function(betas) { 
    1e6 * sum( (y - polyc(x, betas))^2 ) 
  }
  
  degree <- 0
  error <- 1
  while ( error > 1e-10 && degree <= 25 ) { 
    degree <- degree + 1
    theta0 <- rep(0, degree)
    nlmp <- stats::optim(theta0, ss, method = "BFGS", 
                         control = list(trace = FALSE, 
                                        maxit = 1000, 
                                        abstol = .Machine$double.eps))
    error <- nlmp$value
  }
  
  np <- nlmp[["par"]]
  # Plot prediction
  ypred <- polyc(x, np)
  
  if ( error > 1e-10 ) { 
    stop("Could not approximate the given function with a polynomial - your model definition may be too complex for chouca")
  }
  
  return(np)
}

#'@export
print.ca_model <- function(x, ...) { 
  
  force(x) # force the evaluation of x before printing
  cat0 <- function(...) cat(paste0(...), "\n")
  
  cat0("Stochastic Cellular Automaton")
  cat0("")
  cat0("States: ", paste(x[["states"]], collapse = " "))
  cat0("")
  
  for ( tr in x[["transitions"]] ) { 
    cat0("Transition: ", tr[["from"]], " -> ", tr[["to"]])
    cat0(paste0("  ", as.character(tr[["prob"]]) ))
  }
  
  cat0("")
  cat0("Neighborhood: ", x[["neighbors"]], "x", x[["neighbors"]])
  cat0("Wrap: ", x[["wrap"]])
  
  return(invisible(x))
}

#'@export
print.camodel_transition <- function(x, ...) { 
  cat(sprintf("Transition from %s to %s\n", x[["from"]], x[["to"]]))
  prob <- paste0(as.character(x[["prob"]]), collapse = " ")
  cat(sprintf("  %s\n", prob))
  
  invisible(x)
}
