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
#' @param parms a named list of parameters, which should be all numeric, single values
#' 
#' @param all_states the complete list of states (a character vector). If unspecified,
#'   it will be guessed from the transition rules, but it is a good idea pass it here 
#'   to make sure the model definition is correct 
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
#' It is important to remember that \code{chouca} only supports model where the 
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
#' To optimize compiled code used to run the model, very small coefficients in the formula
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
#'   all_states = c("EMPTY", "TREE")
#' )
#' 
#' # A funny plant model 
#' 
#' mod <- camodel(transition("plant", "empty", ~ death * ( 1 - (2*q["plant"]-1)^2) ) , 
#'                transition("empty", "plant", ~ q["plant"]^2), 
#'                all_states = c("empty", "plant"), 
#'                parms = list(death = 0.2496))
#' 
#'@export
camodel <- function(..., 
                    parms = list(), 
                    all_states = NULL, 
                    verbose = TRUE, 
                    check_model = TRUE, 
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
  
  # Compute transition probabilities
  transitions <- lapply(transitions, parse_transition, states, parms, epsilon, 
                        check_model)
  
  msg(sprintf("Using %s cell states\n", nstates))
  msg(sprintf("Using %s transitions\n", length(transitions)))
  
  # Check that there is no NA in transition rates 
  allrates <- plyr::laply(transitions, function(o) { 
    c(o[["X0"]], o[["XP"]], o[["XQ"]], o[["XPSQ"]], o[["XQSQ"]], use.names = FALSE)
  })
  
  if ( any(is.na(allrates)) ) { 
    stop("NAs in computed coefficients, please make sure your model definition is correct")
  }
  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, states))
  
  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

unform <- function(form) { 
  as.expression( as.list(form)[[2]] )
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

parse_transition <- function(tr, state_names, parms, epsilon, check_model) { 
  
  if ( ! inherits(tr, "camodel_transition") ) { 
    m <- paste("The transition definition has not been defined using transition().", 
               "Please do using transition(from = ..., to = ..., prob = ~ ...)")
    stop(m)
  }
  ns <- length(state_names)
  
  # Make sure the environment for the formula is set to the empyt environment. Otherwise
  # it trips up hash computing of the model, which is used to determine whether we have 
  # to recompile or not when using the 'compiled' engine.
  environment(tr[["prob"]]) <- emptyenv()

  pexpr <- unform(tr[["prob"]]) 
  zero <- rep(0, ns)
  names(zero) <- state_names
  
  # Constant probability component (when all the p and q are zero)
  prob_with <- function(p, q) { 
    eval(pexpr, envir = c(parms, list(p = p, q = q)))
  }
  
  X0 <- prob_with(p = zero, q = zero)
  names(X0) <- NULL 
  
  # Global and local density probability component 
  XP <- XPSQ <- XQ <- XQSQ <- numeric(ns)
  
  # See my notes for this
  for ( i in seq.int(ns) ) { 
    
    onevec <- zero
    onevec[i] <- 1 # density of focal state = full
    
    p_one_zero  <- prob_with(p = onevec,   q = zero)
    p_mone_zero <- prob_with(p = - onevec, q = zero)
    XP[i] <- (p_one_zero - p_mone_zero) / 2
    XPSQ[i] <- p_one_zero - X0 - XP[i] 
    
    p_one_zero  <- prob_with(p = zero,  q = onevec)
    p_mone_zero <- prob_with(p = zero,  q = - onevec)

    XQ[i] <- (p_one_zero - p_mone_zero) / 2
    XQSQ[i] <- p_one_zero - X0 - XQ[i]
  }
  
  # Compute probabilities and make sure there is no residual error 
  if ( check_model ) { 
    state_combs <- as.matrix(do.call(expand.grid, 
                                    rep(list(seq(0, 1, l = 4)), length = 2*ns)))
    colnames(state_combs) <- c(paste0("p", seq.int(ns)), paste0("q", seq.int(ns)))
    y <- apply(state_combs, 1, function(v) { 
      p <- v[1:ns]
      q <- v[seq(ns+1, ns+ns)]
      X0 + sum(XP * p) + sum(XPSQ * p^2) + sum(XQ * q) + sum(XQSQ * q^2)
    })
    yok <- apply(state_combs, 1, function(v) { 
      p <- v[seq(1, ns)]
      names(p) <- state_names
      q <- v[seq(ns+1, ns+ns)]
      names(q) <- state_names
      prob_with(p = p, q = q)
    })
    
    if ( ! all( abs(y - yok) < epsilon ) ) { 
      stop("There is residual error in the computed probabilities. Your model looks like it is not supported by chouca. Please look at ?camodel for more details about the supported classes of models")
    }
    
    if ( any( y > 1 ) ) { 
      maxy <- ceiling(max(y))
      warning(sprintf("Model probabilities can go above one, consider increasing the number of 
        substeps to %s or more", maxy))
    }
    
    if ( any( y < 0 ) ) { 
      warning("Model probabilities can be negative!")
    }
    
  }
  
  # Make sure to zero things that are zero 
  eps <- function(x) ifelse(abs(x)<epsilon, 0, x)
  
  c(tr, list(X0 = eps(X0), XP = eps(XP), XQ = eps(XQ), 
             XPSQ = eps(XPSQ), XQSQ = eps(XQSQ)))
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
    cat0("  p = ", as.character(tr[["prob"]])[2] )
  }
  
  return(invisible(x))
}

#'@export
print.camodel_transition <- function(x, ...) { 
  cat(sprintf("Transition from %s to %s\n", x[["from"]], x[["to"]]))
  prob <- paste0(as.character(x[["prob"]]), collapse = " ")
  cat(sprintf("  %s\n", prob))
  
  invisible(x)
}
