# 
# This file contains the function that compiles a declarative model definition into 
# something that is acceptable for the c++ code
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
#'   it will be guessed from the transition rules, but you should really pass it to 
#'   make sure your model definition is correct 
#' 
#' @param verbose whether information should be printed when parsing the model
#'   definition
#' 
#' @details 
#' 
#' This function allows defining a cellular automaton model by its set of transition 
#' rules. These are defined by a set of calls to the \code{transition()} function. Each 
#' of these calls defines the two states of the transition, and the probability (as a 
#' one-sided formula involving constants and p, q, etc.). 
#' 
#' 
#' 
#' 
#'@export
camodel <- function(..., parms = list(), all_states = NULL, verbose = TRUE, 
                    epsilon = sqrt(.Machine[["double.eps"]])) { 
                      
  
  if ( ( ! identical(parms, list()) ) && 
        ( ! is.list(parms) || is.null(names(parms)) || any(names(parms) == "") ) ) { 
    stop("parms must be a named with all elements named")
  }
  
  if ( any( c("p", "q") %in% names(parms) ) ) { 
    stop("no parameters must be named 'q' or 'p', these are used to design densities")
  }
  # Read all possible states 
  msg <- function(txt) if ( verbose ) message(txt, appendLF = FALSE)
  
  msg("Parsing CA model definition\n")
  
  # Read transition objects
  transitions <- list(...)
  
  states <- plyr::ldply(transitions, function(o) as.data.frame(o[c("from", "to")]))
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
  
  msg(sprintf("Using %s cell states: %s\n", nstates, paste(states, collapse = " ")))
  
  # Compute transition probabilities
  transitions <- lapply(transitions, parse_transition, states, parms, epsilon)
  
  # Check that there is no NA in transition rates 
  allrates <- plyr::laply(transitions, function(o) { 
    c(o[["X0"]], o[["XP"]], o[["XQ"]], use.names = FALSE)
  })
  if ( any(is.na(allrates)) ) { 
    stop("NAs in computed coefficients, please make sure your model definition is correct")
  }
  
  # Compute some model information 
  # TODO: does the model has non-zero XQ/XP/XQSP/XPSQ/XQP ?
  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, states))
  
  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

unform <- function(form) { 
  as.expression( as.list(form)[[2]] )
}

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

parse_transition <- function(tr, state_names, parms, epsilon) { 
  
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
  XP <- XPSQ <- XQ <- XQSQ <- XPQ <- numeric(ns)
  
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
  
  # Make sure to zero things that are zero 
  eps <- function(x) ifelse(x<epsilon, 0, x)
  
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

