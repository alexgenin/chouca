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
#' 
camodel <- function(..., parms, all_states = NULL, verbose = TRUE) { 
  if ( ! is.list(parms) || is.null(names(parms)) || any(names(parms) == "") ) { 
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
  transitions <- lapply(transitions, parse_transition, states, parms)
  
  caobj <- list(transitions = transitions, 
                nstates = nstates, 
                states = factor(states, states))
  
  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

unform <- function(form) { 
  as.expression( as.list(form)[[2]] )
}

transition <- function(from, to, prob) { 
#   if ( ! is.numeric(from) && length(from) == 1 ) { 
#     stop("from is not a single-length numeric vector")
#   }
#   if ( ! is.numeric(to) && length(to) == 1 ) { 
#     stop("from is not a single-length numeric vector")
#   }
  if ( ! inherits(prob, "formula") ) { 
    stop("prob must be a formula")
  }
  
  o <- list(from = from, 
            to   = to, 
            prob = prob)
  class(o) <- c("camodel_transition", "list")
  
  return(o)
}

parse_transition <- function(tr, state_names, parms) { 
  
  if ( ! inherits(tr, "camodel_transition") ) { 
    m <- paste("The transition definition has not been defined using transition().", 
               "Please do using transition(from = ..., to = ..., prob = ~ ...)")
    stop(m)
  }
  ns <- length(state_names)
  
  pexpr <- unform(tr[["prob"]]) 
  zero <- rep(0, ns)
  names(zero) <- state_names
  
  # Constant probability component (when all the p and q are zero)
  prob_with <- function(p, q) { 
    eval(pexpr, envir = c(parms, list(p = p, q = q)))
  }
  
  X0 <- prob_with(p = zero, q = zero)
  
  # Global and local density probability component 
  XP <- numeric(ns)
  XQ <- numeric(ns)
  for ( i in seq.int(ns) ) { 
    prob0 <- prob_with(p = zero, q = zero)
    onevec <- zero
    onevec[i] <- 1 # density of focal state = full
    prob1 <- prob_with(p = onevec, q = zero)
    XP[i] <- prob1 - prob0
    
    prob1 <- prob_with(p = zero, q = onevec)
    XQ[i] <- prob1 - prob0
  }
  
  c(tr, list(X0 = X0, XP = XP, XQ = XQ))
}


