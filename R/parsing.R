# 
# This file contains the function that compiles a declarative model definition into 
# something that is acceptable for the c++ code
# 

camodel <- function(..., parms, verbose = TRUE) { 
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
  
  states <- ldply(transitions, function(o) as.data.frame(o[c("from", "to")]))
  uniqstates <- unique(states[ ,"from"], states[ ,"to"])
  nstates <- length(uniqstates)
  msg(sprintf("Found %s states\n", nstates))
  
  if ( ! all( sort(uniqstates) == seq.int(nstates) ) ) { 
    stop("all CA states should be numbered from one to the maximum, without gaps")
  }
  
  # Compute transition probabilities
  transitions <- llply(transitions, parse_transition, 
                       nstates = nstates, parms = parms)
  
  caobj <- list(transitions = transitions, 
                nstates = nstates)
  
  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

unform <- function(form) { 
  as.expression( as.list(form)[[2]] )
}

transition <- function(from, to, prob) { 
  if ( ! is.numeric(from) && length(from) == 1 ) { 
    stop("from is not a single-length numeric vector")
  }
  if ( ! is.numeric(to) && length(to) == 1 ) { 
    stop("from is not a single-length numeric vector")
  }
  if ( ! inherits(prob, "formula") ) { 
    stop("prob must be a formula")
  }
  
  o <- list(from = from, 
            to   = to, 
            prob = prob)
  class(o) <- c("camodel_transition", "list")
  
  return(o)
}

parse_transition <- function(tr, nstates, parms) { 
  
  if ( ! inherits(tr, "camodel_transition") ) { 
    m <- paste("The transition definition has not been defined using transition().", 
               "Please do using transition(from = ..., to = ..., prob = ~ ...)")
    stop(m)
  }
  
  pexpr <- unform(tr[["prob"]]) 
  zero <- rep(0, nstates)
  
  # Constant probability component (when all the p and q are zero)
  prob_with <- function(p, q) { 
    eval(pexpr, envir = c(parms, list(p = p, q = q)))
  }
  
  X0 <- prob_with(p = zero, q = zero)
  
  # Global and local density probability component 
  XP <- numeric(nstates)
  XQ <- numeric(nstates)
  for ( i in seq.int(nstates) ) { 
    prob0 <- prob_with(p = zero, q = zero)
    onevec <- zero
    onevec[i] <- 1
    prob1 <- prob_with(p = onevec, q = zero)
    XP[i] <- prob1 - prob0
    
    prob1 <- prob_with(p = zero, q = onevec)
    XQ[i] <- prob1 - prob0
  }
  
  c(tr, list(X0 = X0, XP = XP, XQ = XQ))
}


