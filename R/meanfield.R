# 
# Run the meanfield version of a model 
# 

#'@title Mean field model 
#'
#'@description Run the mean field model corresponding to a cellular automaton 
#'
#'@param mod The cellular automaton model (produced by \code{\link{camodel}}
#'
#'@param init The initial landscape, or a vector of covers summing to one, whose length 
#'  is equal to the number of states in the model.
#'
#'@param times The points in time for which output is wanted 
#'
#'@param ... other arguments are passed to \code{\link[deSolve]{ode}}
#'
#'@return This function returns the results of the \code{\link[deSolve]{ode}} 
#'  function, which is a matrix with class 'deSolve'
#'
#'@details 
#'  
#'  The mean field approximation to a cellular automaton simply describes the dynamics 
#'    of the global covers using differential equations, assumming that both global 
#'    and local covers are equal (in \code{chouca} model specifications, this assumes 
#'    p = q). 
#'  
#'  For example, if we consider a model with two states 'a' and 'b' and transitions 
#'    from and to each other, then the following system of equation is used to describe 
#'    the variations of the proportions of cells in each state: 
#'  
#'  \deqn{\frac{da}{dt} = p_b P(b \to a) - p_a P(a \to b)}{ da/dt = p[b] P(b to a) - p[a] P(a to b)}
#'  \deqn{\frac{db}{dt} = p_a P(a \to b) - p_b P(b \to a)}{ db/dt = p[a] P(a to b) - p[b] P(b to a)}
#'  
#'  Running mean-field approximations is useful to understand general dynamics in the 
#'    absence of neighborhood interactions between cells, or simply to obtain an 
#'    fast but approximate simulation of the model. 
#'  
#'  Note that this function uses directly the expressions of the probabilities, so any 
#'    cellular automaton is supported, regardless of whether or not it can be simulated 
#'    exactly by \code{\link{run_camodel}}.
#'  
#'@examples
#' 
#' if ( requireNamespace("deSolve") ) { 
#'   # Get the mean-field approximation to the arid vegetation model 
#'   arid <- ca_library("aridvege") 
#'   mod <- ca_library("aridvege")
#'   init <- generate_initmat(mod, rep(1, 3)/3, nrow = 100, ncol = 100)
#'   times <- seq(0, 128)
#'   out <- run_meanfield(mod, init, times)
#'   # This uses the default plot method in deSolve
#'   plot(out)
#'   
#'   # A different model and way to specifiy initial conditions.
#'   coralmod <- ca_library("coralreef")
#'   init <- c(ALGAE = 0.2, CORAL = 0.5, BARE = 0.3)
#'   times <- 10^seq(0, 4, length.out = 64)
#'   out <- run_meanfield(coralmod, init, times, method = "lsoda")
#'   plot(out, ylim = c(0, 1))
#' }
#' 
#'@export
run_meanfield <- function(mod, init, times, ...) { 
  
  if ( ! requireNamespace("deSolve", quietly = TRUE) ) { 
    stop("Mean-field simulations require package 'deSolve', please install it using  install.packages(\"deSolve\")")
  }
  
  trs <- mod[["transitions"]]
  states <- mod[["states"]]
  
  # Handle init 
  if ( inherits(init, "camodel_initmat") ) { 
    init <- sapply(states, function(st) mean(init == st))
    names(init) <- states
  }
  
  if ( is.null(names(init)) ) { 
    names(init) <- states
  }
  
  if ( ! all( names(init) %in% states ) ) { 
    stop("Names in 'init' do not match the model states")
  }
  init <- init[states] # make sure init is in the right order
  
  # X is the vector of states
  dX <- function(time, X, parms, ...) { 
    names(X) <- states
    
    dX <- X * 0
    for ( i in seq_along(trs) ) { 
      tr <- trs[[i]]
      change <- eval(as.expression(as.list(tr[["prob"]])), 
                     envir = c(mod[["parms"]], list(p = X, q = X)), 
                     enclos = environment(tr[["prob"]]))
      change <- change * X[ tr[["from"]] ]
      dX[ tr[["to"]] ] <- dX[ tr[["to"]] ] + change 
      dX[ tr[["from"]] ] <- dX[ tr[["from"]] ] - change 
    }
    
    # Handle out of bound stuff by clamping things to 0-1
    dX <- ifelse(X + dX < 0, -X, dX)
    dX <- ifelse(X + dX > 1, 1 - X , dX)
    
    list(dX)
  }
  
  ode_out <- deSolve::ode(init, times = times, func = dX, ...)
  return(ode_out)
}
