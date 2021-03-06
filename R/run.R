# 
# This file contains function that will initialize and run a model 
#

#' @title Generate an initial matrix for a chouca model 
#' 
#' @description Helper to create a initial landscape (matrix) with specified covers for 
#'   a cellular automaton 
#' 
#' @param mod A stochastic cellular automaton model defined using \code{\link{camodel}}
#' 
#' @param pvec The vector of covers of each state in the initial configuration. 
#' 
#' @param nr The number of rows of the output matrix 
#' 
#' @param nc The number of columns of the output matrix 
#' 
#' @details 
#'   This function is a helper to build a starting configuration (matrix) for a 
#'     stochastic cellular automaton. Based on the definition of the model and the 
#'     specified starting covers (in \code{pvec}). It will produce a landscape with 
#'     expected global cover of each state equal to the covers in \code{pvec}.
#'   
#'   The length of the \code{pvec} vector must match the number of possible cell states 
#'     in the model. If present, the names of \code{pvec} must match the states
#'     defined in the model. In this case, they will be used to determine which state 
#'     gets which starting cover. 
#'   
#'   The \code{pvec} will be normalized to sum to one, emitting a warning if this 
#'     produces a meaningful change in covers. 
#'
#'@export
generate_initmat <- function(mod, pvec, nr, nc = nr) { 
  
  if ( any(is.na(pvec)) ) { 
    stop("NA in pvec are unsupported")
  }
  
  ns <- mod[["nstates"]]
  if ( length(pvec) != ns ) { 
    stop("the length of state covers does not match the number of states in the model definition")
  }
  
  if ( ! is.null(names(pvec)) ) { 
    diff <- setdiff(names(pvec), mod[["states"]])
    if ( length(diff) > 0 ) { 
      stop("State names in 'pvec' do not match the states defined in the model")
    }
    # Make sure order is right 
    pvec <- pvec[mod[["states"]]]
  }
  
  spvec <- sum(pvec)
  if ( sum(pvec) > 1.000001 | sum(pvec) < 0.99999 ) { 
    warning("The initial covers do not sum to one, they will be rescaled")
  }
  pvec <- pvec / spvec
  
  # Generate the matrix 
  m <- matrix(sample(mod[["states"]],
                     replace = TRUE, prob = pvec, size = nr*nc), 
              nrow = nr, ncol = nc)
  
  m <- factor(m, levels = mod[["states"]])
  dim(m) <- c(nr, nc)
  
  return(m)
}

#' @title Running a cellular automata
#' 
#' @description Run a pre-defined stochastic cellular automaton
#' 
#' @param mod A stochastic cellular automaton model defined using \code{\link{camodel}}
#' 
#' @param initmat An initial matrix to use for the simulation, possibly created using 
#'   \code{\link{generate_initmat}} 
#' 
#' @param niter The number of iterations for which to run the model 
#' 
#' @param control a named list with settings to alter how the simulation is run (see 
#'   full list of settings in secion 'Details')
#' 
#' @seealso camodel
#' 
#' @details 
#' 
#' \code{run_camodel()} is the function to run a pre-defined cellular automaton. It loads 
#'   the model definition, and runs the simulation for a pre-defined number of 
#'   iterations. It will run the simulation for \code{niter} iterations, starting from 
#'   the initial landscape matrix \code{initmat}.
#' 
#' The \code{control} list must have named elements, and allows altering the 
#'   way the simulation is run. The possible options are the following: 
#'  
#'  \enumerate{ 
#'    \item \code{substeps} The number of substeps within each iteration. Probabilities 
#'      defined in the model definition may go above one for some models, which produces 
#'      results that are approximate. To avoid this problem, each time step of the 
#'      model can be divided into several 'substeps'. Transitions between cell 
#'      probabilities may occur between each substep, but with a probability divided 
#'      by the number of substeps. This makes sure that every probability evaluated 
#'      during the model run is below one. Default is 1.
#'    
#'    \item \code{save_covers_every} The period of time between which the global covers 
#'      of each state in the landscape is saved. Set to zero to turn off saving 
#'      them. Default is 1 (save covers at each iteration).
#'    
#'    \item \code{save_snapshots_every} Period of time between which a snapshot of the
#'      2D grid is saved. Set to zero to turn off the saving of snapshots (default
#'      option).
#'     
#'    \item \code{console_callback_every} Sets the number of iterations between which 
#'      progress report is printed on the console. Set to zero to turn off progress 
#'      report. Default is to print progress every ten iterations.
#'    
#'    \item \code{neighbors} The number of neighbors to use. Use the integer value 4 to 
#'      use a four-way (von-Neumann) neighborhood, or 8 for an 8-way (Moore) 
#'      neighborhood. Any other values will produce an error. Default is to use 4
#'      neighbors.
#'    
#'    \item \code{wrap} Set to \code{TRUE} to use a toric space in which edges wrap 
#'      around, i.e. that cells on a side of the landscape are considered neighbors of
#'      cells on the other side. 
#'    
#'    \item \code{engine} The engine to use to run the simulations. Accepted values 
#'      are 'r', to use the pure-R engine, 'cpp' to use the C++ engine, or 'compiled', to
#'      compile the model code on the fly. Default is to use the C++ engine. Note that 
#'      the 'compiled' engine uses its own random number generator, and for this reason
#'      may produce results that are different from the two other engines.
#'    
#'    \item \code{olevel} (Compiled engine only) The optimization level to use when
#'      compiling the model code (default, O2, O3 or Ofast). This requires compiling with 
#'      gcc. By default, \code{\link[Rcpp]{sourceCpp}} options are used 
#'      (option 'default').
#'    
#'    \item \code{unroll_loops} (Compiled engine only) Set to \code{TRUE} to unroll loops 
#'      or not when compiling the model. Default is \code{FALSE}. Requires compiling with
#'      gcc.  
#'    
#'    \item \code{precompute_probas} (Compiled engine only) Set to \code{TRUE} to 
#'      precompute probabilities of transitions for all possible combinations of 
#'      neighborhood. When working with a model with a low number of states, this 
#'      can increase simulation speed dramatically. Default is to do it when the 
#'      possible total number of neighborhood combinations is below the number of cells 
#'      in the landscape (so that it is effectively less work to pre-compute transitions 
#'      than to compute them at each iteration). 
#'    
#'    \item \code{verbose_compilation} (Compiled engine only) Set to \code{TRUE} to print 
#'      Rcpp messages when compiling the model. Default is \code{FALSE}. 
#'    
#' }
#' 
#' @return A \code{ca_model_result} objects, which is a list with the following 
#'   components: 
#' 
#' \enumerate{ 
#'   
#'   \item \code{model} The original model used for the model run 
#'   
#'   \item \code{initmat} The initial landscape (matrix) used for the model run 
#'   
#'   \item \code{niter} The number of iterations used 
#'   
#'   \item \code{control} The control list used for the model run, containing the options
#'     used for the run 
#'   
#'   \item \code{output} A named list containing the model outputs. The 'covers' 
#'     component contains a matrix with the first column containing the time step, and the 
#'     other columns the proportions of cells in a given state. The 'snapshots' componentµ
#'     contains the landscapes recorded as matrices, with a 't' attribute indicating 
#'     the corresponding time step of the model run. 
#' }
#' 
#' @examples 
#' 
#' # Run a model with default parameters
#' mod <- ca_library("forestgap")
#' im  <- generate_initmat(mod, c(0.4, 0.6), nr = 100, nc = 50)
#' run_camodel(mod, im, niter = 100) 
#' 
#' # Set some options and use the compiled engine
#' ctrl <- list(engine = "compiled", save_covers_every = 1, save_snapshots_every = 100)
#' run <- run_camodel(mod, im, niter = 200)
#' 
#' covers <- run[["output"]][["covers"]]
#' matplot(covers[ ,1], covers[ ,-1], type = "l")
#' 
#'@export
run_camodel <- function(mod, initmat, niter, 
                        control = list()) { 
  #TODO: use column names in the output$covers matrix 
  # Pack transition coefficients into 3D array that RcppArmadillo understands
  ns <- mod[["nstates"]]
  states <- mod[["states"]]
  
  if ( ! all(levels(initmat) %in% states) ) { 
    stop("States in the initial matrix do not match the model states")
  }
  # Read parameters
  control <- load_control_list(control)
  
  
  # Check that there is no NA in transition rates 
  #TODO:
#   if ( any(is.na(tests)) ) { 
#     stop("NAs in computed coefficients, please make sure your model definition is correct")
#   }
  
  # NOTE: callbacks defined below will modify things in the currenct environment
  
  # Handle cover-storage callback 
  cover_callback <- function(t, ps, n) { }  
  cover_callback_active <- is_positive(control[["save_covers_every"]])
  if ( cover_callback_active ) { 
    nl <- 1 + niter %/% control[["save_covers_every"]]
    global_covers <- matrix(NA_real_, ncol = 1+ns, nrow = nl)
    colnames(global_covers) <- c("t", as.character(states))
    cur_line <- 1
    
    cover_callback <- function(t, ps, n) { 
      global_covers[cur_line, ] <<- c(t, ps / n)
      cur_line <<- cur_line + 1  
    }
  }
  
  # Handle snapshot-storage callback
  snapshot_callback <- function(t, m) { } 
  snapshot_callback_active <- is_positive(control[["save_snapshots_every"]])
  if ( snapshot_callback_active ) { 
    snapshots <- list()
    
    snapshot_callback <- function(t, m) { 
      d <- dim(m)
      m <- factor(states[m+1], levels = states)
      dim(m) <- d
      attr(m, "t") <- t
      snapshots <<- c(snapshots, list(m))
    }
  }
  
  # Handle console output callback 
  console_callback <- function(t, ps, n) { } 
  console_callback_active <- is_positive(control[["console_output_every"]])
  if ( console_callback_active ) { 
    last_time <- proc.time()["elapsed"]
    first_time <- last_time
    last_iter <- 1 
    
    console_callback <- function(t, ps, n) { 
      new_time <- proc.time()["elapsed"]
      iter_per_s <- (t - last_iter) / (new_time - last_time) 
      
      cover_string <- paste(seq.int(length(ps)), format(ps/n, digits = 2), 
                            sep = ":", collapse = " ")
      perc <- paste0(round(100 * (t / niter)), " %")
      speed <- ifelse(is.nan(iter_per_s) | iter_per_s < 0, "", 
                      paste("[", sprintf("%0.2f", iter_per_s), " iter/s]", sep = ""))
      
      tstring <- sprintf("t = %03i", t)
      cat(paste0(tstring, " (", perc, ") ", cover_string, " ", speed, "\n"))
      last_time <<- new_time
      last_iter <<- t 
    }
  }
  
  
  # Convert initmat to internal representation 
  d <- dim(initmat)
  initmat <- as.integer(initmat) - 1
  dim(initmat) <- d 
  
  # Convert transition matrices to internal representation (and to actual matrices)
  fix <- function(x) as.integer(factor(as.character(x), levels = states)) - 1 
  alpha <- mod[["alpha"]]
  pmat  <- mod[["pmat"]]
  qmat  <- mod[["qmat"]]
  for ( c in c("from", "to") ) { 
    pmat[ ,c]  <- fix(pmat[ ,c])
    qmat[ ,c]  <- fix(qmat[ ,c])
    alpha[ ,c] <- fix(alpha[ ,c])
  }
  pmat[ ,"state"] <- fix(pmat[ ,"state"])
  qmat[ ,"state"] <- fix(qmat[ ,"state"])
  alpha <- as.matrix(alpha)
  pmat <- as.matrix(pmat)
  qmat <- as.matrix(qmat)
  
  # Subset matrices for speed 
  pmat <- pmat[ abs(pmat[ ,"coef"]) > mod[["epsilon"]], , drop = FALSE]
  qmat <- qmat[ abs(qmat[ ,"ys"])   > mod[["epsilon"]], , drop = FALSE]
  
  # Take into account substeps 
  qmat[ ,"ys"] <- qmat[ ,"ys"] / control[["substeps"]]
  pmat[ ,"coef"] <- pmat[ ,"coef"] / control[["substeps"]]
  alpha[ ,"a0"] <- alpha[ ,"a0"] / control[["substeps"]]
  
  # TODO: this is ugly
  control_list <- list(substeps = control[["substeps"]], 
                       wrap     = mod[["wrap"]], 
                       init     = initmat, 
                       niter    = niter, 
                       nstates  = ns, 
                       neighbors = mod[["neighbors"]], 
                       xpoints  = mod[["xpoints"]], 
                       console_callback_active  = console_callback_active, 
                       console_callback_every   = control[["console_output_every"]], 
                       cover_callback_active    = cover_callback_active, 
                       cover_callback_every     = control[["save_covers_every"]], 
                       snapshot_callback_active = snapshot_callback_active, 
                       snapshot_callback_every  = control[["save_snapshots_every"]], 
                       olevel = control[["olevel"]], 
                       precompute_probas = control[["precompute_probas"]], 
                       unroll_loops = control[["unroll_loops"]], 
                       verbose_compilation = control[["verbose_compilation"]])
  
  engine <- control[["engine"]][1]
  if ( tolower(engine) == "r" ) { 
    camodel_r_engine(alpha, pmat, qmat, control_list, 
                     console_callback, cover_callback, snapshot_callback)
  } else if ( tolower(engine) %in% c("cpp", "c++") ) { 
    camodel_cpp_engine_wrap(alpha, pmat, qmat, control_list, 
                            console_callback, cover_callback, snapshot_callback)
  } else if ( tolower(engine) %in% c("compiled") ) { 
    camodel_compiled_engine_wrap(alpha, pmat, qmat, control_list, 
                                 console_callback, cover_callback, snapshot_callback)
  } else { 
    stop(sprintf("%s is an unknown CA engine", engine))
  }
  
  # Store artefacts and return result 
  results <- list(model = mod, 
                  initmat = initmat, 
                  niter = niter, 
                  control = control, 
                  output = list())
  
  if ( cover_callback_active ) { 
    results[["output"]][["covers"]] <- global_covers
  }
  
  if ( snapshot_callback_active ) { 
    results[["output"]][["snapshots"]] <- snapshots
  }
  
  class(results) <- list("ca_model_result", "list")
  return(results)
}


load_control_list <- function(l) { 
  
  control_list <- list(
    substeps = 1, 
    save_covers_every = 1, 
    save_snapshots_every = 0, 
    console_output_every = 10, 
    engine = "cpp", 
    # Compiled engine options
    olevel = "default", 
    unroll_loops = FALSE, 
    precompute_probas = "auto", 
    verbose_compilation = FALSE 
  )
  
  for ( nm in names(l) ) { 
    if ( nm %in% names(control_list) ) { 
      control_list[[nm]] <- l[[nm]] 
    } else { 
      stop(sprintf("%s is not a CA model option", nm))
    }
  }
  
  check_length1_integer(control_list[["substeps"]], "substeps", 1)
  check_length1_integer(control_list[["save_covers_every"]], "save_covers_every", 0)
  check_length1_integer(control_list[["save_snapshots_every"]], "save_snapshots_every", 0)
  check_length1_integer(control_list[["console_output_every"]], "console_output_every", 0)
  
  if ( ! control_list[["engine"]] %in% c("cpp", "compiled", "r") ) { 
    stop(sprintf("Engine must be one of 'cpp', 'compiled' or 'r'"))
  }
  
  if ( ! control_list[["olevel"]] %in% c("O0", "O1", "O2", "O3", "Ofast", "default") ) { 
    stop("'olevel' option must be one of 'default', 'O0', 'O1', 'O2', 'O3' or 'Ofast'")
  }
  
  if ( ! is.logical(control_list[["unroll_loops"]]) ) { 
    stop("'unroll_loops' option must be TRUE or FALSE")
  }
  
  if ( ! is.logical(control_list[["verbose_compilation"]]) ) { 
    stop("'verbose_compilation' option must be TRUE or FALSE")
  }
  
  if ( ! ( is.logical(control_list[["precompute_probas"]]) || 
          control_list[["precompute_probas"]] == "auto") ) { 
    stop("precompute probas must be TRUE, FALSE or 'auto'")
  }
  
  control_list
}

check_length1_integer <- function(x, str, minx = 0) { 
  err <- FALSE
  if ( is.null(x) || is.na(x) || length(x) != 1 ) { 
    err <- TRUE
  } else if ( x < minx ) { 
    err <- TRUE
  }
  if ( err ) { 
    msg <- sprintf("Option '%s' must be an integer >= %s", str, minx)
    stop(msg)
  }
  invisible(TRUE)
}

is_positive <- function(x) { 
  if ( is.null(x) || is.na(x) || length(x) == 0 ) { 
    return(FALSE)
  }
  if ( length(x) == 1 && x >= 1 ) { 
    return(TRUE)
  }
  # zero
  return(FALSE)
}
