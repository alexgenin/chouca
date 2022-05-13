# 
# This file contains function that will run the model 
#

generate_initmat <- function(mod, pvec, nr, nc) { 
  
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

run_camodel <- function(mod, initmat, niter, 
                        control = list()) { 
  
  # Pack transition coefficients into 3D array that RcppArmadillo understands
  ns <- mod[["nstates"]]
  states <- mod[["states"]]
  
  if ( ! all(levels(initmat) %in% states) ) { 
    stop("States in the initial matrix do not match the model states")
  }
  
  # TODO: move this to model definition, so it is not done every time we run the model 
  transitions <- do.call(rbind, lapply(mod[["transitions"]], function(o) { 
    data.frame(from = o[["from"]], to = o[["to"]], 
               vec = c(o[["X0"]], o[["XP"]], o[["XQ"]], o[["XPSQ"]], o[["XQSQ"]]))
  }))
  ncoefs <- 1 + ns + ns + ns + ns # X0+XP+XQ+XPSQ+XQSQ
  transpack <- array(0, 
                     dim = list(ncoefs, ns, ns), 
                     dimnames = list(paste0("coef", seq.int(ncoefs)), 
                                     paste0("to", states), 
                                     paste0("from", states)))
  for ( cfrom in states ) { 
    for ( cto in states ) { 
      dat <- transitions[transitions[ ,"from"] == cfrom & transitions[ ,"to"] == cto, ]
      col_from <- which(states == cfrom)
      col_to   <- which(states == cto)
      if ( nrow(dat) > 0 ) { 
        transpack[ , col_to, col_from] <- dat[ ,"vec"]
      }
    }
  }
  
  # Read parameters
  control <- load_control_list(control)
  
  # NOTE: callbacks defined below will modify things in the currenct environment
  
  # Handle cover-storage callback 
  cover_callback <- function(t, ps, n) { }  
  cover_callback_active <- is_positive(control[["save_covers_every"]])
  if ( cover_callback_active ) { 
    nl <- 1 + niter %/% control[["save_covers_every"]]
    global_covers <- matrix(NA_real_, ncol = 1+ns, nrow = nl)
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
  
  # TODO: this is ugly
  control_list <- list(substeps = control[["substeps"]], 
                       wrap     = control[["wrap"]], 
                       init     = initmat, 
                       niter    = niter, 
                       nstates  = ns, 
                       neighbors = control[["neighbors"]], 
                       console_callback_active  = console_callback_active, 
                       console_callback_every   = control[["console_output_every"]], 
                       cover_callback_active    = cover_callback_active, 
                       cover_callback_every     = control[["save_covers_every"]], 
                       snapshot_callback_active = snapshot_callback_active, 
                       snapshot_callback_every  = control[["save_snapshots_every"]], 
                       olevel = control[["olevel"]], 
                       unroll_loops = control[["unroll_loops"]], 
                       verbose_compilation = control[["verbose_compilation"]])
  
  engine <- control[["engine"]][1]
  if ( tolower(engine) == "r" ) { 
    camodel_r_engine(transpack, control_list, 
                     console_callback, cover_callback, snapshot_callback)
  } else if ( tolower(engine) %in% c("cpp", "c++") ) { 
    camodel_cpp_engine(transpack, control_list, 
                       console_callback, cover_callback, snapshot_callback)
  } else if ( tolower(engine) %in% c("compiled") ) { 
    camodel_compiled_engine(transpack, control_list, 
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
    neighbors = 4, 
    wrap = TRUE, 
    engine = "cpp", 
    # Compiled engine options
    olevel = "default", 
    unroll_loops = FALSE, 
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
  
  
  if ( ! ( identical(control_list[["neighbors"]], 8) || 
           identical(control_list[["neighbors"]], 4) ) ) { 
    stop("The 'neighbors' control option must be 8 or 4.")
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
