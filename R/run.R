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
  
  spvec <- sum(pvec)
  if ( sum(pvec) > 1 ) { 
    warning("The sum of initial covers are above one, they will be rescaled")
  }
  pvec <- pvec / spvec
  
  # Generate the matrix 
  m <- matrix(sample(seq.int(ns), replace = TRUE, prob = pvec, size = nr*nc), 
              nrow = nr, ncol = nc)
  
  return(m)
}

run_camodel <- function(mod, initmat, niter, 
                        control = list()) { 
  
  # Pack transition coefficients into 3D array that RcppArmadillo understands
  ns <- mod[["nstates"]]
  transitions <- ldply(mod[["transitions"]], function(o) { 
    data.frame(from = o[["from"]], to = o[["to"]], 
               vec = c(o[["X0"]], o[["XP"]], o[["XQ"]]))
  })
  ncoefs <- 1 + ns + ns
  transpack <- array(0, dim = list(ncoefs, ns, ns), 
                     dimnames = list(paste0("coef", seq.int(ncoefs)), 
                                     paste0("from", seq.int(ns)), 
                                     paste0("to",   seq.int(ns))))
  for ( cfrom in seq.int(ns) ) { 
    for ( cto in seq.int(ns) ) { 
      dat <- subset(transitions, from == cfrom & to == cto) 
      if ( nrow(dat) > 0 ) { 
        transpack[ , cfrom, cto] <- dat[ ,"vec"]
      }
    }
  }
  
  # Read parameters
  control <- load_control_list(control)
  
  # Handle cover-storage callback 
  cover_callback <- function(t, ...) { }  
  cover_callback_active <- control[["save_covers"]] 
  if ( cover_callback_active ) { 
    nl <- niter %/% control[["save_covers_every"]]
    global_covers <- matrix(NA_real_, ncol = 1+ns, nrow = nl)
    cur_line <- 1
    cover_callback <- function(t, ps) { 
      global_covers[cur_line, ] <<- c(t, ps)
      cur_line <<- cur_line + 1  
    }
  }
  
  # Handle snapshot-storage callback
  snapshot_callback <- function(t, ...) { } 
  snapshot_callback_active <- control[["save_snapshots"]]
  if ( snapshot_callback_active ) { 
    snapshots <- list()
    snapshot_callback <- function(t, m) { 
      attr(m, "t") <- t
      snapshots <<- c(snapshots, list(m))
    }
  }
  
  # Handle console output callback 
  console_callback <- function(t, ps) { } 
  console_callback_active <- control[["console_output"]]
  if ( console_callback_active ) { 
    last_time <- proc.time()["elapsed"]
    first_time <- last_time
    last_iter <- 1 
    
    console_callback <- function(t, ps) { 
      new_time <- proc.time()["elapsed"]
      iter_per_s <- (t - last_iter) / (new_time - last_time) 
      
      cover_string <- paste(seq.int(length(ps)), format(ps, digits = 2), 
                            sep = ":", collapse = " ")
      perc <- paste0(round(100 * (t / niter)), " %")
      speed <- paste(format(iter_per_s, digits = 2), " iter/s", sep = "")
      tstring <- sprintf("t = %s", t - 1)
      cat(paste0(tstring, " (", perc, ") ", cover_string, " [", speed, "]", "\n"))
      last_time <<- new_time
      last_iter <<- t 
    }
  }
  
  cover_callback_every <- control[["save_covers_every"]]
  snapshot_callback_every <- control[["save_snapshots_every"]]
  console_callback_every <- control[["console_output_every"]]
  
  control_list <- list(substeps = control[["substeps"]], 
                       wrap     = control[["wrap"]], 
                       console_callback_active  = console_callback_active, 
                       console_callback_every   = console_callback_every, 
                       cover_callback_active    = cover_callback_active, 
                       cover_callback_every     = cover_callback_every, 
                       snapshot_callback_active = snapshot_callback_active, 
                       snapshot_callback_every  = snapshot_callback_every)
                       
  # NOTE: this function can modify some objects in the current environment... 
  camodel_r_engine(transpack, initmat, niter, ns, control_list, 
                   console_callback, cover_callback, snapshot_callback)
  
  # Store artefacts and return result 
  results <- list(model = mod, 
                  initmat = initmat, 
                  niter = niter, 
                  control = control, 
                  output = list())
  
  if ( cover_callback_active ) { 
    results[["output"]][["global_covers"]] <- global_covers
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
    save_covers = TRUE, 
    save_covers_every = 1, 
    save_snapshots = FALSE, 
    save_snapshots_every = 1, 
    console_output = TRUE, 
    console_output_every = 10, 
    wrap = TRUE
  )
  
  for ( nm in names(l) ) { 
    if ( nm %in% names(control_list) ) { 
      control_list[[nm]] <- l[[nm]] 
    } else { 
      stop(sprintf("%s is not a CA model option", nm))
    }
  }
  
  control_list
}
