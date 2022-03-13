# 
# PCA engine in pure R. Mainly here for testing/prototyping purposes
# 

camodel_r_engine <- function(trans, init, niter, nstates, ctrl, 
                             console_callback, cover_callback, snapshot_callback) { 
  
  omat <- nmat <- init 
  
  nr <- nrow(init)
  nc <- ncol(init)
  ns <- nstates 
  substeps <- ctrl[["substeps"]]
  wrap <- ctrl[["wrap"]]
  
  # Compute global densities
  ps <- laply(seq.int(nstates), function(s) { 
    mean(omat == s)
  })
  
  for ( t in seq.int(niter) ) { 
    
    # Make callback to progress display 
    if ( ctrl[["console_callback_active"]] && 
         t %% ctrl[["console_callback_every"]] == 0 ) { 
      console_callback(t, ps)
    }
    
    # Make callback to store global densities 
    if ( ctrl[["cover_callback_active"]] && 
         t %% ctrl[["cover_callback_every"]] == 0 ) { 
      cover_callback(t, ps)
    }
    
    # Make callback to store snapshots 
    if ( ctrl[["snapshot_callback_active"]] && 
         t %% ctrl[["snapshot_callback_every"]] == 0 ) { 
      snapshot_callback(t, omat)
    }
    
    for ( s in seq.int(substeps) ) { 
      
      omat <- nmat 
      
      ps <- laply(seq.int(nstates), function(s) { 
        mean(omat == s)
      })
      
      for ( i in seq.int(nr) ) { 
        for ( j in seq.int(nc) ) { 
          # cat(i, " ", j, "\n")
            
          this_cell_state <- omat[i, j] 
          tprobs <- trans[ ,this_cell_state, ]
          
          # Compute local densities 
          qs <- local_dens(omat, ns, i, j, wrap, nb = 4)
          
          # Compute rates of transitions (probabilities) to other states
          trates <- numeric(ncol(tprobs))
          for ( c in seq.int(ncol(tprobs)) ) { 
            trates[c] <- ( tprobs[1, c] + 
                            sum(tprobs[2:(2+ns-1), c] * ps) + 
                            sum(tprobs[(2+ns):(2+ns+ns-1), c] * qs) ) / substeps
          }
          
          ctrates <- cumsum(trates)
          
          if ( max(ctrates) > 1 ) { 
            warning("Computed probabilities were above one, results will be approximate. Consider increasing control parameter 'substeps'")
          }
          
          # Flip a coin to see if the transition occurs
          makes_transition <- runif(1) < ctrates
          if ( any(makes_transition) ) { 
            new_state <- which(makes_transition)[1]
            nmat[i, j] <- new_state
          }
          
        }
          
      }
      
    }
  
  }
  
  return(nmat)
}

local_dens <- function(m, ns, i, j, wrap, nb) { 
  
  nr <- nrow(m)
  nc <- ncol(m)
  
  # Magic wrapping function 
  north <- south <- west <- east <- NA
  if ( wrap ) { 
    wf <- function(i, n) 1 + (n + i - 1) %% n
  }
  
  if ( wrap || i > 1) { 
    north <- m[wf(i-1, nr), j]
  }
  
  if ( wrap || i < nr) { 
    south <- m[wf(i+1, nr), j]
  }
  
  if ( wrap || j > 1) { 
    west <- m[i, wf(j-1, nc)]
  }
  
  if ( wrap || j < nc) { 
    east <- m[i, wf(j+1, nc)]
  }
  
  #TODO: use 8 neighbors
  
  # Counts 
  counts <- na.omit(c(south, north, east, west))
  
  local_dens <- laply(seq.int(ns), function(s) { 
    mean(counts == s)
  })
  
  return(local_dens)
}

