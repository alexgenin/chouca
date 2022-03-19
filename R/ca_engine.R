# 
# PCA engine in pure R. Mainly here for testing/prototyping purposes
# 

camodel_r_engine <- function(trans, ctrl, 
                             console_callback, cover_callback, snapshot_callback) { 
  
  # Unwrap elements of the ctrl list 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["use_8_neighbors"]] 
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  # Initialize some elements
  # NOTE: we work in integer representation minus one internally, because it makes 
  # way more sense than handling R factors. This means that states (a,b,c) are 
  # represented as (0,1,2)
  nr <- nrow(init)
  nc <- ncol(init)
  
  omat <- nmat <- init
  n <- nr * nc 
  
  # Compute global densities
  ps <- get_global_counts(omat, ns) 
  
  t <- 0 
  while ( t <= niter ) { 
    
    
    # Make callback to progress display 
    if ( ctrl[["console_callback_active"]] && 
         t %% ctrl[["console_callback_every"]] == 0 ) { 
      console_callback(t, ps, n)
    }
    
    # Make callback to store global densities 
    if ( ctrl[["cover_callback_active"]] && 
         t %% ctrl[["cover_callback_every"]] == 0 ) { 
      cover_callback(t, ps, n)
    }
    
    # Make callback to store snapshots 
    if ( ctrl[["snapshot_callback_active"]] && 
         t %% ctrl[["snapshot_callback_every"]] == 0 ) { 
      snapshot_callback(t, omat)
    }
    
    for ( s in seq.int(substeps) ) { 
      
      omat <- nmat 
      
      for ( i in seq.int(nr) ) { 
        for ( j in seq.int(nc) ) { 
          
          this_cell_state <- omat[i, j] 
          tprobs <- trans[ , , this_cell_state + 1] 
          
          # Compute local densities 
          qs <- local_dens(omat, ns, i, j, wrap, use_8_nb)
          qs <- qs / sum(qs)
          
          # Compute rates of transitions (probabilities) to other states
          trates <- numeric(ncol(tprobs))
          for ( col in seq.int(ncol(tprobs)) ) { 
            trates[col] <- ( tprobs[1, col] + 
                            sum(tprobs[2:(2+ns-1), col] * ps/n) + 
                            sum(tprobs[(2+ns):(2+ns+ns-1), col] * qs) ) / substeps
          }
          ctrates <- cumsum(trates)
          
          if ( max(ctrates) > 1 ) { 
            warning("Computed probabilities were above one, results will be approximate. Consider increasing control parameter 'substeps'")
          }
          
          # Flip a coin to see if the transition occurs
          makes_transition <- runif(1) < ctrates
          if ( any(makes_transition) ) { 
            new_state <- which(makes_transition)[1] - 1 # adjust indexing for states
            old_state <- nmat[i, j]
            nmat[i, j] <- new_state
            
            # Adjust counts of cell states
            ps[old_state+1] <- ps[old_state+1] - 1 
            ps[new_state+1] <- ps[new_state+1] + 1 
          }
          
        }
          
      }
      
    }
    
    t <- t + 1 
  }
  
  return(nmat)
}

get_global_counts <- function(m, ns) { 
  do.call(c, lapply(seq.int(ns), function(s) { 
    sum(m == s-1)
  }))
}
