# 
# PCA engine in pure R. Mainly here for testing/prototyping purposes
# 

camodel_r_engine <- function(alpha, pmat, qmat, ctrl, 
                             console_callback, cover_callback, snapshot_callback) { 
  
  # Unwrap elements of the ctrl list 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["neighbors"]] == 8 # whether we use 8 neighbors
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  
  xpoints <- ctrl[["xpoints"]] 
  
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
  delta_ps <- rep(0, ns)
  
  # Discard useless coefficients 
  # TODO: move that in run() or parse() and use epsilon
  # NOTE: this may remove all lines -> it's ok, the for loop will be skipped down there
  qmat <- qmat[ qmat[ ,"ys"] > 1e-8, , drop = FALSE]
  pmat <- pmat[ pmat[ ,"coef"] > 1e-8, , drop = FALSE]
  
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
      
      for ( i in seq.int(nr) ) { 
        for ( j in seq.int(nc) ) { 
          
          this_cell_state <- omat[i, j] 
          
          # Compute local densities 
          qs <- local_dens(omat, ns, i, j, wrap, use_8_nb)
          
          # qpointn corresponds to the line in the lookup table at which we pick up the 
          # value for the local probability
          this_total_nb <- sum(qs)
          qpointn <- as.integer(qs / this_total_nb * (xpoints-1))
#           qs <- qs / sum(qs)
          # Compute probability of transition to other states
          trates <- rep(0, ns)
          
          for ( k in seq.int(nrow(pmat)) ) { 
            trates[ 1+pmat[k, "to"] ] <- trates[ 1+pmat[k, "to"] ] + 
              (pmat[k, "from"] == this_cell_state) * 
                pmat[k, "coef"] * (ps[ 1+pmat[k, "state"] ] / n)^pmat[k, "expo"]
          }
          
          for ( k in seq.int(nrow(qmat)) ) { 
            trates[ 1+qmat[k, "to"] ] <- trates[ 1+qmat[k, "to"] ] + 
              (qmat[k, "from"] == this_cell_state) * 
                ( qmat[k, "qs"] == qpointn[ 1+qmat[k, "state"] ] ) * 
                qmat[k, "ys"]
          }
          
          # This works because the number of states match
          trates[ 1+alpha[ , "to"] ] <- trates[ 1+alpha[ , "to"] ] + 
            (alpha[ , "from"] == this_cell_state) * alpha[ , "a0"] 
          
          # Take substep into account 
          trates <- trates / substeps
          
          if ( this_cell_state == 0 ) { # empty 
            to_tree <- with(kubo_mod$parms, alpha * (ps[2] / n)^1)
            if ( abs(to_tree - trates[2]) > 1e-8 ) browser()
          } else { # tree
            to_empty <- with(kubo_mod$parms, d + delta * (qs[1] / sum(qs)))
            if ( abs(to_empty - trates[1]) > 1e-8 ) browser()
          }
          
          # Compute rates of transitions (probabilities) to other states
          ctrates <- cumsum(trates)
          
          if ( max(ctrates) > 1 + sqrt(.Machine$double.eps) ) { 
            warning("Computed probabilities were above one, results will be approximate. Consider increasing control parameter 'substeps'")
          }
          
          # Flip a coin to see if the transition occurs
          makes_transition <- stats::runif(1) < ctrates
          if ( any(makes_transition) ) { 
            new_state <- which(makes_transition)[1] - 1 # adjust indexing for states
            old_state <- nmat[i, j]
            nmat[i, j] <- new_state
            
            # Adjust counts of cell states
            delta_ps[old_state+1] <- delta_ps[old_state+1] - 1 
            delta_ps[new_state+1] <- delta_ps[new_state+1] + 1 
          }
          
        }
          
      }
      
      # Apply changes
      omat <- nmat 
      ps <- ps + delta_ps
      delta_ps <- delta_ps * 0
      
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
