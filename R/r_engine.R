# 
# PCA engine in pure R. Mainly here for testing/prototyping purposes
# 

camodel_r_engine <- function(ctrl) { 
  
  # Unwrap elements of the ctrl list 
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["neighbors"]] == 8 # whether we use 8 neighbors
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  xpoints <- ctrl[["xpoints"]] 
  substeps <- ctrl[["substeps"]] 
  
  beta_0 <- ctrl[["beta_0"]]
  beta_q <- ctrl[["beta_q"]]
  beta_pp <- ctrl[["beta_pp"]]
  beta_pq <- ctrl[["beta_pq"]]
  beta_qq <- ctrl[["beta_qq"]]
  
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
    
  t <- 0 
  while ( t <= niter ) { 
    
    # Make callback to progress display 
    if ( ctrl[["console_output_every"]] > 0 && 
         t %% ctrl[["console_output_every"]] == 0 ) { 
      ctrl[["console_callback"]](t, ps, n)
    }
    
    # Make callback to store global densities 
    if ( ctrl[["save_covers_every"]] > 0 && 
         t %% ctrl[["save_covers_every"]] == 0 ) { 
      ctrl[["cover_callback"]](t, ps, n)
    }
    
    # Make callback to store snapshots 
    if ( ctrl[["save_snapshots_every"]] > 0 && 
         t %% ctrl[["save_snapshots_every"]] == 0 ) { 
      ctrl[["snapshot_callback"]](t, omat)
    }
    
    # Make custom callback
    if ( ctrl[["custom_output_every"]] > 0 && 
         t %% ctrl[["custom_output_every"]] == 0 ) { 
      ctrl[["custom_callback"]](t, omat)
    }
    
    for ( s in seq.int(substeps) ) { 
      for ( i in seq.int(nr) ) { 
        for ( j in seq.int(nc) ) { 
          
          this_cell_state <- omat[i, j] 
          
          # Compute local densities 
          qs <- local_dens(omat, ns, i, j, wrap, use_8_nb)
          
          # qpointn corresponds to the line in the lookup table at which we pick up the 
          # value for the local probability
          nbs <- sum(qs)
          qpointn <- as.integer(qs / nbs * (xpoints-1))
  #           qs <- qs / sum(qs)
          # Compute probability of transition to other states
          trates <- rep(0, ns)
          
          for ( k in unique(beta_0[ ,"to"]) ) { 
            sub <- which(beta_0[ ,"to"] == k)
            trates[k+1] <- trates[k+1] + sum( 
              (beta_0[sub, "from"] == this_cell_state) * beta_0[sub, "coef"] 
            )
          }
          
          for ( k in unique(beta_q[ ,"to"]) ) { 
            sub <- which(beta_q[ ,"to"] == k)
            trates[k+1] <- trates[k+1] + sum( 
              (beta_q[sub, "from"] == this_cell_state) * 
                ( beta_q[sub, "qs"] == qpointn[ 1+beta_q[sub, "state_1"] ] ) * 
                beta_q[sub, "coef"]
            )
          }
          
          for ( k in unique(beta_pp[ ,"to"]) ) { 
            sub <- which(beta_pp[ ,"to"] == k)
            trates[k+1] <- trates[k+1] + sum(
              (beta_pp[sub, "from"] == this_cell_state) * 
                beta_pp[sub, "coef"] * 
                ( ps[ 1 + beta_pp[sub, "state_1"] ] / n )^beta_pp[sub, "expo_1"] *
                ( ps[ 1 + beta_pp[sub, "state_2"] ] / n ) ^ beta_pp[sub, "expo_2"]
            )
          }
          
          for ( k in unique(beta_pq[ ,"to"]) ) { 
            sub <- which(beta_pq[ ,"to"] == k)
            trates[k+1] <- trates[k+1] + sum(
              (beta_pq[sub, "from"] == this_cell_state) * 
                beta_pq[sub, "coef"] * 
                ( ps[ 1 + beta_pq[sub, "state_1"] ] / n ) ^ beta_pq[sub, "expo_1"] *
                ( qs[ 1 + beta_pq[sub, "state_2"] ] / nbs ) ^ beta_pq[sub, "expo_2"]
            )
          }
          
          
          for ( k in unique(beta_qq[ ,"to"]) ) { 
            sub <- which(beta_qq[ ,"to"] == k)
            trates[k+1] <- trates[k+1] + sum(
              (beta_qq[sub, "from"] == this_cell_state) * 
                beta_qq[sub, "coef"] * 
                ( qs[ 1 + beta_qq[sub, "state_1"] ] / nbs ) ^ beta_qq[sub, "expo_1"] *
                ( qs[ 1 + beta_qq[sub, "state_2"] ] / nbs ) ^ beta_qq[sub, "expo_2"]
            )
          }
          
          # Compute rates of transitions (probabilities) to other states
          ctrates <- cumsum(trates)
          
          if ( max(ctrates) > 1 + sqrt(.Machine$double.eps) ) { 
            warning("Computed probabilities were above one, results will be approximate.")
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
