# 
# This file will make sure models produce spatial patterns
# 

# We need raw_moran from spatialwarnings
if ( requireNamespace("spatialwarnings") ) { 
  
  nr <- 50 
  nc <- 50 
  
  # This model should have spatial autocorrelation 
  mod_pos_autocor <- camodel(transition(from = "a", to = "b", ~ 0.04 * q["b"]), 
                             transition(from = "b", to = "a", ~ 0.01), 
                             all_states = c("a", "b"), 
                             neighbors = 8, 
                             wrap = TRUE)
  im <- generate_initmat(mod_pos_autocor, c(0.5, 0.5), nr, nc)
  
  # Run the thing, then make sure there is autocorrelation
  ctrl <- list(save_snapshots_every = 100, 
               console_output_every = 0)
  run <- run_camodel(mod_pos_autocor, im, 1000, control = ctrl)
  # spatialwarnings::display_matrix(run[["output"]][["snapshots"]])

  a <- lapply(run[["output"]][["snapshots"]], function(m) { 
    a <- matrix(m == "b", nrow = nrow(m), ncol = ncol(m))
    spatialwarnings::raw_moran(a)
  })
  
  # Extract three last matrices to have positive autocorrelation
  expect_true({ 
    all( as.numeric(tail(a, 3)) > 0 ) 
  })
  
  
  # Same model, but with negative autocorrelation
  mod_neg_autocor <- camodel(transition(from = "a", to = "b", ~ 0.01), 
                             transition(from = "b", to = "a", ~ 0.01 * q["b"]), 
                             all_states = c("a", "b"), 
                             neighbors = 8, 
                             wrap = TRUE)
  
  run <- run_camodel(mod_neg_autocor, im, 1000, control = ctrl)
  # spatialwarnings::display_matrix(run[["output"]][["snapshots"]])

  a <- lapply(run[["output"]][["snapshots"]], function(m) { 
    a <- matrix(m == "b", nrow = nrow(m), ncol = ncol(m))
    spatialwarnings::raw_moran(a) 
  })

  expect_true({ 
    all( as.numeric(tail(a, 3)) < 0 ) 
  })
  
  
  # Same model, but with no autocorrelation
  if ( exists("EXTENDED_TESTS") && EXTENDED_TESTS ) { 
    mod_neg_autocor <- camodel(transition(from = "a", to = "b", ~ 0.01), 
                              transition(from = "b", to = "a", ~ 0.01), 
                              all_states = c("a", "b"))
    all_autocors <- replicate(99, { 
      run <- run_camodel(mod_neg_autocor, im, 1000, control = ctrl)
      # spatialwarnings::display_matrix(run[["output"]][["snapshots"]])
      
      a <- lapply(run[["output"]][["snapshots"]], function(m) { 
        a <- matrix(m == "b", nrow = nrow(m), ncol = ncol(m))
        spatialwarnings::raw_moran(a)
      })
      mean(as.numeric(tail(a, 3)))
    })
    
    # We expect very weak autocorrelation on average, positive or negative
    expect_true({ 
      abs( mean(all_autocors) ) < 1e-1
    })
    
  }
  
  
  
}
