# 
# This file will make sure models produce spatial patterns
# 

test_that("autocorrelation in spatial patterns make sense", { 
  # We need raw_moran from spatialwarnings
  if ( requireNamespace("spatialwarnings", quietly = TRUE) ) { 
    
    nr <- 256
    nc <- 256
    ts <- seq(0, 1000)
    
    # This model should have spatial autocorrelation 
    mod_pos_autocor <- camodel(transition(from = "a", to = "b", ~ 0.08 * q["b"]), 
                               transition(from = "b", to = "a", ~ 0.04), 
                               all_states = c("a", "b"), 
                               neighbors = 8, 
                               wrap = TRUE)
    im <- generate_initmat(mod_pos_autocor, c(0.5, 0.5), nr, nc)
    
    # Run the thing, then make sure there is autocorrelation
    ctrl <- list(save_snapshots_every = 100, 
                 console_output_every = 0, 
                 engine = "compiled")
    run <- run_camodel(mod_pos_autocor, im, ts, control = ctrl)
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
    
    run <- run_camodel(mod_neg_autocor, im, ts, control = ctrl)
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
                                 all_states = c("a", "b"), 
                                 neighbors = 8, 
                                 wrap = TRUE)
      
      all_autocors <- replicate(99, { 
        run <- run_camodel(mod_neg_autocor, im, ts, control = ctrl)
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
})
