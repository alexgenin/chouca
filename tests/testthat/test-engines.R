# 
# This tests engines against each other to make sure results are reproduced
# 


# Generate datasets 


if ( ! ( exists("EXTENDED_TESTS") && EXTENDED_TESTS ) ) { 
  # Quick testing. 
  # We use a large tolerance, to make sure the test is not too conservative. We just 
  # catch blatant deviations in results, and do not take the risk of getting unlucky 
  # during a test and getting a false positive. These tests are primarily 
  # here to make sure there is no bad interaction on CRAN with other packages and on 
  # other systems.
  nr <- 20
  nc <- 30
  tolerance <- 5e-2
} else { 
  # Long testing. 
  # We use a large matrix size and a smaller tolerance, to make tests more conservative. 
  nr <- 48
  nc <- 64
  tolerance <- 1e-2
}


models <- list()

# Kubo's forestgap model 
model <- forestgap()
imat <- generate_initmat(model, c(0.5, 0.5), nr, nc)
models[[1]] <- list(model, imat)

# Model with non-zero XPSQ/XQSQ
model <- camodel(transition(from = "A", to = "B", ~ r * ( 1 + p["A"]^2 + q["A"]^2 )), 
                 transition(from = "B", to = "A", ~ r * ( 0.1 + p["B"]^2 + q["B"]^2 )), 
                 parms = list(r = 0.1))
imat <- generate_initmat(model, c(0.2, 0.8), nr, nc)
models[[2]] <- list(model, imat)

# If we do not do extended tests, just do the musselbed model and the non-zero XPSQ/XQSQ
# model
if ( ! exists("EXTENDED_TESTS") || ( ! EXTENDED_TESTS ) ) { 
  models <- models[1]
}

plyr::llply(models, function(modinfo) { 
  
  mod <- modinfo[[1]]
  initmat <- modinfo[[2]]
  
  control <- list(substeps = 1, 
                  console_output_every = 0, 
                  save_covers_every = 1, 
                  save_snapshots_every = 1000, 
                  engine = "cpp")
  
  # Check that we reproduce well the variance and mean of time series between the two 
  # engines. Somehow setti the seed does not 
  engines_ts <- replicate(99, { 
    modcompiled <- run_camodel(mod, initmat, 20, control = { 
      control[["engine"]] <- "compiled" ; control
    })
    modcpp <- run_camodel(mod, initmat, 20, control = { 
      control[["engine"]] <- "cpp" ; control
    })
    modr <- run_camodel(mod, initmat, 20, control = { 
      control[["engine"]] <- "r" ; control
    })
    
    # Time series 
    cbind(modcompiled[["output"]][["covers"]][ ,2:3], 
          modcpp[["output"]][["covers"]][ ,2:3], 
          modr[["output"]][["covers"]][ ,2:3])
  })
  
  emeans <- apply(engines_ts, c(1, 2), mean)
  evars <- apply(engines_ts, c(1, 2), var)
  
  par(mfrow = c(1, 2))
  # Display results 
  plot(emeans[ ,5], type = "n")
  lines(emeans[ ,1], col = "red")
  lines(emeans[ ,3], col = "black")
  lines(emeans[ ,5], col = "green")
  
  # Display variance results
  plot(evars[ ,5], type = "n")
  lines(evars[ ,1], col = "red")
  lines(evars[ ,3], col = "black")
  lines(evars[ ,5], col = "green")
  
  expect_true({ 
    all( abs( emeans[ ,1] - emeans[ ,3] ) < tolerance )
  })

  expect_true({ 
    all( abs( emeans[ ,2] - emeans[ ,4] ) < tolerance )
  })

  expect_true({ 
    all( abs( evars[ ,1] - evars[ ,3] ) < tolerance )
  })

  expect_true({ 
    all( abs( evars[ ,2] - evars[ ,4] ) < tolerance )
  })
  
})

