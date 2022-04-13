# 
# This tests engines against each other to make sure results are reproduced
# 


# Generate datasets 

# Kubo's forest gap model 
nr <- 20
nc <- 30

models <- list()

# Kubo's forestgap model 
model <- forestgap()
imat <- generate_initmat(model, c(0.5, 0.5), nr, nc)
models[[1]] <- list(model, imat)

# Guichard's musselbed model
model <- musselbed()
imat <- generate_initmat(model, c(0.5, 0.1, 0.4), nr, nc)
models[[2]] <- list(model, imat)

# Model with non-zero XPSQ/XQSQ
model <- camodel(transition(from = "A", to = "B", ~ r * ( 1 + p["A"]^2 + q["A"]^2 )), 
                 transition(from = "B", to = "A", ~ r * ( 0.1 + p["B"]^2 + q["B"]^2 )), 
                 parms = list(r = 0.1))
imat <- generate_initmat(model, c(0.4, 0.6), nr, nc)
models[[3]] <- list(model, imat)

# If we do not do extended tests, just do the musselbed model and the non-zero XPSQ/XQSQ
# model
if ( ! exists("EXTENDED_TESTS") || ( ! EXTENDED_TESTS ) ) { 
  models <- models[c(2, 3)]
}

plyr::llply(models, function(modinfo) { 
  mod <- modinfo[[1]]
  initmat <- modinfo[[2]]
  
  control <- list(substeps = 1, 
                  console_output = FALSE, 
                  console_output_every = 20, 
                  save_covers = TRUE, 
                  save_covers_every = 1, 
                  save_snapshots = TRUE, 
                  save_snapshots_every = 1000, 
                  ca_engine = "cpp")
  
  # Check that we reproduce well the variance and mean of time series between the two 
  # engines. Somehow setting the seed does not 
  engines_ts <- replicate(199, { 
    modcompiled <- run_camodel(mod, initmat, 20, control = { 
      control[["ca_engine"]] <- "compiled" ; control
    })
    modcpp <- run_camodel(mod, initmat, 20, control = { 
      control[["ca_engine"]] <- "cpp" ; control
    })
    modr <- run_camodel(mod, initmat, 20, control = { 
      control[["ca_engine"]] <- "r" ; control
    })
    
    # Time series 
    cbind(modcompiled[["output"]][["covers"]][ ,2:3], 
          modcpp[["output"]][["covers"]][ ,2:3], 
          modr[["output"]][["covers"]][ ,2:3])
  })
  
  emeans <- apply(engines_ts, c(1, 2), mean)
  evars <- apply(engines_ts, c(1, 2), var)

  # Display results 
  # plot(emeans[ ,5], type = "n")
  # lines(emeans[ ,1], col = "red")
  # lines(emeans[ ,3], col = "black")
  # lines(emeans[ ,5], col = "green")

  
  expect_true({ 
    all( abs( emeans[ ,1] - emeans[ ,3] ) < 1e-2 )
  })

  expect_true({ 
    all( abs( emeans[ ,2] - emeans[ ,4] ) < 1e-2 )
  })

  expect_true({ 
    all( abs( evars[ ,1] - evars[ ,3] ) < 1e-2 )
  })

  expect_true({ 
    all( abs( evars[ ,2] - evars[ ,4] ) < 1e-2 )
  })

})

