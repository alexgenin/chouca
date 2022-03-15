# 
# This tests engines against each other to make sure results are reproduced
# 

# Kubo's forest gap model 
mod <- forestgap()

nr <- 10
nc <- 10
initmat <- generate_initmat(mod, c(0.5, 0.5), nr, nc)

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
init <- generate_initmat(mod, c(0.1, 0.9), 10, 10)
engines_ts <- replicate(499, { 
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

expect_true({ 
  mean( emeans[ ,1] - emeans[ ,3] ) < 1e-2
})

expect_true({ 
  mean( emeans[ ,2] - emeans[ ,4] ) < 1e-2
})

expect_true({ 
  mean( evars[ ,1] - evars[ ,3] ) < 1e-2
})

expect_true({ 
  mean( evars[ ,2] - evars[ ,4] ) < 1e-2
})


