# 
# Make sure the compiled engine always returns the same results regardless of 
# whether we precompute the probabilities or not, and respects the R seed. 
# 

ncols <- 32
nrows <- ncols * (9/16)

# RPC model. We use it because it is deterministic
mod <- ca_library("rockpaperscissor")
initmm <- generate_initmat(mod, rep(1/3, 3), nrows, ncols)


test_that("Compiled model produces OK results regardless of proba precomputation", { 
  
  iters <- seq(0, 32)
  
  control <- list(save_covers_every = 1, 
                  console_output_every = 0, 
                  engine = "compiled", 
                  force_compilation = TRUE, 
                  precompute_probas = TRUE)
  
  o <- run_camodel(mod, initmm, iters, control = control)
  
  o2 <- run_camodel(mod, initmm, iters, 
                    control = { control[["precompute_probas"]] <- FALSE; control })
  
  par(mfrow = c(1, 2))
  ts1 <- o[["output"]][["covers"]]
  matplot(ts1[ ,1], ts1[ ,-1], type = "l")
  ts2 <- o2[["output"]][["covers"]]
  matplot(ts2[ ,1], ts2[ ,-1], type = "l")
  
  # This should give exactly the same result 
  expect_true({ 
    all( abs(ts1 - ts2) < 1e-8 )
  })
  
  
  
  # Make sure it works also if we do not wraparound
  mod <- update(mod, wrap = FALSE)

  o <- run_camodel(mod, initmm, iters, control = control)
  o2 <- run_camodel(mod, initmm, iters, 
                    control = { control[["precompute_probas"]] <- FALSE; control })
  
  par(mfrow = c(1, 2))
  ts1 <- o[["output"]][["covers"]]
  matplot(ts1[ ,1], ts1[ ,-1], type = "l")
  ts2 <- o2[["output"]][["covers"]]
  matplot(ts2[ ,1], ts2[ ,-1], type = "l")
  
  # This should give exactly the same result 
  expect_true({ 
    all( abs(ts1 - ts2) < 1e-8 )
  })

})


# Make sure things are printed when we use verbose compilation 
test_that("Verbose compilation prints something", { 
  
  control <- list(save_covers_every = 1, 
                  console_output_every = 0, 
                  engine = "compiled", 
                  precompute_probas = TRUE)
  
  mod <- update(mod, wrap = FALSE, parms = list(prob = 0.001))
  o <- capture.output({ 
    a <- run_camodel(mod, initmm, 2, 
                     control = c(control, list(verbose_compilation = TRUE)))
  })
  
  expect_true({ 
    any(grepl("Setting __SUBSTEPS__ to", o))
  })
  
})



