# 
# This tests whether custom callbacks are properly supported 
# 


test_that("Custom callbacks work", { 
  
  nrows <- ncols <- 16
  niters <- 64
  
  
  mod <- ca_library("rockpaperscissor", wrap = TRUE)
  initmm <- generate_initmat(mod, rep(1/3, 3), nrows, ncols)
  
  ccb <- function(t, mat) { 
    bmat <- matrix(mat == "VEGE", nrow = nrow(mat), ncol = ncol(mat))
    data.frame(t = t, 
               cover = mean(mat == "VEGE"), 
               sd  = sd(mat == "VEGE"))
  }
  
  ctrl <- list(engine = "cpp", 
               console_output_every = 0, 
               custom_output_every = 1, 
               custom_output_fun = ccb)
  out <- run_camodel(mod, initmm, niters, ctrl)
  custom_results_cpp <- plyr::rbind.fill(out[["output"]][["custom_output"]])
  
  ctrl[["engine"]] <- "r"
  out <- run_camodel(mod, initmm, niters, ctrl)
  custom_results_r <- plyr::rbind.fill(out[["output"]][["custom_output"]])
  
  expect_true({ 
    all( abs(custom_results_r - custom_results_cpp) < 0.1 )
  })
  
  
  
  # Test compiled engine against cpp engine. Here it works because the model is deterministic
  # so the values are exactly equal. We only check that spatial autocorrelations are 
  # equal. 
  mod <- ca_library("rockpaperscissor", neighbors = 8, wrap = TRUE)
  initmm <- generate_initmat(mod, rep(1/3, 3), nrows, ncols)
  
  ccb <- function(t, mat) { 
    bmat <- matrix(mat == "r", nrow = nrow(mat), ncol = ncol(mat))
    data.frame(t = t, 
               cover = mean(mat == "r"), 
               sd  = sd(mat == "r"))
  }
  
  ctrl <- list(engine = "cpp", 
               console_output_every = 0, 
               custom_output_every = 1, 
               custom_output_fun = ccb)
  
  ctrl[["engine"]] <- "cpp"
  out_cpp <- run_camodel(mod, initmm, niters, ctrl)[["output"]][["custom_output"]]
  out_cpp <- plyr::rbind.fill(out_cpp)
  
  ctrl[["engine"]] <- "compiled"
  out_compiled <- run_camodel(mod, initmm, niters, ctrl)[["output"]][["custom_output"]]
  out_compiled <- plyr::rbind.fill(out_compiled)
  
  expect_true( all( abs(out_cpp - out_compiled) < 1e-10 ) )
  
})
