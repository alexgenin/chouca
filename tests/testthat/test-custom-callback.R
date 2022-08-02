# 
# This tests whether custom callbacks are properly supported 
# 


test_that("Custom callbacks work", { 
  
  # Kubo's forest gap model 
  nrows <- ncols <- 16
  mod <- ca_library("rockpaperscissor", neighbors = 8, wrap = TRUE)
  initmm <- generate_initmat(mod, rep(1/3, 3), nrows, ncols)
  niters <- 64
  
  ccb <- function(t, mat) { 
    data.frame(t = t, 
               cover = mean(mat == "r"), 
               sd  = sd(mat))
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
    all(custom_results_r == custom_results_cpp)
  })
  
  # Test compiled engine against cpp engine. It works because the model is deterministic
  # so the values are exactly equal. We only check that spatial autocorrelations are 
  # equal. 
  ctrl <- list(engine = "cpp", 
                console_output_every = 0, 
                custom_output_every = 1, 
                custom_output_fun = ccb)
  
  ctrl[["engine"]] <- "cpp"
  out_cpp <- run_camodel(mod, initmm, niters, ctrl)[["output"]][["custom_output"]]
  out_cpp <- out_cpp[[length(out_cpp)]][ ,"sd"]
  
  ctrl[["engine"]] <- "compiled"
  out_compiled <- run_camodel(mod, initmm, niters, ctrl)[["output"]][["custom_output"]]
  out_compiled <- out_compiled[[length(out_compiled)]][ ,"sd"]
  
  expect_true( abs(out_cpp - out_compiled) < 1e-10 ) 
  
})
