# 
# Test the custom callbacks
# 

mod <- ca_library("aridvege")
initmm <- generate_initmat(mod, c(0, 0, 1), 10, 10)

ctrl <- list(console_output_every = 0, 
             custom_output_every = 1, 
             custom_output_fun = function(t, mat) { 
  1
})

for ( engine in c("cpp", "r", "compiled") ) { 
  
  ctrl[["engine"]] <- engine
  out <- run_camodel(mod, initmm, 2, ctrl)
  
  expect_true({ 
    all( unlist(out[["output"]][["custom_output"]]) == 1 )
  })
  
  # Test that the plotting of stuff runs (just execute the code)
  ctrl[["custom_output_fun"]] <- landscape_plotter(mod)
  run_camodel(mod, initmm, 2, ctrl)
  
  ctrl[["custom_output_fun"]] <- trace_plotter(mod, initmm)
  run_camodel(mod, initmm, 2, ctrl)
}

