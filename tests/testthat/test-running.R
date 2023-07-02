# 
# Test some runtime options 
# 

test_that("Test console output", { 
  
  for ( engine in c("cpp", "compiled") ) { 
    
    model <- ca_library("aridvege")
    im <- generate_initmat(model, rep(1/3, 3), 21, 55)
    
    ctrl <- list(console_output_every = 1, engine = engine)
    string <- capture.output({ 
      run_camodel(model, im, seq(0, 6), ctrl)
    })
    
    expect_true({ 
      grepl("iter = 3 ( 50 %", string[4], fixed = TRUE)
    })
  }
  
})


test_that("Parameter errors are handled", { 
  
  model <- ca_library("coralreef")
  im <- generate_initmat(model, rep(1/3, 3), 21, 55)
  
  ctrl <- list(console_output_every = 0, engine = "cpp")
  
  expect_error(local({ 
    ctrl[["engine"]] <- "fdsqfdsq" 
    run_camodel(model, im, seq(0, 6), ctrl)
  }))
  
  expect_error(local({ 
    ctrl[["engine"]] <- "compiled" 
    ctrl[["precompute_probas"]] <- "dsqfdsq" 
    run_camodel(model, im, seq(0, 6), ctrl)
  }))
  
  expect_error(local({ 
    ctrl[["custom_output_every"]] <- 1 # no function associated
    run_camodel(model, im, seq(0, 6), ctrl)
  }))
  
  expect_error(local({ 
    ctrl[["console_output_every"]] <- NA
    run_camodel(model, im, seq(0, 6), ctrl)
  }))
  
  
  
  
})
