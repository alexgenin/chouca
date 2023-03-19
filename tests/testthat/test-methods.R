# 
# 
# 

mod <- ca_library("aridvege") 
init <- generate_initmat(mod, rep(1/3, 3), nr = 24)
run <- run_camodel(mod, init, times = seq(1, 64, by = 3), 
                   control = list(console_output_every = 0))

test_that("camodel methods are OK", { 

  # Just run the plotting to make sure the code runs
  expect_true({ 
    plot(mod)
    TRUE
  })
  
  expect_true({ 
    summary(init)
    image(init)
    TRUE
  })
  
})

test_that("camodel run result methods are OK", { 
  
  # Just run the plotting to make sure the code runs
  expect_true({ 
    plot(run)
    image(run)
    TRUE
  })
  
  # Handle errors 
  expect_error({ 
    image(run, snapshot = "fdsq")
  })
  
  # 
  expect_warning({ 
    image(run, snapshot_time = 56)
  })
  
  # Handle the no-saved landscape case
  run_noimg <- run_camodel(mod, init, times = seq(1, 64), 
                     control = list(save_snapshots_every = 0, 
                                    console_output_every = 0))
  expect_error({ 
    image(run_noimg)
  })
  
  
  
})
