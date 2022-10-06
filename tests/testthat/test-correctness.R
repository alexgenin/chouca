# 
# Various tests of model correctness
# 


test_that("pq coefficients work", { 
  
  # Debug model for pq
  mod <- camodel(transition(from = "a", to = "b", ~ 4 * (q["b"] * p["b"])), 
                wrap = TRUE, 
                neighbors = 4)
  initmm <- generate_initmat(mod, c(.25, .75), 2, 2)
  initmm[] <- c("b", "b", "a", "a")
  
  ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = "cpp")
  a <- run_camodel(mod, initmm, 1, control = ctrl)
  mat_final <- a[["output"]][["snapshots"]][[2]]
  expect_true(mat_final[1, 2] == "b")
  
  ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = "compiled")
  a <- run_camodel(mod, initmm, 1, control = ctrl)
  mat_final <- a[["output"]][["snapshots"]][[2]]
  expect_true(mat_final[1, 2] == "b")
  
  ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = "r")
  a <- run_camodel(mod, initmm, 1, control = ctrl)
  mat_final <- a[["output"]][["snapshots"]][[2]]
  expect_true(mat_final[1, 2] == "b")
  
})


test_that("fixed neighborhood works", { 
  
  mod <- camodel(transition(from = "a", to = "b", ~ q["b"]), 
                 wrap = FALSE, 
                 neighbors = 4)
  initmm <- generate_initmat(mod, c(.25, .75), 2, 2)
  initmm[] <- c("a", "b", "b", "a")
  ctrl <- list(console_output_every = 0, 
               save_snapshots_every = 1, 
               engine = "compiled")

  results <- plyr::ldply(seq.int(1999), function(i) { 
    a <- run_camodel(mod, initmm, 1, control = ctrl)
    mat_final <- a[["output"]][["snapshots"]][[2]]
    is_b_varnb <- mat_final[1,1] == "b" # should always be b
    
    ctrl[["fixed_neighborhood"]] <- TRUE
    a <- run_camodel(mod, initmm, 1, control = ctrl)
    mat_final <- a[["output"]][["snapshots"]][[2]]
    is_b_fixednb <- mat_final[1,1] == "b"
    
    data.frame(is_b_varnb = is_b_varnb, 
               is_b_fixednb = is_b_fixednb)
  })
  
  expect_true({ 
    all(results[ ,"is_b_varnb"]) 
  })
  
  # Should be close to ~0.5 because we divide by four instead of two
  expect_true({ 
    abs(mean(results[ ,"is_b_fixednb"]) - 0.5 ) < 0.1
  })
  
})
