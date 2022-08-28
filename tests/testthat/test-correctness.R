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
  
  ctrl <- list(console_output_every = 0, save_snapshots_every = 1, substeps = 2)
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
