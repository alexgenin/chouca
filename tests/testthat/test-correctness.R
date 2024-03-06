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
  mat_final <- a[["output"]][["snapshots"]][[1]]
  expect_true(mat_final[1, 2] == "b")

  ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = "compiled")
  a <- run_camodel(mod, initmm, 1, control = ctrl)
  mat_final <- a[["output"]][["snapshots"]][[1]]
  expect_true(mat_final[1, 2] == "b")

})

test_that("qq coefficients work", {

  # Debug model for pq
  mod <- camodel(transition(from = "a", to = "b", ~ 4 * (q["b"] * q["a"])),
                 wrap = TRUE,
                 neighbors = 4)
  initmm <- generate_initmat(mod, c(.25, .75), 2, 2)
  initmm[] <- c("b", "b", "a", "a")

  for ( engine in c("cpp", "compiled") ) {
    ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = engine)
    a <- run_camodel(mod, initmm, 1, control = ctrl)
    mat_final <- a[["output"]][["snapshots"]][[1]]
    expect_true(mat_final[1, 2] == "b")
  }
})

test_that("Multiple states in the f(q) functions work", {
  # At some point chouca was not computing right the transition probabilities
  # of the form P = f(q1) + f(q2) and would only consider the first state q1,
  # so we added a test here
  mod <- camodel(
    transition(from = "a", to = "b", ~ 0.01 * q["b"] + 0.01 * q["c"] + 3*q["d"]),
    wrap = FALSE,
    all_states = c("a", "b", "c", "d"),
    neighbors = 8
  )

  initmm <- matrix(c("a", "b",
                     "c", "d"),
                   nrow = 2, ncol = 2, byrow = TRUE)
  initmm <- as.camodel_initmat(initmm)
  results <- replicate(19, {
    ctrl <- list(console_output_every = 0, save_snapshots_every = 1, engine = "cpp")
    a <- run_camodel(mod, initmm, 1, control = ctrl)
    o1 <- a[["output"]][["snapshots"]][[1]][1, 1]
    ctrl[["engine"]] <- "compiled"
    b <- run_camodel(mod, initmm, 1, control = ctrl)
    o2 <- b[["output"]][["snapshots"]][[1]][1, 1]
    c(o1, o2)
  })

  # Results should always agree because if we take the second state into
  # account, then proba is always > 1
  expect_true({
    all(results[1, ] == results[2, ])
  })

})

test_that("fixed neighborhood works", {

  mod <- camodel(transition(from = "a", to = "b", ~ q["b"]),
                 wrap = FALSE,
                 neighbors = 4)
  modfnb <- camodel(transition(from = "a", to = "b", ~ q["b"]),
                 wrap = FALSE,
                 neighbors = 4,
                 fixed_neighborhood = TRUE)

  initmm <- generate_initmat(mod, c(.25, .75), 2, 2)
  initmm[] <- c("a", "b", "b", "a")
  ctrl <- list(console_output_every = 0,
               save_snapshots_every = 1,
               engine = "compiled")

  results <- plyr::ldply(seq.int(1999), function(i) {
    a <- run_camodel(mod, initmm, 1, control = ctrl)
    mat_final <- a[["output"]][["snapshots"]][[1]]
    is_b_varnb <- mat_final[1,1] == "b" # should always be b

    a <- run_camodel(modfnb, initmm, 1, control = ctrl)
    mat_final <- a[["output"]][["snapshots"]][[1]]
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


test_that("We can run a model with only one state", {

  mod <- camodel(
    transition(from = "a", to = "a", ~ 1),
    wrap = TRUE,
    neighbors = 4
  )
  init <- generate_initmat(mod, 1, nr = 100)
  expect_true({
    # Just try to run the model like this
    run_camodel(mod, init, seq(0, 2))
    TRUE
  })

})

test_that("The coral reef model is consistent across engines", { 
  # The coral reef model had issues at some point with different engines, so make 
  # sure that is fixed
  mod <- ca_library("coralreef")
  mod <- update(mod, parms = list(m_c = 10^(-3.5), d_0 = -0.5, h_u = 10), 
                wrap = FALSE, neighbors = 8)
  
  ctrl <- list(engine = "cpp", substeps = 2)
  init <- generate_initmat(mod, rep(1/3, 3), nrow = 64, ncol = 64)
  out_cpp <- run_camodel(mod, init, seq(1, 1024), control = ctrl)
  
  ctrl <- list(engine = "compiled", force_compilation = TRUE, substeps = 4)
  out_compiled <- run_camodel(mod, init, seq(1, 1024), control = ctrl)
  
  # Compare output 
  out_cpp_t <- out_cpp[["output"]][["covers"]]
  out_compiled_t <- out_compiled[["output"]][["covers"]]
  
  plot(out_cpp_t[ ,"t"], out_cpp_t[ ,"CORAL"], type = "l")
  lines(out_compiled_t[ ,"t"], out_compiled_t[ ,"CORAL"], type = "l", col = "red")
  
  expect_true({ 
    all(abs(out_cpp_t[ ,"CORAL"] - out_compiled_t[ ,"CORAL"]) < 0.05)
  })
  
})
