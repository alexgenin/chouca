#
# Test the use of normalization functions
#

# TODO: do not warn if using norm function



test_that("softmax runs correctly", {

  mod <- camodel(
    # 0 means as probable as self-interaction
    transition(from = "0", to = "a", ~ 0),
    transition(from = "0", to = "b", ~ 0),
    all_states = c("0", "a", "b"),
    wrap = TRUE,
    neighbors = 4,
    normfun = "softmax"
  )
  initmm <- generate_initmat(mod, c(1, 0, 0), nrow = 256)

  # cpp engine
  run <- run_camodel(mod, initmm, times = c(0, 1), control = list(engine = "cpp"))
  expect_true({
    test <- run[["output"]][["covers"]][2, ]
    all(abs(test - c(1, 1/3, 1/3, 1/3)) < 0.02)
  })

  # compiled engine
  run <- run_camodel(mod, initmm, times = c(0, 1), control = list(engine = "compiled"))
  expect_true({
    test <- run[["output"]][["covers"]][2, ]
    all(abs(test - c(1, 1/3, 1/3, 1/3)) < 0.02)
  })

  expect_true({
    test <- run[["output"]][["covers"]][2, ]
    all(abs(test - c(1, 1/3, 1/3, 1/3)) < 0.02)
  })

  # This takes a while so bail now on CRAN
  skip_on_cran()

  # Non-identical probas
  p1 <- rnorm(1)
  p2 <- rnorm(1)
  Z <- exp(0) + exp(p1) + exp(p2)
  ps <- c(exp(0), exp(p1), exp(p2)) / Z # expected proportions

  mod <- camodel(
    # 0 means as probable as self-interaction
    transition(from = "0", to = "a", ~ p1),
    transition(from = "0", to = "b", ~ p2),
    all_states = c("0", "a", "b"),
    wrap = TRUE,
    neighbors = 4,
    normfun = "softmax"
  )
  initmm <- generate_initmat(mod, c(1, 0, 0), nrow = 256)

  run <- run_camodel(mod, initmm, times = c(0, 1), control = list(engine = "compiled"))
  expect_true({
    test <- run[["output"]][["covers"]][2, ]
    all(abs(test - c(1, ps)) < 0.02)
  })

  run <- run_camodel(mod, initmm, times = c(0, 1), control = list(engine = "cpp"))
  expect_true({
    test <- run[["output"]][["covers"]][2, ]
    all(abs(test - c(1, ps)) < 0.02)
  })

})

test_that("there is no false transition", {

  mod <- camodel(
    transition(from = "0", to = "a", ~ 0.1 + 0.2 * q["a"] + 0.3 * q["b"]),
    transition(from = "0", to = "b", ~ 0.1 - 0.2 * q["a"] + 0.1 * q["b"]),
    wrap = TRUE,
    neighbors = 4,
    normfun = "softmax"
  )

  initmm <- generate_initmat(mod, c(0.2, 0.4, 0.4), nrow = 100)
  times <- c(0, 1)
  run <- run_camodel(mod, initmm, times, control = list(engine = "cpp"))

  a <- run[["output"]][["snapshots"]]
  expect_true({
    ! any( a[[1]] == "a" & a[[2]] == "b" )
  })


  initmm <- generate_initmat(mod, c(0.2, 0.4, 0.4), nrow = 100)
  times <- c(0, 1)
  run <- run_camodel(mod, initmm, times, control = list(engine = "compiled"))

  a <- run[["output"]][["snapshots"]]
  expect_true({
    ! any( a[[1]] == "a" & a[[2]] == "b" )
  })

})


test_that("parsing is fine with broader range of rate values", {

  expect_warning({
    mod <- camodel(
      # 0 means as probable as self-interaction
      transition(from = "0", to = "a", ~ -1),
      transition(from = "0", to = "b", ~ 0),
      all_states = c("0", "a", "b"),
      wrap = TRUE,
      neighbors = 4,
      normfun = "identity"
    )
  }, regexp = "below 0")

  # No warning
  mod <- camodel(
    # 0 means as probable as self-interaction
    transition(from = "0", to = "a", ~ -1),
    transition(from = "0", to = "b", ~ 0),
    all_states = c("0", "a", "b"),
    wrap = TRUE,
    neighbors = 4,
    normfun = "softmax"
  )
})

test_that("parsing is correct with normfun", {
  # No warning
  expect_error({
    mod <- camodel(
      # 0 means as probable as self-interaction
      transition(from = "0", to = "a", ~ -1),
      transition(from = "0", to = "b", ~ 0),
      all_states = c("0", "a", "b"),
      wrap = TRUE,
      neighbors = 4,
      normfun = "fsdqfeqreza"
    )
  }, regexp = "normalization function")
})
