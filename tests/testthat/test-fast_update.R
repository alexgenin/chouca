#
# Update a fast model
#

# devtools::document()

mod <- camodel(
  transition(from = "0", to = "p",
             ~ ( q["p"] > 0.5 ) * a +
                 b * p["p"] +
                 c * p["p"] * p["p"] +
                 d * p["p"] * q["p"] +
                 c*d * q["p"] * q["0"]
            ),
  transition(from = "p", to = "0", ~ 0.1),
  wrap = TRUE,
  neighbors = nb_kernel("moore", 2),
  parms = list(a = 0.1, b = 0.2, c = 0.3, d = 0.4)
)

mod_new <- camodel(
  transition(from = "0", to = "p",
             ~ ( q["p"] > 0.5 ) * a +
                 b * p["p"] +
                 c * p["p"] * p["p"] +
                 d * p["p"] * q["p"] +
                 c*d * q["p"] * q["0"]
            ),
  transition(from = "p", to = "0", ~ 0.1),
  wrap = TRUE,
  neighbors = nb_kernel("moore", 2),
  parms = list(a = 0.09, b = 0.31, c = 0.53, d = 0.12)
)


test_that("Updating produces correct results", {
  mod2 <- update(mod, parms = mod$parms)
  expect_true( all.equal(mod, mod2) )

  mod2 <- update(mod, parms = mod_new$parms)
  expect_true({
    all.equal(mod_new, mod2) # this takes into account small differences in floats
  })

})

test_that("We warn on model update if the fast update fails", {

  mod <- camodel(
    transition(from = "0", to = "a", ~ r),
    transition(from = "0", to = "b", ~ r),
    transition(from = "a", to = "0", ~ r - m[1,1] * q["a"] - m[1,2] * q["b"]),
    transition(from = "b", to = "0", ~ r - m[2,1] * q["a"] - m[2,2] * q["b"]),
    parms = list(r = 0.1, m = diag(2) * 0),
    neighbors = 4,
    wrap = TRUE,
    normfun = "identity"
  )

  new_pars <- list(
    r = 0.2,
    m = {
      m <- diag(2) * 0
      m[1, 1] <- 0.1
      m[2, 1] <- - 0.1
      m
  })

  mod_new_ref <- camodel(
    transition(from = "0", to = "a", ~ r),
    transition(from = "0", to = "b", ~ r),
    transition(from = "a", to = "0", ~ r - m[1,1] * q["a"] - m[1,2] * q["b"]),
    transition(from = "b", to = "0", ~ r - m[2,1] * q["a"] - m[2,2] * q["b"]),
    parms = new_pars,
    neighbors = 4,
    wrap = TRUE,
    normfun = "identity"
  )

  # We show a warning that we can produce wrong results
  expect_warning({
    mod_new <- update(mod, parms = new_pars, check_model = "none")
  })

  # We actually produce the wrong results
  suppressWarnings({
    mod_new <- update(mod, parms = new_pars, check_model = "none")
  })
  expect_true({
    ! identical(mod_new_ref, mod_new)
  })

  # If we check the model, then this is caught and we produce good results
  mod_new <- update(mod, parms = new_pars, check_model = "quick")
  expect_true({
    all.equal(mod_new_ref, mod_new)
  })

})
