#
# Test the computation of individual transition rates
#

test_that("Transition rates are well-computed by the helper function", {

  for ( nf in c("identity", "softmax") ) {

    if ( nf == "identity" ) {
      cs <- runif(5) / 10
    } else {
      cs <- rnorm(5)
    }

    mod <- camodel(
      transition(from = "0", to = "a",
                 ~ cs[1]  +
                    cs[2] * q["a"] +
                    cs[3] * q["a"] * p["0"] +
                    cs[4] * p["0"] * p["0"] +
                    cs[5] * q["a"] * q["0"]),
      all_states = c("0", "a"),
      wrap = TRUE,
      neighbors = 4,
      normfun = nf
    )

    n <- 16
    initm <- generate_initmat(mod, c(.5, .5), nrow = n)

    ps <- next_state_probs(mod, initm, log = FALSE)

    # Correct dims
    expect_true({
      all( dim(ps) == c(n, n, length(mod[["states"]])) )
    })

    p <- purrr::map_dbl(mod[["states"]], function(s) mean(initm == s))
    names(p) <- mod[["states"]]

    # Get a cell in the 0 state and not on the sides of the matrix
    ok_cell <- FALSE
    while ( ! ok_cell ) {
      xy <- sample(2:15, size = 2) # not on a side
      if ( initm[xy[1], xy[2]] == "0" ) {
        ok_cell <- TRUE
      }
    }

    qa <- c(initm[xy[1] - 1, xy[2]],
            initm[xy[1] + 1, xy[2]],
            initm[xy[1], xy[2] - 1],
            initm[xy[1], xy[2] + 1])
    q <- purrr::map_dbl(mod[["states"]], function(s) mean(qa == s))
    names(q) <- mod[["states"]]

    p_to_a <- eval(as.expression(as.list(mod[["transitions"]][[1]][["prob"]])[[2]]),
                   envir = list(cs = cs, p = p, q = q))
    p_to_a <- as.numeric(p_to_a)

    if ( mod[["normfun"]] == "identity" ) {
      expected <- c(1 - p_to_a, p_to_a)
    } else if ( mod[["normfun"]] == "softmax" ) {
      expected <- exp( c(0, p_to_a) )
      expected <- expected / sum(expected)
    }

    returned <- ps[xy[1], xy[2], ]
#     print(expected)
#     print(
#       ps[xy[1], xy[2], ]
#     )
#     if ( ! all( abs(expected - ps[xy[1], xy[2], ] ) < 1e-2 ) ) { browser() }

    expect_true({
      all( abs(expected - returned ) < 1e-8 )
    })

    if ( ! all( abs(returned - expected) < 1e-8 ) ) {
      returned
      expected
      browser()
    }

    #   for the last one, do it as log too
    ps <- next_state_probs(mod, initm, log = TRUE)
    returned <- ps[xy[1], xy[2], ]

    # The test is not as stringent because there is numerical error here
    expect_true({
      all( abs(exp(returned) - expected) < 1e-4 )
    })
  }


})

# stop(")

test_that("Simulation and computation of transition rates match (identity norm)", {

  pars <- rnorm(3) / 10
  mod <- camodel(
    transition(from = "b", to = "c", ~ 0.4 + pars[1] * q["0"] + pars[2] * p["a"] +
                                          pars[3] * p["b"]),
    all_states = c("0", "a", "b", "c"),
    wrap = FALSE,
    neighbors = 8
  )

  initm <- as.camodel_initmat(matrix(c("0", "a", "a",
                                       "0", "b", "a",
                                       "0", "0", "0"),
                                     nrow = 3, ncol = 3, byrow = TRUE),
                              levels = mod[["states"]])
  p_b_to_c <- replicate(999, {
    a <- run_camodel(mod, initm, times = c(1))
    a[["output"]][["covers"]][ ,"c"] > 0
  })

  est_p_b_to_c <- mean(p_b_to_c)

  q <- c(5, 3, 0, 0) / 8
  p <- c(5, 3, 1, 0) / 8
  names(q) <- names(p) <- mod[["states"]]
  true_p_b_to_c <- eval(as.expression(as.list(mod[["transitions"]][[1]][["prob"]])[[2]]),
                        envir = list(pars = pars, q = q, p = p))

  # Returned value by get_proba
  ps <- next_state_probs(mod, initm, log = FALSE)
  calc_p_b_to_c <- ps[2, 2, 4]

  expect_true({
    abs(est_p_b_to_c - true_p_b_to_c) < 0.2
  })

  expect_true({
    abs(calc_p_b_to_c - true_p_b_to_c) < 0.01
  })

})


test_that("Simulation and computation of transition rates match (softmax norm)", {

  pars <- rnorm(3) / 10
  mod <- camodel(
    transition(from = "b", to = "c", ~ 0.4 + pars[1] * q["0"] + pars[2] * p["a"] +
                                          pars[3] * p["b"]),
    all_states = c("0", "a", "b", "c"),
    wrap = FALSE,
    normfun = "softmax",
    neighbors = 8
  )

  initm <- as.camodel_initmat(matrix(c("0", "a", "a",
                                       "0", "b", "a",
                                       "0", "0", "0"),
                                     nrow = 3, ncol = 3, byrow = TRUE),
                              levels = mod[["states"]])
  p_b_to_c <- replicate(999, {
    a <- run_camodel(mod, initm, times = 1, control = list(engine = "compiled"))
    a[["output"]][["covers"]][ ,"c"] > 0
  })

  est_p_b_to_c <- mean(p_b_to_c)

  q <- c(5, 3, 0, 0) / 8
  p <- c(5, 3, 1, 0) / 8
  names(q) <- names(p) <- mod[["states"]]
  true_p_b_to_c <- eval(as.expression(as.list(mod[["transitions"]][[1]][["prob"]])[[2]]),
                        envir = list(pars = pars, q = q, p = p))
  true_p_b_to_c <- exp(true_p_b_to_c) / ( 1.0 + exp(true_p_b_to_c) )

  # Returned value by get_proba
  ps <- next_state_probs(mod, initm, log = FALSE)
  calc_p_b_to_c <- ps[2, 2, 4]

  expect_true({
    abs(est_p_b_to_c - true_p_b_to_c) < 0.2
  })

  expect_true({
    abs(calc_p_b_to_c - true_p_b_to_c) < 0.01
  })

})
