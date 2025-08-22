
# Only one chouca dll loaded
options(chouca.maxdll = 1)

# Make 20 models with compilation
mods <- lapply(1+sample.int(50, size = 2), function(i) {
  # Voter model
  gamma <- 0.5
  kappa <- i
  coefs <- array(0, dim = rep(kappa, 3))
  for ( state in seq.int(kappa) ) {
    coefs[ , state, state] <- gamma
    coefs[state, state, ] <- 0
  }

  mod <- camodel_mat(
    beta_q = coefs,
    wrap = TRUE,
    all_states = as.character(seq.int(kappa)),
    neighbors = 8,
    normfun = "identity",
    build_transitions = TRUE
  )
})

# This will fail if maxdll = 3 and we do not handle proper removal of generated
# Rcpp functions, as the third time we run the lapply() call then we will reuse
# a model that has its DLL unloaded.
test_that("We cleanup things properly when DLL unloading", {
  options(chouca.maxdll.warn = FALSE)
  lapply(c(mods, mods), function(mod) {

    initmm <- generate_initmat(mod, rep(1/mod[["nstates"]], mod[["nstates"]]),
                               nrow = 32)

    # Time sequence
    iters <- seq(0, 2)

    # Simulation options
    control <- list(console_output_every = floor(max(iters)/5),
                    engine = "compiled")

    # Simulate the model
    mod_out <- run_camodel(mod, initmm, iters, control = control)

    return(TRUE)
  })

  expect_true(TRUE)
})

# Reset options
options(chouca.maxdll.warn = NULL)
options(chouca.maxdll = NULL)
