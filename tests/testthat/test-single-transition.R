# 
# This file contains tests to make sure individual transitions occur correctly. We test 
# the first transition happening in all simulations. Conveniently, this also makes 
# sure the first random number we generate is unbiased. 
# 


# Debug
mod <- camodel(
  transition(from = "0", to = "a", ~ a + b * q["a"] + c * p["a"]), 
  parms = list(a = 0, b = 0, c = 0), 
  wrap = FALSE, 
  neighbors = 8
)
initmat <- as.camodel_initmat(
  matrix(c("0", "a", 
           "a", "a"), 
         nrow = 2, ncol = 2)
)

# mod <- ca_library("rock-paper-scissor", parms = list(prob = 0.5))
# initmat <- generate_initmat(mod, rep(1/3, 3), nrow = nr, ncol = nc)

control <- list(console_output_every = 0, 
                save_covers_every = 1, 
                save_snapshots_every = 0, 
                precompute_probas = FALSE)

repmod <- function(n, this_mod, init) { 
  # Check that we reproduce well the variance and mean of time series between the two 
  # engines. 
  engines_ts <- replicate(n, { 
    niter <- seq(0, 1)
    modcompiled <- run_camodel(this_mod, init, niter, control = { 
      control[["engine"]] <- "compiled" ; control
    })
    modcpp <- run_camodel(this_mod, init, niter, control = { 
      control[["engine"]] <- "cpp" ; control
    })
    
    switch_compiled <- diff(modcompiled[["output"]][["covers"]][ ,2]) != 0
    switch_cpp <- diff(modcpp[["output"]][["covers"]][ ,2]) != 0
    cbind(switch_compiled, switch_cpp)
  })
  
  apply(engines_ts, c(1, 2), mean)
}


# Transition with no dependence on neighborhood
ptrans <- runif(1, 0.7, 0.8)
mod <- update(mod, list(a = ptrans, b = 0, c = 0))
rhos <- repmod(999, mod, initmat)
expect_true( abs(rhos[ ,"switch_compiled"] - ptrans) < 0.05 )
expect_true( abs(rhos[ ,"switch_cpp"] - ptrans) < 0.05 )

# Transition with dependence on local neighborhood
ptrans <- runif(1, 0.7, 0.8)
mod <- update(mod, list(a = 0, b = ptrans, c = 0))
q_a <- mean(initmat[-1] == "a")
rhos <- repmod(999, mod, initmat)
expect_true( abs(rhos[ ,"switch_compiled"] - (ptrans*q_a)) < 0.05 )
expect_true( abs(rhos[ ,"switch_cpp"] - (ptrans*q_a)) < 0.05 )

# Transition with dependence on global neighborhood
ptrans <- .75 # runif(1, 0.7, 0.8)
mod <- update(mod, list(a = 0, b = 0, c = ptrans))
rhos <- repmod(999, mod, initmat)
p_a <- mean(initmat == "a")
expect_true( abs(rhos[ ,"switch_compiled"] - (ptrans*p_a)) < 0.05 )
expect_true( abs(rhos[ ,"switch_cpp"] - (ptrans*p_a)) < 0.05 )


