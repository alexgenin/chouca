# 
# Test the parsing of models from matrices
# 

genmod <- function(ns) { 
  stopifnot(ns == 4) 
  
  beta0 <- matrix(0, nrow = ns, ncol = ns)
  betap <- array(0, dim = c(ns, ns, ns))
  betaq <- array(0, dim = c(ns, ns, ns))

  beta0[1,  ] <- 0.2
  for ( i in seq.int(ns) ) { 
    betap[1, i, i] <- 0.1
    betaq[1, i, i] <- 0.01
  }
  betaq[1, 4, 2] <- 0.04
  
  beta0[ , 1] <- .4

  # Set zeros 
  for ( i in seq.int(ns) ) { 
    beta0[i, i] <- 0
    betap[i, i, ] <- 0
    betaq[i, i, ] <- 0
  }
  
  camodel_mat(beta_0 = beta0, 
              beta_p = betap, 
              beta_q = betaq, 
              all_states = c("0", paste0("sp", seq.int(ns-1))), 
              neighbors = 4, 
              wrap = TRUE, 
              build_transitions = TRUE)
}

NS <- 4 # 0 + 3 species
mod_adv <- genmod(NS)
mod <- camodel(
  transition(from = "0", to = "sp1", ~ 0.2 + 0.1 * p["sp1"] + 0.01 * q["sp1"]), 
  transition(from = "0", to = "sp2", ~ 0.2 + 0.1 * p["sp2"] + 0.01 * q["sp2"]), 
  transition(from = "0", to = "sp3", ~ 0.2 + 0.1 * p["sp3"] + 0.04 * q["sp1"] + 
                                         0.01 * q["sp3"]), 
  transition(from = "sp1", to = "0", ~ 0.4), 
  transition(from = "sp2", to = "0", ~ 0.4), 
  transition(from = "sp3", to = "0", ~ 0.4), 
  wrap = TRUE, 
  neighbors = 4, 
  all_states = c("0", "sp1", "sp2", "sp3")
)

test_that("Advanced model-building builds correct things", { 
  expect_true({ 
    all(names(mod) == names(mod_adv))
  })
  
  comp_test <- names(mod)[ ! names(mod) %in% c("max_error", "max_rel_error", 
                                               "transitions") ]
  for ( i in comp_test ) { 
    expect_true({ 
      all.equal(mod[[i]], mod_adv[[i]])
    })
  }
  
# We do not test transitions because their formula expression often do not match, but if
# internal tables match, then this means transitions match by definition
#   capture.output({ 
#     expect_true({ 
#       all.equal(dput(mod[["transitions"]]), 
#                 dput(mod_adv[["transitions"]]))
#     })
#   })
})

test_that("Advanced model-building runs identically", { 
  ctrl <- list(engine = "compiled", force_compilation = FALSE, 
               console_output_every = 0)
  init <- generate_initmat(mod, rep(1/NS, NS), nrow = 256, ncol = 256)
  
  run_adv <- run_camodel(mod_adv, init, seq.int(128), control = ctrl)
  run_old <- run_camodel(mod, init, seq.int(128), control = ctrl)
  
  adv <- run_adv[["output"]][["covers"]]
  old <- run_old[["output"]][["covers"]]
  
  expect_true({ 
    max(abs(adv - old)) < 0.05
  })
})

