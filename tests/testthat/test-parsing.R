# Test the parsing of models 


# Test that a warning is emitted when probabilities can be above one
expect_warning({ 
  forestgap(parms = list(d = 10, 
                         delta = 0.5, 
                         alpha = 0.2))
})

# Test that a warning is emitted when probabilities can be below zero
expect_warning({ 
  forestgap(parms = list(d = -1, 
                         delta = 0.5, 
                         alpha = 0.2))
})


tr <- forestgap(parms = list(d = 0.125, 
                             delta = 0.5, 
                             alpha = 0.2))[["transitions"]]

expect_true({ 
  all(c(
    all( tr[[1]][["X0"]] == 0.125 ), 
    all( tr[[1]][["XP"]] == c(0, 0) ), 
    all( tr[[1]][["XQ"]] == c(0.5, 0) ), 
    all( tr[[1]][["XPSQ"]] == c(0, 0) ), 
    all( tr[[1]][["XQSQ"]] == c(0, 0) ), 
    
    all( tr[[2]][["X0"]] == 0 ), 
    all( tr[[2]][["XP"]] == c(0, 0.2) ), 
    all( tr[[2]][["XQ"]] == c(0, 0) ), 
    all( tr[[2]][["XPSQ"]] == c(0, 0) ), 
    all( tr[[2]][["XQSQ"]] == c(0, 0) )
    
  ))
})



modsq <- camodel(transition("a", "b", 
                            ~ 10 + q["a"] + 0.1 * p["b"] + 
                                4 * p["b"]^2 + 0.2 * q["a"]^2), 
                 check_model = FALSE)
tr <- modsq[["transitions"]]

expect_true({ 
  all(c(
    all( abs(tr[[1]][["X0"]] - 10 )          < 1e-8), 
    all( abs(tr[[1]][["XP"]] - c(0, 0.1) )   < 1e-8), 
    all( abs(tr[[1]][["XQ"]] - c(1, 0) )     < 1e-8), 
    all( abs(tr[[1]][["XPSQ"]] - c(0, 4) )   < 1e-8), 
    all( abs(tr[[1]][["XQSQ"]] - c(0.2, 0) ) < 1e-8)
  ))
})



# This model has trouble with XPSQ
mod <- camodel(transition(from = "A", to = "E", 
                          ~ r * ( 1 - p["A"]^2 + p["A"] + q["A"] + q["A"]^2 )), 
               transition(from = "E", to = "A", 
                          ~ r * (a * q["A"])), 
               all_states = c("A", "E"), 
               parms = list(r = 1, a = 1), 
               check_model = FALSE)
tr <- mod[["transitions"]]

expect_true({ 
  all(c(
    all( abs(tr[[1]][["X0"]] - 1 )          < 1e-8), 
    all( abs(tr[[1]][["XP"]] - c(1, 0) )   < 1e-8), 
    all( abs(tr[[1]][["XQ"]] - c(1, 0) )     < 1e-8), 
    all( abs(tr[[1]][["XPSQ"]] - c(-1, 0) )   < 1e-8), 
    all( abs(tr[[1]][["XQSQ"]] - c(1, 0) ) < 1e-8)
  ))
})

# This model is unsupported by chouca (^3 coefficient)
expect_error({ 
  mod <- camodel(transition(from = "A", to = "E", 
                            ~ r * ( 1 - p["A"]^2 + p["A"] + q["A"] + 0.2*q["A"]^3 )), 
                  transition(from = "E", to = "A", 
                            ~ r * (a * q["A"])), 
                 all_states = c("A", "E"), 
                 parms = list(r = 1, a = 1), 
                 check_model = TRUE)
})

# This model is supported by chouca because the ^3 coefficients cancel out
mod <- camodel(transition(from = "A", to = "E", 
                          ~ r * ( 1 - p["A"]^2 - 0.2*q["A"] * (q["A"]^2 + q["A"]) 
                                  + p["A"] + 0.2*q["A"]^3 )), 
                transition(from = "E", to = "A", 
                          ~ r * (a * q["A"])), 
                all_states = c("A", "E"), 
                parms = list(r = 1, a = 1), 
               check_model = FALSE)


# Test updating of models. Internal coefficient matrix should have changed 
mod <- forestgap()
mod2 <- update(mod, parms = list(d = 0, delta = 0.5, alpha = 0.01))
expect_true({ 
  all(mod[["parms"]] != mod[["mod2"]])
})
expect_true({ 
  ! all(mod[["transmatrix"]] == mod2[["transmatrix"]])
})



# Test generate_initmat()
mod <- forestgap()

# Wrong number of states -> error 
expect_error({ 
  generate_initmat(mod, c(1, 0, 1), 2, 2)
})
# Named vector, but not the right names -> error
expect_error({ 
  generate_initmat(mod, c(A = 1, B = 1), 2, 2)
})
# NAs -> error
expect_error({ 
  generate_initmat(mod, c(A = NA, B = 1), 2, 2)
})
# Sum of states above/below 1 -> warning
expect_warning({ 
  generate_initmat(mod, c(1, 2), 2, 2)
})
expect_warning({ 
  generate_initmat(mod, c(0, 0.1), 2, 2)
})


