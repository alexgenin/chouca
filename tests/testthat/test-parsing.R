# Test the parsing of models 




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
    
    all( tr[[2]][["X0"]] == 0.2 ), 
    all( tr[[2]][["XP"]] == c(0, 0) ), 
    all( tr[[2]][["XQ"]] == c(0, 0) ), 
    all( tr[[2]][["XPSQ"]] == c(0, 0) ), 
    all( tr[[2]][["XQSQ"]] == c(0, 0) )
    
  ))
})



modsq <- camodel(transition("a", "b", 
                            ~ 10 + q["a"] + 0.1 * p["b"] + 
                                4 * p["b"]^2 + 0.2 * q["a"]^2))
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
                parms = list(r = 1, a = 1))
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
                  parms = list(r = 1, a = 1))
})

# This model is supported by chouca because the ^3 coefficients cancel out
mod <- camodel(transition(from = "A", to = "E", 
                          ~ r * ( 1 - p["A"]^2 - 0.2*q["A"] * (q["A"]^2 + q["A"]) 
                                  + p["A"] + 0.2*q["A"]^3 )), 
                transition(from = "E", to = "A", 
                          ~ r * (a * q["A"])), 
                all_states = c("A", "E"), 
                parms = list(r = 1, a = 1))
