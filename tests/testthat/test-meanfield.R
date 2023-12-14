# 
# 
# 

if ( requireNamespace("deSolve", quietly = TRUE) ) { 
  mod <- camodel(transition(from = "0", to = "+", ~ r * p["+"]), 
                  wrap = TRUE, 
                  parms = list(r = 0.1), 
                  neighbors = 8)

  init <- generate_initmat(mod, c(0.9, 0.1), nr = 1024)

  ctrl <- list(console_output_every = 0)

  test_that("Mean field should be a very close approximation in logistic growth", { 
    out <- run_camodel(mod, init, times = seq(0, 100), control = ctrl)
    ts <- out[["output"]][["covers"]]
    plot(out)
    
    o <- run_meanfield(mod, init, times = seq(0, 100))
    plot(o[ ,"time"], o[ ,"+"], type = "l")
    lines(ts[ ,"t"], ts[ ,"+"], col = "red")
    
    expect_true({ 
      all(abs(o[ ,"+"] - ts[ ,"+"]) < 0.05)
    })
    
  })

  test_that("Mean field can be used with proportion of states", { 
    # Check that 
    o1 <- run_meanfield(mod, c(0.9, 0.1), times = seq(0, 100))
    o2 <- run_meanfield(mod, init, times = seq(0, 100))
    
    expect_true({ 
      all(abs(o1 - o2) < 0.05)
    })
    
  })
}
