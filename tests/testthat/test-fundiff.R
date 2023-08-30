# 
# This file makes sure that chouca has the same behavior as differential equations when 
# it should be. 
# 

test_that("Behavior is OK compared to differential equation", { 
  
  dy <- function(x) { 
    0.05 * x * ( 1 - x )
  }

  times <- seq(0, 512)
  y <- times * 0
  y[1] <- 0.001
  for ( i in seq(2, length(times)) ) { 
    y[i] <- y[i-1] + dy(y[i-1]) * ( times[i] - times[i-1] )
  }

  # Logistic growth
  # plot(times, y)

  mod <- camodel(
    transition("0", "+", ~ 0.05 * p["+"] ), 
    wrap = TRUE, 
    neighbors = 4
  )
  initmm <- generate_initmat(mod, c(`0` = 1 - 0.001, `+` = 0.001), nr = 1024)
  run <- run_camodel(mod, initmm, times, 
                     control = list(engine = "compiled", 
                                    precompute_probas = FALSE, 
                                    force_compilation = TRUE, 
                                    console_output_every = 0))
  
  covs <- as.data.frame(run[["output"]][["covers"]])
  covs[ ,"diff+"] <- y
  
  # matplot(covs[ ,1], covs[ ,-1], type = "l")
  
  expect_true({ 
    all( abs(covs[ ,"+"] - covs[ ,"diff+"]) < 0.1 )
  })

  slope <- coef(lm(covs[ ,"+"] ~ covs[ ,"diff+"]))[2]
  expect_true({ 
    abs(slope - 1) < 0.1
  })
})
