

test_that("all_qs is generated correctly", { 
  
  max_nb <- 8
  ns <- 5 
  wrap <- TRUE 
  
  all_qs <- rep( list(seq(0, max_nb) ), each = ns)
  all_qs <- as.matrix(do.call(expand.grid, all_qs)) 
  
  # revert because the math expects the states in that order
  all_qs <- all_qs[ ,seq(ncol(all_qs), 1)]
  colnames(all_qs) <- rownames(all_qs) <- NULL

  # If the number of neighbors is constant, we can discard the data that is present 
  # every 4 or 8 neighbors. 
  if ( wrap ) { 
    all_qs <- all_qs[seq(0, nrow(all_qs)-1) %% max_nb == 0, ]
  }

  all_qs <- all_qs[-1, ] # this has sum zero which produces division by zero. Discard it.
  
  all_qs <- cbind(all_qs, rowSums(all_qs))
  
  all_qs2 <- generate_all_qs(max_nb, ns, filter = wrap, line_cap = 0)
  
  expect_true({ 
    all( all_qs == all_qs2 ) 
  })
})

