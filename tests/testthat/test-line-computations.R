# 
# This file contains test to make sure we get the line and line changes right in 
# 




# Test that the function getline() works and produces correct results 
lapply(c(2, 4), function(ns) { 
  lapply(c(4, 8), function(nb) { 
    all_qs_example <- as.matrix( do.call(expand.grid, rep(list(seq(0, nb)), ns)) )
    all_qs_example <- all_qs_example[apply(all_qs_example, 1, sum) %% nb == 0, ]
    
    for ( i in seq.int(nrow(all_qs_example)) ) { 
      # we need a +1 here because getline() returns the index for c++ arrays
      expect_true({ 
        ( getline(all_qs_example[i, ], nb, ns) + 1 ) == i 
      })
    }
    return(TRUE)
  })
})



# 
nb <- 4
ns <- 3

lapply(c(4, 8), function(nb) { 
  all_qs_example <- as.matrix( do.call(expand.grid, rep(list(seq(0, nb)), ns)) )
  all_qs_example <- all_qs_example[apply(all_qs_example, 1, sum) == nb, ]
  all_qs_example <- all_qs_example[ ,seq(ncol(all_qs_example), 1)]

  expect_true({ 
    all( generate_all_qs(nb, ns, filter = TRUE) == all_qs_example )
  })
})


if ( FALSE ) { 
  
  a <- plyr::ldply(seq.int(7), function(ns) { 
    plyr::ldply(c(8, 4), function(nb) { 
      times <- system.time({ 
        n <- nrow(generate_all_qs(nb, ns, filter = TRUE))
      })
      times2 <- system.time({ 
        all_qs <- rep( list(seq(0, nb) ), each = ns)
        all_qs <- as.matrix(do.call(expand.grid, all_qs)) 
        # revert because the math expects the states in that order
        all_qs <- all_qs[ ,seq(ncol(all_qs), 1), drop = FALSE]
        colnames(all_qs) <- rownames(all_qs) <- NULL
        all_qs <- cbind(all_qs, apply(all_qs, 1, sum))
      })
      
      nbase <- (nb+1)^ns / nb
      data.frame(nb = nb, ns = ns, nrows = n, 
                 as.data.frame(as.list(times)), 
                 as.data.frame(as.list(times2)), 
                 nbase = nbase)
    })
  })
  
  library(ggplot2)
  ggplot(a, aes(x = ns, y = nrows, 
                color = as.factor(nb))) + 
    geom_point() + 
    geom_line() + 
    geom_line(aes(y = nbase)) + 
    scale_y_continuous(trans = "log10") 
  
}


