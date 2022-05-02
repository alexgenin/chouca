# 
# This file tests that the grabbing of neighbors is OK 
# 



mtest <- matrix(c(0, 1, 2, 3, 
                  0, 0, 1, 1, 
                  0, 1, 1, 1, 
                  2, 1, 1, 1), 
                byrow = TRUE, ncol = 4, nrow = 4)
nstates <- 4

testnb <- function(i, j, wrap, use_8_nb) { 
  local_dens(mtest, nstates, i, j, wrap, use_8_nb)
}
testnbcol <- function(j, wrap, use_8_nb) { 
  local_dens_col(mtest, nstates, j, wrap, use_8_nb)
}



expect_true(all(
  testnb(1, 1, wrap = FALSE, use_8_nb = FALSE) == c(1, 1, 0, 0)
))

expect_true(all(
  testnb(1, 1, wrap = TRUE, use_8_nb = FALSE) == c(1, 1, 1, 1)
))
  
expect_true(all(
  testnb(1, 1, wrap = TRUE, use_8_nb = TRUE) == c(2, 4, 1, 1)
))
  
expect_true(all(
  testnb(1, 1, wrap = FALSE, use_8_nb = TRUE) == c(2, 1, 0, 0)
))
  
expect_true(all(
  testnb(1, 4, wrap = FALSE, use_8_nb = FALSE) == c(0, 1, 1, 0)
))
  
expect_true(all(
  testnb(1, 4, wrap = TRUE, use_8_nb = FALSE) == c(1, 2, 1, 0)
))
  
expect_true(all(
  testnb(1, 4, wrap = TRUE, use_8_nb = TRUE) == c(2, 4, 2, 0)
))
  
expect_true(all(
  testnb(1, 4, wrap = FALSE, use_8_nb = TRUE) == c(0, 2, 1, 0)
))



# Test the computation of neighbors by columns 
for ( wrap in c(TRUE, FALSE) ) { 
  for ( use_8_nb in c(TRUE, FALSE) ) { 
      
  cols <- testnbcol(1, wrap, use_8_nb)
  for ( row in seq(1, nrow(mtest)) ) { 
    expect_true( all(cols[row, ] == testnb(row, 1, wrap, use_8_nb)) )
  }
  }
}
