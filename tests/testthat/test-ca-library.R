# 
# 
# TODO: 
#   coralreef 
#   
# 

test_that("Example models produce correct results", { 

  # 
  # GAME OF LIFE
  # ------------
  gol <- ca_library("gameoflife")
  imat <- matrix(0, nrow = 6, ncol = 6)
  glider <- matrix(c(0, 0, 1, 
                     1, 0, 1, 
                     0, 1, 1), byrow = TRUE, ncol = 3)
  imat[1:3, 1:3] <- glider
  fmat <- factor(ifelse(imat > 0, "LIVE", "DEAD"), levels = c("DEAD", "LIVE"))
  dim(fmat) <- dim(imat)
  
  control <- list(console_output_every = 0, 
                  save_covers_every = 1, 
                  save_snapshots_every = 1, 
                  engine = "cpp")

  a <- run_camodel(gol, fmat, 10, control = control)

  sn <- a[["output"]][["snapshots"]]
#   spatialwarnings::display_matrix(sn)

  m <- tail(sn, 1)[[1]]
  expect_true({ 
    all( c(m[3, 5] == "LIVE", 
           m[4, 6] == "LIVE", 
           m[5, 4] == "LIVE", 
           m[5, 5] == "LIVE", 
           m[5, 6] == "LIVE", 
           m[4, 5] == "DEAD") )
  })

})


test_that("Errors in ca_library are handled", { 
  expect_error({ 
    ca_library("gsqgcxwbwfsdcxwvcwxvefaezaq")
  })
})
