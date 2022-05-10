# 
# Test results against caspr 
# 

# Test that results match those of caspr for models that are in common
if ( dir.exists("./caspr_output") ) { 
  
  kubo_caspr <- readRDS("./caspr_output/caspr_ts_forestgap.rds")
  rownames(kubo_caspr) <- NULL
  
  # Make chouca run 
  p <- list(alpha = 0.2, d = 0.01, delta = 0.17)   # set parameters
  kubo_mod <- forestgap(p)
  initmat <- generate_initmat(kubo_mod, c(0.4, 0.6), 
                              nr = 100, nc = 100)
  ctrl <- list(substeps = 10, ca_engine = "compiled", console_output = FALSE)
  run <- run_camodel(kubo_mod, initmat, 200, control = ctrl)
  kubo_chouca <- as.data.frame(run[["output"]][["covers"]])
  names(kubo_chouca) <- c("time", "0", "+")
  
  # Display graphs
  # with(kubo_caspr, 
  #      plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
  # with(kubo_caspr,  lines(time, `+`, col = "red"))
  # with(kubo_chouca, lines(time, `+`, col = "red", lty = 2))
  # with(kubo_caspr,  lines(time, `0`, col = "blue"))
  # with(kubo_chouca, lines(time, `0`, col = "blue", lty = 2))
  
  # Make sure things are close 
  tolerance <- 0.1
  expect_true({ 
    all( abs(kubo_chouca[ ,"+"] - kubo_caspr[ ,"+"]) < tolerance )
  })
  
}

