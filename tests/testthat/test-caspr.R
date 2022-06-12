# 
# Test results against caspr 
# 

# Test that results match those of caspr for models that are in common
if ( dir.exists("./caspr_output") ) { 
  
  kubo_caspr <- readRDS("./caspr_output/caspr_ts_forestgap.rds")
  rownames(kubo_caspr) <- NULL
  
  # Make chouca run 
  p <- list(alpha = 0.2, d = 0.01, delta = 0.17)   # set parameters
  kubo_mod <- forestgap(parms = p)
  initmat <- generate_initmat(kubo_mod, c(0.4, 0.6), 
                              nr = 100, nc = 100)
  ctrl <- list(substeps = 1, engine = "cpp", console_output_every = 1)
  run <- run_camodel(kubo_mod, initmat, 100, control = ctrl)
  kubo_chouca <- as.data.frame(run[["output"]][["covers"]])
  names(kubo_chouca) <- c("time", "0", "+")
  
  # Display graphs
  with(kubo_caspr, 
       plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
  with(kubo_caspr,  lines(time, `+`, col = "red"))
  with(kubo_chouca, lines(time, `+`, col = "red", lty = 2))
  with(kubo_caspr,  lines(time, `0`, col = "blue"))
  with(kubo_chouca, lines(time, `0`, col = "blue", lty = 2))
  
  # Make sure things are close 
  tolerance <- 0.1
  expect_true({ 
    all( abs(kubo_chouca[ ,"+"] - kubo_caspr[ ,"+"]) < tolerance )
  })
  
  # NOTE: guichard's musselbed model cannot be implemented in caspr because it has 
  # rules like "if one neighbor in a given state then", 
  # 
  # guichard_caspr <- readRDS("./caspr_output/caspr_ts_musselbed.rds")
  # rownames(guichard_caspr) <- NULL
  # 
  # # Make chouca run 
  # p <- list(alpha = 0.4, d = 0.01, delta = 0.9)   # set parameters
  # guichard_mod <- musselbed(p)
  # initmat <- generate_initmat(guichard_mod, c(0.6, 0.2, 0.2), 
  #                             nr = 100, nc = 100)
  # ctrl <- list(substeps = 1, engine = "compiled", console_output_every = FALSE, 
  #              neighbors = 8)
  # run <- run_camodel(guichard_mod, initmat, 200, control = ctrl)
  # guichard_chouca <- as.data.frame(run[["output"]][["covers"]])
  # names(guichard_chouca) <- c("time", "m", "e", "d")
  # 
  # # Display graphs
  # with(guichard_caspr, 
  #      plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
  # with(guichard_caspr,  lines(time, `+`, col = "red"))
  # with(guichard_chouca, lines(time, `m`, col = "red", lty = 2))
  # with(guichard_caspr,  lines(time, `0`, col = "blue"))
  # with(guichard_chouca, lines(time, `e`, col = "blue", lty = 2))
  # with(guichard_caspr,  lines(time, `-`, col = "darkgreen"))
  # with(guichard_chouca, lines(time, `d`, col = "darkgreen", lty = 2))
  # 
}

