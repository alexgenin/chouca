# 
# Test results against caspr 
# 

# Test that results match those of caspr for models that are in common
if ( dir.exists("./caspr_output") ) { 
  
  kubo_caspr <- readRDS("./caspr_output/caspr_ts_forestgap.rds")
  rownames(kubo_caspr) <- NULL
  ts <- seq(0, 200)
  
  # Make chouca run 
  p <- list(alpha = 0.2, d = 0.01, delta = 0.17)   # set parameters
  kubo_mod <- ca_library("forestgap", parms = p)
  initmat <- generate_initmat(kubo_mod, c(0.4, 0.6), 
                              nr = 100, nc = 100)
  ctrl <- list(engine = "cpp", console_output_every = 0, substeps = 10)
  run <- run_camodel(kubo_mod, initmat, ts, control = ctrl)
  kubo_chouca <- as.data.frame(run[["output"]][["covers"]])
  names(kubo_chouca) <- c("time", "0", "+")
  
  # Display graphs
  if ( FALSE ) { 
    with(kubo_caspr, 
        plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
    with(kubo_caspr,  lines(time, `+`, col = "red"))
    with(kubo_chouca, lines(time, `+`, col = "red", lty = 2))
    with(kubo_caspr,  lines(time, `0`, col = "blue"))
    with(kubo_chouca, lines(time, `0`, col = "blue", lty = 2))
  }
  
  # Make sure things are close 
  tolerance <- 0.1
  expect_true({ 
    all( abs(kubo_chouca[ ,"+"] - kubo_caspr[ ,"+"]) < tolerance )
  })
  
  
  # Guichard's Mussel bed
  # 
  guichard_caspr <- readRDS("./caspr_output/caspr_ts_musselbed.rds")
  rownames(guichard_caspr) <- NULL
  
  # Make chouca run 
  p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
  guichard_mod <- update(ca_library("musselbed"), parms = p)
  
  initmat <- generate_initmat(guichard_mod, c(0.6, 0.2, 0.2), 
                              nr = 100, nc = 100)
  ctrl <- list(engine = "cpp", console_output_every = 0, substeps = 10)
  run <- run_camodel(guichard_mod, initmat, ts, control = ctrl)
  guichard_chouca <- as.data.frame(run[["output"]][["covers"]])
  names(guichard_chouca) <- c("time", "m", "e", "d")
  
  # Display graphs
  if ( FALSE ) { 
    with(guichard_caspr, 
         plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
    with(guichard_caspr,  lines(time, `+`, col = "red"))
    with(guichard_chouca, lines(time, `m`, col = "red", lty = 2))
    with(guichard_caspr,  lines(time, `0`, col = "blue"))
    with(guichard_chouca, lines(time, `e`, col = "blue", lty = 2))
    with(guichard_caspr,  lines(time, `-`, col = "darkgreen"))
    with(guichard_chouca, lines(time, `d`, col = "darkgreen", lty = 2))
  }
  
  expect_true({ 
    all( c(abs(guichard_caspr[ ,"+"] - guichard_chouca[ ,"m"]) < 0.1, 
           abs(guichard_caspr[ ,"0"] - guichard_chouca[ ,"e"]) < 0.1, 
           abs(guichard_caspr[ ,"-"] - guichard_chouca[ ,"d"]) < 0.1) )
  })
  
  
  
  
  
  # Kéfi arid vegetation model 
  # 
  aridvege_caspr <- readRDS("./caspr_output/caspr_ts_grazing.rds")
  p <- list(delta = 0.9, b = 0.5, c = 0.2, m0 = 0.05, g0 = 0.2, 
            r = 0.01, f = 0.9, d = 0.1, pr = 1)
  
  aridvege_mod <- ca_library("aridvege", parms = p)
  # Init configuration in caspr:
  # l <- init_landscape(c("+","0", "-"),  c(0.6, 0.2, 0.2), width = 100) 
  initmat <- generate_initmat(aridvege_mod, c(VEGE = 0.6, EMPTY = 0.2, DEGR = 0.2), 
                              nr = 100, nc = 100)
  ctrl <- list(engine = "cpp", console_output_every = 0, substeps = 10)
  run <- run_camodel(aridvege_mod, initmat, ts, control = ctrl)
  aridvege_chouca <- as.data.frame(run[["output"]][["covers"]])
  names(aridvege_chouca) <- c("time", "-", "0", "+")
  
  # Display graphs
  if ( FALSE ) { 
    with(aridvege_caspr, 
         plot(time, sample(c(0, 1), size = length(time), replace = TRUE), type = "n"))
    with(aridvege_caspr,  lines(time, `+`, col = "red"))
    with(aridvege_chouca, lines(time, `+`, col = "red", lty = 2))
    with(aridvege_caspr,  lines(time, `0`, col = "blue"))
    with(aridvege_chouca, lines(time, `0`, col = "blue", lty = 2))
    with(aridvege_caspr,  lines(time, `-`, col = "darkgreen"))
    with(aridvege_chouca,  lines(time, `-`, col = "darkgreen", lty = 2))
  }
  
  expect_true({ 
    all( c(abs(aridvege_caspr[ ,"+"] - aridvege_chouca[ ,"+"]) < 0.1, 
           abs(aridvege_caspr[ ,"0"] - aridvege_chouca[ ,"0"]) < 0.1, 
           abs(aridvege_caspr[ ,"-"] - aridvege_chouca[ ,"-"]) < 0.1) )
  })
  
}

