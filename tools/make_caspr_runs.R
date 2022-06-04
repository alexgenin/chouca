# 
# Make the caspr runs 

# Load caspr 
devtools::load_all("/home/alex/work/2022/chouca/tools/caspr-test/caspr")

# 
if ( FALSE ) { 
  
  # Run simulation of forestgap with caspr 
  l <- init_landscape(c("+","0"), c(0.6, 0.4), width = 100) 
  p <- list(alpha = 0.2, d = 0.01, delta = 0.17)   # set parameters
  r <- ca(l, model = forestgap, parms = p, t_max = 200)    # run simulation
  plot(r) 
  
  covers <- as.data.frame(with(r, cbind(time, cover)), 
                          row.names = NULL)
  saveRDS(covers, 
          file = file.path("./tests/testthat/caspr_output/", 
                           "caspr_ts_forestgap.rds"))
  
}

# Make caspr benchmarks 
if ( FALSE ) { 
  
  # Run simulation of musselbed with caspr 
  sizes <- c(16, 32, 64, 128, 256, 512)
  timings <- ldply(sizes, function(size) { 
    a <- system.time({ 
      l <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = size) 
      p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
      r <- ca(l, model = musselbed, parms = p, t_max = 512)    # run simulation
    })
    as.data.frame(as.list(a))
  }, .progress = "time")
  
}
