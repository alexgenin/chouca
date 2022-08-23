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

# 
if ( FALSE ) { 
  
  # Run simulation of musselbed with caspr 
  p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
  l <- init_landscape(c("+","0", "-"),  c(0.6, 0.2, 0.2), width = 100) 
  
  r <- ca(l, model = musselbed, parms = p, t_max = 200, subs = 10)    # run simulation
  plot(r) 
  
  covers <- as.data.frame(with(r, cbind(time, cover)), 
                          row.names = NULL)
  saveRDS(covers, 
          file = file.path("./tests/testthat/caspr_output/", 
                           "caspr_ts_musselbed.rds"))
  
}

# 
if ( FALSE ) { 
  
  # Run simulation of grazing with caspr 
  p <- list(del = 0.9, b = 0.5, c_ = 0.2, m0 = 0.05, g = 0.2, r = 0.01, 
            f = 0.9, d = 0.1, p = 1)
  l <- init_landscape(c("+","0", "-"),  c(0.6, 0.2, 0.2), width = 100) 
  r <- ca(l, model = grazing, parms = p, t_max = 200, subs = 10)    # run simulation
  plot(r) 
  
  covers <- as.data.frame(with(r, cbind(time, cover)), 
                          row.names = NULL)
  
  saveRDS(covers, 
          file = file.path("./tests/testthat/caspr_output/", 
                           "caspr_ts_grazing.rds"))
}
