
#'@export
show_landscape <- function(mod, col = NULL, fps_cap = 24, transpose = TRUE, ...) { 
  
  # Check that palette makes sense 
  if ( length(col) != mod[["nstates"]] ) { 
    stop(paste0("The number of colors in col does not correspond to the ", 
                "number of states in the model."))
  }
  
  function(t, m) { 
    
    this_call_time <- Sys.time()
    dtime <- as.numeric(difftime(this_call_time, last_call_time, units = "secs"))
    
    # If {1/fps_cap} has not passed since last time, we wait a little bit
    if ( dtime < (1/fps_cap) ) { 
      Sys.sleep( (1/fps_cap) - dtime )
    }
    
    # Hold plot update, and flush it when exiting the function and everything is drawn
    dev.hold() 
    on.exit({ dev.flush() })
    
    if ( transpose ) { 
      m <- t(m)
    }
    
    # m is a factor, so limits should be 1:nstates
    image.camodel_initmat(m, col = col, zlim = c(1, mod[["nstates"]]), ...)
    last_call_time <<- Sys.time()
    return(TRUE)
  }
  
}
