# 
# 
# These functions are functions factory for custom callbacks
# 
# 

#' @title Plot simulation landscapes 
#' 
#' @description This function creates a function to plot the model landscapes during 
#'   a \code{chouca} simulation, as it is being run. 
#' 
#' @param mod The model being run 
#' 
#' @param col a set of colors (character vector) of the same length as the number of 
#'   states in the model. 
#' 
#' @param fps_cap The maximum number of frame per seconds at which to plot the results 
#' 
#' @param burn_in A number of iterations to skip before plotting
#' 
#' @param transpose Setting this to \code{TRUE} will transpose the landscape matrix 
#'   before passing it to  \code{\link[graphics]{image}}
#' 
#' @param ... other arguments are passed to \code{\link[graphics]{image}}
#' 
#' @details
#' 
#'   This function creates another function that is suitable for use with \code{chouca}. 
#'   It allows plotting the landscape as it is being simulated, using the base function
#'   \code{\link[graphics]{image}}. You can set the colors using the argument \code{col}, 
#'   or tranpose the landscape before plotting using \code{transpose}. 
#'   
#'   \code{\link[graphics]{image}} is quite slow at displaying matrices, especially if 
#'   it is reasonably large, but if it is still to fast for your test, you can cap the 
#'   number of landscape displayed per seconds by setting the argument \code{fps_cap}. 
#'   
#'   It is important to note that this function will probably massively slow down a 
#'   simulation, so this function is mostly here for exploratory analyses, or just to 
#'   have a good look of what is happening in a model. 
#'   
#' @examples 
#' 
#' \dontrun{ 
#' 
#' # Display the psychedelic spirals of the rock-paper-scissor model as the model is 
#' # being run 
#' mod <- ca_library("rock-paper-scissor")
#' colors <- c("#fb8500", "#023047", "#8ecae6")
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = landscape_plotter(mod, col = colors))
#' init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)
#' run_camodel(mod, init, niter = 128, control = ctrl)
#' 
#' # Arid vegetation model 
#' mod <- ca_library("aridvege")
#' init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)
#' colors <- c("gray80", "white", "darkgreen")
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = landscape_plotter(mod, col = colors, xaxt = "n", 
#'                                                    yaxt = "n"))
#' run_camodel(mod, init, niter = 128, control = ctrl)
#' 
#' # Game of life 
#' mod <- ca_library("gameoflife") 
#' init <- generate_initmat(mod, c(0.5, 0.5), 100, 178) 
#' colors <- c("white", "black") 
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = landscape_plotter(mod, col = colors, 
#'                                                    xaxt = "n", 
#'                                                    yaxt = "n"))
#' run_camodel(mod, init, niter = 128, control = ctrl)
#' 
#' } 
#' 
#'@export
landscape_plotter <- function(mod, 
                              col = NULL, 
                              fps_cap = 24, 
                              burn_in = 0, 
                              transpose = TRUE, 
                              ...) { 
  
  # col is nothing, set default
  if ( is.null(col) ) { 
    col <- grDevices::hcl.colors(mod[["nstates"]], "viridis")
  }
  
  # col is a palette
  if ( is.character(col) && length(col) == 1 ) { 
    col <- grDevices::hcl.colors(mod[["nstates"]], col)
  }
  
  # Check that palette makes sense 
  if ( length(col) != mod[["nstates"]] ) { 
    stop(paste0("The number of colors in col does not correspond to the ", 
                "number of states in the model."))
  }
  
  # We open a new graphic
  
  last_call_time <- Sys.time()
  
  function(t, mat) { 
    
    # If we are before the burn_in iterations, return without doing anything
    if ( t < burn_in ) { 
      return( NULL ) 
    }

    this_call_time <- Sys.time()
    dtime <- as.numeric(difftime(this_call_time, last_call_time, units = "secs"))
    
    # If {1/fps_cap} has not passed since last time, we wait a little bit
    if ( dtime < (1/fps_cap) ) { 
      Sys.sleep( (1/fps_cap) - dtime )
    }
    
    if ( is.null(grDevices::dev.list()) ) { 
      grDevices::dev.new()
    }
    
    # Hold plot update, and flush it when exiting the function and everything is drawn
    grDevices::dev.hold() 
    on.exit({ 
      grDevices::dev.flush() 
    })
    
    if ( transpose ) { 
      mat <- t(mat)
    }
    
    # mat is a factor, so limits should be 1:nstates
    image.camodel_initmat(mat, col = col, zlim = c(1, mod[["nstates"]]), ...)
    last_call_time <<- Sys.time()
    
    return(NULL)
  }
  
}

#'@export
trace_plotter <- function(mod, initmat, 
                          fun = function(m) { 
                            sapply(mod[["states"]], function(s) mean(m == s))
                          }, 
                          col = NULL, 
                          max_samples = 256, 
                          burn_in = 0, 
                          ...) { 
  
  # col is nothing, set default
  if ( is.null(col) ) { 
    col <- grDevices::hcl.colors(mod[["nstates"]], "viridis")
  }
  
  ex_res <- fun(initmat)
  if ( length(ex_res) == 0 || ! is.atomic(ex_res) || ! is.numeric(ex_res) ) { 
    stop("The 'tracer_fun' function must return a vector of numeric values.") 
  }
  
  backlog <- matrix(NA_real_, ncol = length(ex_res) + 1, nrow = max_samples)
  backlog_line <- 1
  states <- mod[["states"]]
  
  function(t, mat) { 
    
    # backlog persists across runs, so we need to zero it out when re-running a 
    # simulation with the same control list
    if ( t == 0 ) { 
      backlog[] <<- NA_real_
    }
    
    # If we are before the burn_in iterations, return without doing anything
    if ( t < burn_in ) { 
      return( NULL ) 
    }
    
    # Compute covers and store them in backlog 
    backlog[backlog_line, ] <<- c(t, fun(mat))
    
    # Use only non-NA values and sort them by time
    backlog_sorted <- backlog[ ! is.na(backlog[ ,1]), , drop = FALSE]
    backlog_sorted <- backlog_sorted[order(backlog_sorted[ ,1]), , drop = FALSE]
    
    if ( nrow(backlog_sorted) > 1 ) { 
      
      if ( is.null(grDevices::dev.list()) ) { 
        grDevices::dev.new()
      }
      
      grDevices::dev.hold() 
      on.exit({ 
        grDevices::dev.flush() 
      })
      graphics::matplot(backlog_sorted[ ,1], 
                        backlog_sorted[ ,-1], 
                        col = col, 
                        type = "l", 
                        xlab = "time", 
                        ylab = "covers", 
                        ...)
    }
    
    # Rollback to one if we go above the max number of lines. We need to store it in 
    # the enclosing environment so it is persistent across calls. 
    backlog_line <<- backlog_line + 1
    if ( backlog_line > max_samples ) { 
      backlog_line <<- 1 
    }
    
    return(NULL)
  }
  
}
