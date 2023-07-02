# 
# 
# These functions are function factories for custom callbacks
# 
# 

#' @title Plot simulation landscapes 
#' 
#' @description This function creates another function that plot the model 
#'   landscape during the of simulation a stochastic cellular automaton. 
#' 
#' @param mod The model being used (created by \code{link{camodel}})
#' 
#' @param col a set of colors (character vector) of length equal to the number of 
#'   states in the model. 
#' 
#' @param fps_cap The maximum number of frame per seconds at which to plot the results 
#' 
#' @param burn_in A number of iterations to skip before plotting
#' 
#' @param transpose Set to \code{TRUE} to transpose the landscape matrix 
#'   before displaying it \code{\link[graphics]{image}}
#' 
#' @param new_window Controls whether the plots are displayed in a new window, 
#'   or in the default device (typically the plot panel in Rstudio)
#' 
#' @param ... other arguments are passed to \code{\link[graphics]{image}}
#' 
#' @details
#' 
#'   This function creates another function that is suitable for use with 
#'   \code{\link{run_camodel}}. It allows plotting the landscape as it is being 
#'   simulated, using the base function \code{\link[graphics]{image}}. You can set 
#'   colors using the argument \code{col}, or tranpose the landscape before 
#'   plotting using \code{transpose}. The resulting function must by passed to 
#'   \code{\link{run_camodel}} as the control argument \code{custom_output_fun}. 
#'   
#'   \code{\link[graphics]{image}} tends to be quite slow at displaying matrices, 
#'   especially if the matrix is large, but if it is still too fast for your visualize 
#    your results, you can cap the number of landscapes displayed per seconds by setting 
#    the argument \code{fps_cap}. 
#'   
#'   It is important to note that this function will probably massively slow down a 
#'   simulation, so this is mostly for exploratory analyses, or just to 
#'   have a good look of what is happening in a model. 
#' 
#' @seealso trace_plotter, run_camodel
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
#' init <- generate_initmat(mod, rep(1, 3)/3, nr = 100, nc = 178)
#' run_camodel(mod, init, times = seq(0, 128), control = ctrl)
#' 
#' # Arid vegetation model 
#' mod <- ca_library("aridvege")
#' init <- generate_initmat(mod, rep(1, 3)/3, nr = 100, nc = 178)
#' colors <- c("gray80", "white", "darkgreen")
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = landscape_plotter(mod, col = colors, xaxt = "n", 
#'                                                    yaxt = "n"))
#' run_camodel(mod, init, times = seq(0, 128), control = ctrl)
#' 
#' # Game of life, turning off plot decoration using extra arguments read by matplot()
#' mod <- ca_library("gameoflife") 
#' init <- generate_initmat(mod, c(0.5, 0.5), nr = 100, nc = 178)
#' colors <- c("white", "black") 
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = landscape_plotter(mod, col = colors, 
#'                                                    xaxt = "n", 
#'                                                    yaxt = "n"))
#' run_camodel(mod, init, times = seq(0, 128), control = ctrl)
#' 
#' } 
#' 
#'@export
landscape_plotter <- function(mod, 
                              col = NULL, 
                              fps_cap = 24, 
                              burn_in = 0, 
                              transpose = TRUE, 
                              new_window = TRUE, 
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
  
  last_call_time <- Sys.time()
  has_set_par <- FALSE
  
  function(t, mat) { 
    
    # If we are before the burn_in iterations, return without doing anything
    if ( t < burn_in ) { 
      return( NULL ) 
    }
    
    if ( ! has_set_par ) { 
      setup_par(new_window)
      has_set_par <<- TRUE
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

#' @title Plot simulation covers 
#' 
#' @description This function returns another helper function to display covers or any 
#'   other metric as a stochastic cellular automaton simulation is being run by 
#'   \code{\link{run_camodel}}. 
#' 
#' @param mod The model being used (created by \code{link{camodel}})
#' 
#' @param initmat The initial landscape given to \code{\link{run_camodel}} 
#' 
#' @param fun The function used to summarise the landscape into summary metrics. By 
#'   default, gloal covers of each state are computed. It must return a vector of 
#'   numeric values. 
#' 
#' @param col a set of colors (character vector) of length equal to the number of 
#'   values returned by fun.
#' 
#' @param fps_cap The maximum number of frame per seconds at which to plot the results 
#' 
#' @param max_samples The maximum number of samples to display in the plot 
#' 
#' @param burn_in A number of iterations to skip before plotting
#' 
#' @param new_window Controls whether the plots are displayed in a new window, 
#'   or in the default device (typically the plot panel in Rstudio)
#' 
#' @param ... other arguments passed to \code{\link[graphics]{matplot}}, which is used to 
#'   display the trends. 
#' 
#' @details 
#' 
#'   This function creates another function that is suitable for use with 
#'   \code{\link{run_camodel}}. It allows plotting any matric computed on the landscape 
#'   as it is being simulated, using the base function \code{\link[graphics]{matplot}}. 
#'   The resulting function must by passed to \code{\link{run_camodel}} as the control
#'   argument \code{custom_output_fun}.
#'   
#'   By default, the global covers of each state in the landscape will be displayed, but 
#'   you can pass any function as argument \code{fun} to compute something else, as long 
#'   as \code{fun} returns a numeric vector of length at least 1. 
#'   
#'   \code{\link[graphics]{matplot}} tends to be quite slow at displaying results, but 
#'   if it is still too fast for your taste, you can cap the refresh rate by setting 
#    the argument \code{fps_cap}. 
#'   
#'   It is important to note that this function will probably massively slow down a 
#'   simulation, so it is mostly here for exploratory analyses, or just to 
#'   have a good look of what is happening in a model. 
#' 
#' 
#' @seealso landscape_plotter, run_camodel
#' 
#' # Display covers of the rock/paper/scissor model as it is running 
#' mod <- ca_library("rock-paper-scissor")
#' init <- generate_initmat(mod, rep(1, 3)/3, nr = 100, nc = 178)
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = trace_plotter(mod, init))
#' run_camodel(mod, init, times = seq(0, 256), control = ctrl)
#' 
#' # Adjust colors of the previous example and increase speed 
#' colors <- c("#fb8500", "#023047", "#8ecae6")
#' ctrl <- list(custom_output_every = 1, 
#'              custom_output_fun = trace_plotter(mod, init, fps_cap = 60, col = colors))
#' run_camodel(mod, init, times = seq(0, 300), control = ctrl)
#' 
#' # Display statistics on autocorrelation on the fly for the arid 
#' # vegetation model 
#' if ( requireNamespace("spatialwarnings") ) { 
#'   mod <- ca_library("aridvege")
#'   init <- generate_initmat(mod, rep(1, 3)/3, nr = 100, nc = 178)
#'   moran <- function(mat) {
#'     m <- matrix(mat == "VEGE", nrow = nrow(mat), ncol = ncol(mat))
#'     spatialwarnings::raw_moran(m)
#'   }
#'   ctrl <- list(custom_output_every = 1, 
#'                custom_output_fun = trace_plotter(mod, init, 
#'                                                  fun = moran))
#'   run_camodel(mod, init, times = seq(0, 128), control = ctrl)
#' }
# 
#'@export
trace_plotter <- function(mod, initmat, 
                          fun = NULL, 
                          col = NULL, 
                          max_samples = 256, 
                          fps_cap = 24, 
                          burn_in = 0, 
                          new_window = TRUE, 
                          ...) { 
  
  if ( is.null(fun) ) { 
    fun <- function(m) { 
      sapply(mod[["states"]], function(s) mean(m == s))
    }
  }
  
  ex_res <- fun(initmat)
  if ( length(ex_res) == 0 || ! is.atomic(ex_res) || ! is.numeric(ex_res) ) { 
    stop("The function 'fun' must return a vector of numeric values.") 
  }
  
  # col is nothing, set default
  if ( is.null(col) ) { 
    col <- grDevices::hcl.colors(length(ex_res), "viridis")
  }
  
  if ( length(col) != length(ex_res) ) { 
    stop("The number of colors does not match the number of values returned by 'fun'")
  }
  
  backlog <- matrix(NA_real_, ncol = length(ex_res) + 1, nrow = max_samples)
  backlog_line <- 1
  states <- mod[["states"]]
  has_set_par <- FALSE
  
  last_call_time <- Sys.time()

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
    
    if ( ! has_set_par ) { 
      setup_par(new_window)
      has_set_par <<- TRUE
    }

    this_call_time <- Sys.time()
    dtime <- as.numeric(difftime(this_call_time, last_call_time, units = "secs"))
    
    # If {1/fps_cap} has not passed since last time, we wait a little bit
    if ( dtime < (1/fps_cap) ) { 
      Sys.sleep( (1/fps_cap) - dtime )
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
    
    last_call_time <<- Sys.time()
    
    return(NULL)
  }
  
}


setup_par <- function(external) { 
  devtype <- NA
  ostype <- tolower(Sys.info()["sysname"])
  if ( grepl("^linux", ostype) || grepl("^solaris", ostype) ) { 
    devtype <- "x11"
  } else if ( grepl("^darwin", ostype) ) { 
    devtype <- "quartz"
  } else if (grepl("^win", ostype) ) { 
    devtype <- "windows"
  }
  
  if ( external ) { 
    grDevices::dev.new(devtype)
  }
  
  if ( is.na(devtype) ) { 
    stop("Could not detect OS type to set up graphics (please report a bug!)")
  }
  
  graphics::par(mfrow = c(1, 1))
}
