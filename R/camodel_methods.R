# 
# Standard methods for ca_model objects. 
# 

ERROR_ABS_MAX <- 0.001 # ERROR_REL_MAX defined in transition.R

#'@export
print.ca_model <- function(x, ...) { 
  
  force(x) # force the evaluation of x before printing
  cat0 <- function(...) cat(paste0(...), "\n")
  
  cat0("Stochastic Cellular Automaton")
  cat0("")
  cat0("States: ", paste(x[["states"]], collapse = " "))
  cat0("")
  
  for ( tr in x[["transitions"]] ) { 
    cat0("Transition: ", tr[["from"]], " -> ", tr[["to"]])
    cat0(paste0("  ", as.character(tr[["prob"]]) ))
  }
  
  cat0("")
  cat0("Neighborhood: ", x[["neighbors"]], "x", x[["neighbors"]])
  cat0("Wrap: ", x[["wrap"]])
  cat0("Max error: ", format(max(x[["max_error"]])), " (", 
       ifelse(max(x[["max_error"]]) < ERROR_ABS_MAX, 
              "OK", "WARNING"), ")" )
  
  cat0("Max rel error: ", format(max(x[["max_rel_error"]])), " (", 
       ifelse(max(x[["max_rel_error"]]) < ERROR_REL_MAX, 
              "OK", "WARNING"), ")" )
  
  return(invisible(x))
}

#'@export
print.camodel_transition <- function(x, ...) { 
  cat(sprintf("Transition from %s to %s\n", x[["from"]], x[["to"]]))
  prob <- paste0(as.character(x[["prob"]]), collapse = " ")
  cat(sprintf("  %s\n", prob))
  
  invisible(x)
}




#'@export 
summary.ca_model_result <- function(object, ...) { 
  
  cat("Stochastic Cellular Automaton simulation\n") 
  cat("\n")
  cat(sprintf("Iterations: %s\n", object[["niter"]]))
  cat(sprintf("Landscape size: %sx%s\n", 
              nrow(object[["initmat"]]), ncol(object[["initmat"]])))
  
  covers <- object[["output"]][["covers"]]
  if ( ! is.null(covers) ) { 
    cat("Global covers:\n")
    print(utils::tail(covers, n = 6))
  } else { 
    cat("Covers were not saved during simulation\n")
  }
  
  snapshots <- object[["output"]][["snapshots"]]
  if ( ! is.null(covers) ) { 
    cat(sprintf("Landscapes saved: %s\n", length(snapshots)))
  } else { 
    cat("Landscapes were not saved during simulation\n")
  }
  
  cat("\n")
  cat("The following methods are available: \n")
  cat(list_methods(class(object)[1]), "\n")

  return(invisible(object))
}

#'@export 
plot.ca_model_result <- function(x, 
                                 legend.x = NULL, 
                                 legend.y = NULL, 
                                 ...) { 
  
  covers <- x[["output"]][["covers"]]
  
  if ( is.null(covers) ) { 
    stop("Covers are not present in the model result object, and cannot be displayed. Adjust your simulation parameters to make sure they are saved")
  }
  
  
  colors <- seq.int(ncol(covers)-1)
  graphics::matplot(covers[ ,1], 
                    covers[ ,-1],
                    type = "l", 
                    col = colors, 
                    xlab = "Time", 
                    ylab = "Covers", 
                    ...)
  
  legend.x <- ifelse(is.null(legend.x), min(covers[ ,1]), legend.x) # max t 
  legend.y <- ifelse(is.null(legend.y), max(covers[ ,-1]), legend.y) # max covers 
  
  graphics::legend(legend.x, legend.y, 
                   legend = colnames(covers)[-1], 
                   lty = 1,
                   col = colors)
  
}



#'@export 
image.ca_model_result <- function(x, snapshot_time = "last", ...) { 
  
  snapshots <- x[["output"]][["snapshots"]]
  if ( is.null(snapshots) ) { 
    stop("Snapshots are not present in the model result object, and cannot be displayed. Adjust your simulation parameters to make sure they are saved")
  }
  
  
  if ( snapshot_time == "last" ) { 
    snap_to_show <- snapshots[[length(snapshots)]]
  } else { 
    if ( ! is.numeric(snapshot_time) ) { 
      stop("Argument 'snapshot_time' must be a numeric value, or the string 'last'")
    }
    
    times_available <- plyr::laply(snapshots, function(o) attr(o, "t"))
    time_closest <- abs(times_available - snapshot_time) 
    time_closest <- times_available[time_closest == min(time_closest)]
    
    if ( time_closest != snapshot_time ) { 
      warning(sprintf("Exact snapshot not found for t=%s in the saved data, returning closest one at t=%s", snapshot_time, time_closest))
    }
    
    snap_to_show <- snapshots[[which(times_available == time_closest)[1]]]
  }
  
  image.camodel_initmat(snap_to_show, ...)
  graphics::title(sprintf("Landscape at t=%s", attr(snap_to_show, "t")))
}

# image() does not work with matrices of factors, so we define a method here for the 
# camodel_initmat class. 
# 
#'@export
image.camodel_initmat <- function(x, ...) { 
  x <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
  graphics::image(x, ...)
}



# List the methods available for an S3 class
list_methods <- function(class, 
                         exclude = c("print", "summary")) { 
  
  all_methods <- lapply(class, function(class) { 
    tab <- attr(utils::methods(class = class), "info")
#     tab[tab[ ,"from"] == "spatialwarnings", "generic"]
    tab[ ,"generic"]
  })
  all_methods <- sort(unique(unlist(all_methods)))
  
  # Exclude some reported methods 
  all_methods <- all_methods[! all_methods %in% exclude]
  
  return(all_methods)
}
