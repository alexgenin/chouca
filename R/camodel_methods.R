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
plot.ca_model <- function(x, y, ...) { 
  if ( requireNamespace("igraph", quietly = TRUE) ) { 
    gr <- igraph::as.igraph(x)
    plot(gr, edge.label = igraph::E(gr)$name, 
#          layout = igraph::layout_with_kk, 
         edge.label.family = "sans", 
         vertex.label.family = "sans", 
         edge.curved = 0.4,
         vertex.color = grDevices::rgb(0, 0, 0, 0),
         vertex.size = 30, 
         ...)
  }
}

as.igraph.ca_model <- function(x, ...) { 
  trans_df <- plyr::ldply(x[["transitions"]], function(tr) { 
    data.frame(from = tr[["from"]], 
               to   = tr[["to"]], 
               name = as.character(tr[["prob"]])[-1])
  })
  igraph::graph_from_data_frame(trans_df, ...)
}

#'@export
print.camodel_transition <- function(x, ...) { 
  cat(sprintf("Transition from %s to %s\n", x[["from"]], x[["to"]]))
  prob <- paste0(as.character(x[["prob"]]), collapse = " ")
  cat(sprintf("  %s\n", prob))
  
  invisible(x)
}



#'@export 
summary.ca_model_result <- function(object, .name = NULL, ...) { 
  
  cat("Stochastic Cellular Automaton simulation\n") 
  cat("\n")
  cat(sprintf("Times: %s to %s\n", min(object[["times"]]), max(object[["times"]])))
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
  if ( ! is.null(snapshots) ) { 
    cat(sprintf("Landscapes saved: %s\n", length(snapshots)))
  } else { 
    cat("Landscapes were not saved during simulation\n")
  }
  
  cat("\n")
  cat("The following methods are available: \n")
  cat("  ", list_methods(class(object)[1]), "\n")
  
  # Display methods to extract output 
  objname <- ifelse(is.null(.name), as.character(substitute(x)), .name)
  cat("\n")
  cat("Extract simulation results using one of:\n")
  for ( name in names(object[["output"]]) ) { 
    cat(sprintf("  %s[[\"output\"]][[\"%s\"]]\n", 
                as.character(substitute(object)), name))
  }
  
  return( invisible(object) )
}

#'@export 
print.ca_model_result <- function(x, ...) { 
  summary.ca_model_result(x)
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
                    lty = 1, 
                    ...)
  
  legend.x <- ifelse(is.null(legend.x), min(covers[ ,1]), legend.x) # max t 
  legend.y <- ifelse(is.null(legend.y), max(covers[ ,-1]), legend.y) # max covers 
  
  graphics::legend(legend.x, legend.y, 
                   legend = colnames(covers)[-1], 
                   lty = 1,
                   col = colors)
  
}


# Define levels() methods to access states
#'@export
levels.ca_model <- function(x) { 
  x[["states"]]
}

#'@export
levels.ca_model_result <- function(x) { 
  levels.ca_model(x[["mod"]][["states"]])
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
    # This can have several values, always take the first one
    time_closest <- times_available[time_closest == min(time_closest)][1]
    
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

#'@export
summary.camodel_initmat <- function(object, ...) { 
  cat("Cellular Automaton landscape\n") 
  cat("\n")
  cat(sprintf("Size: %sx%s\n", nrow(object), ncol(object)))
  
  cat("\n")
  xfmt <- format(object, nmax = 5)
  cat(xfmt)
  
  rhos <- as.numeric(table(object))/(nrow(object)*ncol(object))
  names(rhos) <- levels(object)
  cat("\n")
  cat("Global covers\n")
  print(rhos)
  
  cat("\n")
  cat("The following methods are available: \n")
  cat("  ", list_methods(class(object)[1]), "\n")
}

format.camodel_initmat <- function(x, nmax, ...) { 
  
  xfmt <- x[1:min(nmax, nrow(x)), 1:min(nmax, ncol(x))]
  xfmt <- matrix(as.character(levels(x)[xfmt]), nrow(xfmt), ncol(xfmt))
  
  # Add dots for extra cols
  if ( ncol(x) > nmax ) { 
    xfmt <- cbind(xfmt, rep("...", nrow(xfmt)))
  }
  
  # Add dots for extra rows
  if ( nrow(x) > nmax ) { 
    xfmt <- rbind(xfmt, rep("...", ncol(xfmt)))
  }
  
  # Collapse into string representations 
  xfmt <- apply(format(xfmt), 1, paste, collapse = " ")
  xfmt <- paste(xfmt, collapse = "\n")
  xfmt <- paste(xfmt, "\n")
  
  return(xfmt)
}


# List the methods available for an S3 class
list_methods <- function(class, 
                         exclude = c("print", "summary", "format")) { 
  
  all_methods <- lapply(class, function(class) { 
    tab <- attr(utils::methods(class = class), "info")
    tab[ ,"generic"]
  })
  all_methods <- sort(unique(unlist(all_methods)))
  
  # Exclude some reported methods 
  all_methods <- all_methods[! all_methods %in% exclude]
  
  return(all_methods)
}

