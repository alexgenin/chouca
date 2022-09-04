# 
# Standard methods for ca_model objects. 
# 

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
  
  return(invisible(x))
}

#'@export
print.camodel_transition <- function(x, ...) { 
  cat(sprintf("Transition from %s to %s\n", x[["from"]], x[["to"]]))
  prob <- paste0(as.character(x[["prob"]]), collapse = " ")
  cat(sprintf("  %s\n", prob))
  
  invisible(x)
}
