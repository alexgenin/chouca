# 
# 
# 


# Plot a CAmodel: represent the transitions as functions of q/p[k]
#'@export
plot.ca_model <- function(x, ...) { 
  
  zero <- rep(0, length(x[["states"]]))
  
  # For all transitions, compute the curve that the model uses
  mod_trans <- expand.grid(from = x[["states"]], to = x[["states"]], 
                           covstate = x[["states"]])
  
  # Dependency of probability to the local cover q
  mod_transq <- ddply(mod_trans, ~ from + to + covstate, function(tinfo) { 
    x0s <- seq(0, 1, length.out = 25)
    transa0 <- subset(x[["transa0"]], from == tinfo[ ,"from"] & to == tinfo[ ,"to"])
    transq <- subset(x[["transq"]], { 
      from == tinfo[ ,"from"] & to == tinfo[ ,"to"] & covstate == tinfo[ ,"covstate"]
    })
    
    # No such transition 
    if ( nrow(transa0) == 0 && nrow(transq) == 0 ) { 
      return( NULL ) # ddply will discard the results
    }
    
    ys <- sapply(x0s, function(x0) { 
        transa0[ ,"a0"] + 
              sum(transq[ ,"coef"] * x0 ^ transq[ ,"expo"])
    })
    
    data.frame(xs = x0s, ys = ys)
  })
  
  ggplot(NULL, aes(x = xs, y = ys)) + 
    geom_step(aes(color = covstate), data = mod_transq) + 
    geom_point(aes(color = covstate), data = mod_transq) + 
    facet_grid( from ~ to ) 
  
}
