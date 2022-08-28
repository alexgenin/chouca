# Test the parsing of models 


# Test that a warning is emitted when probabilities can be below zero
expect_warning({ 
  ca_library("forestgap", parms = list(d = -1, 
                                       delta = 0.5, 
                                       alpha = 0.2))
})

# Error in reserved names for parameters
expect_error({ 
  camodel(transition(from = "a", to = "b", ~ 1), 
          parms = list(p = 1), 
          neighbors = 4, 
          wrap = TRUE)
}, "no parameters must be named")

# Error when things not defined using transition()
expect_error({ 
  camodel(list(a = "fdsqfdqs", ~ 1), 
          neighbors = 4, 
          wrap = TRUE)
}, "model transition was not processed by the transition()")

# Test updating of models. Internal coefficient matrix should have changed 
mod <- ca_library("forestgap")
mod2 <- update(mod, parms = list(d = 0, delta = 0.5, alpha = 0.01))
expect_true({ 
  any(mod[["beta_p"]] != mod2[["beta_p"]])
})


# Test printing of camodel() objects 
mod <- ca_library("forestgap")
expect_true({ 
  any( grepl("States: EMPTY TREE", capture.output( print(mod) ) ) )
})

tr <- transition(from = "a", to = "b", ~ 1 )
expect_true({ 
  any( grepl("Transition from a to b", capture.output( print(tr) ) ) )
})



# Test generate_initmat()
mod <- ca_library("forestgap")

# Wrong number of states -> error 
expect_error({ 
  generate_initmat(mod, c(1, 0, 1), 2, 2)
})
# Named vector, but not the right names -> error
expect_error({ 
  generate_initmat(mod, c(A = 1, B = 1), 2, 2)
})
# NAs -> error
expect_error({ 
  generate_initmat(mod, c(A = NA, B = 1), 2, 2)
})
# Sum of states above/below 1 -> warning
expect_warning({ 
  generate_initmat(mod, c(1, 2), 2, 2)
})
expect_warning({ 
  generate_initmat(mod, c(0, 0.1), 2, 2)
})


# Make a quick run, extract the first step and make sure it matches what we expect 
im <- generate_initmat(mod, c(TREE = 0.2, EMPTY = 0.8), 100, 100)
props <- sapply(c("EMPTY", "TREE"), function(s) mean(im == s))
expect_true( abs(props["TREE"] - 0.2) < 0.05 )
expect_true( abs(props["EMPTY"] - 0.8) < 0.05 )




