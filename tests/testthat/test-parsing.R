# Test the parsing of models 


# Test that a warning is emitted when probabilities can be above one
expect_warning({ 
  forestgap(parms = list(d = 10, 
                         delta = 0.5, 
                         alpha = 0.2))
})

# Test that a warning is emitted when probabilities can be below zero
expect_warning({ 
  forestgap(parms = list(d = -1, 
                         delta = 0.5, 
                         alpha = 0.2))
})


# Test updating of models. Internal coefficient matrix should have changed 
#TODO
mod <- ca_library("forestgap")
mod2 <- update(mod, parms = list(d = 0, delta = 0.5, alpha = 0.01))



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


