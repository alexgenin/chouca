# 
# This file makes sure that the compiled model produces the same simulations given 
# the same random seed. In can be different from the cpp simulations because it uses 
# a different PRNG, but the seed should be respected. 
# 

mod <- ca_library("forestgap") 
init <- generate_initmat(mod, c(0.5, 0.5), nr = 64)

set.seed(123)
run1 <- run_camodel(mod, init, seq(0, 8), control = list(engine = "compiled", 
                                                         console_output_every = 0))

set.seed(123)
run2 <- run_camodel(mod, init, seq(0, 8), control = list(engine = "compiled", 
                                                         console_output_every = 0))

expect_true({
  all( abs( run1[["output"]][["covers"]] - run2[["output"]][["covers"]] ) < 1e-8 )
})
