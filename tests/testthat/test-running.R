# 
# This file contains tests to check 
# 

# Bad init matrix (states do not match) -> error
expect_error({ 
  initmat <- generate_initmat(aridvege(), c(0.5, 0.5), 10, 10)
  mod <- forestgap()
  run_camodel(mod, initmat, 100)
})

