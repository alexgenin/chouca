# 
# Make sure the number of substeps does not change the output of the simulation 
# 

test_that("The number of substeps is correctly implemented", { 
  
  p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
  mod <- update(ca_library("musselbed"), parms = p, 
                ctrl = list(console_output_every = 0))
  init <- generate_initmat(mod, c(0.5, 0.4, 0.1), 64)
  
  ctrl <- list(substeps = 10)
  out1 <- run_camodel(mod, init, 100, ctrl)
  ctrl <- list(substeps = 20)
  out2 <- run_camodel(mod, init, 100, ctrl)

  with(out1[["output"]], plot(covers[ ,"t"],  covers[ , "MUSSEL"], type = "n"))
  with(out2[["output"]], lines(covers[ ,"t"], covers[ , "MUSSEL"], col = "red"))
  with(out1[["output"]], lines(covers[ ,"t"], covers[ , "MUSSEL"], col = "blue"))
  
  expect_true({ 
    all( abs(out1[["output"]][["covers"]] - out2[["output"]][["covers"]]) < 0.2 )
  })
  
})

