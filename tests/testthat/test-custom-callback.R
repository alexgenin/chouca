#
# This tests whether custom callbacks are properly supported
#


test_that("Trace plotting works", {

  nrows <- ncols <- 16
  niters <- 64

  # Test compiled engine against cpp engine. Here it works because the model is
  # deterministic so the values are exactly equal. We only check that spatial
  # autocorrelations are equal.
  mod <- ca_library("rockpaperscissor", neighbors = 8, wrap = TRUE)
  initmm <- generate_initmat(mod, rep(1/3, 3), nrows, ncols)

  ccb <- function(t, mat) {
    bmat <- matrix(mat == "r", nrow = nrow(mat), ncol = ncol(mat))
    data.frame(t = t,
               cover = mean(mat == "r"),
               sd  = sd(mat == "r"))
  }

  ctrl <- list(engine = "cpp",
               console_output_every = 0,
               custom_output_every = 1,
               custom_output_fun = ccb)

  ctrl[["engine"]] <- "cpp"
  out_cpp <- run_camodel(mod, initmm, seq(0, niters), ctrl)[["output"]][["custom"]]
  out_cpp <- plyr::rbind.fill(out_cpp)

  ctrl[["engine"]] <- "compiled"
  out_compiled <- run_camodel(mod, initmm, seq(0, niters), ctrl)[["output"]][["custom"]]
  out_compiled <- plyr::rbind.fill(out_compiled)

  expect_true( all( abs(out_cpp - out_compiled) < 1e-10 ) )

})

test_that("Landscape plotting works", {


  mod <- ca_library("aridvege", neighbors = 4, wrap = TRUE)
  initmm <- generate_initmat(mod, c(0, 0.5, 0.5), 32, 32)



  control <- list(console_output_every = 0,
                  custom_output_every = 1,
                  custom_output_fun = trace_plotter(mod, initmm,
                                                    lty = 1,
                                                    mar = rep(6, 4), # par() arg
                                                    fps_cap = 5,
                                                    max_samples = 4),
  #                 custom_output_fun = landscape_plotter(mod, fps_cap = 5),
                  engine = "cpp")

  tmax <- 16
  aa <- run_camodel(mod, initmm, times = seq(1, tmax), control = control)

  # Make sure custom output is printed
  expect_true({ any(
    grepl("aa[[\"output\"]][[\"custom\"]]",
          capture.output(summary(aa)),
          fixed = TRUE)
  )})

  expect_true({
    length(aa[["output"]][["custom"]]) == tmax
  })

  expect_true({
    all(sapply(aa[["output"]][["custom"]], is.null))
  })



  control <- list(console_output_every = 0,
                  custom_output_every = 1,
                  custom_output_fun = landscape_plotter(mod,
                                                        col = c("red", "blue", "green"),
                                                        transpose = TRUE,
                                                        mar = rep(0, 4),
                                                        fps_cap = 5),
                  engine = "compiled")

  aa <- run_camodel(mod, initmm, times = seq(1, 4), control = control)

  # Make sure custom output is printed
  expect_true({ any(
    grepl("aa[[\"output\"]][[\"custom\"]]",
          capture.output(summary(aa)),
          fixed = TRUE)
  )})

  expect_true({
    length(aa[["output"]][["custom"]]) == 4
  })

  expect_true({
    all(sapply(aa[["output"]][["custom"]], is.null))
  })

})

test_that("We can re-run plotting without redefining control list", { 
  # This could produce errors when landscape_plotter or trace_plotter did not maintain
  # their internal state correctly, so make sure it runs 
  
  # Close all devices 
  while ( ! is.null(plyr::tryNULL(dev.off())) ) { 
    1
  }
  
  mod <- ca_library("rock-paper-scissor")
  init <- generate_initmat(mod, rep(1, 3)/3, nrow = 100, ncol = 178)
  
  # Trace plotter
  ctrl <- list(custom_output_every = 1,
                custom_output_fun = trace_plotter(mod, init, new_window = FALSE))
  run_camodel(mod, init, times = seq(0, 2), control = ctrl)
  dev.off()
  run_camodel(mod, init, times = seq(0, 2), control = ctrl)
  
  
  # Landscape plotter
  ctrl <- list(custom_output_every = 1,
               custom_output_fun = landscape_plotter(mod, new_window = FALSE))
  run_camodel(mod, init, times = seq(0, 2), control = ctrl)
  dev.off()
  run_camodel(mod, init, times = seq(0, 2), control = ctrl)
  
  expect_true(TRUE)
})

test_that("Graphical parameters are unchanged", { 
  
  par(mar = rep(2, 4))
  plot(1:10, 1:10 + rnorm(10)*2)
  opar <- par(no.readonly = TRUE)
  odev <- dev.cur()
  
  mod <- ca_library("rock-paper-scissor")
  init <- generate_initmat(mod, rep(1, 3)/3, nrow = 100, ncol = 178)
  
  # Display covers of the rock/paper/scissor model as it is running. Here we assume that 
  # regardless of whether we use a new window or an old one, at the end of the
  # simulation, the active device should be the old one with its old pars
  lapply(c(TRUE, FALSE), function(new_win) { 
    
    # Trace plotter
    ctrl <- list(custom_output_every = 1,
                 custom_output_fun = trace_plotter(mod, init, new_window = new_win))
    run_camodel(mod, init, times = seq(0, 2), control = ctrl)
    
    npar <- par(no.readonly = TRUE)
    ndev <- dev.cur()
    expect_true({ 
      all(c(all.equal(opar, npar), 
            all.equal(odev, ndev)))
    })
    
    ctrl <- list(custom_output_every = 1,
                 custom_output_fun = landscape_plotter(mod, new_window = FALSE))
    run_camodel(mod, init, times = seq(0, 2), control = ctrl)
    
    npar <- par(no.readonly = TRUE)
    ndev <- dev.cur()
    expect_true({ 
      all(c(all.equal(opar, npar), 
            all.equal(odev, ndev)))
    })
    
  })
  
  expect_true({ 
    all.equal(opar, par(no.readonly = TRUE))
  })
  
  expect_true({ 
    all.equal(odev, dev.cur())
  })
  
})
