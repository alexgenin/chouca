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
