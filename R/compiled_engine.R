# 
# Interface to the compiled engine
# 

camodel_compiled_engine_wrap <- function(alpha, pmat, qmat, control_list, 
                                    console_callback, cover_callback, snapshot_callback) { 
  
  
  # Split alpha 
  alpha_index <- intmat(alpha[ ,c("from", "to"), drop = FALSE])
  alpha_vals  <- as.numeric(alpha[ ,c("a0")]) # vector of length nstates
  
  # Split pmat 
  pmat_index <- intmat(pmat[ ,c("from", "to", "state"), drop = FALSE])
  pmat_vals  <- pmat[ ,c("coef", "expo"), drop = FALSE]
  
  # Split qmat 
  qmat_index <- intmat(qmat[ ,c("from", "to", "state", "qs"), drop = FALSE])
  qmat_vals  <- qmat[ ,"ys"] # vector
  
  # Reduce pmat/qmat sizes to non-zero coefficients
  non_zero_pmat <- which(pmat_vals[ ,"coef"] > 1e-8)
  pmat_vals <- pmat_vals[non_zero_pmat, , drop = FALSE]
  pmat_index <- pmat_index[non_zero_pmat, , drop = FALSE]
  
  non_zero_qmat <- which(qmat_vals > 1e-8)
  qmat_vals <- qmat_vals[non_zero_qmat]
  qmat_index <- qmat_index[non_zero_qmat, , drop = FALSE]
  
  camodel_compiled_engine(alpha_index, 
                          alpha_vals, 
                          pmat_index, 
                          pmat_vals, 
                          qmat_index, 
                          qmat_vals, 
                          control_list, console_callback, cover_callback, snapshot_callback)
  
}


camodel_compiled_engine <- function(alpha_index, 
                                    alpha_vals, 
                                    pmat_index, 
                                    pmat_vals, 
                                    qmat_index, 
                                    qmat_vals, 
                                    ctrl, 
                                    console_callback, 
                                    cover_callback, 
                                    snapshot_callback) { 
  
  # Unwrap elements of the ctrl list 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["neighbors"]] == 8
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  
  # Read file
  cmaxfile <- system.file("compiled_engine.cpp", package = "chouca")
  cmaxlines <- readLines(cmaxfile) 
  
  if ( ctrl[["verbose_compilation"]] ) { 
  # 
    cat("Compilation options:\n")
    gsubf <- function(a, b, lines) { 
      cat(sprintf("Setting %s to %s\n", a, b))
      lines <- gsub(a, b, lines)
    }
  } else { 
    gsubf <- gsub
  }
  
  # Replace template values
  cmaxlines <- gsubf("__NR__", format(nrow(init)), cmaxlines)
  cmaxlines <- gsubf("__NC__", format(ncol(init)), cmaxlines)
  cmaxlines <- gsubf("__NS__", format(ns), cmaxlines)
  cmaxlines <- gsubf("__WRAP__", ifelse(wrap, "true", "false"), cmaxlines)
  cmaxlines <- gsubf("__USE_8_NB__", ifelse(use_8_nb, "true", "false"), cmaxlines)
  cmaxlines <- gsubf("__SUBSTEPS__", format(substeps), cmaxlines)
  cmaxlines <- gsubf("__XPOINTS__", format(ctrl[["xpoints"]]), cmaxlines)
  cmaxlines <- gsubf("__ALPHA_NROW__", format(nrow(alpha_index)), cmaxlines)
  cmaxlines <- gsubf("__PMAT_NROW__", format(nrow(pmat_index)), cmaxlines)
  cmaxlines <- gsubf("__QMAT_NROW__", format(nrow(qmat_index)), cmaxlines)
  cmaxlines <- gsubf("__COMMON_HEADER__", 
                     system.file("common.h", package = "chouca"), cmaxlines)
  
  # Set #define on whether to precompute transition probabilities 
  precompute_probas <- ctrl[["precompute_probas"]]
  if ( precompute_probas == "auto" ) { 
    precompute_probas <- ns^ifelse(use_8_nb, 8, 4) < prod(dim(init))
  }
  cmaxlines <- gsubf("__PRECOMPUTE_TRANS_PROBAS_DEFINE__", 
                    ifelse(precompute_probas, "#define PRECOMPUTE_TRANS_PROBAS", ""), 
                    cmaxlines)
  
  # Set GCC options on command line 
  # Emit optimize pragma 
  olvl <- gsub(" ", "", ctrl[["olevel"]])
  olvl_str <- ifelse(olvl == "default", "", 
                     sprintf('#pragma GCC optimize("%s")', olvl))
  cmaxlines <- gsub("__OLEVEL__", olvl_str, cmaxlines)
  
  # Emit loop unrolling pragma
  loop_unroll <- isTRUE(ctrl[["unroll_loops"]])
  loop_unroll_str <- ifelse(loop_unroll, '#pragma GCC optimize("unroll-loops")', "")
  cmaxlines <- gsub("__OUNROLL__", loop_unroll_str, cmaxlines)
  
  # Make hash of file and replace function name 
  # We make the hash depend on the model too, just in case the user changes models, but
  # the rest is different. Unlikely, but who knows.
  hash <- digest::digest(list(cmaxlines), algo = "md5")
  cmaxlines <- gsub("__FPREFIX__", hash, cmaxlines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")
  
  # Make the table with all combinations of qs 
  if ( precompute_probas ) { 
    max_nb <- ifelse(use_8_nb, 8, 4)
    all_qs <- rep( list(seq(0, max_nb) ), each = ns)
    # !!! very large matrix -> TODO: find a way to make this smaller directly, not 
    # after the fact when we subset it. 
    all_qs <- as.matrix(do.call(expand.grid, all_qs)) 
    # revert because the math expects the states in that order
    all_qs <- all_qs[ ,seq(ncol(all_qs), 1)]
    colnames(all_qs) <- rownames(all_qs) <- NULL
    all_qs <- cbind(all_qs, apply(all_qs, 1, sum))
    
    # If the number of neighbors is constant, we can discard the data that is present 
    # every 4 or 8 neighbors. 
    if ( wrap ) { 
      all_qs <- all_qs[seq(0, nrow(all_qs)-1) %% max_nb == 0, ]
    }
    
    all_qs <- all_qs[-1, ] # this has sum zero which produces division by zero. Discard it.
  } else { 
    # This is a dummy matrix just to make sure we pass something to the c++ function.
    all_qs <- matrix(0, nrow = 1, ncol = ns)
  }
  
  # Replace size in compiled code
  cmaxlines <- gsubf("__ALL_QS_NROW__", format(nrow(all_qs)), cmaxlines)
  
  # Source cpp if needed 
  if ( ! exists(fname) ) { 
    # We compile from the file, so that lines can be put in a debug run
    funs <- sourceCpp(code = paste(cmaxlines, collapse = "\n"), 
                      verbose = ctrl[["verbose_compilation"]], 
                      cleanupCacheDir = FALSE)
  }
  runf <- get(fname)
  runf(alpha_index, alpha_vals, 
       pmat_index, 
       pmat_vals, 
       qmat_index, 
       qmat_vals, 
       all_qs, 
       ctrl, 
       console_callback, 
       cover_callback, 
       snapshot_callback)
}


#'@title Benchmark a chouca model 
#'
#'@description Benchmark compilation options for a chouca model 
#'
#'@param mod A chouca model, as defined by \code{\link{camodel}}
#'
#'@param init An initial matrix, as produced by \code{\link{generate_initmat}}
#'
#'@param niter The number of iterations to use for the benchmark 
#'
#'@param olevel The optimizations levels to try (one or several of 'O0', 'O1', 'O2', 'O3', 
#' 'Ofast' or 'default'). 
#'
#'@param unroll_loops The loop unrolling options to try (one or both of TRUE and FALSE)
#'
#'@param control The options to use for the simulations. See the full list of options 
#'  documented in \code{\link{run_camodel}}. 
#'
#'@param nrepeats The number of samples to take when measuring run time. 
#'
#'@details 
#'  
#'  This function will take a chouca model, try the different compilation optimisation
#'    options of the 'compiled' engine (see \code{\link{run_camodel}}) and measure 
#'    run times. This allows deciding on which optimization options to use when running 
#'    simulations. 
#'    
#'@return
#' 
#' A data.frame whose lines are ordered by model run time with the following columns:
#' \enumerate{
#'   \item \code{olevel} The optimization level 
#'   
#'   \item \code{unroll_loops} The type of loop unrolling 
#'   
#'   \item \code{iter_per_s} The measured number of iteration per second
#'   
#'   \item \code{speedup} The speedup of a set of options compared to the slowest result. 
#' }
#' 
#'@examples 
#' 
#' mod <- ca_library("forestgap")
#' inimat <- generate_initmat(mod, c(0.5, 0.5), 1024) 
#' \dontrun{ 
#'   benchmark_compiled_model(mod, inimat, niter = 100) 
#' }
#'@export
benchmark_compiled_model <- function(mod, init, niter = 100, 
                                     olevel = c("O2", "O3", "Ofast"), 
                                     unroll_loops = c(TRUE, FALSE), 
                                     control = list(), 
                                     nrepeats = 1) { 
  
  all_combs <- expand.grid(olevel = olevel, 
                           unroll_loops = unroll_loops, 
                           rep = seq.int(nrepeats), 
                           stringsAsFactors = FALSE)
  
  timings <- plyr::ldply(seq.int(nrow(all_combs)), function(i) { 
    control[["olevel"]]       <- all_combs[i, "olevel"]
    control[["unroll_loops"]] <- all_combs[i, "unroll_loops"]
    control[["engine"]]       <- "compiled"
    # warmup for compilation
    system.time( run_camodel(mod, init, niter, control = control) )
    # benchmark running time 
    timings <- system.time( run_camodel(mod, init, niter, control = control) )
    data.frame(all_combs[i, ], iter_per_s = niter / timings["elapsed"] )
  })
  
  # Make average 
  timings <- plyr::ddply(timings, ~ olevel + unroll_loops, function(df) { 
    data.frame(iter_per_s = mean(df[ ,"iter_per_s"]))
  })
  
  # Order timings 
  timings <- timings[order(timings[ ,"iter_per_s"]), ]
  
  # Compute speedups 
  timings[ ,"speedup"] <- timings[ ,"iter_per_s"] / min(timings[ ,"iter_per_s"])
  
  return(timings)
}

