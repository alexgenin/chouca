# 
# Interface to the compiled engine
# 


# Here we use the local() call so that the function can keep some internal state, which 
# is here the set of functions that have been compiled by the engine. 
camodel_compiled_engine_wrap <- local({
  
  # This environment is here to store the functions that will be created here 
  function_envir <- environment()
  
  function(ctrl, 
           console_callback, cover_callback, snapshot_callback) { 
  
  # Split betas into floats/integers matrices for c++ code, and augment the control 
  # list with them. 
  for ( tab in c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq") ) { 
    tabix <- ctrl[[tab]] 
    tabix <- tabix[ ,intersect(colnames(tabix), c("from", "to", "state_1", "state_2",
                                                  "qs", "expo_1", "expo_2")), 
                   drop = FALSE]
    ctrl[[paste0(tab, "_index")]] <- intmat(tabix)
    
    tabfl <- ctrl[[tab]] 
    tabfl <- tabfl[ ,intersect(colnames(tabfl), "coef"), drop = FALSE]
    ctrl[[paste0(tab, "_vals")]] <- tabfl
  }
  
  # Unwrap elements of the ctrl list that we need here 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["neighbors"]] == 8
  fixed_neighborhood <- ctrl[["fixed_neighborhood"]]
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  
  
  
  
  # Read file
  cmaxfile <- system.file("compiled_engine.cpp", package = "chouca")
  cmaxlines <- readLines(cmaxfile) 
  
  if ( ctrl[["verbose_compilation"]] ) { 
    cat("Compilation options:\n")
    gsubf <- function(a, b, lines) { 
      cat(sprintf("Setting %s to '%s'\n", a, b))
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
  cmaxlines <- gsubf("__BETA_0_NROW__",  format(nrow(ctrl[["beta_0"]])),  cmaxlines)
  cmaxlines <- gsubf("__BETA_Q_NROW__",  format(nrow(ctrl[["beta_q"]])), cmaxlines)
  cmaxlines <- gsubf("__BETA_PP_NROW__", format(nrow(ctrl[["beta_pp"]])), cmaxlines)
  cmaxlines <- gsubf("__BETA_QQ_NROW__", format(nrow(ctrl[["beta_qq"]])), cmaxlines)
  cmaxlines <- gsubf("__BETA_PQ_NROW__", format(nrow(ctrl[["beta_pq"]])), cmaxlines)
  cmaxlines <- gsubf("__COMMON_HEADER__", 
                     system.file("common.h", package = "chouca"), cmaxlines)
  
  # Decide whether we fixed the number of neighbors or not 
  cmaxlines <- gsubf("__FIXED_NEIGHBOR_NUMBER__", 
                     ifelse(fixed_neighborhood, "true", "false"), 
                     cmaxlines)
  
  # Set code lines that control the number of cores 
  cores <- ctrl[["cores"]]
  cmaxlines <- gsubf("__CORES__", format(cores), cmaxlines)
  cmaxlines <- gsubf("__USE_OMP__", ifelse(cores > 1, 1, 0), cmaxlines)
  
  # Set #define on whether to precompute transition probabilities 
  precompute_probas <- ctrl[["precompute_probas"]]
  if ( precompute_probas == "auto" ) { 
    precompute_probas <- ns^ifelse(use_8_nb, 8, 4) < prod(dim(init))
  }
  cmaxlines <- gsubf("__PRECOMP_PROBA_VALUE__", ifelse(precompute_probas, 1, 0), 
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
  hash <- digest::digest(list(cmaxlines, ctrl[["cores"]]), algo = "md5")
  cmaxlines <- gsub("__FPREFIX__", hash, cmaxlines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")
  
  # Make the table with all combinations of qs 
  if ( precompute_probas ) { 
    max_nb <- ifelse(use_8_nb, 8, 4)
    all_qs <- generate_all_qs(max_nb, ns, filter = wrap)
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
                      cleanupCacheDir = FALSE, 
                      env = function_envir)
  }
  
  runf <- get(fname, envir = function_envir)
  runf(all_qs, ctrl)
}
  
}) # end of local() block


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
#'@param precompute_probas The different values to be taken for precompute_probas (TRUE, 
#'  FALSE, or the vector c(TRUE, FALSE)).
#'
#'@param control Other options to use for the simulations. See the full list of options 
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
#' 
#'@export
benchmark_compiled_model <- function(mod, init, 
                                     control = list(), 
                                     niter = 100, 
                                     olevel = c("default", "Ofast"), 
                                     unroll_loops = c(TRUE, FALSE), 
                                     precompute_probas = c(TRUE, FALSE), 
                                     nrepeats = 1) { 
  
  all_combs <- expand.grid(olevel = olevel, 
                           unroll_loops = unroll_loops, 
                           rep = seq.int(nrepeats), 
                           precompute_probas = precompute_probas, 
                           stringsAsFactors = FALSE)
  
  timings <- plyr::ldply(seq.int(nrow(all_combs)), function(i) { 
    control[["olevel"]]       <- all_combs[i, "olevel"]
    control[["unroll_loops"]] <- all_combs[i, "unroll_loops"]
    control[["precompute_probas"]]       <- all_combs[i, "precompute_probas"]
    control[["engine"]]       <- "compiled"
    
    # warmup for compilation
    system.time( run_camodel(mod, init, niter, control = control) )
    
    # benchmark running time 
    timings <- system.time( run_camodel(mod, init, niter, control = control) )
    data.frame(all_combs[i, ], iter_per_s = niter / timings["elapsed"] )
  })
  
  # Make average 
  timings <- plyr::ddply(timings, ~ olevel + unroll_loops + precompute_probas, 
                         function(df) { 
    data.frame(iter_per_s = mean(df[ ,"iter_per_s"]))
  })
  
  # Order timings 
  timings <- timings[order(timings[ ,"iter_per_s"]), ]
  
  # Compute speedups 
  timings[ ,"speedup"] <- timings[ ,"iter_per_s"] / min(timings[ ,"iter_per_s"])
  
  return(timings)
}
