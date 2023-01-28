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
  wrap       <- ctrl[["wrap"]]
  continuous <- ctrl[["continuous"]]
  use_8_nb   <- ctrl[["neighbors"]] == 8
  fixed_neighborhood <- ctrl[["fixed_neighborhood"]]
  init       <- ctrl[["init"]]
  niter      <- ctrl[["niter"]]
  ns         <- ctrl[["nstates"]]
  
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
  boolstr <- function(x) ifelse(x, "true", "false")
  
  # Replace template values
  cmaxlines <- gsubf("__NR__", format(nrow(init)), cmaxlines)
  cmaxlines <- gsubf("__NC__", format(ncol(init)), cmaxlines)
  cmaxlines <- gsubf("__NS__", format(ns), cmaxlines)
  cmaxlines <- gsubf("__WRAP__", boolstr(wrap), cmaxlines)
  cmaxlines <- gsubf("__CONTINUOUS_SCA__", boolstr(continuous), cmaxlines)
  cmaxlines <- gsubf("__DELTA_T__", ctrl[["delta_t"]], cmaxlines)
  cmaxlines <- gsubf("__USE_8_NB__", boolstr(use_8_nb), cmaxlines)
  cmaxlines <- gsubf("__SUBSTEPS__", format(ctrl[["substeps"]]), cmaxlines)
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
  cmaxlines <- gsubf("__PRECOMP_PROBA_VALUE__", boolstr(precompute_probas), cmaxlines)
  
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

