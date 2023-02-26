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
  ns         <- ctrl[["nstates"]]
  
  # Read file
  cppfile <- system.file("compiled_engine.cpp", package = "chouca")
  clines <- readLines(cppfile) 
  
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
  clines <- gsubf("__NR__", format(nrow(init)), clines)
  clines <- gsubf("__NC__", format(ncol(init)), clines)
  clines <- gsubf("__NS__", format(ns), clines)
  clines <- gsubf("__WRAP__", boolstr(wrap), clines)
  clines <- gsubf("__CONTINUOUS_SCA__", boolstr(continuous), clines)
  clines <- gsubf("__DELTA_T__", ctrl[["delta_t"]], clines)
  clines <- gsubf("__USE_8_NB__", boolstr(use_8_nb), clines)
  clines <- gsubf("__SUBSTEPS__", format(ctrl[["substeps"]]), clines)
  clines <- gsubf("__XPOINTS__", format(ctrl[["xpoints"]]), clines)
  clines <- gsubf("__BETA_0_NROW__",  format(nrow(ctrl[["beta_0"]])),  clines)
  clines <- gsubf("__BETA_Q_NROW__",  format(nrow(ctrl[["beta_q"]])), clines)
  clines <- gsubf("__BETA_PP_NROW__", format(nrow(ctrl[["beta_pp"]])), clines)
  clines <- gsubf("__BETA_QQ_NROW__", format(nrow(ctrl[["beta_qq"]])), clines)
  clines <- gsubf("__BETA_PQ_NROW__", format(nrow(ctrl[["beta_pq"]])), clines)
  clines <- gsubf("__COMMON_HEADER__", 
                     system.file("common.h", package = "chouca"), clines)
  
  # Handle absorbing states. We get the states we go to but never get out from 
  if ( length(ctrl[["absorb_states"]]) > 0 ) { 
    clines <- gsubf("__HAS_ABSORB_STATES__", boolstr(TRUE), clines)
    clines <- gsubf("__N_ABSORB_STATES__", length(ctrl[["absorb_states"]]), clines)
    
    abs_states_array <- write_cpp_array_1d(ctrl[["absorb_states"]])
    clines <- gsubf("__ABSORB_STATES_ARRAY__", abs_states_array, clines)
    
  } else { 
    clines <- gsubf("__HAS_ABSORB_STATES__", boolstr(FALSE), clines)
  }
  
  # Write transition matrix 
  tmat_array <- write_cpp_array_2d(ctrl[["transition_mat"]])
  clines <- gsubf("__TMATRIX_ARRAY__", tmat_array, clines)
  
  # Decide whether we fixed the number of neighbors or not 
  clines <- gsubf("__FIXED_NEIGHBOR_NUMBER__", 
                     ifelse(fixed_neighborhood, "true", "false"), 
                     clines)
  
  # Set code lines that control the number of cores 
  cores <- ctrl[["cores"]]
  clines <- gsubf("__CORES__", format(cores), clines)
  clines <- gsubf("__USE_OMP__", ifelse(cores > 1, 1, 0), clines)
  
  # Set #define on whether to precompute transition probabilities 
  precompute_probas <- ctrl[["precompute_probas"]]
  if ( precompute_probas == "auto" ) { 
    precompute_probas <- ns^ifelse(use_8_nb, 8, 4) < prod(dim(init))
  }
  clines <- gsubf("__PRECOMP_PROBA_VALUE__", boolstr(precompute_probas), clines)
  
  # Make hash of file and replace function name 
  # We make the hash depend on the model too, just in case the user changes models, but
  # the rest is different. Unlikely, but who knows.
  hash <- digest::digest(list(clines, ctrl[["cores"]]), algo = "md5")
  clines <- gsub("__FPREFIX__", hash, clines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")
  
  # Make the table with all combinations of qs 
  if ( precompute_probas ) { 
    max_nb <- ifelse(use_8_nb, 8, 4)
    # filter = 1 when wrap = TRUE -> keep multiple of neighbors
    all_qs <- generate_all_qs(max_nb, ns, filter = wrap, line_cap = 0)
  } else { 
    # This is a dummy matrix just to make sure we pass something to the c++ function.
    all_qs <- matrix(0, nrow = 1, ncol = ns)
  }
  
  # Replace size in compiled code
  clines <- gsubf("__ALL_QS_NROW__", format(nrow(all_qs)), clines)
  
  # If we need to write the file somewhere, do it 
  if ( ! is.null(ctrl[["write_source"]]) ) { 
    writeLines(clines, ctrl[["write_source"]])
  }
  
  # Source cpp if needed 
  if ( ! exists(fname) || ctrl[["force_compilation"]] ) { 
    
    # We compile from the file, so that lines can be put in a debug run
    funs <- sourceCpp(code = paste(clines, collapse = "\n"), 
                      verbose = ctrl[["verbose_compilation"]], 
                      cleanupCacheDir = FALSE, 
                      rebuild = TRUE, # always force rebuild has we have our own cache
                      env = function_envir)
  }
  
  runf <- get(fname, envir = function_envir)
  runf(all_qs, ctrl)
}
  
}) # end of local() block


# Functions to write to c++ arrays
write_cpp_array_1d <- function(v) { 
  paste0("{", paste(v, collapse = ", "), "}")
}

# Functions to write to c++ arrays
write_cpp_array_2d <- function(m) { 
  m <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  write_cpp_array_1d(apply(m, 1, write_cpp_array_1d))
}

