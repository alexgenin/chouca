# 
# Interface to the compiled engine
# 

camodel_compiled_engine <- function(trans, ctrl, 
                                    console_callback, cover_callback, snapshot_callback) { 
  
  # Unwrap elements of the ctrl list 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["neighbors"]] == 8
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  coefs    <- nrow(trans) 
  
  # Read file
  cmaxfile <- system.file("compiled_engine_precomputed.cpp", package = "chouca")
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
  cmaxlines <- gsubf("__NCOEFS__", format(coefs), cmaxlines)
  cmaxlines <- gsubf("__WRAP__", ifelse(wrap, "true", "false"), cmaxlines)
  cmaxlines <- gsubf("__USE_8_NB__", ifelse(use_8_nb, "true", "false"), cmaxlines)
  cmaxlines <- gsubf("__SUBSTEPS__", format(substeps), cmaxlines)
  cmaxlines <- gsubf("__COMMON_HEADER__", 
                     system.file("common.h", package = "chouca"), cmaxlines)
  
  # Probability components: turn on or off in compiled code as needed
  totX0 <- sum(trans[1, , ]) # X0
  totXP <- sum(trans[2:(2+ns-1), , ]) # XP
  totXQ <- sum(trans[(2+ns):(2+2*ns-1), , ]) # XQ
  totXPSQ <- sum(trans[(2+2*ns):(2+3*ns-1), , ]) # XPSQ
  totXQSQ <- sum(trans[(2+3*ns):(2+4*ns-1), , ]) # XQSQ
  
  boolstr <- function(x) ifelse(x>0, "true", "false")
  cmaxlines <- gsubf("__HAS_X0__", boolstr(totX0), cmaxlines)
  cmaxlines <- gsubf("__HAS_XP__", boolstr(totXP), cmaxlines)
  cmaxlines <- gsubf("__HAS_XQ__", boolstr(totXQ), cmaxlines)
  cmaxlines <- gsubf("__HAS_XPSQ__", boolstr(totXPSQ), cmaxlines)
  cmaxlines <- gsubf("__HAS_XQSQ__", boolstr(totXQSQ), cmaxlines)
  
  # Make the table with all combinations of qs 
  max_nb <- ifelse(use_8_nb, 8, 4)
  all_qs <- rep( list(seq(0, max_nb) ), each = ns)
  all_qs <- as.matrix(do.call(expand.grid, all_qs)) # !!! very large matrix
  all_qs <- all_qs[ ,seq(ncol(all_qs), 1)]
  colnames(all_qs) <- rownames(all_qs) <- NULL
  
  # Precompute neighbors sum and add it as the last column
  all_qs <- cbind(all_qs, rowSums(all_qs))
  
  # Set the number of lines of all_qs in c file
  cmaxlines <- gsubf("__TPROB_SIZE__", nrow(all_qs), cmaxlines)
  
  
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
  
  # Source cpp if needed 
  if ( ! exists(fname) ) { 
      
    # We compile from the file, so that lines can be put in a debug run
    funs <- sourceCpp(code = paste(cmaxlines, collapse = "\n"), 
                      verbose = TRUE, cacheDir = ".", cleanupCacheDir = TRUE)
  }
  
  runf <- get(fname)
  runf(trans, ctrl, console_callback, cover_callback, snapshot_callback, 
       all_qs)
}



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

