# 
# Interface to the compiled engine
# 

camodel_compiled_engine <- function(trans, ctrl, 
                                    console_callback, cover_callback, snapshot_callback) { 
  
  # Unwrap elements of the ctrl list 
  substeps <- ctrl[["substeps"]]
  wrap     <- ctrl[["wrap"]]
  use_8_nb <- ctrl[["use_8_neighbors"]] 
  init     <- ctrl[["init"]]
  niter    <- ctrl[["niter"]]
  ns       <- ctrl[["nstates"]]
  coefs    <- nrow(trans) 
  
  # Read file
  cmaxfile <- system.file("cmax_engine.cpp", package = "chouca")
  cmaxlines <- readLines(cmaxfile) 
  
  # Replace template values
  cmaxlines <- gsub("__NR__", format(nrow(init)), cmaxlines)
  cmaxlines <- gsub("__NC__", format(ncol(init)), cmaxlines)
  cmaxlines <- gsub("__NS__", format(ns), cmaxlines)
  cmaxlines <- gsub("__NCOEFS__", format(coefs), cmaxlines)
  cmaxlines <- gsub("__WRAP__", ifelse(wrap, "true", "false"), cmaxlines)
  cmaxlines <- gsub("__USE_8_NB__", ifelse(use_8_nb, "true", "false"), cmaxlines)
  cmaxlines <- gsub("__SUBSTEPS__", format(substeps), cmaxlines)
  
  # Make hash of file and replace function name 
  # We make the hash depend on the model too, just in case the user changes models, but
  # the rest is different. Unlikely, but who knows.
  hash <- digest::digest(list(cmaxlines, trans), algo = "md5")
  cmaxlines <- gsub("__FPREFIX__", hash, cmaxlines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")
  
  # Source cpp if needed 
  if ( ! exists(fname) ) { 
    # We compile from the file, so that lines can be put in a debug run
    funs <- sourceCpp(code = paste(cmaxlines, collapse = "\n"), 
                      verbose = TRUE, cacheDir = ".")
  }
  
  runf <- get(fname)
  runf(trans, ctrl, console_callback, cover_callback, snapshot_callback)
}

