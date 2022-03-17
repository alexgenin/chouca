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
  
  # Compiling options 
  write_to_file <- ctrl[["code_file"]]
  
  # Read file, compile, and write
  cmaxfile <- "./inst/cmax_engine.cpp"
  cmaxlines <- readLines(cmaxfile) 
  cmaxlines <- gsub("__NR__", format(nrow(init)), cmaxlines)
  cmaxlines <- gsub("__NC__", format(ncol(init)), cmaxlines)
  cmaxlines <- gsub("__NS__", format(ns), cmaxlines)
  cmaxlines <- gsub("__NCOEFS__", format(coefs), cmaxlines)
  cmaxlines <- gsub("__WRAP__", ifelse(wrap, "true", "false"), cmaxlines)
  cmaxlines <- gsub("__USE_8_NB__", ifelse(use_8_nb, "true", "false"), cmaxlines)
  cmaxlines <- gsub("__SUBSTEPS__", format(substeps), cmaxlines)
  
  # Make hash of file and replace function name 
  hash <- digest::digest(cmaxlines, algo = "md5")
  cmaxlines <- gsub("__FPREFIX__", hash, cmaxlines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")
  
  
  # Source cpp 
  if ( ! is.null(write_to_file) ) { 
    # Write code to temporary file 
    tfile <- paste0(write_to_file, ".cpp")
    writeLines(cmaxlines, tfile)
    # We compile from the file, so that lines can be put in a debug run
    funs <- sourceCpp(tfile, verbose = TRUE)
  } else { 
    # We compile directly from the read code, letting Rcpp choose its temporary file 
    funs <- sourceCpp(code = paste(cmaxlines, collapse = "\n"), verbose = TRUE)
  }
  
  runf <- get(funs[["functions"]][1])
  runf(trans, ctrl, console_callback, cover_callback, snapshot_callback)
}
