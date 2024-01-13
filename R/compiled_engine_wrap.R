#
# Interface to the compiled engine
#


# Here we use the local() call so that the function can keep some internal state.
# Otherwise the parent environment of the function is the package's internal
# environment, which cannot be written into.
camodel_compiled_engine_wrap <- local({

# This environment is here to store the functions that will be created here
function_envir <- environment()

function(ctrl, console_callback, cover_callback, snapshot_callback) {

  # Unwrap elements of the ctrl list that we need here
  wrap       <- ctrl[["wrap"]]
  use_8_nb   <- ctrl[["neighbors"]] == 8
  fixed_neighborhood <- ctrl[["fixed_neighborhood"]]
  init       <- ctrl[["init"]]
  ns         <- ctrl[["nstates"]]
  max_nb <- ifelse(use_8_nb, 8, 4)

  # Read cpp file that we will compile
  cppfile <- system.file("compiled_engine.cpp", package = "chouca")
  clines <- readLines(cppfile)

  # gsubf() replaces lines in the read cpp file, and print to the console if that
  # was asked for
  if ( ctrl[["verbose_compilation"]] ) {
    cat("Compilation options:\n")
    gsubf <- function(a, b, lines) { #
      cat(sprintf("Setting %s to '%s'\n", a, b))
      lines <- gsub(a, b, lines)
    }
  } else {
    gsubf <- gsub
  }
  # Convert a logical to a cpp bool string
  boolstr <- function(x) ifelse(x, "true", "false")

  # Replace values in cpp template
  clines <- gsubf("__NR__", format(nrow(init)), clines)
  clines <- gsubf("__NC__", format(ncol(init)), clines)
  clines <- gsubf("__NS__", format(ns), clines)
  clines <- gsubf("__WRAP__", boolstr(wrap), clines)
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


  # Write transition matrix. This contains 1/0 depending on whether the transition
  # from one state to another exists.
  tmat_array <- write_cpp_array_2d(ctrl[["transition_mat"]])
  clines <- gsubf("__TMATRIX_ARRAY__", tmat_array, clines)

  # Write indices over which to iterate for each transition. These small matrices
  # contains where to jump in the big table for each transition

  # We pack all coefficients needed to compute all transition probabilities into a big
  # table. Recall that all transition probabilities in chouca have the form:
  #  beta0 +
  #    sum( f(q) ) +
  #    sum( c p_a^b p_c^d ) for various a,b,c,d
  #    sum( c q_a^b q_c^d ) for various a,b,c,d
  #    sum( c p_a^b q_c^d ) for various a,b,c,d
  #
  # All coefficients (beta0, a, b, c, d, etc.) are stored in different tables at this
  # stage, depending on which component above they correspond to:
  #    - 'c' (beta_0),
  #    - 'f(q)' (beta_q),
  #    - 'c*p[x]*p[y]' (beta_pp),
  #    - 'c*q[x]*q[y]' (beta_qq)
  #    - 'c*p[x]*q[y]' (beta_pq)
  #
  # Here we merge all that into a single big table. This seems to help with cache
  # locality when computing transition probabilities
  tables <- c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq")
  coef_table <- plyr::ldply(tables, function(tab) {
    if ( nrow(ctrl[[tab]]) > 0 ) {
      data.frame(table = tab, ctrl[[tab]])
    } else {
      data.frame(table = character(0), ctrl[[tab]])
    }
  })

  # 'betas_index' contains the starting and ending index to pick relevant coefficients in
  # coef_table for each transition. It is a 3d array where the first columns correspond
  # to the state being transitioned from and to, and the third dimension contains the
  # starting and ending indices in coef_table. 
  # This is a bit cumbersome, but it allows having all the different tables beta_0, 
  # beta_q, etc. in one big table, which helps with cache locality
  # 
  betas_index <- array(nrow(coef_table)+1, dim = list(ns, ns, 2 * length(tables)))
  for ( i in seq_along(tables) ) {
    for ( from in seq(0, ns-1) ) {
      for ( to in seq(0, ns-1) ) {
        this_tab <- tables[i]
        which_rows <- which(this_tab == coef_table[ ,"table"] &
                            coef_table[ ,"from"] == from &
                            coef_table[ ,"to"] == to)

        if ( length(which_rows) > 0 ) {
          # Note the +/- 1 to take into account c++/R indexing differences
          betas_index[from+1, to+1, 1 + 2*(i-1) ] <- min(which_rows) - 1
          betas_index[from+1, to+1, 1 + 2*(i-1)+1 ] <- max(which_rows) - 1
        } else {
          betas_index[from+1, to+1, 1 + 2*(i-1) ] <- -1
          betas_index[from+1, to+1, 1 + 2*(i-1)+1 ] <- -2
        }
      }
    }
  }

  # Trim coef table and write it to c++
  coef_tab_ints <- as.matrix(coef_table[ ,c("from", "to", "state_1", "state_2",
                                            "expo_1", "expo_2", "qs")])
  coef_tab_ints[is.na(coef_tab_ints)] <- 99
  coef_tab_dbls <- as.matrix(coef_table[ ,"coef", drop = FALSE])
  ctrl[["coef_tab_dbls"]] <- coef_tab_dbls
  ctrl[["coef_tab_ints"]] <- coef_tab_ints
  clines <- gsubf("__COEF_TAB_NROW__", nrow(coef_tab_ints), clines)

  # Write it as c++ array in the file
  betas_index_str <- write_cpp_array_3d(betas_index)
  clines <- gsubf("__BETAS_INDEX__", betas_index_str, clines)

  # Decide whether we fixed the number of neighbors or not
  clines <- gsubf("__FIXED_NEIGHBOR_NUMBER__",
                  ifelse(fixed_neighborhood, "true", "false"),
                  clines)

  # Set code lines that control the number of cores
  cores <- ctrl[["cores"]]
  clines <- gsubf("__CORES__", format(cores), clines)
  clines <- gsubf("__USE_OMP__", ifelse(cores > 1, 1, 0), clines)

  # Set #define on whether to precompute transition probabilities or not
  precompute_probas <- ctrl[["precompute_probas"]]

  if ( precompute_probas == "auto" ) {
    precompute_probas <- ns^ifelse(use_8_nb, 8, 4) < prod(dim(init)) &
                           ( max_nb^ns < 1e8 )
  }
  if ( precompute_probas && ( max_nb^ns > 1e8 ) ) {
    warning("Model too complex to precompute probabilities, not doing it.")
    precompute_probas <- FALSE
  }

  clines <- gsubf("__PRECOMP_PROBA_VALUE__", boolstr(precompute_probas), clines)

  # Make hash of file and replace function name
  # We make the hash depend on the model too, just in case the user changes models, but
  # the rest is different. Unlikely, but who knows.
  hash <- digest::digest(list(clines, ctrl[["cores"]]), algo = "md5")
  clines <- gsub("__FPREFIX__", hash, clines)
  fname <- paste0("aaa", hash, "camodel_compiled_engine")

  # Make the table with all combinations of qs. If we wrap, then we can discard
  # combinations that are not multiples of the number of neighbors (this is what the
  # 'filter' argument does below)
  if ( precompute_probas ) {
    all_qs <- generate_all_qs(max_nb, ns, 
                              filter = wrap, 
                              line_cap = 0)
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


# Functions to write a c++ arrays from an R matrix, 1d, 2d and 3d
write_cpp_array_1d <- function(v) {
  paste0("{", paste(v, collapse = ", "), "}")
}
write_cpp_array_2d <- function(m) {
  m <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  write_cpp_array_1d(apply(m, 1, write_cpp_array_1d))
}
write_cpp_array_3d <- function(m) {
  write_cpp_array_1d(apply(m, 3, write_cpp_array_2d))
}

