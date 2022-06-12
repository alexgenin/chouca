# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(alpha, pmat, qmat, control_list, 
                                    console_callback, cover_callback, snapshot_callback) { 
  
  
  # Split alpha 
  alpha_index <- intmat(alpha[ ,c("from", "to")])
  alpha_vals  <- as.numeric(alpha[ ,c("a0")]) # vector of length nstates
  
  # Split pmat 
  pmat_index <- intmat(pmat[ ,c("from", "to", "state")])
  pmat_vals  <- pmat[ ,c("coef", "expo")]
  
  # Split qmat 
  qmat_index <- intmat(qmat[ ,c("from", "to", "state", "qs")])
  qmat_vals  <- qmat[ ,"ys"] # vector
  
  # Divide values by number of substeps
  substep <- control_list[["substeps"]]
  pmat_vals[ ,"coef"] <- pmat_vals[ ,"coef"] / substep 
  qmat_vals <- qmat_vals / substep # all of it, it is a vector
  
  # Reduce pmat/qmat sizes to non-zeo coefficients
  non_zero_pmat <- which(pmat_vals[ ,"coef"] > 1e-8)
  pmat_vals <- pmat_vals[non_zero_pmat, , drop = FALSE]
  pmat_index <- pmat_index[non_zero_pmat, , drop = FALSE]
  
  non_zero_qmat <- which(qmat_vals > 1e-8)
  qmat_vals <- qmat_vals[non_zero_qmat]
  qmat_index <- qmat_index[non_zero_qmat, , drop = FALSE]
  
  camodel_cpp_engine(alpha_index, 
                     alpha_vals, 
                     pmat_index, 
                     pmat_vals, 
                     qmat_index, 
                     qmat_vals, 
                     control_list, console_callback, cover_callback, snapshot_callback)
  
}

intmat <- function(m) { 
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}
