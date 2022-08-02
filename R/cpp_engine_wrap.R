# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(alpha, pmat, qmat, ctrl) { 
  
  
  # Split alpha 
  alpha_index <- intmat(alpha[ ,c("from", "to"), drop = FALSE])
  alpha_vals  <- as.numeric(alpha[ ,c("a0")]) # vector of length nstates
  
  # Split pmat 
  pmat_index <- intmat(pmat[ ,c("from", "to", "state"), drop = FALSE])
  pmat_vals  <- pmat[ ,c("coef", "expo"), drop = FALSE]
  
  # Split qmat 
  qmat_index <- intmat(qmat[ ,c("from", "to", "state", "qs"), drop = FALSE])
  qmat_vals  <- qmat[ ,"ys"] # vector
  
  camodel_cpp_engine(alpha_index, 
                     alpha_vals, 
                     pmat_index, 
                     pmat_vals, 
                     qmat_index, 
                     qmat_vals, 
                     ctrl)
  
}


intmat <- function(m) { 
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}
