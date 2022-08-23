# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(ctrl) { 
  
  # Split coefficient tables
  alpha <- ctrl[["alpha"]]
  pmat <- ctrl[["pmat"]]
  qmat <- ctrl[["qmat"]]
  pqmat <- ctrl[["pqmat"]]
  
  ctrl <- c(ctrl, 
            list(alpha_index = intmat(alpha[ ,c("from", "to"), drop = FALSE]), 
                 alpha_vals = as.numeric(alpha[ ,c("a0")]), # vector
                 pmat_index = intmat(pmat[ ,c("from", "to", "state"), drop = FALSE]), 
                 pmat_vals = pmat[ ,c("coef", "expo"), drop = FALSE], 
                 qmat_index = intmat(qmat[ ,c("from", "to", "state", "qs"), drop = FALSE]), 
                 qmat_vals = qmat[ ,"ys"],
                 pqmat_index = intmat(pqmat[ ,c("from", "to", "state"), drop = FALSE]), 
                 pqmat_vals = pqmat[ ,c("coef", "expo"), drop = FALSE]))
  
  camodel_cpp_engine(ctrl)
}


intmat <- function(m) { 
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}
