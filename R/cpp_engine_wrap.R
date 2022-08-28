# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(ctrl) { 
  
  # Split coefficient tables
  beta_0 <- ctrl[["beta_0"]]
  beta_p <- ctrl[["beta_p"]]
  beta_q <- ctrl[["beta_q"]]
  beta_pq <- ctrl[["beta_pq"]]
  
  ctrl <- c(ctrl, 
    list(beta_0_index  = intmat(beta_0[ ,c("from", "to"), drop = FALSE]), 
         beta_0_vals   = as.numeric(beta_0[ ,"ys"]), # vector
         beta_p_index  = intmat(beta_p[ ,c("from", "to", "state"), drop = FALSE]), 
         beta_p_vals   = beta_p[ ,c("coef", "expo"), drop = FALSE], 
         beta_q_index  = intmat(beta_q[ ,c("from", "to", "state", "qs"), drop = FALSE]), 
         beta_q_vals   = beta_q[ ,"ys"],
         beta_pq_index = intmat(beta_pq[ ,c("from", "to", "state"), drop = FALSE]), 
         beta_pq_vals  = beta_pq[ ,c("coef", "expo"), drop = FALSE])
    )
  
  camodel_cpp_engine(ctrl)
}


intmat <- function(m) { 
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}
