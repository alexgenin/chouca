# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(ctrl) { 
  
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
  
  camodel_cpp_engine(ctrl)
}


intmat <- function(m) { 
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}

