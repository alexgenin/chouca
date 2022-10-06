// 
// A few helper functions
// 

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

#include <RcppArmadillo.h>

#define _state_p 0
#define _state_q 1
#define _expo_p  2
#define _expo_q  3

//[[Rcpp::export]]
arma::vec quick_pred_cpp(const arma::vec coefs, 
                         const arma::mat ps, 
                         const arma::mat qs, 
                         const arma::mat vals) { 
  const arma::uword n = ps.n_rows; 
  const arma::uword nvals = vals.n_rows; 
  arma::vec yp(n); 
  
  for ( arma::uword i=0; i<n; i++ ) { 
    for ( arma::uword k=0; k<nvals; k++ ) { 
      // Note that we adjust the index 
      yp(i) += coefs(k) * 
        pow( ps(i, vals(k, _state_p) - 1), vals(k, _expo_p) ) * 
        pow( qs(i, vals(k, _state_q) - 1), vals(k, _expo_q) ); 
    }
  }
  
  return yp; 
}
