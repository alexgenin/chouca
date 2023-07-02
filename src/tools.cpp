// 
// A few helper functions
// 

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

#include <RcppArmadillo.h>
#include "../inst/include/chouca.h"

#define _state_p 0
#define _state_q 1
#define _expo_p  2
#define _expo_q  3

static inline double fintpow(double a, arma::uword k) { 
  switch ( k ) { 
    case 0: return 1.0; 
    case 1: return a; 
    case 2: return a*a; 
    case 3: return a*a*a; 
    case 4: return a*a*a*a; 
    case 5: return a*a*a*a*a; 
    case 6: return a*a*a*a*a*a; 
    case 7: return a*a*a*a*a*a*a; 
  }
  // Should never happen
  return 1.0; 
}

//[[Rcpp::export]]
arma::vec quick_pred_cpp(const arma::vec coefs, 
                         const arma::mat ps, 
                         const arma::mat qs, 
                         const arma::mat vals) { 
  const arma::uword n = ps.n_rows; 
  const arma::uword nvals = vals.n_rows; 
  arma::vec yp(n); 
  yp.fill(0.0); 
  
  for ( arma::uword i=0; i<n; i++ ) { 
    for ( arma::uword k=0; k<nvals; k++ ) { 
      // Note that we adjust the index 
      yp(i) += coefs(k) * 
        fintpow( ps(i, vals(k, _state_p) - 1), vals(k, _expo_p) ) * 
        fintpow( qs(i, vals(k, _state_q) - 1), vals(k, _expo_q) ); 
    }
  }
  
  return yp; 
}
