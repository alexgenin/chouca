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

//[[Rcpp::export]]
arma::vec get_cell_qvec(const arma::Mat<unsigned short>& m,
                        const arma::uword& i_R,
                        const arma::uword& j_R,
                        const arma::uword ns, 
                        const bool& wrap,
                        const arma::Mat<unsigned short>& kernel) {
  
  arma::vec qs(ns); 
  qs.fill(0);
  arma::uword nc = m.n_cols;
  arma::uword nr = m.n_rows;
  
  arma::uword i = i_R - 1; 
  arma::uword j = j_R - 1; 
  
  const arma::sword kernel_semiheight = ( kernel.n_rows - 1 ) / 2;
  const arma::sword kernel_semiwidth  = ( kernel.n_cols - 1 ) / 2;
  
  // We loop over the required offsets
  for ( arma::sword o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) { 
    for ( arma::sword o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) { 
      // Rcpp::Rcout << "o_r: " << o_r << " o_c: " << o_c << "\n"; 
      
      if ( wrap ) { 
        const ushort state = m( (nr + i + o_r) % nr, 
                                (nc + j + o_c) % nc ); 
        qs(state) += kernel(o_r + kernel_semiheight,
                               o_c + kernel_semiwidth); 
        
      } else { 
        // If we don't wrap, then we need to add a bound check to make sure 
        // we are in the matrix
        const arma::uword i_target = i + o_r; 
        const arma::uword j_target = j + o_c; 
        if ( i_target >= 0 && j_target >= 0 && i_target < nr && j_target < nc ) { 
          const ushort state = m(i_target, j_target); 
          qs(state) += kernel(o_r + kernel_semiheight, 
                              o_c + kernel_semiwidth); 
        }
      }
    }
  }
  
  return qs; 
}
