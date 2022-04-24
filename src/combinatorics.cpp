
#include <RcppArmadillo.h>
using namespace arma; 

// Return the line in the table of neighbor combinations to get the probability of 
// transition
// 
// [[Rcpp::export]]
arma::uword getline(arma::uvec qs, 
                    arma::uword nb, 
                    arma::uword ns) { 
  
  uword lines_before = 0; 
  uword remaining = nb; 
  
  for ( uword k=1; k<ns; k++ ) { 
    uword choose_low = ns - 1 - k; 
    // Compute choose(remaining - qi + choose_low, choose_low)
    uword prod_low = 1; 
    for ( uword l=1; l<=choose_low; l++ ) { 
      prod_low *= l; 
    }
    uword curqs = qs(k-1); 
    
    for ( uword qi=0; qi<curqs; qi++ ) { 
      uword choose_high = remaining - qi + ns - 1 - k; 
      uword prod_high = 1; 
      for ( uword l=choose_high; l>=choose_high - choose_low + 1; l-- ) { 
        prod_high *= l;
      }
      // Add to the total of lines seen 
      lines_before += prod_high / prod_low; 
    }
    
    remaining -= qs(k-1); 
  }
  
  return( lines_before ); 
}
