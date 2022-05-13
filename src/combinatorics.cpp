
#include <RcppArmadillo.h>
using namespace arma; 

// Return the line in the table of neighbor combinations to get the probability of 
// transition
// 
// 
arma::uword factorial(arma::uword n) { 
  return n > 0 ? n * factorial( n - 1 ) : 1;  
}

// [[Rcpp::export]]
arma::uword getline(const arma::uvec & qs, 
                    const arma::uword & nb, 
                    const arma::uword & ns) { 
  
  uword lines_before = 0; 
  uword remaining = nb; 
  
  // Init choose low 
  uword choose_low = ns - 1; 
  // Max prod_low = (ns-1-0)!
  uword prod_low = factorial(ns-1); 
  
  for ( uword k=1; k<ns; k++ ) { 
    
    // Compute choose(remaining - qi + choose_low, choose_low)
    // Adjust choose low and prod low
    choose_low--; 
    prod_low /= (choose_low + 1); 
    
    uword curqs = qs(k-1); 
    
    uword choose_high = remaining + ns - 1 - k + 1; 
    
    for ( uword qi=0; qi<curqs; qi++ ) { 
      choose_high--; //  = remaining - qi + ns - 1 - k; 
//       Rcpp::Rcout << "clow: " << choose_low << " chigh: " << choose_high << "\n"; 
//       if ( choose_low > choose_high ) return( 1 ) ; 
      uword prod_high = 1; 
      for ( uword l=remaining - qi + 1; l<=choose_high; l++ ) { 
        prod_high *= l;
      }
      // Add to the total of lines seen 
      // TODO: investigate possible recursion relationship?
      lines_before += prod_high / prod_low; 
    }
    
    remaining -= curqs; 
  }
  
  return( lines_before ); 
}

// Simple sum
// [[Rcpp::export]]
arma::uword simple_sum(const arma::uvec &qs, 
                       const arma::uword &nb, 
                       const arma::uword &ns) { 
  uword line = 0; 
  for ( uword k = 0; k<ns; k++ ) { 
    line = line * (1+nb) + qs(k); 
  }
  return line; 
}
