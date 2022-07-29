
#include <RcppArmadillo.h>
using namespace arma; 

// Return the line in the table of neighbor combinations to get the probability of 
// transition
// TODO: expand the rationale behind this
// 
// When working with pre-computed probabilities, we do not need to store the local 
// densities: we just need to know at which line to go grab the probability we are 
// interested in. 
// 
// When a cell changes states, we need to update the line. However, computing the line 
// itself is quite expensive. Instead we precompute how the line needs to be adjusted 
// in a matrix, so on state change we can just go get the adjustment value. There is 
// probability a mathematical formula to compute the change in lines, but this it is 
// too hard for my small ecologist brain.
// 
// 
// [[Rcpp::export]]
arma::uword getline(arma::urowvec qs, 
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




inline uword intpow(arma::uword a, 
                    arma::uword b) { 
  uword p = 1; 
  for ( uword k=0; k<b; k++ ) { 
    p *= a; 
  }
  return p; 
}

// 
// Generate the matrix all_qs. This function allows creating it directly instead of 
// using expand.grid(), then filtering only values that sum to the desired number of 
// neighbors.
// The all_qs matrix contains all the valid combinations of neighbors. 
// 
//[[Rcpp::export]]
arma::umat generate_all_qs(arma::uword nb, 
                           arma::uword ns, 
                           bool filter) { 
  
  // The last qs vector has all of its neighbors in the highest state 
  urowvec last_qs(ns); 
  last_qs(0) = nb; 
  uword out_nrows = intpow(nb+1, ns) / ( filter ? nb : 1 ); 
  
  umat all_qs(out_nrows, ns); 
  
  urowvec this_qs(ns); 
  uword line=0; 
  uword total = accu(this_qs); 
  
  while ( line < out_nrows ) { 
    
    // Update this_qs 
    this_qs(ns-1)++; 
    total++; // we added one to the total of neighbors
    for ( sword k=(ns-1); k>0; k-- ) { 
      // Carry over what is above nb to the next state
      if ( this_qs(k) > nb ) { 
        this_qs(k) = 0; 
        this_qs(k-1)++; 
        // we removed (nb+1) to the total (the value of this_qs(k) before setting it to
        // zero), and added one 
        total += 1 - (nb + 1); 
      }
    }
    
    // If we sum to nb, save the value 
    if ( filter && total % nb == 0 ) { 
      all_qs.row(line) = this_qs; 
      line++; 
    }
    
    R_CheckUserInterrupt();
  }
  
  return all_qs; 
}

