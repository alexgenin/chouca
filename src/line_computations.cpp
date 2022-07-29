
#include <RcppArmadillo.h>
using namespace arma; 



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
  
  umat all_qs(out_nrows, ns+1); 
  
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
    
    // If we every 4/8 neighbor values, save the value. Or if we do not filter, then 
    // keep everything. 
    if ( (! filter) || total % nb == 0 ) { 
      all_qs.submat(line, 0, line, ns-1) = this_qs; 
      all_qs(line, ns) = total; 
      line++; 
    }
    
    R_CheckUserInterrupt();
  }
  
  return all_qs; 
}

