
#include <RcppArmadillo.h>
#include "../inst/include/chouca.h"

using namespace arma;


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif


inline unsigned long long intpow(arma::uword a,
                                 arma::uword b) {
  unsigned long long p = 1;
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
                           arma::uword filter,
                           arma::uword line_cap) {

  // keep only multiple of neighbor number (0, 8, 16, 32, etc.)
  constexpr uword FILTER_MULTIPLE_NEIGHB = 1;
  // keep only when neighbor number is obtained (4 or 8)
  constexpr uword FILTER_OK_NEIGHBOR_SUM = 2;

  // Max number of lines to store by default. This can be a big number
  unsigned long long out_nrows = intpow(nb+1, ns);

  out_nrows = filter == FILTER_MULTIPLE_NEIGHB ? out_nrows / nb : out_nrows;
  // https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)
  out_nrows = filter == FILTER_OK_NEIGHBOR_SUM ?
                (uword) Rf_choose(nb + ns - 1, ns - 1) :
                out_nrows;

  // Handle the line cap
  uword save_every = 1;
  if ( line_cap > 0 && out_nrows > line_cap ) {
    uword save_every = out_nrows / line_cap;
    out_nrows /= save_every;
  }

  umat all_qs(out_nrows, ns+1);

  // This follows out_nrows and can be very long
  unsigned long long line=0;

  uvec this_qs(ns);
  uword total = accu(this_qs);
  uword this_save = 0;

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
      } else {
        break;
      }
    }

    // If we every 4/8 neighbor values, save the value. Or if we do not filter, then
    // keep everything.
    bool store = true;
    if ( filter == FILTER_MULTIPLE_NEIGHB ) {
      store = total % nb == 0;
    } else if ( filter == FILTER_OK_NEIGHBOR_SUM ) {
      store = total == nb;
    }

    if ( store && ( this_save % save_every == 0 ) ) {

      // Store into matrix. Do not use submat which does not seem to handle
      // very big matrices
      for ( arma::uword st=0; st<ns; st++ ) {
        all_qs(line, st) = this_qs(st);
      }
      all_qs(line, ns) = total;

      line++;
      this_save++;
    }

    R_CheckUserInterrupt();
  }

  return all_qs;
}

