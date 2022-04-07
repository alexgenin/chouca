


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

// We tell gcc to unroll loops, as we have many small loops. This can double 
// performance (!!!)
__OLEVEL__
__OUNROLL__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// These strings will be replaced by their values 
constexpr arma::uword nr = __NR__; 
constexpr arma::uword nc = __NC__; 
constexpr char ns = __NS__; 
constexpr arma::uword ncoefs = __NCOEFS__; 
constexpr bool wrap = __WRAP__; 
constexpr bool use_8_nb = __USE_8_NB__; 
constexpr arma::uword substeps = __SUBSTEPS__; 
constexpr double fsubsteps = __SUBSTEPS__; 
constexpr double n = nr * nc; 

using namespace arma;



/* This is xoshiro256+ 
 * https://prng.di.unimi.it/xoshiro256plus.c 
 */ 
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static uint64_t s[4];

static uint64_t nextr(void) {
	const uint64_t result = s[0] + s[3];
	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

static inline double randunif() { 
  uint64_t x = nextr(); 
  double xf = (x >> 11) * 0x1.0p-53; 
  return xf; 
}

inline double number_of_neighbors(const arma::uword i, 
                                  const arma::uword j) { 
  // When we do not wrap, the total number of neighbors depends on where we are 
  // in the matrix, e.g. first column or last row has missing neighbors. In that 
  // case, we need to correct for that. When wrapping is on, then the 
  // number of neighbors is constant.
  // 
  // We use a double because this number is then used in math using doubles.
  // 
  // The compiler will remove the if{} statements below at compile time if we are 
  // wrapping, making this function return a constant, defined on the line below:
  double nnb = use_8_nb ? 8.0 : 4.0; 
  
  // If we do not wrap and we use 8 neighbors, then we just need to substract from the 
  // total (maximum) counts. 
  if ( ! wrap && ! use_8_nb ) { 
    if ( i == 0 || i == (nr-1) ) { nnb--; }  // First or last row
    if ( j == 0 || j == (nc-1) ) { nnb--; }  // First or last column
  }
  
  // If we do not wrap and we use 8 neighbors, then it is a bit more subtle
  if ( ! wrap && use_8_nb ) { 
    if ( i == 0 || i == (nr-1) ) { nnb = nnb-3; }  // First or last row
    if ( j == 0 || j == (nc-1) ) { nnb = nnb-3; }  // First or last column
    
    // If we are in the corners, then one removed neighbor is being counted twice, so we
    // need to add it back
    if ( (i == 0 && j == 0) || 
         (i == 0 && j == (nc-1)) || 
         (i == (nr-1) && j == 0) || 
         (i == (nr-1) && j == (nc-1)) ) { 
      nnb++;
    }
    
  }
  
  return nnb; 
}

inline void get_local_counts(char qs[ns], 
                             const char m[nr][nc], 
                             arma::uword i, 
                             arma::uword j) { 
  
  // Set all counts to zero
  for ( char k=0; k<ns; k++ ) { 
    qs[k] = 0; 
  }
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    char state_left = m[i][(nc + j - 1) % nc];
    qs[state_left]++; // left 
    
    char state_right = m[i][(nc + j + 1) % nc];
    qs[state_right]++; // right
    
    char state_up = m[(nr + i - 1) % nr][j];
    qs[state_up]++; // up
    
    char state_down = m[(nr + i + 1) % nr][j];
    qs[state_down]++; // down
    
    if ( use_8_nb ) { 
      
      char state_upleft = m[(nr + i - 1) % nr][(nc + j - 1) % nc]; 
      qs[state_upleft]++; // upleft
      
      char state_upright = m[(nr + i - 1) % nr][(nc + j + 1) % nc]; 
      qs[state_upright]++; // upright
      
      char state_downleft = m[(nr + i + 1) % nr][(nc + j - 1) % nc]; 
      qs[state_downleft]++; // downleft
      
      char state_downright = m[(nr + i + 1) % nr][(nc + j + 1) % nc]; 
      qs[state_downright]++; // downright
    }
    
  } else { 
    
    if ( i > 0 ) { 
      char state_up = m[i-1][j];
      qs[state_up]++; // up
    }
    if ( i < (nr-1) ) { 
      char state_down = m[i+1][j]; 
      qs[state_down]++; // down
    }
    if ( j > 0 ) { 
      char state_left = m[i][j-1]; 
      qs[state_left]++; // left
    }
    if ( j < (nc-1) ) { 
      char state_right = m[i][j+1]; 
      qs[state_right]++; // right
    }
    
    if ( use_8_nb ) { 
      if ( i > 0 && j > 0 ) { 
        char state_upleft = m[i-1][j-1]; 
        qs[state_upleft]++; // upleft
      }
      if ( i > 0 && j < (nc-1) ) { 
        char state_upright = m[i-1][j+1]; 
        qs[state_upright]++; // upright
      }
      if ( i < (nr-1) && j > 0 ) { 
        char state_downleft = m[i+1][j-1]; 
        qs[state_downleft]++; // downleft
      }
      if ( i < (nr-1) && j < (nc-1) ) { 
        char state_downright = m[i+1][j+1]; 
        qs[state_downright]++; // downright
      }
    }
    
  }
  
}

void adjust_XQproba_at_ij(double XQproba[nr][nc][ns], 
                          const char m[nr][nc], 
                          const arma::uword i, 
                          const arma::uword j, 
                          const double ctrans[ns][ns][ncoefs]) { 
  
  // Compute local densities, etc. 
  char qs[ns]; 
  get_local_counts(qs, m, i, j); 
  double qs_total = number_of_neighbors(i, j); 
  
  // Make sum of local density components for transition to state 'to'
  char cstate = m[i][j]; 
//   Rcpp::Rcout << "current state at " << i << "/" << j << ": " << cstate << "\n"; 
  
  for ( char to=0; to<ns; to++ ) { 
    XQproba[i][j][to] = 0; 
    // Loop over coefficients to be added
    for ( char k=0; k<ns; k++ ) { 
//       Rcpp::Rcout << "neighbors of " << (int) k << ": " << (int) qs[k] << "\n"; 
//       Rcpp::Rcout << "coef of " << (int)k << ": " << ctrans[to][cstate][1+k+ns] << "\n"; 
      XQproba[i][j][to] += ctrans[to][cstate][1+k+ns] * 
                             ( (double) qs[k] / qs_total ) / fsubsteps; 
    }
//     Rcpp::Rcout << "proba to " << (int) to << ": " << XQproba[i][j][to] << "\n"; 
  }
  
}

void adjust_XQproba_around_ij(double XQproba[nr][nc][ns], 
                              const char m[nr][nc], 
                              const arma::uword i, 
                              const arma::uword j, 
                              const double ctrans[ns][ns][ncoefs]) { 
  
  // Adjust probability on the cell itself, which just changed state
  adjust_XQproba_at_ij(XQproba, m, i, j, ctrans); 
  
  if ( wrap ) { 
    uword col_left = (nc + j - 1) % nc; 
    adjust_XQproba_at_ij(XQproba, m, i, col_left, ctrans); 
    
    uword col_right = (nc + j + 1) % nc; 
    adjust_XQproba_at_ij(XQproba, m, i, col_right, ctrans); 
    
    uword row_up = (nr + i - 1) % nr; 
    adjust_XQproba_at_ij(XQproba, m, row_up, j, ctrans); 
    
    uword row_down = (nr + i + 1) % nr; 
    adjust_XQproba_at_ij(XQproba, m, row_down, j, ctrans); 
    
    if ( use_8_nb ) { 
      adjust_XQproba_at_ij(XQproba, m, row_up,   col_left,  ctrans); 
      adjust_XQproba_at_ij(XQproba, m, row_up,   col_right, ctrans); 
      adjust_XQproba_at_ij(XQproba, m, row_down, col_left,  ctrans); 
      adjust_XQproba_at_ij(XQproba, m, row_down, col_right, ctrans); 
    }
  
  } else { 
    Rcpp::Rcout << "not implemented\n";
  }
  
}

void console_callback_wrap(const arma::uword iter, 
                           const arma::uword ps[ns], 
                           const double n, 
                           Rcpp::Function console_callback) { 
  uvec ps_arma(ns); 
  for ( char k=0; k<ns; k++ ) { 
    ps_arma(k) = ps[k]; 
  }
  console_callback(iter, ps_arma, n); 
}

void snapshot_callback_wrap(const arma::uword iter, 
                            const char omat[nr][nc],
                            Rcpp::Function snapshot_callback) { 
  
  // Make arma array to give back to R
  Mat<ushort> m(nr, nc);
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      m(i,j) = omat[i][j]; 
    }
  }
  
  snapshot_callback(iter, m); 
}

void cover_callback_wrap(const arma::uword iter, 
                         const arma::uword ps[ns], 
                         const double n, 
                         Rcpp::Function cover_callback) { 
  uvec ps_arma(ns); 
  for ( char k=0; k<ns; k++ ) { 
    ps_arma(k) = ps[k]; 
  }
  cover_callback(iter, ps_arma, n); 
}


// 
// [[Rcpp::export]]
void aaa__FPREFIX__camodel_compiled_engine(const arma::cube trans, 
                                           const Rcpp::List ctrl, 
                                           const Rcpp::Function console_callback, 
                                           const Rcpp::Function cover_callback, 
                                           const Rcpp::Function snapshot_callback) { 
  
  // Unpack control list
  Mat<ushort> init = ctrl["init"]; // this is ushort because init is an arma mat
  uword niter      = ctrl["niter"]; // TODO: think about overflow in those values
  
  bool console_callback_active = ctrl["console_callback_active"]; 
  uword console_callback_every = ctrl["console_callback_every"]; 
  
  bool cover_callback_active = ctrl["cover_callback_active"]; 
  uword cover_callback_every = ctrl["cover_callback_every"]; 
  
  bool snapshot_callback_active = ctrl["snapshot_callback_active"]; 
  uword snapshot_callback_every = ctrl["snapshot_callback_every"]; 
  
  // Initialize some things as c arrays
  // Note: we allocate omat/nmat on the heap since they can be big matrices and blow up 
  // the size of the C stack beyond what is acceptable.
  auto omat = new char[nr][nc];
  auto nmat = new char[nr][nc];
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      omat[i][j] = (char) init(i, j);
      nmat[i][j] = (char) init(i, j);
    }
  }
  
  // Convert armadillo array to c-style array. Note that we change the order things are 
  // stored so that coefficients are contiguous in memory.
  double ctrans[ns][ns][ncoefs];
  for ( uword i=0; i<ns; i++ ) { 
    for ( uword j=0; j<ns; j++ ) { 
      for ( uword coef=0; coef<ncoefs; coef++ ) { 
        ctrans[i][j][coef] = (double) trans(coef, i, j);
      }
    }
  }
  
  // Initialize vector with counts of cells in each state in the landscape (used to 
  // compute global densities)
  uword ps[ns];
  int delta_ps[ns]; // Needs to be a signed integer
  memset(ps, 0, sizeof(ps));
  memset(delta_ps, 0, sizeof(ps));
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      ps[ omat[i][j] ]++; 
    }
  }
  
  // Initialize random number generator 
  // Initialize rng state using armadillo random integer function
  for ( uword i=0; i<4; i++ ) { 
    s[i] = randi<long>(); 
  }
  
  // The first number returned by the RNG is garbage, probably because randi returns 
  // something too short for s (64 bits long)
  // TODO: find a way to feed 64 bits integers to xoshiro
  randunif(); 
  
  // Initialize XQ component of probability 
  auto XQproba_old = new double [nr][nc][ns]; 
  auto XQproba_new = new double [nr][nc][ns]; 
  
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      adjust_XQproba_at_ij(XQproba_old, omat, i, j, ctrans); 
      adjust_XQproba_at_ij(XQproba_new, omat, i, j, ctrans); 
//       for ( char k=0; k<ns; k++ ) { 
//         Rcpp::Rcout << "XQproba_old: " << i << " " << j << " " << (int)k << ": " << 
//           XQproba_old[i][j][k] << "\n"; 
//       }
    }
  }
  
  // Allocate some things we will reuse later 
  double ptrans[ns];
  
  uword iter = 0; 
  
  while ( iter <= niter ) { 
    
    // Call callbacks 
    if ( console_callback_active && iter % console_callback_every == 0 ) { 
      console_callback_wrap(iter, ps, n, console_callback); 
    }
    
    if ( cover_callback_active && iter % cover_callback_every == 0 ) { 
      cover_callback_wrap(iter, ps, n, cover_callback); 
    }
    
    if ( snapshot_callback_active && iter % snapshot_callback_every == 0 ) { 
      snapshot_callback_wrap(iter, omat, snapshot_callback); 
    }
    
    for ( uword substep=0; substep < substeps; substep++ ) { 
      
      for ( uword i=0; i<nr; i++ ) { 
        
        for ( uword j=0; j<nc; j++ ) { 
          
          double rn = randunif();
          
          char cstate = omat[i][j]; 
          
          // Compute probabilities of transition
          for ( char col=0; col<ns; col++ ) { 
            ptrans[col] = ctrans[col][cstate][0] / fsubsteps; 
            // We already computed the sum of the c'(k) * q(k), so just add it here
            ptrans[col] += XQproba_old[i][j][col]; 
            // Add the coefficients linked to global densities c(k) * p(k). Here we need
            // the sum as there is no precomputed value. 
            for ( char k=0; k<ns; k++ ) { 
              double pcoef = ctrans[col][cstate][1+k]; 
              ptrans[col] += pcoef * ( ps[k] / n ) / fsubsteps; 
            }
//             Rcpp::Rcout << i << "/" << j << ": " << 
//               "to " << (int)col << ": " << ptrans[col] << "\n"; 
          }
          
          // Check if cell switches. 
          char new_cell_state = cstate; 
          double left = 0; 
          double right = 0; 
          for ( char k=0; k<ns; k++ ) { 
            right += ptrans[k]; 
            new_cell_state = ( left < rn && rn < right ) ? k : new_cell_state; 
            left = right; 
          }
          
          
          // TODO: Forget precomputing XQproba, too costly, but think about adjusting qs only. 
          // Copying is cheap.
          
          if ( new_cell_state != cstate ) { 
            delta_ps[new_cell_state]++;
            delta_ps[cstate]--;
            nmat[i][j] = new_cell_state; 
            
            // Adjust changes in local component of probabilities around the cell that
            // changes
            adjust_XQproba_around_ij(XQproba_new, nmat, i, j, ctrans); 
          }
        }
      }
      
      // Apply changes in global densities 
      for ( char k=0; k<ns; k++ ) { 
        ps[k] += delta_ps[k]; 
        delta_ps[k] = 0; 
      }
       
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          omat[i][j] = nmat[i][j];
        }
      }
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          for ( char k=0; k<ns; k++ ) { 
            XQproba_old[i][j][k] = XQproba_new[i][j][k]; // no memcpy???  
          }
        }
      }
      
    } // end of substep loop
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] omat; 
  delete [] nmat; 
  delete [] XQproba_old; 
  delete [] XQproba_new; 
  
}

  
