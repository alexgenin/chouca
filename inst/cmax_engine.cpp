


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

// We tell gcc to unroll loops, as we have many small loops. This can double 
// performance (!!!)
#pragma GCC optimize("unroll-loops")

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

inline void get_local_densities(char qs[nr][nc][ns], 
                                const char m[nr][nc]) { 
  
  // Set all counts to zero
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      for ( char k=0; k<ns; k++ ) { 
        qs[i][j][k] = 0; 
      }
    }
  }
  
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      // Get neighbors to the left 
      if ( wrap ) { 
        
        char state_left = m[i][(nc + j - 1) % nc];
        qs[i][j][state_left]++; // left 
        
        char state_right = m[i][(nc + j + 1) % nc];
        qs[i][j][state_right]++; // right
        
        char state_up = m[(nr + i - 1) % nr][j];
        qs[i][j][state_up]++; // up
        
        char state_down = m[(nr + i + 1) % nr][j];
        qs[i][j][state_down]++; // down
        
        if ( use_8_nb ) { 
          
          char state_upleft = m[(nr + i - 1) % nr][(nc + j - 1) % nc]; 
          qs[i][j][state_upleft]++; // upleft
          
          char state_upright = m[(nr + i - 1) % nr][(nc + j + 1) % nc]; 
          qs[i][j][state_upright]++; // upright
          
          char state_downleft = m[(nr + i + 1) % nr][(nc + j - 1) % nc]; 
          qs[i][j][state_downleft]++; // downleft
          
          char state_downright = m[(nr + i + 1) % nr][(nc + j + 1) % nc]; 
          qs[i][j][state_downright]++; // downright
        }
        
      } else { 
        
        if ( i > 0 ) { 
          char state_up = m[i-1][j];
          qs[i][j][state_up]++; // up
        }
        if ( i < (nr-1) ) { 
          char state_down = m[i+1][j]; 
          qs[i][j][state_down]++; // down
        }
        if ( j > 0 ) { 
          char state_left = m[i][j-1]; 
          qs[i][j][state_left]++; // left
        }
        if ( j < (nc-1) ) { 
          char state_right = m[i][j+1]; 
          qs[i][j][state_right]++; // right
        }
        
        if ( use_8_nb ) { 
          if ( i > 0 && j > 0 ) { 
            char state_upleft = m[i-1][j-1]; 
            qs[i][j][state_upleft]++; // upleft
          }
          if ( i > 0 && j < (nc-1) ) { 
            char state_upright = m[i-1][j+1]; 
            qs[i][j][state_upright]++; // upright
          }
          if ( i < (nr-1) && j > 0 ) { 
            char state_downleft = m[i+1][j-1]; 
            qs[i][j][state_downleft]++; // downleft
          }
          if ( i < (nr-1) && j < (nc-1) ) { 
            char state_downright = m[i+1][j+1]; 
            qs[i][j][state_downright]++; // downright
          }
        }
        
      }
    }
  }
  
}

inline void adjust_local_densities(char qs[nr][nc][ns], 
                                   const uword i, 
                                   const uword j, 
                                   const char from, 
                                   const char to) { 
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    // left 
    qs[i][(nc + j - 1) % nc][to]++; 
    qs[i][(nc + j - 1) % nc][from]--; 
    
    // right
    qs[i][(nc + j + 1) % nc][to]++; 
    qs[i][(nc + j + 1) % nc][from]--; 
    
    // up
    qs[(nr + i - 1) % nr][j][to]++; 
    qs[(nr + i - 1) % nr][j][from]--; 
    
    // down
    qs[(nr + i + 1) % nr][j][to]++;
    qs[(nr + i + 1) % nr][j][from]--;
    
    if ( use_8_nb ) { 
      
      qs[(nr + i - 1) % nr][(nc + j - 1) % nc][to]++; // upleft
      qs[(nr + i - 1) % nr][(nc + j - 1) % nc][from]--; 
      
      qs[(nr + i - 1) % nr][(nc + j + 1) % nc][to]++; // upright
      qs[(nr + i - 1) % nr][(nc + j + 1) % nc][from]--; 
      
      qs[(nr + i + 1) % nr][(nc + j - 1) % nc][to]++; // downleft
      qs[(nr + i + 1) % nr][(nc + j - 1) % nc][from]--; 
      
      qs[(nr + i + 1) % nr][(nc + j + 1) % nc][to]++; // downright
      qs[(nr + i + 1) % nr][(nc + j + 1) % nc][from]--;
    }
    
  } else { 
    
    if ( i > 0 ) { 
      qs[i-1][j][to]++; // up
      qs[i-1][j][from]--;
    }
    if ( i < (nr-1) ) { 
      qs[i+1][j][to]++; // down
      qs[i+1][j][from]--; 
    }
    if ( j > 0 ) { 
      qs[i][j-1][to]++; // left
      qs[i][j-1][from]--; 
    }
    if ( j < (nc-1) ) { 
      qs[i][j+1][to]++; // right
      qs[i][j+1][from]--; 
    }
    
    if ( use_8_nb ) { 
      if ( i > 0 && j > 0 ) { 
        qs[i-1][j-1][to]++; // upleft
        qs[i-1][j-1][from]--; 
      }
      if ( i > 0 && j < (nc-1) ) { 
        qs[i-1][j+1][to]++; // upright
        qs[i-1][j+1][from]++; 
      }
      if ( i < (nr-1) && j > 0 ) { 
        qs[i+1][j-1][from]++; // downleft
        qs[i+1][j-1][to]--; 
      }
      if ( i < (nr-1) && j < (nc-1) ) { 
        qs[i+1][j+1][to]++; // downright
        qs[i+1][j+1][from]--; 
      }
    }
    
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
  
  // Compute local densities 
  char qs[nr][nc][ns]; 
  
  // Initialize random number generator 
  // Initialize rng state using armadillo random integer function
  for ( uword i=0; i<4; i++ ) { 
    s[i] = randi<long>(); 
  }
  
  // The first number returned by the RNG is garbage, probably because randi returns 
  // something too short for s (64 bits long)
  // TODO: find a way to feed 64 bits integers to xoshiro
  randunif(); 
  
  // Allocate some things we will reuse later 
  double ptrans[ns];
  double qs_norm = use_8_nb ? 8.0 : 4.0; 
  
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
    
    for ( uword substep = 0; substep < substeps; substep++ ) { 
      
      // Compute local densities
      get_local_densities(qs, omat); 
      
      for ( uword i=0; i<nr; i++ ) { 
        
        for ( uword j=0; j<nc; j++ ) { 
          
          double rn = randunif();
          
          char cstate = omat[i][j]; 
          
          // Get total of neighbors
          uword total_qs = 0;
          for ( char k=0; k<ns; k++ ) { 
            total_qs += qs[i][j][k]; 
          }
          
          // Compute probabilities of transition
          for ( char col=0; col<ns; col++ ) { 
            
            ptrans[col] = (double) ctrans[col][cstate][0] / substeps; 
            for ( char k=0; k<ns; k++ ) { 
              // 40% time spent here
              // c(k) * p(k)
              double pcoef = (double) ctrans[col][cstate][1+k]; 
              ptrans[col] += pcoef * ( ps[k] / (double) n ) / substeps; 
              // c'(k) * q(k)
              double qcoef = (double) ctrans[col][cstate][1+k+ns]; 
              ptrans[col] += qcoef * ( qs[i][j][k] /  (double) total_qs ) / substeps; 
            }
            
          }
          
          // Compute cumsum
          for ( char k=1; k<ns; k++ ) { 
            ptrans[k] = ptrans[k-1] + ptrans[k]; 
          }
          
          // Check if cell switches 
          char new_cell_state = 0; 
          bool cell_switch = false; 
          if ( rn < ptrans[0] ) { 
            new_cell_state = 0;
            cell_switch = true; 
          } else { 
            for ( char newstate=1; newstate<ns; newstate++ ) { 
              // Very bad branch prediction, consider predicating (6% time spent)
              if ( rn > ptrans[newstate-1] && rn < ptrans[newstate] ) { 
                new_cell_state = newstate;
                cell_switch = true; 
              }
            }
          }
          
          if ( cell_switch ) {
            delta_ps[new_cell_state]++;
            delta_ps[cstate]--;
            nmat[i][j] = new_cell_state; 
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
      
      
    }
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] omat; 
  delete [] nmat; 
  
}

  
