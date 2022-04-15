


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

// We tell gcc to unroll loops, as we have many small loops. This can double 
// performance (!!!)
__OLEVEL__
__OUNROLL__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// These strings will be replaced by their values 
constexpr arma::uword nr = __NR__; 
constexpr arma::uword nc = __NC__; 
constexpr char ns = __NS__; 
constexpr arma::uword ncoefs = __NCOEFS__; 
constexpr bool wrap = __WRAP__; 
constexpr bool use_8_nb = __USE_8_NB__; 
constexpr arma::uword substeps = __SUBSTEPS__; 
constexpr double ndbl = nr * nc; 

// Components of model 
constexpr bool has_X0 = __HAS_X0__; 
constexpr bool has_XP = __HAS_XP__; 
constexpr bool has_XQ = __HAS_XQ__; 
constexpr bool has_XPSQ = __HAS_XPSQ__; 
constexpr bool has_XQSQ = __HAS_XQSQ__; 

// Some things 
constexpr arma::uword max_nb = use_8_nb ? 8 : 4; 
// Fancy expression for max_nb^ns. The number of permutations in neighbors.
constexpr arma::uword tprob_size = (uword) exp( log( (double) (1+max_nb) ) * (double) ns ); 

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

inline void adjust_local_densities(signed char delta_qs[nr][nc][ns], 
                                   const uword i, 
                                   const uword j, 
                                   const char from, 
                                   const char to) { 
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    // left 
    delta_qs[i][(nc + j - 1) % nc][to]++; 
    delta_qs[i][(nc + j - 1) % nc][from]--; 
    
    // right
    delta_qs[i][(nc + j + 1) % nc][to]++; 
    delta_qs[i][(nc + j + 1) % nc][from]--; 
    
    // up
    delta_qs[(nr + i - 1) % nr][j][to]++; 
    delta_qs[(nr + i - 1) % nr][j][from]--; 
    
    // down
    delta_qs[(nr + i + 1) % nr][j][to]++;
    delta_qs[(nr + i + 1) % nr][j][from]--;
    
    if ( use_8_nb ) { 
      
      delta_qs[(nr + i - 1) % nr][(nc + j - 1) % nc][to]++; // upleft
      delta_qs[(nr + i - 1) % nr][(nc + j - 1) % nc][from]--; 
      
      delta_qs[(nr + i - 1) % nr][(nc + j + 1) % nc][to]++; // upright
      delta_qs[(nr + i - 1) % nr][(nc + j + 1) % nc][from]--; 
      
      delta_qs[(nr + i + 1) % nr][(nc + j - 1) % nc][to]++; // downleft
      delta_qs[(nr + i + 1) % nr][(nc + j - 1) % nc][from]--; 
      
      delta_qs[(nr + i + 1) % nr][(nc + j + 1) % nc][to]++; // downright
      delta_qs[(nr + i + 1) % nr][(nc + j + 1) % nc][from]--;
    }
    
  } else { 
    
    if ( i > 0 ) { 
      delta_qs[i-1][j][to]++; // up
      delta_qs[i-1][j][from]--;
    }
    if ( i < (nr-1) ) { 
      delta_qs[i+1][j][to]++; // down
      delta_qs[i+1][j][from]--; 
    }
    if ( j > 0 ) { 
      delta_qs[i][j-1][to]++; // left
      delta_qs[i][j-1][from]--; 
    }
    if ( j < (nc-1) ) { 
      delta_qs[i][j+1][to]++; // right
      delta_qs[i][j+1][from]--; 
    }
    
    if ( use_8_nb ) { 
      if ( i > 0 && j > 0 ) { 
        delta_qs[i-1][j-1][to]++; // upleft
        delta_qs[i-1][j-1][from]--; 
      }
      if ( i > 0 && j < (nc-1) ) { 
        delta_qs[i-1][j+1][to]++; // upright
        delta_qs[i-1][j+1][from]++; 
      }
      if ( i < (nr-1) && j > 0 ) { 
        delta_qs[i+1][j-1][from]++; // downleft
        delta_qs[i+1][j-1][to]--; 
      }
      if ( i < (nr-1) && j < (nc-1) ) { 
        delta_qs[i+1][j+1][to]++; // downright
        delta_qs[i+1][j+1][from]--; 
      }
    }
    
  }
  
}



void console_callback_wrap(const arma::uword iter, 
                           const arma::uword ps[ns], 
                           const double n, 
                           Rcpp::Function console_callback) { 
  uvec ps_arma(ns); // double
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
      m(i,j) = (ushort) omat[i][j]; 
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


// Compute transition probabilities between all possible qs states 
void compute_transition_probabilites(double tprobs[tprob_size][ns][ns], 
                                     const double ctrans[ns][ns][ncoefs], 
                                     const arma::uword ps[ns], 
                                     const char all_qs[tprob_size][ns+1]) { 
  
  for ( uword l=0; l<tprob_size; l++ ) { 
    // Get total of neighbors = last column of all_qs
    double total_qs = all_qs[l][ns]; 
    
    for ( char from=0; from<ns; from++ ) { 
      
      for ( char to=0; to<ns; to++ ) { 
        
        if ( has_X0 ) { 
          tprobs[l][from][to] = ctrans[to][from][0]; 
        } else { 
          tprobs[l][from][to] = 0.0; 
        }
        
        // Loop over coefs
        for ( char k=0; k<ns; k++ ) { 
          
          // XP
          if ( has_XP ) { 
            double coef = ctrans[to][from][1+k]; 
            tprobs[l][from][to] += coef * ( ps[k] / ndbl ); 
          }
          
          // XQ
          if ( has_XQ ) { 
            double coef = ctrans[to][from][1+k+ns]; 
            tprobs[l][from][to] += coef * ( all_qs[l][k] / total_qs ); 
          }
          
          // XPSQ
          if ( has_XPSQ ) { 
            double coef = ctrans[to][from][1+k+2*ns]; 
            tprobs[l][from][to] += coef * ( ps[k] / ndbl ) * ( ps[k] / ndbl ); 
          }
          
          // XQSQ
          if ( has_XQSQ ) { 
            double coef = ctrans[to][from][1+k+3*ns]; 
            tprobs[l][from][to] += coef * ( all_qs[l][k] / total_qs ) * 
                                     ( all_qs[l][k] / total_qs ); 
          }
        }
        
      }
      
      // Compute cumsum 
      for ( char to=1; to<ns; to++ ) { 
        tprobs[l][from][to] += tprobs[l][from][to-1];
      }
      
    }
  }
  
  
}

// 
// [[Rcpp::export]]
void aaa__FPREFIX__camodel_compiled_engine(const arma::cube trans, 
                                           const Rcpp::List ctrl, 
                                           const Rcpp::Function console_callback, 
                                           const Rcpp::Function cover_callback, 
                                           const Rcpp::Function snapshot_callback,
                                           const arma::Mat<ushort> all_qs_combinations) { 
  
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
  
  // Convert all_qs_combinations to a c array
  char all_qs[tprob_size][ns+1]; 
  for ( uword l=0; l<tprob_size; l++ ) { 
    for ( uword k=0; k<(ns+1); k++ ) { 
      all_qs[l][k] = all_qs_combinations(l, k); 
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
  
  // We divide all coefficients by substep to take that into account 
  for ( uword i=0; i<ns; i++ ) { 
    for ( uword j=0; j<ns; j++ ) { 
      for ( uword coef=0; coef<ncoefs; coef++ ) { 
        ctrans[i][j][coef] /= (double) substeps; 
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
  auto qs = new char[nr][nc][ns]; 
  get_local_densities(qs, omat); 
  
  // Array that will store changes in local densities
  auto delta_qs = new signed char[nr][nc][ns]; 
  
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
  double trans_probs[tprob_size][ns][ns]; 
  
  uword iter = 0; 
  
  while ( iter <= niter ) { 
    
    // Call callbacks 
    if ( console_callback_active && iter % console_callback_every == 0 ) { 
      console_callback_wrap(iter, ps, ndbl, console_callback); 
    }
    
    if ( cover_callback_active && iter % cover_callback_every == 0 ) { 
      cover_callback_wrap(iter, ps, ndbl, cover_callback); 
    }
    
    if ( snapshot_callback_active && iter % snapshot_callback_every == 0 ) { 
      snapshot_callback_wrap(iter, omat, snapshot_callback); 
    }
    
    for ( uword substep=0; substep < substeps; substep++ ) { 
      
      // Compute transition probabilities 
      compute_transition_probabilites(trans_probs, ctrans, ps, all_qs); 
      
      // Adjust changes in local densities to zero 
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          for ( char k=0; k<ns; k++ ) { 
            delta_qs[i][j][k] = 0; 
          }
        }
      }
      
      for ( uword i=0; i<nr; i++ ) { 
        
        for ( uword j=0; j<nc; j++ ) { 
          
          double rn = randunif();
          
          char cstate = omat[i][j]; 
          
          // Get total of neighbors
          double total_qs = number_of_neighbors(i, j);
          
          // Compute probabilities of transition
          // TODO: find a way to not compute the transition to self. We could consider 
          // trans = 0 the self-transition, then go from there?
          // TODO: consider pre-computing probability transitions, this solves the above 
          // problem
          // for ( char col=0; col<ns; col++ ) { 
          //   
          //   if ( has_X0 ) { 
          //     ptrans[col] = ctrans[col][cstate][0]; 
          //   } else { 
          //     ptrans[col] = 0.0; 
          //   }
          //   
          //   for ( char k=0; k<ns; k++ ) { 
          //     
          //     // XP
          //     if ( has_XP ) { 
          //       double coef = ctrans[col][cstate][1+k]; 
          //       ptrans[col] += coef * ( ps[k] / ndbl ); 
          //     }
          //     
          //     // XQ
          //     if ( has_XQ ) { 
          //       double coef = ctrans[col][cstate][1+k+ns]; 
          //       ptrans[col] += coef * ( qs[i][j][k] / total_qs ); 
          //     }
          //     
          //     // XPSQ
          //     if ( has_XPSQ ) { 
          //       double coef = ctrans[col][cstate][1+k+2*ns]; 
          //       ptrans[col] += coef * ( ps[k] / ndbl ) * ( ps[k] / ndbl ); 
          //     }
          //     
          //     // XQSQ
          //     if ( has_XQSQ ) { 
          //       double coef = ctrans[col][cstate][1+k+3*ns]; 
          //       ptrans[col] += coef * ( qs[i][j][k] / total_qs ) * 
          //                               ( qs[i][j][k] / total_qs ); 
          //     }
          //   }
          //   
          //   // Get line in precomputed proba for ptrans[col];  
          //   // uword line = 0; 
          //   // for ( char k = 0; k<ns; k++ ) { 
          //   //   line = line * (1+max_nb) + qs[i][j][k]; 
          //   // }
          //   // 
          //   // double ptrans2 = trans_probs[line][cstate][col]; 
          //   // Rcpp::Rcout << "i/j:" << i << "/" << j << ": "; 
          //   // for ( char k = 0; k<ns; k++ ) { 
          //   //   Rcpp::Rcout << (int) qs[i][j][k] << "(" << (int) all_qs[line][k] << "), " ; 
          //   // }
          //   // 
          //   // Rcpp::Rcout << "\n"; 
          //   // Rcpp::Rcout << (int) cstate << " to " << (int) col << " old: " << ptrans[col] << 
          //   //   " new: " << ptrans2 << " (line: " << line << ")\n"; 
          // }
          
          // Read from pre-computed probability
          uword line = 0; 
          for ( char k = 0; k<ns; k++ ) { 
            line = line * (1+max_nb) + qs[i][j][k]; 
          }
          
          // Check if we actually transition.  
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation since the random number is always below one. 
          char new_cell_state = cstate; 
          // for ( char k=1; k<ns; k++ ) { // cumsum
          //   ptrans[k] += ptrans[k-1]; 
          // }
          // for ( char k=0; k<ns; k++ ) { 
          //   Rcpp::Rcout << "c: " << trans_probs[line][cstate][k] << " o: " << 
          //     ptrans[k] << "\n"; 
          // }
          for ( char k=(ns-1); k>=0; k-- ) { 
            new_cell_state = rn < trans_probs[line][cstate][k] ? k : new_cell_state; 
          }
          
          if ( new_cell_state != cstate ) { 
            delta_ps[new_cell_state]++;
            delta_ps[cstate]--;
            nmat[i][j] = new_cell_state; 
            // Adjust delta_qs of neighbors
            adjust_local_densities(delta_qs, i, j, cstate, new_cell_state); 
          }
        }
      }
      
      // Apply changes in global densities 
      for ( char k=0; k<ns; k++ ) { 
        ps[k] += delta_ps[k]; 
        delta_ps[k] = 0; 
      }
      
      // Apply changes in local densities
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          for ( char k=0; k<ns; k++ ) { 
            qs[i][j][k] += delta_qs[i][j][k]; 
//             Rcpp::Rcout << "i: " << i << " j: " << j << " k: "  << (int) k << 
//               "qs: " << (int) qs[i][j][k] << 
//               " delta_qs: " << (int) delta_qs[i][j][k] << "\n"; 
          }
        }
      }
      
      // Apply changes in matrices (TODO: use memcpy)
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          omat[i][j] = nmat[i][j]; 
        }
      }
      
    } // end of substep loop
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] omat; 
  delete [] nmat; 
  delete [] qs; 
  delete [] delta_qs; 
  
}

  
