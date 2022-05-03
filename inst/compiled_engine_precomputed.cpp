


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

// Define an unsigned char as uchar for more legibility 
typedef unsigned char uchar; 

// These strings will be replaced by their values 
constexpr arma::uword nr = __NR__; 
constexpr arma::uword nc = __NC__; 
constexpr uchar ns = __NS__; 
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

// Include functions and type declarations
#include "__COMMON_HEADER__"

// Some things 
constexpr arma::uword max_nb = use_8_nb ? 8 : 4; 
// The number of permutations in neighbors.
constexpr arma::uword tprob_size = __TPROB_SIZE__;
constexpr uchar tprob_interval = wrap ? ( use_8_nb ? 8 : 4 ) : 1; 

// Compute transition probabilities between all possible qs states 
void precompute_transition_probabilites(double tprobs[tprob_size][ns][ns], 
                                        const double ctrans[ns][ns][ncoefs], 
                                        const arma::uword ps[ns], 
                                        const uchar all_qs[tprob_size][ns+1]) { 
  
  // Note tprob_interval here. In all combinations of neighbors, only some of them 
  // can be observed in the wild. If we wraparound, then the number of neighbors is 
  // constant, it is either 8 or 4. So we can compute the values in tprobs only at 
  // indices every 4 or 8 values. 
  
  for ( uword l=0; l<tprob_size; l+=tprob_interval ) { 
    // Get total of neighbors = last column of all_qs
    double total_qs = all_qs[l][ns]; 
    
    for ( uchar from=0; from<ns; from++ ) { 
      
      for ( uchar to=0; to<ns; to++ ) { 
        
        // Useless
        // if ( from == to ) { 
        //   continue; 
        // }
        
        if ( has_X0 ) { 
          tprobs[l][from][to] = ctrans[from][to][0]; 
        } else { 
          tprobs[l][from][to] = 0.0; 
        }
        
        // Loop over coefs
        for ( uchar k=0; k<ns; k++ ) { 
          
          // All the if{} blocks below will be removed appropriately by the compiler.
          
          // XP
          if ( has_XP ) { 
            double coef = ctrans[from][to][1+k]; 
            tprobs[l][from][to] += coef * ( ps[k] / ndbl ); 
          }
          
          // XQ
          // TODO: we don't have to divide by total_qs here, we can do it beforehand, 
          // but this is not the hottest loop in town when precomputing transitions
          if ( has_XQ ) { 
            double coef = ctrans[from][to][1+k+ns]; 
            tprobs[l][from][to] += coef * ( all_qs[l][k] / total_qs ); 
          }
          
          // XPSQ
          if ( has_XPSQ ) { 
            double coef = ctrans[from][to][1+k+2*ns]; 
            tprobs[l][from][to] += coef * ( ps[k] / ndbl ) * ( ps[k] / ndbl ); 
          }
          
          // XQSQ
          if ( has_XQSQ ) { 
            double coef = ctrans[from][to][1+k+3*ns]; 
            tprobs[l][from][to] += coef * ( all_qs[l][k] / total_qs ) * 
                                     ( all_qs[l][k] / total_qs ); 
          }
          
        }
        
      }
      
      // Compute cumsum 
      for ( uchar to=1; to<ns; to++ ) { 
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
  
  // Convert all_qs_combinations to a c array
  uchar all_qs[tprob_size][ns+1]; 
  for ( uword l=0; l<tprob_size; l++ ) { 
    for ( uword k=0; k<(ns+1); k++ ) { 
      all_qs[l][k] = all_qs_combinations(l, k); 
    }
  }
  
  // Convert armadillo array to c-style array. Note that we change the order things are 
  // stored so that coefficients are contiguous in memory. We also divide by substeps 
  // here so that rates are always under one. 
  double ctrans[ns][ns][ncoefs];
  for ( uword from=0; from<ns; from++ ) { 
    for ( uword to=0; to<ns; to++ ) { 
      for ( uword coef=0; coef<ncoefs; coef++ ) { 
        ctrans[from][to][coef] = trans(coef, to, from) / (double) substeps; 
      }
    }
  }
  
  // Initialize some things as c arrays
  // Note: we allocate omat/nmat on the heap since they can be big matrices and blow up 
  // the size of the C stack beyond what is acceptable.
  auto mat_a = new uchar[nr][nc];
  auto mat_b = new uchar[nr][nc];
  uchar (*old_mat)[nc] = mat_a; // pointer to first element of array wich is char[nc]
  uchar (*new_mat)[nc] = mat_b; 
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      old_mat[i][j] = (uchar) init(i, j);
      new_mat[i][j] = (uchar) init(i, j);
    }
  }
  
  // Initialize vector with counts of cells in each state in the landscape (used to 
  // compute global densities)
  uword ps_a[ns];
  uword ps_b[ns];
  memset(ps_a, 0, sizeof(ps_a));
  memset(ps_b, 0, sizeof(ps_b));
  
  // Start with ps_a as old_ps, and ps_b as new_ps
  uword *old_ps = ps_a; 
  uword *new_ps = ps_b; 
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      old_ps[ old_mat[i][j] ]++; 
      new_ps[ old_mat[i][j] ]++; 
    }
  }
  
  // Compute local densities 
  auto qs_a = new uchar[nr][nc][ns]; 
  auto qs_b = new uchar[nr][nc][ns]; 
  
  uchar (*old_qs)[nc][ns] = qs_a;  
  uchar (*new_qs)[nc][ns] = qs_b;  
  
  get_local_densities(old_qs, old_mat); 
  get_local_densities(new_qs, old_mat); 
  
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
  double trans_probs[tprob_size][ns][ns]; 
  
  uword iter = 0; 
  
  while ( iter <= niter ) { 
    
    // Call callbacks 
    if ( console_callback_active && iter % console_callback_every == 0 ) { 
      console_callback_wrap(iter, old_ps, console_callback); 
    }
    
    if ( cover_callback_active && iter % cover_callback_every == 0 ) { 
      cover_callback_wrap(iter, old_ps, cover_callback); 
    }
    
    if ( snapshot_callback_active && iter % snapshot_callback_every == 0 ) { 
      snapshot_callback_wrap(iter, old_mat, snapshot_callback); 
    }
    
    for ( uword substep=0; substep < substeps; substep++ ) { 
      // Compute transition probabilities 
      precompute_transition_probabilites(trans_probs, ctrans, old_ps, all_qs); 
      
      for ( uword i=0; i<nr; i++ ) { 
        
        for ( uword j=0; j<nc; j++ ) { 
          
          double rn = randunif();
          
          uchar cstate = old_mat[i][j]; 
          
          // Get line in pre-computed transition probability table 
          uword line = 0; 
          for ( uchar k = 0; k<ns; k++ ) { 
            line = line * (1+max_nb) + old_qs[i][j][k]; 
          }
          
          // Check if we actually transition.  
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation since the random number is always below one. 
          // TODO: the right number of substeps can be determined beforehand !!! 
          uchar new_cell_state = cstate; 
          for ( signed char k=(ns-1); k>=0; k-- ) { 
            new_cell_state = rn < trans_probs[line][cstate][k] ? k : new_cell_state; 
          } 
          
          if ( new_cell_state != cstate ) { 
            new_ps[new_cell_state]++; 
            new_ps[cstate]--; 
            new_mat[i][j] = new_cell_state; 
            adjust_local_densities(new_qs, i, j, cstate, new_cell_state); 
          }
        }
      }
      
      // Switch pointers to old/new ps/qs arrays, etc
      old_ps = new_ps; 
      old_qs = new_qs; 
      old_mat = new_mat; 
      
    } // end of substep loop
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] mat_a; 
  delete [] mat_b; 
  delete [] qs_a; 
  delete [] qs_b; 
  
}

  
