


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
// Fancy expression for max_nb^ns converted to unsigned integer. 
// The number of permutations in neighbors.
constexpr arma::uword tprob_size = (uword) exp( log( (double) (1+max_nb) ) * (double) ns ); 
constexpr char tprob_interval = wrap ? ( use_8_nb ? 8 : 4 ) : 1; 

// Include functions 
#include "__COMMON_HEADER__"

// Compute transition probabilities between all possible qs states 
void compute_transition_probabilites(double tprobs[tprob_size][ns][ns], 
                                     const double ctrans[ns][ns][ncoefs], 
                                     const arma::uword ps[ns], 
                                     const char all_qs[tprob_size][ns+1]) { 
  
  // Note tprob_interval here. In all combinations of neighbors, only some of them 
  // can be observed in the wild. If we wraparound, then the number of neighbors is 
  // constant, it is either 8 or 4. So we can compute the values in tprobs only at 
  // indices every 4 or 8 values. 
  
  for ( uword l=0; l<tprob_size; l+=tprob_interval ) { 
    // Get total of neighbors = last column of all_qs
    double total_qs = all_qs[l][ns]; 
    
    for ( char from=0; from<ns; from++ ) { 
      
      for ( char to=0; to<ns; to++ ) { 
        
        if ( has_X0 ) { 
          tprobs[l][from][to] = ctrans[from][to][0]; 
        } else { 
          tprobs[l][from][to] = 0.0; 
        }
        
        // Loop over coefs
        for ( char k=0; k<ns; k++ ) { 
          
          // All the if{} blocks below will be removed appropriately by the compiler.
          
          // XP
          if ( has_XP ) { 
            double coef = ctrans[from][to][1+k]; 
            tprobs[l][from][to] += coef * ( ps[k] / ndbl ); 
          }
          
          // XQ
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
      for ( char to=1; to<ns; to++ ) { 
        tprobs[l][from][to] += tprobs[l][from][to-1];
      }
      
    }
  }
  
  
}

constexpr arma::uword factorial(arma::uword n) { 
  return n > 0 ? n * factorial( n - 1 ) : 1;  
}

constexpr arma::uword prod_low_base = factorial(ns - 1); 

// Get line in the table of transition probabilities 
// This was more sweat that I would like to admit 
inline arma::uword getline(const char qs[nr][nc][ns], 
                           const arma::uword i, 
                           const arma::uword j) { 
  
  uword lines_before = 0; 
  uword remaining = max_nb; 
  
  // Init choose low 
  uword choose_low = ns - 1; 
  // Max prod_low = (ns-1-0)!
  uword prod_low = prod_low_base; 
  
  for ( uword k=1; k<ns; k++ ) { 
    
    // Compute choose(remaining - qi + choose_low, choose_low)
    // Adjust choose low and prod low
    choose_low--; 
    prod_low /= (choose_low + 1); 
    
    uword curqs = qs[i][j][k-1]; 
    
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
  
  // This is the line at which we want to pick up the probability vector
  return( lines_before ); 
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
  double trans_probs[tprob_size][ns][ns]; 
  
  uword iter = 0; 
  
  while ( iter <= niter ) { 
    
    // Call callbacks 
    if ( console_callback_active && iter % console_callback_every == 0 ) { 
      console_callback_wrap(iter, ps, console_callback); 
    }
    
    if ( cover_callback_active && iter % cover_callback_every == 0 ) { 
      cover_callback_wrap(iter, ps, cover_callback); 
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
          
          // Get line in pre-computed table
          uword line = getline(qs, i, j); 
//           for ( char k = 0; k<ns; k++ ) { 
//             line = line * (1+max_nb) + qs[i][j][k]; 
//           }
          
//           Rcpp::Rcout << "qs :" << 
//             (int) qs[i][j][0] << " " << (int) qs[i][j][1] << " " << (int) qs[i][j][2] << " -> line " << line << "\n"; 
          
          // Check if we actually transition.  
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation since the random number is always below one. 
          char new_cell_state = cstate; 
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

  
