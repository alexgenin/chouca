// 
// 
// 
// 


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

// Define some column names for clarity
#define _from 0
#define _to 1
#define _state 2
#define _qs 3
#define _coef 0
#define _expo 1

// We tell gcc to unroll loops, as we have many small loops. This can double 
// performance in some cases (!)
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
constexpr bool wrap = __WRAP__; 
constexpr bool use_8_nb = __USE_8_NB__; 
constexpr uword xpoints = __XPOINTS__; 
constexpr arma::uword substeps = __SUBSTEPS__; 
constexpr double ncells = nr * nc; 
constexpr uword alpha_nrow = __ALPHA_NROW__; 
constexpr uword pmat_nrow = __PMAT_NROW__; 
constexpr uword qmat_nrow = __QMAT_NROW__; 
constexpr uword all_qs_nrow = __ALL_QS_NROW__; 

// Whether we want to precompute probabilities or not 
__PRECOMPUTE_TRANS_PROBAS_DEFINE__

// The maximum number of neighbors
constexpr arma::uword max_nb = use_8_nb ? 8 : 4; 

// Include functions and type declarations
#include "__COMMON_HEADER__"

// Compute transition probabilities between all possible qs states 
inline void precompute_transition_probabilites(double tprobs[all_qs_nrow][ns][ns], 
                                               const uchar all_qs[all_qs_nrow][ns+1], 
                                               const arma::Mat<ushort>& alpha_index, 
                                               const arma::Col<double>& alpha_vals, 
                                               const arma::Mat<ushort>& pmat_index, 
                                               const arma::Mat<double>& pmat_vals, 
                                               const arma::Mat<ushort>& qmat_index, 
                                               const arma::Col<double>& qmat_vals, 
                                               const arma::uword ps[ns]) { 
  
  // Note tprob_interval here. In all combinations of neighbors, only some of them 
  // can be observed in the wild. If we wraparound, then the number of neighbors is 
  // constant, it is either 8 or 4. So we can compute the values in tprobs only at 
  // indices every 4 or 8 values. 
  
  for ( uword l=0; l<all_qs_nrow; l++ ) { 
    
    for ( uchar from=0; from<ns; from++ ) { 
      
      for ( uchar to=0; to<ns; to++ ) { 
        
        // Factor to convert the number of neighbors into the point at which the 
        // dependency on q is sampled.
        // all_qs[ ][ns] holds the total number of neighbors
        uword qpointn_factorf = (xpoints - 1) / all_qs[l][ns]; 
        
        // Init probability
        tprobs[l][from][to] = 0; 
        
        // Scan the table of alphas 
        for ( uword k=0; k<alpha_nrow; k++ ) { 
          tprobs[l][from][to] += 
            ( alpha_index(k, _from) == from ) * 
            ( alpha_index(k, _to) == to) * 
            alpha_vals(k); 
        }
        
        // Scan the table of pmat to reconstruct probabilities -> where is ps?
        for ( uword k=0; k<pmat_nrow; k++ ) { 
          tprobs[l][from][to] += 
            ( pmat_index(k, _from) == from ) * 
            ( pmat_index(k, _to) == to) * 
            pmat_vals(k, _coef) * pow( ps[pmat_index(k, _state)] / ncells, 
                                        pmat_vals(k, _expo) );
        }
        
        // Scan the table of qmat to reconstruct probabilities 
        for ( uword k=0; k<qmat_nrow; k++ ) { 
          
          // Lookup which point in the qs function we need to use for the 
          // current neighbor situation.
          uword qthis = all_qs[l][qmat_index(k, _state)] * qpointn_factorf;
          
          tprobs[l][from][to] += 
            ( qmat_index(k, _from) == from ) * 
            ( qmat_index(k, _to) == to) * 
            // Given the observed local abundance of this state, which line in 
            // qmat should be retained ? 
            ( qmat_index(k, _qs) == qthis ) * 
            qmat_vals(k); 
        }
      }
      
      // Compute cumsum 
      for ( ushort k=1; k<ns; k++ ) { 
        tprobs[l][from][k] += tprobs[l][from][k-1];
      }
      
    }
  }
  
}


// [[Rcpp::export]]
void aaa__FPREFIX__camodel_compiled_engine(const arma::Mat<ushort> alpha_index, 
                                           const arma::Col<double> alpha_vals, 
                                           const arma::Mat<ushort> pmat_index, 
                                           const arma::Mat<double> pmat_vals, 
                                           const arma::Mat<ushort> qmat_index, 
                                           const arma::Col<double> qmat_vals, 
                                           const arma::Mat<ushort> all_qs_arma, 
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
  auto old_mat = new uchar[nr][nc];
  auto new_mat = new uchar[nr][nc];
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      old_mat[i][j] = (uchar) init(i, j);
      new_mat[i][j] = (uchar) init(i, j);
    }
  }
  
#ifdef PRECOMPUTE_TRANS_PROBAS
  // Convert all_qs to char array 
  auto all_qs = new uchar[all_qs_nrow][ns+1]; 
  for ( uword i=0; i<all_qs_nrow; i++ ) { 
    for ( uword k=0; k<(ns+1); k++ ) { 
      all_qs[i][k] = (uchar) all_qs_arma(i, k); 
    }
  }
#endif
  
  // Initialize vector with counts of cells in each state in the landscape (used to 
  // compute global densities)
  uword old_ps[ns];
  uword new_ps[ns];
  memset(old_ps, 0, sizeof(old_ps));
  memset(new_ps, 0, sizeof(new_ps));
  
  // Compute local densities
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      old_ps[ old_mat[i][j] ]++; 
    }
  }
  memcpy(new_ps, old_ps, sizeof(uword)*ns); 
  
#ifdef PRECOMPUTE_TRANS_PROBAS
  // Matrix holding probability line in the precomputed table 
  auto old_pline = new uword[nr][nc]; 
  auto new_pline = new uword[nr][nc]; 
  initialize_prob_line(old_pline, old_mat); 
  memcpy(new_pline, old_pline, sizeof(uword)*nr*nc); 
  
  // Initialize table with precomputed probabilities 
  auto tprobs = new double[all_qs_nrow][ns][ns]; 
#else
  // Create tables that hold local densities 
  auto old_qs = new uchar[nr][nc][ns]; 
  auto new_qs = new uchar[nr][nc][ns]; 
  get_local_densities(old_qs, old_mat); 
  memcpy(new_qs, old_qs, sizeof(uchar)*nr*nc*ns); 
#endif
  
  // Initialize random number generator 
  // Initialize rng state using armadillo random integer function
  for ( uword i=0; i<4; i++ ) { 
    s[i] = randi<long>(); 
  }
  
  // The first number returned by the RNG is garbage, probably because randi returns 
  // something too short for s (i.e. not 64 bits long if I got it right).
  // TODO: find a way to feed 64 bits integers to xoshiro
  //   -> we could feed it 64/8 = 8 chars taken from randi(). 
  randunif(); 
  
  // Allocate some things we will reuse later
  double ptrans[ns]; 
  
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
    
#ifdef PRECOMPUTE_TRANS_PROBAS
    precompute_transition_probabilites(tprobs, 
                                       all_qs, 
                                       alpha_index, 
                                       alpha_vals, 
                                       pmat_index, 
                                       pmat_vals, 
                                       qmat_index, 
                                       qmat_vals, 
                                       old_ps); 
    // for ( uword l=0; l<all_qs_nrow; l++ ) { 
    //   for ( ushort j=0; j<ns; j++ ) {
    //     for ( ushort k=0; k<ns; k++ ) {
    //       Rcpp::Rcout << "l: " << l << " j: " << j << " k: " << k << " " << 
    //         tprobs[l][j][k] << " \n"; 
    //     }
    //   }
    // }
#endif 
    
    for ( uword i=0; i<nr; i++ ) { 
      
      for ( uword j=0; j<nc; j++ ) { 
          
        uchar cstate = old_mat[i][j]; 
        
#ifndef PRECOMPUTE_TRANS_PROBAS
        // Normalized local densities to proportions
        uword qs_total = number_of_neighbors(i, j); 
        
        // Factor to convert the number of neighbors into the point at which the 
        // dependency on q is sampled.
        uword qpointn_factorf = (xpoints - 1) / qs_total; 
        
        // Compute probability transitions 
        for ( ushort to=0; to<ns; to++ ) { 
          // Init probability
          ptrans[to] = 0; 
          
          // Scan the table of alphas 
          for ( uword k=0; k<alpha_nrow; k++ ) { 
            ptrans[to] += 
              ( alpha_index(k, _from) == cstate ) * 
              ( alpha_index(k, _to) == to) * 
              alpha_vals(k); 
          }
          
          // Scan the table of pmat to reconstruct probabilities -> where is ps?
          for ( uword k=0; k<pmat_nrow; k++ ) { 
            ptrans[to] += 
              ( pmat_index(k, _from) == cstate ) * 
              ( pmat_index(k, _to) == to) * 
              pmat_vals(k, _coef) * pow( old_ps[pmat_index(k, _state)] / ncells, 
                                          pmat_vals(k, _expo) );
          }
          
          // Scan the table of qmat to reconstruct probabilities 
          for ( uword k=0; k<qmat_nrow; k++ ) { 
            
            // Lookup which point in the qs function we need to use for the 
            // current neighbor situation.
            uword qthis = old_qs[i][j][qmat_index(k, _state)] * qpointn_factorf;
            
            ptrans[to] += 
              ( qmat_index(k, _from) == cstate ) * 
              ( qmat_index(k, _to) == to) * 
              // Given the observed local abundance of this state, which line in 
              // qmat should be retained ? 
              ( qmat_index(k, _qs) == qthis ) * 
              qmat_vals(k); 
          }
        }
        
        // Compute cumsum 
        for ( uchar k=1; k<ns; k++ ) { 
          ptrans[k] += ptrans[k-1];
        }
#endif
        // Check if we actually transition.  
        // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
        //               ^ p0 < rn < (p0+p1) => p1 wins
        // Of course the sum of probabilities must be lower than one, otherwise we are 
        // making an approximation since the random number is always below one. 
        for ( uword substep=0; substep < substeps; substep++ ) { 
          
          uchar new_cell_state = cstate; 
          double rn = randunif(); // get random number
          for ( signed char k=(ns-1); k>=0; k-- ) { 
#ifdef PRECOMPUTE_TRANS_PROBAS
            uword line = old_pline[i][j]; 
            new_cell_state = rn < tprobs[line][cstate][k] ? k : new_cell_state; 
#else 
            new_cell_state = rn < ptrans[k] ? k : new_cell_state; 
#endif
          } 
          
          if ( new_cell_state != cstate ) { 
            new_ps[new_cell_state]++; 
            new_ps[cstate]--; 
            new_mat[i][j] = new_cell_state; 
#ifdef PRECOMPUTE_TRANS_PROBAS
            // if ( j < (nc-1) ) { 
            //   Rcpp::Rcout << "from: " << (int) cstate << " to: " << (int) new_cell_state << 
            //     " at " << i << ", " << j << "\n"; 
            //   Rcpp::Rcout << "old_nbs at i/j+1: "; 
            //   
            //   for ( uword k=0; k<ns; k++ ) { 
            //     Rcpp::Rcout << (int) old_qs[(nr + i + 1) % nr][j][k] << " "; 
            //   }
            //   Rcpp::Rcout << "\n"; 
            //   
            //   Rcpp::Rcout << "new_nbs at i/j+1: "; 
            //   for ( uword k=0; k<ns; k++ ) { 
            //     Rcpp::Rcout << (int) new_qs[(nr + i + 1) % nr][j][k] << " "; 
            //   }
            //   Rcpp::Rcout << "\n"; 
            //   
            // }
            // TODO: this will update every neighbors of the target cell, doing work 
            // multiple times when two neighbors change state. Consider doing a 
            // dirty/clean matrix, and update all required lines at once. 
            adjust_nb_plines(new_pline, i, j, cstate, new_cell_state); 
            
            // Consider the neighbor below
            // if ( j < nc ) { 
            //   Rcpp::Rcout << "oline: " << old_pline[i][j+1] << 
            //   " nline: " << new_pline[i][j+1] << "\n"; 
            // }
            
            
#else 
            adjust_local_density(new_qs, i, j, cstate, new_cell_state); 
#endif
            // If the transition occurs in one of the substeps, we break the for loop
            break; 
          } 
        } // end of substep loop
      }
    } 
    
    // Copy old matrix to new, etc. 
    memcpy(old_ps,  new_ps,  sizeof(uword)*ns); 
    memcpy(old_mat, new_mat, sizeof(uchar)*nr*nc); 
#ifdef PRECOMPUTE_TRANS_PROBAS
    memcpy(old_pline, new_pline, sizeof(uword)*nr*nc); 
#else 
    memcpy(old_qs,  new_qs,  sizeof(uchar)*nr*nc*ns); // TODO: get rid of this
#endif
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] old_mat; 
  delete [] new_mat; 
#ifdef PRECOMPUTE_TRANS_PROBAS
  delete [] old_pline; 
  delete [] new_pline; 
#else
  delete [] old_qs; 
  delete [] new_qs; 
#endif
}

  
