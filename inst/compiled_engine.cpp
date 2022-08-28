// 
// This is the source code for the 'compiled' chouca engine. Its fields will be replaced and 
// compiled to run a PCA model. 
// 
// 


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

#define USE_OMP __USE_OMP__

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
constexpr arma::uword nr     = __NR__; 
constexpr arma::uword nc     = __NC__; 
constexpr uchar ns           = __NS__; 
constexpr bool wrap          = __WRAP__; 
constexpr bool use_8_nb      = __USE_8_NB__; 
constexpr bool fixed_nb      = __FIXED_NEIGHBOR_NUMBER__; 
constexpr uword xpoints      = __XPOINTS__; 
constexpr arma::uword substeps = __SUBSTEPS__; 
constexpr double ncells      = nr * nc; 
constexpr uword beta_0_nrow  = __BETA_0_NROW__; 
constexpr uword beta_p_nrow  = __BETA_P_NROW__; 
constexpr uword beta_q_nrow  = __BETA_Q_NROW__; 
constexpr uword beta_pq_nrow = __BETA_PQ_NROW__; 
constexpr uword all_qs_nrow  = __ALL_QS_NROW__; 
constexpr uword cores        = __CORES__; 

// Whether we want to precompute probabilities or not 
// TODO: convert to value def for clarity
__PRECOMPUTE_TRANS_PROBAS_DEFINE__

// The maximum number of neighbors
constexpr arma::uword max_nb = use_8_nb ? 8 : 4; 

// Include functions and type declarations
#include "__COMMON_HEADER__"

// Compute transition probabilities between all possible qs states 
inline void precompute_transition_probabilites(double tprobs[all_qs_nrow][ns][ns], 
                                               const uchar all_qs[all_qs_nrow][ns+1], 
                                               const arma::Mat<ushort>& beta_0_index, 
                                               const arma::Col<double>& beta_0_vals, 
                                               const arma::Mat<ushort>& beta_p_index, 
                                               const arma::Mat<double>& beta_p_vals, 
                                               const arma::Mat<ushort>& beta_q_index, 
                                               const arma::Col<double>& beta_q_vals, 
                                               const arma::Mat<ushort>& beta_pq_index, 
                                               const arma::Mat<double>& beta_pq_vals, 
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
        for ( uword k=0; k<beta_0_nrow; k++ ) { 
          tprobs[l][from][to] += 
            ( beta_0_index(k, _from) == from ) * 
            ( beta_0_index(k, _to) == to) * 
            beta_0_vals(k); 
        }
        
        // Scan the table of pmat to reconstruct probabilities -> where is ps?
        for ( uword k=0; k<beta_p_nrow; k++ ) { 
          double p = ps[beta_p_index(k, _state)] / ncells; 
          
          tprobs[l][from][to] += 
            ( beta_p_index(k, _from) == from ) * 
            ( beta_p_index(k, _to) == to) * 
            beta_p_vals(k, _coef) * pow(p, beta_p_vals(k, _expo));
        }
        
        // Scan the table of beta_q to reconstruct probabilities 
        for ( uword k=0; k<beta_q_nrow; k++ ) { 
          
          // Lookup which point in the qs function we need to use for the 
          // current neighbor situation.
          uword qthis = all_qs[l][beta_q_index(k, _state)] * qpointn_factorf;
          
          tprobs[l][from][to] += 
            ( beta_q_index(k, _from) == from ) * 
            ( beta_q_index(k, _to) == to) * 
            // Given the observed local abundance of this state, which line in 
            // beta_q should be retained ? 
            ( beta_q_index(k, _qs) == qthis ) * 
            beta_q_vals(k); 
        }
        
        // Scan the table of beta_pq to reconstruct probabilities 
        for ( uword k=0; k<beta_pq_nrow; k++ ) { 
          
          // Lookup which point in the qs function we need to use for the 
          // current neighbor situation.
          // all_qs[ ][ns] holds the total number of neighbors
          double q = (double) all_qs[l][beta_pq_index(k, _state)] / all_qs[l][ns];
          double p = (double) ps[beta_pq_index(k, _state)] / ncells; 
          
          tprobs[l][from][to] += 
            ( beta_pq_index(k, _from) == from ) * 
            ( beta_pq_index(k, _to) == to) * 
            beta_pq_vals(k, _coef) * 
            pow(q * p, beta_pq_vals(k, _expo));
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
void aaa__FPREFIX__camodel_compiled_engine(const arma::Mat<ushort> all_qs_arma, 
                                           const Rcpp::List ctrl) { 
  
  // Unpack control list
  const Mat<ushort> init = ctrl["init"]; // this is ushort because init is an arma mat
  const uword niter      = ctrl["niter"]; // TODO: think about overflow in those values
  
  const uword console_callback_every = ctrl["console_output_every"]; 
  const bool console_callback_active = console_callback_every > 0; 
  Rcpp::Function console_callback = ctrl["console_callback"]; 
  
  const uword cover_callback_every = ctrl["save_covers_every"]; 
  const bool cover_callback_active = cover_callback_every > 0; 
  Rcpp::Function cover_callback = ctrl["cover_callback"]; 
  
  const uword snapshot_callback_every = ctrl["save_snapshots_every"]; 
  const bool snapshot_callback_active = snapshot_callback_every > 0; 
  Rcpp::Function snapshot_callback = ctrl["snapshot_callback"]; 
  
  const uword custom_callback_every = ctrl["custom_output_every"]; 
  const bool custom_callback_active = custom_callback_every > 0; 
  Rcpp::Function custom_callback = ctrl["custom_callback"]; 
  
  // Extract things from list 
  const arma::Mat<ushort> beta_0_index  = ctrl["beta_0_index"];
  const arma::Col<double> beta_0_vals   = ctrl["beta_0_vals"];
  const arma::Mat<ushort> beta_p_index  = ctrl["beta_p_index"];
  const arma::Mat<double> beta_p_vals   = ctrl["beta_p_vals"];
  const arma::Mat<ushort> beta_q_index  = ctrl["beta_q_index"];
  const arma::Col<double> beta_q_vals   = ctrl["beta_q_vals"];
  const arma::Mat<ushort> beta_pq_index = ctrl["beta_pq_index"];
  const arma::Mat<double> beta_pq_vals  = ctrl["beta_pq_vals"];
  
  // Copy some things as c arrays. Convert
  // Note: we allocate omat/nmat on the heap since they can be big matrices and blow up 
  // the size of the C stack beyond what is acceptable.
  auto old_mat = new uchar[nr][nc];
  auto new_mat = new uchar[nr][nc];
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      old_mat[i][j] = (uchar) init(i, j);
    }
  }
  memcpy(new_mat, old_mat, sizeof(uchar)*nr*nc); 
  
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
  init_local_densities(old_qs, old_mat); 
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
  
#if USE_OMP
  // This matrix holds all random numbers when we use openmp. 
  auto randnums = new double[nr][nc]; 
#endif
  
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
    
    if ( custom_callback_active && iter % custom_callback_every == 0 ) { 
      custom_callback_wrap(iter, old_mat, custom_callback); 
    }
    
#ifdef PRECOMPUTE_TRANS_PROBAS
    precompute_transition_probabilites(tprobs, 
                                       all_qs, 
                                       beta_0_index, 
                                       beta_0_vals, 
                                       beta_p_index, 
                                       beta_p_vals, 
                                       beta_q_index, 
                                       beta_q_vals, 
                                       beta_pq_index, 
                                       beta_pq_vals, 
                                       old_ps); 
#endif 
    
    for ( uword substep=0; substep < substeps; substep++ ) { 
      
#if USE_OMP
      // Fill array with random numbers
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          randnums[i][j] = randunif(); 
        }
      }
#pragma omp parallel for num_threads(cores) schedule(static) default(shared) private(ptrans) 
#endif
      for ( uword i=0; i<nr; i++ ) { 
       // TODO: when parallel, we could walk the matrix with an offset: this guarantees 
       // that we will not update two neighboring cells, and thus that there will not 
       // be race conditions when updating.   
        for ( uword j=0; j<nc; j++ ) { 
          
          uchar cstate = old_mat[i][j]; 
          
#ifdef PRECOMPUTE_TRANS_PROBAS
#else
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
            for ( uword k=0; k<beta_0_nrow; k++ ) { 
              ptrans[to] += 
                ( beta_0_index(k, _from) == cstate ) * 
                ( beta_0_index(k, _to) == to) * 
                beta_0_vals(k); 
            }
            
            // Scan the table of pmat to reconstruct probabilities -> where is ps?
            for ( uword k=0; k<beta_p_nrow; k++ ) { 
              ptrans[to] += 
                ( beta_p_index(k, _from) == cstate ) * 
                ( beta_p_index(k, _to) == to) * 
                beta_p_vals(k, _coef) * pow(old_ps[beta_p_index(k, _state)] / ncells, 
                                            beta_p_vals(k, _expo) );
            }
            
            // Scan the table of qmat to reconstruct probabilities 
            for ( uword k=0; k<beta_q_nrow; k++ ) { 
              
              // Lookup which point in the qs function we need to use for the 
              // current neighbor situation.
              uword qthis = old_qs[i][j][beta_q_index(k, _state)] * qpointn_factorf;
              
              ptrans[to] += 
                ( beta_q_index(k, _from) == cstate ) * 
                ( beta_q_index(k, _to) == to) * 
                // Given the observed local abundance of this state, which line in 
                // qmat should be retained ? 
                ( beta_q_index(k, _qs) == qthis ) * 
                beta_q_vals(k); 
            }
            
            // Scan the table of pqmat to reconstruct probabilities 
            for ( uword k=0; k<beta_pq_nrow; k++ ) { 
              
              // Lookup which point in the qs function we need to use for the 
              // current neighbor situation.
              // all_qs[ ][ns] holds the total number of neighbors
              double pq = (double) old_qs[i][j][beta_pq_index(k, _state)] / 
                            (double) qs_total;
              pq *= (double) old_ps[beta_pq_index(k, _state)] / ncells; 
              
              ptrans[to] += 
                ( beta_pq_index(k, _from) == cstate ) * 
                ( beta_pq_index(k, _to) == to) * 
                beta_pq_vals(k, _coef) * 
                pow(pq, beta_pq_vals(k, _expo));
            }
            
          }
          
          // Compute cumsum 
          for ( uchar k=1; k<ns; k++ ) { 
            ptrans[k] += ptrans[k-1];
          }
#endif
          
#if USE_OMP
          double rn = randnums[i][j]; // get random number
#else 
          double rn = randunif(); 
#endif
          
          // Check if we actually transition.  
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation since the random number is always below one. 
          uchar new_cell_state = cstate; 
          for ( signed char k=(ns-1); k>=0; k-- ) { 
#ifdef PRECOMPUTE_TRANS_PROBAS
            uword line = old_pline[i][j]; 
            new_cell_state = rn < tprobs[line][cstate][k] ? k : new_cell_state; 
#else 
            new_cell_state = rn < ptrans[k] ? k : new_cell_state; 
#endif
          } 
          
          if ( new_cell_state != cstate ) { 
#if USE_OMP
            // Only one thread can call this block at any time. It is not thread safe as 
            // shared matrices are used. This tanks the performance unfortunately when
            // probabilities are precomputed, because each for loop iteration has trivial
            // work to do, and in model with lots of cell switches, the lock is hit 
            // very often.
#pragma omp critical (line_adjustment) hint(omp_sync_hint_uncontended)
            {
#endif
            new_ps[new_cell_state]++; 
            new_ps[cstate]--; 
            new_mat[i][j] = new_cell_state; 
#ifdef PRECOMPUTE_TRANS_PROBAS
            adjust_nb_plines(new_pline, i, j, cstate, new_cell_state); 
#else 
            adjust_local_density(new_qs, i, j, cstate, new_cell_state); 
#endif
#if USE_OMP
            }
#endif
            
          } 
        } 
      }
      
      // Copy old matrix to new, etc. 
      memcpy(old_ps,  new_ps,  sizeof(uword)*ns); 
      memcpy(old_mat, new_mat, sizeof(uchar)*nr*nc); 
#ifdef PRECOMPUTE_TRANS_PROBAS
      memcpy(old_pline, new_pline, sizeof(uword)*nr*nc); 
#else 
      memcpy(old_qs,  new_qs,  sizeof(uchar)*nr*nc*ns); 
#endif
    } // end of substep loop
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] old_mat; 
  delete [] new_mat; 
#if USE_OMP
  delete [] randnums; 
#endif
#ifdef PRECOMPUTE_TRANS_PROBAS
  delete [] old_pline; 
  delete [] new_pline; 
#else
  delete [] old_qs; 
  delete [] new_qs; 
#endif
}

  
