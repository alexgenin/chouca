
// #ifndef ARMA_NO_DEBUG
// #define ARMA_NO_DEBUG
// #endif 

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// using namespace Rcpp;

// Define some column names
#define _from 0
#define _to 1
#define _state 2
#define _qs 3
#define _coef 0
#define _expo 1

using namespace arma;

inline void get_local_densities(arma::Col<uword>& qs, 
                                const arma::Mat<ushort>& m, 
                                const arma::uword i, 
                                const arma::uword j, 
                                const bool wrap, 
                                const bool use_8_nb) { 
  qs.fill(0); 
  uword nc = m.n_cols; 
  uword nr = m.n_rows; 
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    ushort state_left = m( i, (nc + j - 1) % nc);
    qs( state_left )++; // left 
    
    ushort state_right = m( i, (nc + j + 1) % nc);
    qs( state_right )++; // right
    
    ushort state_up = m( (nr + i - 1) % nr, j);
    qs( state_up )++; // up
    
    ushort state_down = m( (nr + i + 1) % nr, j);
    qs( state_down )++; // down
    
    if ( use_8_nb ) { 
      
      ushort state_upleft = m( (nr + i - 1) % nr, (nc + j - 1) % nc); 
      qs( state_upleft )++; // upleft
      
      ushort state_upright = m( (nr + i - 1) % nr, (nc + j + 1) % nc); 
      qs( state_upright )++; // upright
      
      ushort state_downleft = m( (nr + i + 1) % nr, (nc + j - 1) % nc); 
      qs( state_downleft )++; // downleft
      
      ushort state_downright = m( (nr + i + 1) % nr, (nc + j + 1) % nc); 
      qs( state_downright )++; // downright
    }
    
  } else { 
    
    if ( j > 0 ) { 
      ushort state_left = m(i, j-1); 
      qs( state_left )++; // left
    }
    if ( j < (nc-1) ) { 
      ushort state_right = m(i, j+1); 
      qs( state_right )++; // right
    }
    if ( i > 0 ) { 
      ushort state_up = m(i-1, j);
      qs( state_up )++; // up
    }
    if ( i < (nr-1) ) { 
      ushort state_down = m(i+1, j); 
      qs( state_down )++; // down
    }
    
    if ( use_8_nb ) { 
      if ( i > 0 && j > 0 ) { 
        ushort state_upleft = m(i-1, j-1); 
        qs( state_upleft )++; // upleft
      }
      if ( i > 0 && j < (nc-1) ) { 
        ushort state_upright = m(i-1, j+1); 
        qs( state_upright )++; // upright
      }
      if ( i < (nr-1) && j > 0 ) { 
        ushort state_downleft = m(i+1, j-1); 
        qs( state_downleft )++; // downleft
      }
      if ( i < (nr-1) && j < (nc-1) ) { 
        ushort state_downright = m(i+1, j+1); 
        qs( state_downright )++; // downright
      }
    }
    
  }
  
}

// This is a function that returns the local state counts to R. i and j must be indexed 
// the R-way (1-indexing)
// [[Rcpp::export]]
arma::Col<arma::uword> local_dens(const arma::Mat<ushort> m, 
                                  const arma::uword nstates, 
                                  const arma::uword i, 
                                  const arma::uword j, 
                                  const bool wrap, 
                                  const bool use_8_nb) { 
  arma::Col<uword> newq(nstates);
  newq.fill(0); 
  // m, i and j must be adjusted to take into account the way things are stored in R
  get_local_densities(newq, m, i-1, j-1, wrap, use_8_nb); 
  
  return(newq);
}

inline void get_local_densities_column(arma::Mat<uword>& qs, 
                                       const arma::Mat<ushort>& m, 
                                       const arma::uword j, 
                                       const bool wrap, 
                                       const bool use_8_nb) { 
  
  qs.fill(0); 
  uword nc = m.n_cols; 
  uword nr = m.n_rows; 
  
  if ( wrap ) { 
    for ( uword i=0; i<nr; i++ ) { 
      // left column 
      ushort state_left = m( i, (nc + j - 1) % nc);
      qs(i, state_left)++; 
      
      // right column 
      ushort state_right = m( i, (nc + j + 1) % nc);
      qs(i, state_right)++; 
    
      ushort state_up = m( (nr + i - 1) % nr, j);
      qs(i, state_up)++; // up
      
      ushort state_down = m( (nr + i + 1) % nr, j);
      qs(i, state_down)++; // down
    }
    
    if ( use_8_nb ) { 
      for ( uword i=0; i<nr; i++ ) { 
        ushort state_upleft = m( (nr + i - 1) % nr, (nc + j - 1) % nc); 
        qs(i, state_upleft)++; // upleft
        
        ushort state_upright = m( (nr + i - 1) % nr, (nc + j + 1) % nc); 
        qs(i, state_upright)++; // upright
        
        ushort state_downleft = m( (nr + i + 1) % nr, (nc + j - 1) % nc); 
        qs(i, state_downleft)++; // downleft
        
        ushort state_downright = m( (nr + i + 1) % nr, (nc + j + 1) % nc); 
        qs(i, state_downright)++; // downright
      }
    }
  
  } else { 
    for ( uword i=0; i<nr; i++ ) { 
      if ( j > 0 ) { 
        ushort state_left = m(i, j-1); 
        qs(i, state_left)++; // left
      }
      if ( j < (nc-1) ) { 
        ushort state_right = m(i, j+1); 
        qs(i, state_right)++; // right
      }
      if ( i > 0 ) { 
        ushort state_up = m(i-1, j);
        qs(i, state_up)++; // up
      }
      if ( i < (nr-1) ) { 
        ushort state_down = m(i+1, j); 
        qs(i, state_down)++; // down
      }
      
      if ( use_8_nb ) { 
        if ( i > 0 && j > 0 ) { 
          ushort state_upleft = m(i-1, j-1); 
          qs(i, state_upleft)++; // upleft
        }
        if ( i > 0 && j < (nc-1) ) { 
          ushort state_upright = m(i-1, j+1); 
          qs(i, state_upright)++; // upright
        }
        if ( i < (nr-1) && j > 0 ) { 
          ushort state_downleft = m(i+1, j-1); 
          qs(i, state_downleft)++; // downleft
        }
        if ( i < (nr-1) && j < (nc-1) ) { 
          ushort state_downright = m(i+1, j+1); 
          qs(i, state_downright)++; // downright
        }
      }
    }
    
  }
  
  
}

// This is a function that returns the local state counts to R, for the full column of 
// a matrix. j must be indexed the R-way (1-indexing)
//[[Rcpp::export]]
arma::Mat<arma::uword> local_dens_col(const arma::Mat<ushort> m, 
                                      const arma::uword nstates, 
                                      const arma::uword j, 
                                      const bool wrap, 
                                      const bool use_8_nb) { 
  arma::Mat<uword> newq(m.n_rows, nstates);
  newq.fill(0); 
  // m, i and j must be adjusted to take into account the way things are stored in R
  get_local_densities_column(newq, m, j-1, wrap, use_8_nb); 
  
  return(newq);
}

// 
// This file contains the main c++ engine to run CAs 
// 
// [[Rcpp::export]]
void camodel_cpp_engine(const arma::Mat<ushort> alpha_index, 
                        const arma::Col<double> alpha_vals, 
                        const arma::Mat<ushort> pmat_index, 
                        const arma::Mat<double> pmat_vals, 
                        const arma::Mat<ushort> qmat_index, 
                        const arma::Col<double> qmat_vals, 
                        const Rcpp::List ctrl, 
                        const Rcpp::Function console_callback, 
                        const Rcpp::Function cover_callback, 
                        const Rcpp::Function snapshot_callback) { 
  
  // Unpack control list
  const uword substeps     = ctrl["substeps"]; 
  const bool  wrap         = ctrl["wrap"]; 
  const Mat<ushort> init   = ctrl["init"];
  const uword niter        = ctrl["niter"]; // TODO: think about overflow in those values
  const ushort ns          = ctrl["nstates"]; 
  const ushort nbs = ctrl["neighbors"]; // Needed to make sure conversion from list is OK
  const bool use_8_nb      = nbs == 8 ? true : false; 
  
  // Number of samples for qs
  const ushort xpoints = ctrl["xpoints"]; 
  
  bool console_callback_active = ctrl["console_callback_active"]; 
  uword console_callback_every = ctrl["console_callback_every"]; 
  
  bool cover_callback_active = ctrl["cover_callback_active"]; 
  uword cover_callback_every = ctrl["cover_callback_every"]; 
  
  bool snapshot_callback_active = ctrl["snapshot_callback_active"]; 
  uword snapshot_callback_every = ctrl["snapshot_callback_every"]; 
  uword nr = init.n_rows; 
  uword nc = init.n_cols; 
  double n = (double) nr * (double) nc;
  
  // Initialize some things 
  Mat<ushort> omat = init; 
  Mat<ushort> nmat = init; 
  
  // Initialize vector with counts of cells in each state in the landscape (used to 
  // compute global densities)
  Col<int> ps(ns); 
  Col<int> delta_ps(ns); 
  ps.fill(0); 
  delta_ps.fill(0); 
  
  for ( uword j=0; j<nc; j++ ) { 
    for ( uword i=0; i<nr; i++ ) { 
      ps(init(i, j))++; 
    }
  }
  
  uword iter = 0; 
  
  // Allocate some things we will reuse later 
  Mat<uword> qs(nr, ns);
  Col<double> qs_prop(ns);
  Col<double> ptrans(ns); 
  
  while ( iter <= niter ) { 
    
    // Call callbacks 
    if ( console_callback_active && iter % console_callback_every == 0 ) { 
      console_callback(iter, ps, n); 
    }
    
    if ( cover_callback_active && iter % cover_callback_every == 0 ) { 
      cover_callback(iter, ps, n); 
    }
    if ( snapshot_callback_active && iter % snapshot_callback_every == 0 ) { 
      snapshot_callback(iter, omat); 
    }
    
    for ( uword substep = 0; substep < substeps; substep++ ) { 
      
      for ( uword j=0; j<nc; j++ ) { 
        
        // Get local state counts for the column. Getting it for a column at once is 
        // slightly more cache-friendly.
        get_local_densities_column(qs, omat, j, wrap, use_8_nb); 
        
        for ( uword i=0; i<nr; i++ ) { 
          
          
          // Normalized local densities to proportions
          uword qs_total = accu(qs.row(i)); 
          
          // Factor to convert the number of neighbors into the point at which the 
          // dependency on q is sampled.
          uword qpointn_factorf = (xpoints - 1) / qs_total; 
          
          // Get current state 
          ushort cstate = omat(i, j); 
          
          // Compute probability transitions 
          for ( ushort to=0; to<ns; to++ ) { 
            
            // Init probability
            ptrans(to) = 0; 
            
            // Scan the table of alphas 
            for ( uword k=0; k<alpha_index.n_rows; k++ ) { 
              ptrans(to) += 
                ( alpha_index(k, _from) == cstate ) * 
                ( alpha_index(k, _to) == to) * 
                alpha_vals(k); 
            }
            
            // Scan the table of pmat to reconstruct probabilities -> where is ps?
            for ( uword k=0; k<pmat_index.n_rows; k++ ) { 
              ptrans(to) += 
                ( pmat_index(k, _from) == cstate ) * 
                ( pmat_index(k, _to) == to) * 
                pmat_vals(k, _coef) * pow( ps(pmat_index(k, _state)) / (double) n, 
                                           pmat_vals(k, _expo) );
            }
            
            // Scan the table of qmat to reconstruct probabilities 
            for ( uword k=0; k<qmat_index.n_rows; k++ ) { 
              
              // Lookup which point in the qs function we need to use for the 
              // current neighbor situation.
              uword qthis = qs(i, qmat_index(k, _state)) * qpointn_factorf;
              
              ptrans(to) += 
                ( qmat_index(k, _from) == cstate ) * 
                ( qmat_index(k, _to) == to) * 
                // Given the observed local abundance of this state, which line in 
                // qmat should be retained ? 
                ( qmat_index(k, _qs) == qthis ) * 
                qmat_vals(k); 
            }
            
//             ptrans(to) /= substeps; 
          }
          
          // Check if we actually transition. We scan all states and switch to the 
          // one with the highest probability. 
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          //                                           ^ rn > everything => no transition
          //       ^ 0 < rn < p0 => p0 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation and may never consider a given transition. 
          // 
          // ptrans = cumsum(ptrans); // alternative code, but slower because it needs
          //                          // a copy
          for ( ushort k=1; k<ns; k++ ) { 
            ptrans(k) += ptrans(k-1);
          }
          
          //TODO: consider substepping only the probability step, maybe that works 
          // better. 
          ushort new_cell_state = cstate; 
          double rn = Rf_runif(0, 1); 
          for ( signed short k=(ns-1); k>=0; k-- ) { 
            new_cell_state = rn < ptrans(k) ? k : new_cell_state; 
          } 
          
          if ( new_cell_state != cstate ) {
            ushort old_cell_state = nmat(i, j);
            delta_ps(new_cell_state)++; 
            delta_ps(old_cell_state)--;
            nmat(i, j) = new_cell_state; 
          }
          
        }
      }
      
      // Apply iteration changes in global densities and matrix state
      for ( ushort k=0; k<ns; k++ ) { 
        ps(k) += delta_ps(k); 
        delta_ps(k) = 0; 
      }
      omat = nmat; 
      
    } // end of substep loop
//     return 1; 
    iter++; 
  }
  
}
