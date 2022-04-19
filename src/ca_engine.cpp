
#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// using namespace Rcpp;

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
arma::Mat<arma::uword>local_dens_col(const arma::Mat<ushort> m, 
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
void camodel_cpp_engine(const arma::cube trans, 
                        const Rcpp::List ctrl, 
                        const Rcpp::Function console_callback, 
                        const Rcpp::Function cover_callback, 
                        const Rcpp::Function snapshot_callback) { 
  
  // Unpack control list
  uword substeps     = ctrl["substeps"]; 
  bool  wrap         = ctrl["wrap"]; 
  Mat<ushort> init   = ctrl["init"];
  uword niter        = ctrl["niter"]; // TODO: think about overflow in those values
  ushort ns          = ctrl["nstates"]; 
  bool use_8_nb      = ctrl["use_8_neighbors"]; 
  
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
  
  // Take into account substeps by dividing rates 
  arma::cube trans_sub = trans; 
  trans_sub /= (double) substeps; 
  
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
          
          // Get a random number 
          double rn = Rf_runif(0, 1); 
          
          // Normalized local densities to proportions
          double qs_total = accu(qs.row(i)); 
          
          // Get current state 
          ushort cstate = omat(i, j); 
          
          // Compute probability transitions 
          for ( ushort col=0; col<ns; col++ ) { 
            // X0
            ptrans(col) = trans_sub(0, col, cstate); 
            // Add global and local density components
            for ( ushort k=0; k<ns; k++ ) { 
              // XP
              ptrans(col) += trans_sub(1+k, col, cstate) * 
                               ( ps(k) / (double) n);  
              // XQ
              ptrans(col) += trans_sub(1+k+ns, col, cstate) * 
                               ( qs(i, k) / (double) qs_total ); 
              // XPSQ
              ptrans(col) += trans_sub(1+k+2*ns, col, cstate) * 
                               ( ps(k) * ps(k) / (double) n / (double) n ); 
              // XQSQ
              ptrans(col) += trans_sub(1+k+3*ns, col, cstate) * 
                        ( qs(i, k) * qs(i, k) / (double) qs_total / (double) qs_total);
            }
          }
          
          // Check if we actually transition. We scan all states and switch to the 
          // one with the highest probability. 
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are 
          // making an approximation. 
          // 
          // ptrans = cumsum(ptrans); // alternative code, but slower because it needs
          //                          // a copy
          for ( ushort k=1; k<ptrans.n_elem; k++ ) { 
            ptrans(k) += ptrans(k-1);
          }
          
          ushort new_cell_state=0; 
          bool cell_switch = false; 
          if ( rn < ptrans(0) ) { 
            new_cell_state = 0;
            cell_switch = true; 
          } else { 
            for ( ushort newstate=1; newstate<ns; newstate++ ) { 
              if ( rn > ptrans(newstate-1) && rn < ptrans(newstate) ) { 
                new_cell_state = newstate;
                cell_switch = true; 
              }
            }
          }
          
          if ( cell_switch ) {
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
    
    iter++; 
  }
  
}
