
#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

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

//BUG: WHEN wrap is false, the total number of possible neighbors is 2 or 3 and not 4, 
// so it is more complicated than just deciding whether to divide by 4 or 8... 

using namespace arma;

inline void get_local_densities_c(arma::uword qs[ns], 
                                  const char m[nr][nc], 
                                  const arma::uword i, 
                                  const arma::uword j) { 
  
  // Set all counts to zero
  for ( char k=0; k<ns; k++ ) { 
    qs[k] = 0; 
  }
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    char state_left = m[i][(nc + j - 1) % nc];
    qs[ state_left ]++; // left 
    
    char state_right = m[i][(nc + j + 1) % nc];
    qs[ state_right ]++; // right
    
    char state_up = m[(nr + i - 1) % nr][j];
    qs[ state_up ]++; // up
    
    char state_down = m[(nr + i + 1) % nr][j];
    qs[ state_down ]++; // down
    
    if ( use_8_nb ) { 
      
      char state_upleft = m[(nr + i - 1) % nr][(nc + j - 1) % nc]; 
      qs[ state_upleft ]++; // upleft
      
      char state_upright = m[(nr + i - 1) % nr][(nc + j + 1) % nc]; 
      qs[ state_upright ]++; // upright
      
      char state_downleft = m[(nr + i + 1) % nr][(nc + j - 1) % nc]; 
      qs[ state_downleft ]++; // downleft
      
      char state_downright = m[(nr + i + 1) % nr][(nc + j + 1) % nc]; 
      qs[ state_downright ]++; // downright
    }
    
  } else { 
    
    if ( i > 0 ) { 
      char state_up = m[i-1][j];
      qs[ state_up ]++; // up
    }
    if ( i < (nr-1) ) { 
      char state_down = m[i+1][j]; 
      qs[ state_down ]++; // down
    }
    if ( j > 0 ) { 
      char state_left = m[i][j-1]; 
      qs[ state_left ]++; // left
    }
    if ( j < (nc-1) ) { 
      char state_right = m[i][j+1]; 
      qs[ state_right ]++; // right
    }
    
    if ( use_8_nb ) { 
      if ( i > 0 && j > 0 ) { 
        char state_upleft = m[i-1][j-1]; 
        qs[ state_upleft ]++; // upleft
      }
      if ( i > 0 && j < (nc-1) ) { 
        char state_upright = m[i-1][j+1]; 
        qs[ state_upright ]++; // upright
      }
      if ( i < (nr-1) && j > 0 ) { 
        char state_downleft = m[i+1][j-1]; 
        qs[ state_downleft ]++; // downleft
      }
      if ( i < (nr-1) && j < (nc-1) ) { 
        char state_downright = m[i+1][j+1]; 
        qs[ state_downright ]++; // downright
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
// This file contains the main c++ engine to run CAs 
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
  // NOTE: we allocate omat/nmat on the heap since they can be big matrices and blow up 
  // the size of the C stack beyond what is acceptable.
  auto omat = new char[nr][nc];
  auto nmat = new char[nr][nc];
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      omat[i][j] = (char) init(i, j);
      nmat[i][j] = (char) init(i, j);
    }
  }
  
  double ctrans[ncoefs][ns][ns];
  for ( uword i=0; i<ns; i++ ) { 
    for ( uword j=0; j<ns; j++ ) { 
      for ( uword k=0; k<ncoefs; k++ ) { 
        ctrans[k][i][j] = (double) trans(k, i, j);
      }
    }
  }
  
  // Initialize vector with counts of cells in each state in the landscape (used to 
  // compute global densities)
  uword ps[ns];
  memset(ps, 0, sizeof(ps));
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      ps[ omat[i][j] ]++; 
    }
  }
  
  uword iter = 0; 
  
  // Allocate some things we will reuse later 
  uword qs[ns];
  double ptrans[ns];
  double qs_norm = use_8_nb ? 8.0 : 4.0; 
  
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
      
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          omat[i][j] = nmat[i][j]; 
        }
      }
      
      for ( uword i=0; i<nr; i++ ) { 
        for ( uword j=0; j<nc; j++ ) { 
          
          double rn = Rf_runif(0, 1);
          
          char cstate = omat[i][j]; 
          
          // Get local densities
          get_local_densities_c(qs, omat, i, j); 
          uword total_qs = 0;
          for ( char k=0; k<ns; k++ ) { 
            total_qs += qs[k]; 
          }
          
          // Compute probabilities of transition
          for ( char col=0; col<ns; col++ ) { 
            
            ptrans[col] = (double) ctrans[0][col][cstate] / substeps; 
            for ( char k=0; k<ns; k++ ) { 
              
              double pcoef = (double) ctrans[1+k][col][cstate]; 
              ptrans[col] += pcoef * ( ps[k] / (double) n ) / substeps; 
              double qcoef = (double) ctrans[1+k+ns][col][cstate]; 
              ptrans[col] += qcoef * ( qs[k] /  (double) total_qs )  / substeps; 
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
              if ( rn > ptrans[newstate-1] && rn < ptrans[newstate] ) { 
                new_cell_state = newstate;
                cell_switch = true; 
              }
            }
          }
          
          if ( cell_switch ) {
            char old_cell_state = nmat[i][j];
            ps[new_cell_state]++; 
            ps[old_cell_state]--;
            nmat[i][j] = new_cell_state; 
          }
        }
      }
      
    }
    
    iter++; 
  }
  
  // Free up heap-allocated arrays
  delete [] omat; 
  delete [] nmat; 
  
}

  
