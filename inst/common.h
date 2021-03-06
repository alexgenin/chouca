
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


inline uword intpow(uchar a, 
                    uchar b) { 
  uword p = 1; 
  for ( uchar k=0; k<b; k++ ) { 
    p *= a; 
  }
  return p; 
}

inline void get_local_densities(uchar qs[nr][nc][ns], 
                                const uchar m[nr][nc]) { 
  
  // Set all counts to zero
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      for ( uchar k=0; k<ns; k++ ) { 
        qs[i][j][k] = 0; 
      }
    }
  }
  
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      // Get neighbors to the left 
      if ( wrap ) { 
        
        uchar state_left = m[i][(nc + j - 1) % nc];
        qs[i][j][state_left]++; // left 
        
        uchar state_right = m[i][(nc + j + 1) % nc];
        qs[i][j][state_right]++; // right
        
        uchar state_up = m[(nr + i - 1) % nr][j];
        qs[i][j][state_up]++; // up
        
        uchar state_down = m[(nr + i + 1) % nr][j];
        qs[i][j][state_down]++; // down
        
        if ( use_8_nb ) { 
          
          uchar state_upleft = m[(nr + i - 1) % nr][(nc + j - 1) % nc]; 
          qs[i][j][state_upleft]++; // upleft
          
          uchar state_upright = m[(nr + i - 1) % nr][(nc + j + 1) % nc]; 
          qs[i][j][state_upright]++; // upright
          
          uchar state_downleft = m[(nr + i + 1) % nr][(nc + j - 1) % nc]; 
          qs[i][j][state_downleft]++; // downleft
          
          uchar state_downright = m[(nr + i + 1) % nr][(nc + j + 1) % nc]; 
          qs[i][j][state_downright]++; // downright
        }
        
      } else { 
        
        if ( i > 0 ) { 
          uchar state_up = m[i-1][j];
          qs[i][j][state_up]++; // up
        }
        if ( i < (nr-1) ) { 
          uchar state_down = m[i+1][j]; 
          qs[i][j][state_down]++; // down
        }
        if ( j > 0 ) { 
          uchar state_left = m[i][j-1]; 
          qs[i][j][state_left]++; // left
        }
        if ( j < (nc-1) ) { 
          uchar state_right = m[i][j+1]; 
          qs[i][j][state_right]++; // right
        }
        
        if ( use_8_nb ) { 
          if ( i > 0 && j > 0 ) { 
            uchar state_upleft = m[i-1][j-1]; 
            qs[i][j][state_upleft]++; // upleft
          }
          if ( i > 0 && j < (nc-1) ) { 
            uchar state_upright = m[i-1][j+1]; 
            qs[i][j][state_upright]++; // upright
          }
          if ( i < (nr-1) && j > 0 ) { 
            uchar state_downleft = m[i+1][j-1]; 
            qs[i][j][state_downleft]++; // downleft
          }
          if ( i < (nr-1) && j < (nc-1) ) { 
            uchar state_downright = m[i+1][j+1]; 
            qs[i][j][state_downright]++; // downright
          }
        }
        
      }
    }
  }
  
}

inline uword number_of_neighbors(const arma::uword i, 
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
  uword nnb = use_8_nb ? 8 : 4; 
  
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

// This function will set the probability lines
inline void initialize_prob_line(arma::uword prob_lines[nr][nc], 
                                 const uchar m[nr][nc]) { 
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      uword this_qs[ns]; 
      memset(this_qs, 0, sizeof(this_qs)); 
      
      // Get neighbors to the left 
      if ( wrap ) { 
        
        uchar state_left = m[i][(nc + j - 1) % nc];
        this_qs[state_left]++; // left 
        
        uchar state_right = m[i][(nc + j + 1) % nc];
        this_qs[state_right]++; // right
        
        uchar state_up = m[(nr + i - 1) % nr][j];
        this_qs[state_up]++; // up
        
        uchar state_down = m[(nr + i + 1) % nr][j];
        this_qs[state_down]++; // down
        
        if ( use_8_nb ) { 
          
          uchar state_upleft = m[(nr + i - 1) % nr][(nc + j - 1) % nc]; 
          this_qs[state_upleft]++; // upleft
          
          uchar state_upright = m[(nr + i - 1) % nr][(nc + j + 1) % nc]; 
          this_qs[state_upright]++; // upright
          
          uchar state_downleft = m[(nr + i + 1) % nr][(nc + j - 1) % nc]; 
          this_qs[state_downleft]++; // downleft
          
          uchar state_downright = m[(nr + i + 1) % nr][(nc + j + 1) % nc]; 
          this_qs[state_downright]++; // downright
        }
        
      } else { 
        
        if ( i > 0 ) { 
          uchar state_up = m[i-1][j];
          this_qs[state_up]++; // up
        }
        if ( i < (nr-1) ) { 
          uchar state_down = m[i+1][j]; 
          this_qs[state_down]++; // down
        }
        if ( j > 0 ) { 
          uchar state_left = m[i][j-1]; 
          this_qs[state_left]++; // left
        }
        if ( j < (nc-1) ) { 
          uchar state_right = m[i][j+1]; 
          this_qs[state_right]++; // right
        }
        
        if ( use_8_nb ) { 
          if ( i > 0 && j > 0 ) { 
            uchar state_upleft = m[i-1][j-1]; 
            this_qs[state_upleft]++; // upleft
          }
          if ( i > 0 && j < (nc-1) ) { 
            uchar state_upright = m[i-1][j+1]; 
            this_qs[state_upright]++; // upright
          }
          if ( i < (nr-1) && j > 0 ) { 
            uchar state_downleft = m[i+1][j-1]; 
            this_qs[state_downleft]++; // downleft
          }
          if ( i < (nr-1) && j < (nc-1) ) { 
            uchar state_downright = m[i+1][j+1]; 
            this_qs[state_downright]++; // downright
          }
        }
        
      }
      
      // Get line in pre-computed transition probability table 
      uword line = 0; 
      for ( uchar k = 0; k<ns; k++ ) { 
        line = line * (1+max_nb) + this_qs[k]; 
      }
      line -= 1; 
      
      // If we have constant number of neighbors, then all_qs only contains the 
      // values at each max_nb values, so we need to divide by max_nb here to fall on the 
      // right line in the table of probabilities. 
      if ( wrap ) { 
        line = (line+1) / max_nb - 1; 
      }
      
      prob_lines[i][j] = line; 
    }
  }
  
}

// Adjust the local densities of the neighboring cells of a cell along with the lines 
// of the precomputations. 
inline void adjust_local_density(uchar qs[nr][nc][ns], 
                                 const uword i, 
                                 const uword j, 
                                 const uchar from, 
                                 const uchar to) { 
  
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

inline void adjust_nb_plines(uword pline[nr][nc], 
                             const uword i, 
                             const uword j, 
                             const uchar from, 
                             const uchar to) { 
  
  // NOTE: we used a signed integer here because the pline adjustment can be negative 
  sword adj = intpow(max_nb+1, (ns-1) - to) - intpow(max_nb+1, (ns-1) - from); 
  
  // If we wrap, then what would go max_nb by max_nb goes instead one by one, so 
  // we need to divide the line jump here. 
  adj = wrap ? adj / ( (sword) max_nb ) : adj; 
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    // left 
//     Rcpp::Rcout << "adjusting_pline" << "\n"; 
//     Rcpp::Rcout << "pline_old: " << pline[i][(nc + j - 1) % nc] << "\n"; 
//     Rcpp::Rcout << "pline_adj: " << pline_adjustment(from, to) << "\n"; 
    pline[i][(nc + j - 1) % nc] += adj; 
//     Rcpp::Rcout << "pline_new: " << pline[i][(nc + j - 1) % nc] << "\n"; 
    
    // right
    pline[i][(nc + j + 1) % nc] += adj; 
    
    // up
    pline[(nr + i - 1) % nr][j] += adj; 
    
    // down
    pline[(nr + i + 1) % nr][j] += adj; 
    
    if ( use_8_nb ) { 
      // upleft
      pline[(nr + i - 1) % nr][(nc + j - 1) % nc] += adj; 
      
      // upright
      pline[(nr + i - 1) % nr][(nc + j + 1) % nc] += adj; 
      
      // downleft
      pline[(nr + i + 1) % nr][(nc + j - 1) % nc] += adj; 
      
      // downright
      pline[(nr + i + 1) % nr][(nc + j + 1) % nc] += adj; 
    }
    
  } else { 
    
    // left
    if ( i > 0 ) { 
      pline[i-1][j] += adj; 
    }
    // right
    if ( i < (nr-1) ) { 
      pline[i+1][j] += adj; 
    }
    // up
    if ( j > 0 ) { 
      pline[i][j-1] += adj; 
    }
    // down
    if ( j < (nc-1) ) { 
      pline[i][j+1] += adj; 
    }
    
    if ( use_8_nb ) { 
      // upleft
      if ( i > 0 && j > 0 ) { 
        pline[i-1][j-1] += adj; 
      }
      // upright
      if ( i > 0 && j < (nc-1) ) { 
        pline[i-1][j+1] += adj; 
      }
      // downleft
      if ( i < (nr-1) && j > 0 ) { 
        pline[i+1][j-1] += adj; 
      }
      // downrighgt
      if ( i < (nr-1) && j < (nc-1) ) { 
        pline[i+1][j+1] += adj; 
      }
    }
    
  }
  
}

void console_callback_wrap(const arma::uword iter, 
                           const arma::uword ps[ns], 
                           const Rcpp::Function console_callback) { 
  uvec ps_arma(ns); // double
  for ( uchar k=0; k<ns; k++ ) { 
    ps_arma(k) = ps[k]; 
  }
  console_callback(iter, ps_arma, ncells); 
}

void snapshot_callback_wrap(const arma::uword iter, 
                            const uchar omat[nr][nc],
                            const Rcpp::Function snapshot_callback) { 
    
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
                         Rcpp::Function cover_callback) { 
  uvec ps_arma(ns); 
  for ( uchar k=0; k<ns; k++ ) { 
    ps_arma(k) = ps[k]; 
  }
  cover_callback(iter, ps_arma, ncells); 
}
