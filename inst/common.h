
/* This is xoshiro256+ 
 * https://prng.di.unimi.it/xoshiro256plus.c 
 */ 
static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

// Seeds 
static uint64_t s[cores][4];

static uint64_t nextr(uchar core) {
  const uint64_t result = s[core][0] + s[core][3];
  const uint64_t t = s[core][1] << 17;
  
  s[core][2] ^= s[core][0];
  s[core][3] ^= s[core][1];
  s[core][1] ^= s[core][2];
  s[core][0] ^= s[core][3];
  
  s[core][2] ^= t;
  
  s[core][3] = rotl(s[core][3], 45);
  
  return result;
}

static inline double randunif(uchar core) { 
  uint64_t x = nextr(core); 
  double xf = (x >> 11) * 0x1.0p-53; // use upper 53 bits only 
  return xf; 
}

inline uword intpow(const uchar a, 
                    const uchar b) { 
  uword p = 1; 
  for ( uchar k=0; k<b; k++ ) { 
    p *= a; 
  }
  return p; 
}

// Fast power when b is integer between zero and five (the maximum degree considered 
// for models)
static inline double fintpow(const double a, 
                             const uchar  b) { 
  switch ( b ) { 
    case 0: return 1.0; 
    case 1: return a; 
    case 2: return a*a; 
    case 3: return a*a*a; 
    case 4: return a*a*a*a; 
    case 5: return a*a*a*a*a; 
    case 6: return a*a*a*a*a*a; 
    case 7: return a*a*a*a*a*a*a; 
  }
  
  return 1.0; 
}

inline void init_local_densities(uchar qs[nr][nc][ns], 
                                 const uchar m[nr][nc]) { 
  
  // Set all counts to zero
  memset(qs, 0, sizeof(uchar)*nr*nc*ns); 
  
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
  
  // If we use a fixed number of neighbors, then return early. In this case, cells on 
  // the edges will have a lower chance of switching, but this may be an approximation 
  // we are willing to make for better performance. 
  if ( fixed_nb ) { 
    return( nnb ); 
  }
  
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
  
  // NOTE: we use a signed integer here because the pline adjustment can be negative 
  sword adj = intpow(max_nb+1, (ns-1) - to) - intpow(max_nb+1, (ns-1) - from); 
  
  // If we wrap, then what would go max_nb by max_nb goes instead one by one, so 
  // we need to divide the line jump here. 
  adj = wrap ? adj / ( (sword) max_nb ) : adj; 
  
  // Get neighbors to the left 
  if ( wrap ) { 
    
    // left 
    pline[i][(nc + j - 1) % nc] += adj; 
    
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

// Compute probability components 
static inline void compute_rate(double tprob_line[ns], 
                                const uchar qs[ns], 
                                const uword ps[ns], 
                                const arma::uword & qpointn, 
                                const arma::uword & total_nb, 
                                const arma::uword & from, 
                                const double delta_t) { 
  
  for ( ushort to=0; to<ns; to++ ) { 
    double total = 0.0; 
    
    // constant component
    for ( uword k=0; k<beta_0_nrow; k++ ) { 
      total += 
        ( beta_0_index(k, _from) == from ) * 
        ( beta_0_index(k, _to) == to) * 
        beta_0_vals(k); 
    }
    
    // f(q) component
    for ( uword k=0; k<beta_q_nrow; k++ ) { 
      
      // Lookup which point in the qs function we need to use for the 
      // current neighbor situation.
      uword qthis = qs[beta_q_index(k, _state_1)] * qpointn;
      
      total += 
        ( beta_q_index(k, _from) == from ) * 
        ( beta_q_index(k, _to) == to) * 
        // Given the observed local abundance of this state, which line in 
        // beta_q should be retained ? 
        ( beta_q_index(k, _qs) == qthis ) * 
        beta_q_vals(k); 
    }

    // pp
    for ( uword k=0; k<beta_pp_nrow; k++ ) { 
      double p1 = ps[beta_pp_index(k, _state_1)] / ncells; 
      double p2 = ps[beta_pp_index(k, _state_2)] / ncells; 
      
      total += 
        ( beta_pp_index(k, _from) == from ) * 
        ( beta_pp_index(k, _to) == to) * 
        beta_pp_vals(k, _coef) * fintpow(p1, beta_pp_index(k, _expo_1)) * 
          fintpow(p2, beta_pp_index(k, _expo_2));
    }
    
    // qq
    for ( uword k=0; k<beta_qq_nrow; k++ ) { 
      double q1 = (double) qs[beta_qq_index(k, _state_1)] / total_nb;
      double q2 = (double) qs[beta_qq_index(k, _state_2)] / total_nb; //TODO state_1 ????
      
      total += 
        ( beta_qq_index(k, _from) == from ) * 
        ( beta_qq_index(k, _to) == to) * 
        beta_qq_vals(k, _coef) * fintpow(q1, beta_qq_index(k, _expo_1)) * 
          fintpow(q2, beta_qq_index(k, _expo_2));
    }
    
    // pq
    // Rcpp::Rcout << beta_pq_vals << "\n";
    // Rcpp::Rcout << beta_pq_index << "\n";
    for ( uword k=0; k<beta_pq_nrow; k++ ) { 
      double p1 = ps[beta_pq_index(k, _state_1)] / ncells; 
      double q1 = (double) qs[beta_pq_index(k, _state_2)] / total_nb;
      
      // Rcpp::Rcout << "k: " << k << "\n"; 
      // Rcpp::Rcout << "beta_pq_nrow: " << beta_pq_nrow << "\n"; 
      // Rcpp::Rcout << "beta_pq_ix_ncol: " << beta_pq_index.n_cols << "\n"; 
      // Rcpp::Rcout << "beta_pq_vals_ncol: " << beta_pq_vals.n_cols << "\n"; 
      // Rcpp::Rcout << "_expo_1: " << _expo_1 << "\n"; 
      // Rcpp::Rcout << "_expo_2: " << _expo_2 << "\n"; 
      // Rcpp::Rcout << "_to: " << _to << "\n"; 
      // Rcpp::Rcout << "_from: " << _from << "\n"; 
      // Rcpp::Rcout << beta_pq_vals << "\n"; 
      
      total += 
        ( beta_pq_index(k, _from) == from ) * 
        ( beta_pq_index(k, _to) == to) * 
        beta_pq_vals(k, _coef) * 
          fintpow(p1, beta_pq_index(k, _expo_1)) * 
          fintpow(q1, beta_pq_index(k, _expo_2));
    }
    
    tprob_line[to] = total; 
  }

  
  // tprobs holds the rate at which things transition, we need to transform it if we 
  // are working with a continuous-time SCA
  // P(i -> j) = ( 1 - exp( - rate * delta_t ) ) * ( rate / sum(rates) ) 
  // TODO: This is "wrong", it only works asymptotically as delta_t -> 0. 
  // We need to compute the probability that after xx time units 
  // have passed, the chain is in a given state. We need to simulate the switches 
  // between time steps. 
  if ( continuous_sca ) { 
    
    double sum_rates = 0; 
    for ( ushort k=0; k<ns; k++ ) { 
      sum_rates += tprob_line[k]; 
    }
    
    for ( ushort k=0; k<ns; k++ ) { 
      double rate = tprob_line[k]; 
      // Rcpp::Rcout << "from: " << from << " to: " << (int) k << ": " << tprob_line[k] << 
        // " sum: " << sum_rates << " delta_t: " << delta_t << "\n"; 
      tprob_line[k] = ( 1 - exp( - rate * delta_t ) ) * ( rate / sum_rates ); 
    }
    
    
    // Now tprobs holds the probabilities of switching for a continuous SCA. If we 
    // were working with a discrete SCA, there would be nothing to do here. 
  }
  
  // Cap probabilities to values above zero 
  // for ( ushort k=0; k<ns; k++ ) { 
    // tprob_line[k] = tprob_line[k] < 0 ? 0 : tprob_line[k]; 
  // }
  
  // Compute cumsum of probabilities
  for ( ushort k=1; k<ns; k++ ) { 
    tprob_line[k] += tprob_line[k-1];
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

void snapshot_callback_wrap(const double t, 
                            const uchar omat[nr][nc],
                            const Rcpp::Function snapshot_callback) { 
    
  // Make arma array to give back to R
  Mat<ushort> m(nr, nc);
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      m(i,j) = (ushort) omat[i][j]; 
    }
  }
  
  snapshot_callback(t, m); 
}

void custom_callback_wrap(const double t, 
                          const uchar omat[nr][nc],
                          const Rcpp::Function custom_callback) { 
    
  // Make arma array to give back to R
  Mat<ushort> m(nr, nc);
  for ( uword i=0; i<nr; i++ ) { 
    for ( uword j=0; j<nc; j++ ) { 
      m(i,j) = (ushort) omat[i][j]; 
    }
  }
  
  custom_callback(t, m); 
}

void cover_callback_wrap(const double t, 
                         const arma::uword ps[ns], 
                         Rcpp::Function cover_callback) { 
  uvec ps_arma(ns); 
  for ( uchar k=0; k<ns; k++ ) { 
    ps_arma(k) = ps[k]; 
  }
  
  cover_callback(t, ps_arma, ncells); 
}
