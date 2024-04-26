
// Depth of the betas_index array at which to pick the start/end for each coefficient
constexpr uword _beta_0_start  = 0;
constexpr uword _beta_0_end    = 1;
constexpr uword _beta_q_start  = 2;
constexpr uword _beta_q_end    = 3;
constexpr uword _beta_pp_start = 4;
constexpr uword _beta_pp_end   = 5;
constexpr uword _beta_qq_start = 6;
constexpr uword _beta_qq_end   = 7;
constexpr uword _beta_pq_start = 8;
constexpr uword _beta_pq_end   = 9;

// Define some column names for clarity
#define _state_1 2
#define _state_2 3
#define _expo_1 4
#define _expo_2 5

/* This is xoshiro256+
 * https://prng.di.unimi.it/xoshiro256plus.c
 */
static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

// Seeds
static uint64_t s[CORES][4];

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

// Not very efficient integer power of another value. We should probably use
// exponentiation by squaring, though since we almost always deal with powers <=5,
// typically 1-2, this may be overkill.
// Note that we know at compile time the maximum value for b, so there may be
// some optimization we are missing, especially when max(b) <= 1.
inline uword intpow(const u_nbcount a,
                    const u_nbcount b) {
  uword p = 1;
  for ( uchar k=0; k<b; k++ ) {
    p *= a;
  }
  return p;
}

// Another try at fast power. Remarks above also apply, this is probably not very
// efficient. for information tests seem to say ~10% cpu time spent in here when
// not memoization probas of transition.
static inline double fintpow(const pfloat a,
                             const u_nbcount b) {
  
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

// Initalize the count of local densities for each cell in the landscape
void init_local_densities(u_nbcount qs[NR][NC][NS],
                          const u_state m[NR][NC]) {

  // Set all counts to zero
  memset(qs, 0, sizeof(u_nbcount)*NR*NC*NS);

  for ( uword i=0; i<NR; i++ ) {
    for ( uword j=0; j<NC; j++ ) {
      // Get neighbors to the left
      if ( WRAP ) {

        u_state state_left = m[i][(NC + j - 1) % NC];
        qs[i][j][state_left]++; // left

        u_state state_right = m[i][(NC + j + 1) % NC];
        qs[i][j][state_right]++; // right

        u_state state_up = m[(NR + i - 1) % NR][j];
        qs[i][j][state_up]++; // up

        u_state state_down = m[(NR + i + 1) % NR][j];
        qs[i][j][state_down]++; // down

        if ( USE_8_NB ) {

          u_state state_upleft = m[(NR + i - 1) % NR][(NC + j - 1) % NC];
          qs[i][j][state_upleft]++; // upleft

          u_state state_upright = m[(NR + i - 1) % NR][(NC + j + 1) % NC];
          qs[i][j][state_upright]++; // upright

          u_state state_downleft = m[(NR + i + 1) % NR][(NC + j - 1) % NC];
          qs[i][j][state_downleft]++; // downleft

          u_state state_downright = m[(NR + i + 1) % NR][(NC + j + 1) % NC];
          qs[i][j][state_downright]++; // downright
        }

      } else {

        if ( i > 0 ) {
          u_state state_up = m[i-1][j];
          qs[i][j][state_up]++; // up
        }
        if ( i < (NR-1) ) {
          u_state state_down = m[i+1][j];
          qs[i][j][state_down]++; // down
        }
        if ( j > 0 ) {
          u_state state_left = m[i][j-1];
          qs[i][j][state_left]++; // left
        }
        if ( j < (NC-1) ) {
          u_state state_right = m[i][j+1];
          qs[i][j][state_right]++; // right
        }

        if ( USE_8_NB ) {
          if ( i > 0 && j > 0 ) {
            u_state state_upleft = m[i-1][j-1];
            qs[i][j][state_upleft]++; // upleft
          }
          if ( i > 0 && j < (NC-1) ) {
            u_state state_upright = m[i-1][j+1];
            qs[i][j][state_upright]++; // upright
          }
          if ( i < (NR-1) && j > 0 ) {
            u_state state_downleft = m[i+1][j-1];
            qs[i][j][state_downleft]++; // downleft
          }
          if ( i < (NR-1) && j < (NC-1) ) {
            u_state state_downright = m[i+1][j+1];
            qs[i][j][state_downright]++; // downright
          }
        }

      }
    }
  }
}

inline u_nbcount number_of_neighbors(const u_xyint i,
                                     const u_xyint j) {
  // When we do not wrap, the total number of neighbors depends on where we are
  // in the matrix, e.g. first column or last row has missing neighbors. In that
  // case, we need to correct for that. When wrapping is on, then the
  // number of neighbors is constant.
  //
  // The compiler will remove the if{} statements below at compile time if we are
  // wrapping, making this function return a constant, defined on the line below:
  u_nbcount nnb = USE_8_NB ? 8 : 4;

  // If we use a fixed number of neighbors, then return early. In this case, cells on
  // the edges will have a lower chance of switching, but this may be an approximation
  // we are willing to make for better performance.
  if ( FIXED_NB ) {
    return( nnb );  
  }

  // If we do not wrap and we use 8 neighbors, then we just need to substract from the
  // total (maximum) counts.
  if ( ! WRAP && ! USE_8_NB ) {
    if ( i == 0 || i == (NR-1) ) { nnb--; }  // First or last row
    if ( j == 0 || j == (NC-1) ) { nnb--; }  // First or last column
  }

  // If we do not wrap and we use 8 neighbors, then it is a bit more subtle because we
  // also need to check we are not in a corner
  if ( ! WRAP && USE_8_NB ) {
    if ( i == 0 || i == (NR-1) ) { nnb = nnb-3; }  // First or last row
    if ( j == 0 || j == (NC-1) ) { nnb = nnb-3; }  // First or last column

    // If we are in the corners, then one removed neighbor is being counted twice, so we
    // need to add it back
    if ( (i == 0 && j == 0) ||
         (i == 0 && j == (NC-1)) ||
         (i == (NR-1) && j == 0) ||
         (i == (NR-1) && j == (NC-1)) ) {
      nnb++;
    }

  }

  return nnb;
}

// This function will set the "probability lines", ie fill in the array that contains
// for each neighborhood configuration, where to pick its probability of transitions
void initialize_prob_line(u_pline prob_lines[NR][NC],
                          const u_state m[NR][NC]) {
  
  for ( u_xyint i=0; i<NR; i++ ) {
    for ( u_xyint j=0; j<NC; j++ ) {
      
      u_nbcount this_qs[NS];
      memset(this_qs, 0, sizeof(this_qs));

      // Get neighbors to the left
      if ( WRAP ) {

        u_state state_left = m[i][(NC + j - 1) % NC];
        this_qs[state_left]++; // left

        u_state state_right = m[i][(NC + j + 1) % NC];
        this_qs[state_right]++; // right

        u_state state_up = m[(NR + i - 1) % NR][j];
        this_qs[state_up]++; // up

        u_state state_down = m[(NR + i + 1) % NR][j];
        this_qs[state_down]++; // down

        if ( USE_8_NB ) {

          u_state state_upleft = m[(NR + i - 1) % NR][(NC + j - 1) % NC];
          this_qs[state_upleft]++; // upleft

          u_state state_upright = m[(NR + i - 1) % NR][(NC + j + 1) % NC];
          this_qs[state_upright]++; // upright

          u_state state_downleft = m[(NR + i + 1) % NR][(NC + j - 1) % NC];
          this_qs[state_downleft]++; // downleft

          u_state state_downright = m[(NR + i + 1) % NR][(NC + j + 1) % NC];
          this_qs[state_downright]++; // downright
        }

      } else {

        if ( i > 0 ) {
          u_state state_up = m[i-1][j];
          this_qs[state_up]++; // up
        }
        if ( i < (NR-1) ) {
          u_state state_down = m[i+1][j];
          this_qs[state_down]++; // down
        }
        if ( j > 0 ) {
          u_state state_left = m[i][j-1];
          this_qs[state_left]++; // left
        }
        if ( j < (NC-1) ) {
          u_state state_right = m[i][j+1];
          this_qs[state_right]++; // right
        }

        if ( USE_8_NB ) {
          if ( i > 0 && j > 0 ) {
            u_state state_upleft = m[i-1][j-1];
            this_qs[state_upleft]++; // upleft
          }
          if ( i > 0 && j < (NC-1) ) {
            u_state state_upright = m[i-1][j+1];
            this_qs[state_upright]++; // upright
          }
          if ( i < (NR-1) && j > 0 ) {
            u_state state_downleft = m[i+1][j-1];
            this_qs[state_downleft]++; // downleft
          }
          if ( i < (NR-1) && j < (NC-1) ) {
            u_state state_downright = m[i+1][j+1];
            this_qs[state_downright]++; // downright
          }
        }

      }

      // Get line in pre-computed transition probability table
      u_pline line = 0;
      for ( u_state k = 0; k<NS; k++ ) {
        line = line * (1+max_nb) + this_qs[k];
      }
      line -= 1;

      // If we have constant number of neighbors, then all_qs only contains the
      // values at each max_nb values, so we need to divide by max_nb here to fall on the
      // right line in the table of probabilities.
      if ( WRAP ) {
        line = (line+1) / max_nb - 1;
      }

      prob_lines[i][j] = line;
    }
  }

}

// Adjust the local densities of the neighboring cells of a cell that changed state
// from state 'from' to state 'to' at position 'i', 'j'.
inline void adjust_local_density(u_nbcount qs[NR][NC][NS],
                                 const u_xyint i,
                                 const u_xyint j,
                                 const u_state from,
                                 const u_state to) {

  // Get neighbors to the left
  if ( WRAP ) {

    // left
    qs[i][(NC + j - 1) % NC][to]++;
    qs[i][(NC + j - 1) % NC][from]--;

    // right
    qs[i][(NC + j + 1) % NC][to]++;
    qs[i][(NC + j + 1) % NC][from]--;

    // up
    qs[(NR + i - 1) % NR][j][to]++;
    qs[(NR + i - 1) % NR][j][from]--;

    // down
    qs[(NR + i + 1) % NR][j][to]++;
    qs[(NR + i + 1) % NR][j][from]--;

    if ( USE_8_NB ) {
      qs[(NR + i - 1) % NR][(NC + j - 1) % NC][to]++; // upleft
      qs[(NR + i - 1) % NR][(NC + j - 1) % NC][from]--;

      qs[(NR + i - 1) % NR][(NC + j + 1) % NC][to]++; // upright
      qs[(NR + i - 1) % NR][(NC + j + 1) % NC][from]--;

      qs[(NR + i + 1) % NR][(NC + j - 1) % NC][to]++; // downleft
      qs[(NR + i + 1) % NR][(NC + j - 1) % NC][from]--;

      qs[(NR + i + 1) % NR][(NC + j + 1) % NC][to]++; // downright
      qs[(NR + i + 1) % NR][(NC + j + 1) % NC][from]--;
    }

  } else {

    if ( i > 0 ) {
      qs[i-1][j][to]++; // up
      qs[i-1][j][from]--;
    }
    if ( i < (NR-1) ) {
      qs[i+1][j][to]++; // down
      qs[i+1][j][from]--;
    }
    if ( j > 0 ) {
      qs[i][j-1][to]++; // left
      qs[i][j-1][from]--;
    }
    if ( j < (NC-1) ) {
      qs[i][j+1][to]++; // right
      qs[i][j+1][from]--;
    }

    if ( USE_8_NB ) {
      if ( i > 0 && j > 0 ) {
        qs[i-1][j-1][to]++; // upleft
        qs[i-1][j-1][from]--;
      }
      if ( i > 0 && j < (NC-1) ) {
        qs[i-1][j+1][to]++; // upright
        qs[i-1][j+1][from]--;
      }
      if ( i < (NR-1) && j > 0 ) {
        qs[i+1][j-1][to]++; // downleft
        qs[i+1][j-1][from]--;
      }
      if ( i < (NR-1) && j < (NC-1) ) {
        qs[i+1][j+1][to]++; // downright
        qs[i+1][j+1][from]--;
      }
    }

  }

}

// When we memoize probabilities, we do not need to adjust the local densities, but
// we need to change the line at which to pick the probability of transition
inline void adjust_nb_plines(u_pline pline[NR][NC],
                             const u_xyint i,
                             const u_xyint j,
                             const u_state from,
                             const u_state to) {

  // Remember that the table of possible neighborhood configurations is organized like 
  // this: 
  //        q_3  q_2  q_1  total
  //       [,1] [,2] [,3] [,4]
  // [1,]    0    0    0    0
  // [2,]    0    0    1    1
  // [3,]    0    0    2    2
  // [4,]    0    0    3    3
  // [5,]    0    0    4    4
  // [6,]    0    1    0    1
  // [7,]    0    1    1    2
  // [8,]    0    1    2    3
  // [9,]    0    1    3    4
  // [10,]   0    1    4    5
  // [11,]   0    2    0    2
  // 
  // where the three first columns are the states, and the last column is the total 
  // number of neighbors. This means that the number of lines before a configuration 
  // is equal to 
  // 
  //  \sum_k=1^NS (nb+1)^(k-1) q_k where 
  //  
  //   NS is the number of states 
  //   nb is the max number of neighbors 
  //   q_k is the number of neighbors in state k 
  // 
  // for example, for line 10 above the number of lines to skip to reach it is equal to 
  //   (nb+1)^0 * 4 + (nb+1)^1*1 + (nb+1)^2 * 0 = 4 + 5 = 9 
  // 
  //  ...and there are indeed 9 lines before line 10. There is some adjustment to do 
  //  to take into account 0-indexing vs 1-indexing, but this is the gist of it. 
  // 
  // Because this is a sum, when one neighbor changes state, this means that one of 
  // those columns changes, so we need to remove the term corresponding 
  // to the old state, and add the term corresponding to the new state, hence the 
  // following adjustment factor: 
  // 
  // NOTE: we use a signed integer here because the pline adjustment can be negative (we 
  // go up in the table)
  s_pline adj = intpow(max_nb+1, (NS-1) - to) - intpow(max_nb+1, (NS-1) - from);
  
  // If we have a fixed number of neighbors, then we discard the combinations that 
  // do not sum to this number of neighbors. In the table above for example, this is 
  // not done so you have configurations with, for example, 3 neighbors, that do not 
  // correspond to something that occurs in the grid. This effectively removes every 
  // line but those that fall every 'nb' (the number of neighbors), so we need to divide
  // the adjustment factor accordingly (i.e. instead of going nb-by-nb, we go one by
  // one). 
  adj = WRAP ? adj / ( (s_pline) max_nb ) : adj;

  // Get neighbors to the left
  if ( WRAP ) {

    // left
    pline[i][(NC + j - 1) % NC] += adj;

    // right
    pline[i][(NC + j + 1) % NC] += adj;

    // up
    pline[(NR + i - 1) % NR][j] += adj;

    // down
    pline[(NR + i + 1) % NR][j] += adj;

    if ( USE_8_NB ) {
      // upleft
      pline[(NR + i - 1) % NR][(NC + j - 1) % NC] += adj;

      // upright
      pline[(NR + i - 1) % NR][(NC + j + 1) % NC] += adj;

      // downleft
      pline[(NR + i + 1) % NR][(NC + j - 1) % NC] += adj;

      // downright
      pline[(NR + i + 1) % NR][(NC + j + 1) % NC] += adj;
    }

  } else {

    // left
    if ( i > 0 ) {
      pline[i-1][j] += adj;
    }
    // right
    if ( i < (NR-1) ) {
      pline[i+1][j] += adj;
    }
    // up
    if ( j > 0 ) {
      pline[i][j-1] += adj;
    }
    // down
    if ( j < (NC-1) ) {
      pline[i][j+1] += adj;
    }

    if ( USE_8_NB ) {
      // upleft
      if ( i > 0 && j > 0 ) {
        pline[i-1][j-1] += adj;
      }
      // upright
      if ( i > 0 && j < (NC-1) ) {
        pline[i-1][j+1] += adj;
      }
      // downleft
      if ( i < (NR-1) && j > 0 ) {
        pline[i+1][j-1] += adj;
      }
      // downrighgt
      if ( i < (NR-1) && j < (NC-1) ) {
        pline[i+1][j+1] += adj;
      }
    }

  }

}

// Compute probability of transitions to other states, and store in tprob_line
inline void compute_rate(pfloat tprob_line[NS],
                         const u_nbcount qs[NS],
                         const u_pscount ps[NS],
                         const u_nbcount & total_nb,
                         const u_state & from) {

  // Set all probs to zero
  memset(tprob_line, 0, sizeof(pfloat)*NS);

  for ( u_state to=0; to<NS; to++ ) {

    // Check if we will ever transition into the state. This assumes
    // probability is set to zero above, otherwise tprob_line contains undefined or
    // wrong values.
    if ( ! transition_matrix[from][to] ) {
      continue;
    }

    // double total = 0.0. 
    pfloat total = 0.0; 
    s_pline kstart, kend;

    // constant component
    if ( BETA_0_NROW > 0 ) {
      kstart = betas_index[_beta_0_start][from][to];
      kend   = betas_index[_beta_0_end][from][to];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        total += coef_tab_flts[k];
      }
    }

    // q component
    if ( BETA_Q_NROW > 0 ) {
      kstart = betas_index[_beta_q_start][from][to];
      kend = betas_index[_beta_q_end][from][to];
      uword n_coef_to_sum = (kend - kstart + 1) / XPOINTS;
      if ( kstart >= 0 ) {
        // qthis holds the current point at which to extract the value of the coefficient
        // for q. When a transition depends on multiple neighbor expressions, i.e.
        // when P(a -> b) = f(q_1) + f(q_2), then we need to sum both the
        // corresponding coefficient for f(q_1), but also the one for f(q_2).
        // We know they are spaced xpoints apart from each other in the
        // coefficient table, so we just sum values every xpoints
        
        // qpointn is a number used to convert the number of neighbors into 
        // the point at which the value of f(q) is stored in the coefficient tables. 
        // all_qs[ ][NS] holds the total number of neighbors, and is passed as 
        // total_nb. However, here fixed_nb is constexpr so this if() should be 
        // optimized away when we use a fixed number of neighbors as qpointn
        // is then known at compile time. 
        const u_nbcount qpointn = ( XPOINTS - 1 ) / ( FIXED_NB ? N_NB : total_nb ); 
        
        for ( uword coef=0; coef<n_coef_to_sum; coef++) {
          uword qthis = qs[ coef_tab_ints[kstart + XPOINTS*coef][_state_1] ] * qpointn;
          total += coef_tab_flts[kstart + qthis + XPOINTS*coef];
        }
      }
    }

    // pp
    if ( BETA_PP_NROW > 0 ) {
      kstart = betas_index[_beta_pp_start][from][to];
      kend   = betas_index[_beta_pp_end][from][to];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat p1 = ps[ coef_tab_ints[k][_state_1] ] / FNCELLS;
        const pfloat p2 = ps[ coef_tab_ints[k][_state_2] ] / FNCELLS;
        total += coef_tab_flts[k] *
          fintpow(p1, coef_tab_ints[k][_expo_1]) *
          fintpow(p2, coef_tab_ints[k][_expo_2]);
      }
    }

    // qq
    if ( BETA_QQ_NROW > 0 ) {
      kstart = betas_index[_beta_qq_start][from][to];
      kend   = betas_index[_beta_qq_end][from][to];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat q1 = (pfloat) qs[ coef_tab_ints[k][_state_1] ] / total_nb;
        const pfloat q2 = (pfloat) qs[ coef_tab_ints[k][_state_2] ] / total_nb;
        total += coef_tab_flts[k] *
          fintpow(q1, coef_tab_ints[k][_expo_1]) *
          fintpow(q2, coef_tab_ints[k][_expo_2]);
      }
    }

    // pq
    if ( BETA_PQ_NROW > 0 ) {
      kstart = betas_index[_beta_pq_start][from][to];
      kend   = betas_index[_beta_pq_end][from][to];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat p1 = ps[ coef_tab_ints[k][_state_1] ] / FNCELLS;
        const pfloat q1 = (pfloat) qs[ coef_tab_ints[k][_state_2] ] / total_nb;
        total += coef_tab_flts[k] *
          fintpow(p1, coef_tab_ints[k][_expo_1]) *
          fintpow(q1, coef_tab_ints[k][_expo_2]);
      }
    }

    // NOTE: it seems like a good idea to cap probabilities to 0-1, but it
    // is really not. The fact that things are above 1 should be handled
    // by increasing the number of substeps, otherwise the rates of
    // probabilities become unbalanced between transitions. For example,
    // if we have P(s -> a) = 2 and P(s -> b) = 0.4, then by setting
    // P(s -> a) to 1.0, we reduce its relative probability compared to
    // P(s -> b) when using substeps >= 2.
    // Cap the proba to 0-1 range
    // total = total > 1.0 ? 1.0 : total;
    total = total < 0.0 ? 0.0 : total;
    
    // Total contains the probability for this transition. We do the cumsum right away, 
    // while we used to do it seperately below. 
    tprob_line[to] = total + ( to == 0 ? 0 : tprob_line[to-1] );
  }

  // Compute cumsum of probabilities. This is already included when saving the 
  // resulting value in the code line above. 
  // for ( u_state k=1; k<NS; k++ ) {
  //  tprob_line[k] += tprob_line[k-1];
  //}

}


// These wrappers are here to convert back c-style arrays to armadillo arrays, which
// Rcpp can understand.
//
void console_callback_wrap(const u_tstep iter,
                           const u_pscount ps[NS],
                           const Rcpp::Function console_callback) {
  uvec ps_arma(NS); 
  for ( u_state k=0; k<NS; k++ ) {
    ps_arma(k) = (arma::uword) ps[k];
  }
  console_callback(iter, ps_arma, FNCELLS);
}

void snapshot_callback_wrap(const u_tstep t,
                            const u_state omat[NR][NC],
                            const Rcpp::Function snapshot_callback) {

  // Make arma array to give back to R
  arma::Mat<ushort> m(NR, NC);
  for ( uword i=0; i<NR; i++ ) {
    for ( uword j=0; j<NC; j++ ) {
      m(i,j) = (ushort) omat[i][j];
    }
  }

  snapshot_callback(t, m);
}

void custom_callback_wrap(const u_tstep t,
                          const u_state omat[NR][NC],
                          const Rcpp::Function custom_callback) {

  // Make arma array to give back to R
  arma::Mat<ushort> m(NR, NC);
  for ( uword i=0; i<NR; i++ ) {
    for ( uword j=0; j<NC; j++ ) {
      m(i,j) = (ushort) omat[i][j];
    }
  }

  custom_callback(t, m);
}

void cover_callback_wrap(const u_tstep t,
                         const u_pscount ps[NS],
                         Rcpp::Function cover_callback) {
  uvec ps_arma(NS);
  for ( u_state k=0; k<NS; k++ ) {
    ps_arma(k) = ps[k];
  }
  
  cover_callback(t, ps_arma, FNCELLS);
}
