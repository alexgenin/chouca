
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

// Define some things known at compile time
const s_xyint kernel_semiheight = (NB_KERNEL_NR - 1) / 2; 
const s_xyint kernel_semiwidth  = (NB_KERNEL_NC - 1) / 2; 

// Define a maximum function
inline s_xyint MAX(s_xyint a, 
                   s_xyint b) { 
  if ( b <= a ) { 
    return a; 
  } else { 
    return b; 
  }
}

/* This is xoshiro256+
 * https://prng.di.unimi.it/xoshiro256plus.c
 */
static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

// Seeds
static uint64_t s[CORES][4];

static uint64_t nextr(arma::uword core) {
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

static inline double randunif(arma::uword core) {
  uint64_t x = nextr(core);
  double xf = (x >> 11) * 0x1.0p-53; // use upper 53 bits only
  return xf;
}

// Not very efficient integer power of another value. We should probably use
// exponentiation by squaring, though since we almost always deal with powers <=5,
// typically 1-2, this may be overkill.
// Note that we know at compile time the maximum value for b, so there may be
// some optimization we are missing, especially when max(b) <= 1.
inline arma::uword intpow(const u_nbcount a,
                          const u_nbcount b) {
  uword p = 1;
  for ( uchar k=0; k<b; k++ ) {
    p *= a;
  }
  return p;
}

// Another try at fast power. Remarks above also apply, this is probably not very
// efficient. for information tests seem to say ~10% cpu time spent in here when
// probas of transition are not memoized.
static inline double fintpow(const pfloat x,
                             const u_nbcount b) {
  // Here we make available to the compiler the maximum value for b, this may 
  // enable some optimizations. 
  const u_nbcount bmax = b <= MAX_POW_DEGREE ? b : MAX_POW_DEGREE; 
  
  pfloat ans = 1.0; 
  for ( u_nbcount i=0; i<bmax; i++ ) { 
    ans *= x; 
  }

  return ans;
}

// Return the product x^a * y^b, while enabling some optimizations when possible
static inline double xy_ab_product(const pfloat x, 
                                   const pfloat y, 
                                   const u_nbcount a,   // aka EXPO_1
                                   const u_nbcount b) { // aka EXPO_2
  
  pfloat ans; 
  
  // If we will never call this function with b > 0, then discard that part of 
  // the multiplication
  // if ( MAX_PP_EXPO_2 == 0 && MAX_PQ_EXPO_2 == 0 && MAX_QQ_EXPO_2 == 0 ) { 
  //   
  //   // If the exponents for a are one or zero, then this function will be always called 
  //   // with a == 1, so we just return x. Note that a is guaranteed to be above zero, 
  //   // because the corresponding coefficients are removed from the tables when expo_1 is 
  //   // zero
  //   if ( MAX_PP_EXPO_1 == 1 && MAX_PQ_EXPO_1 == 1 && MAX_QQ_EXPO_1 == 1 ) {
  //     return x; 
  //   }
  //   
  //   ans = fintpow(x, a); 
  //   return ans; 
  // } 
  
  ans = fintpow(x, a) * fintpow(y, b); 
  
  return ans; 
}

// TODO: what happens when there is no neighbors set in the neighborhood kernel? 

inline u_nbcount number_of_neighbors(const u_xyint i,
                                     const u_xyint j) {
  // When we do not wrap, the total number of neighbors depends on where we are
  // in the matrix, e.g. first column or last row has missing neighbors. In that
  // case, we need to correct for that. When wrapping is on, then the
  // number of neighbors is constant.
  //
  // The compiler will remove the if{} statements below at compile time if we are
  // wrapping, making this function return a constant, defined on the line below:
  u_nbcount nnb = NB_KERNEL_MAXN; 

  // If we use a fixed number of neighbors, then things are fine with just the 
  // compile time-value. We can assume this is what we want as an approximation, even 
  // if we don't wrap. In this case, cells on the edges will have a lower chance of
  // switching, but this may be an approximation we are willing to make for better
  // performance.
  // 
  // If we don't assume fixed number of neighbors, then we need to count them for 
  // real. This is what we do below. Note that if we wraparound, then by definition 
  // the numer of neighbors is fixed (but FIXED_NB is not necessarily true). 
  if ( ! FIXED_NB && ! WRAP ) {
    nnb = 0; 
    for ( s_xyint o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) { 
      for ( s_xyint o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) { 
        // If we don't wrap, then we need to add a bound check to make sure 
        // we are in the matrix
        const s_xyint i_target = i + o_r; 
        const s_xyint j_target = j + o_c; 
        if ( i_target >= 0 && j_target >= 0 && i_target < NR && j_target < NC ) { 
          nnb += NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        }
      }
    }
  }
  
  // Rcpp::Rcout << "i: " << (int) i << " j: " << (int) j << " nbs: " << (int) nnb << "\n"; 
  return nnb;
}

// Build the q vector for a given cell
void inline count_qs(u_state this_qs[NS], 
                     u_xyint i, 
                     u_xyint j, 
                     const u_state m[NR][NC]) { 
  
  memset(this_qs, 0, sizeof(u_state) * NS); 
  
  // We loop over the required offsets
  for ( s_xyint o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) { 
    for ( s_xyint o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) { 
      
      if ( WRAP ) { 
        u_state state = m[(NR + i + o_r) % NR][(NC + j + o_c) % NC]; 
        this_qs[state] += 
          NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
      
      } else { 
        // If we don't wrap, then we need to add a bound check to make sure 
        // we are in the matrix
        const s_xyint i_target = i + o_r; 
        const s_xyint j_target = j + o_c; 
        if ( i_target >= 0 && j_target >= 0 && i_target < NR && j_target < NC ) { 
          u_state state = m[i_target][j_target]; 
          this_qs[state] += 
            NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        }
      }
    }
  }
    // if ( i == 0 && j == 0 ) { 
    //   for ( int k=0; k<NS; k++ ) { 
    //     Rcpp::Rcout << "this_qs[" << (int) k << "]: " << (int) this_qs[k] << "\n";
    //   }
    // }
  
}

// This function will set the "probability lines", ie fill in the array that contains
// for each neighborhood configuration, where to pick its probability of transitions
void initialize_prob_line(u_pline prob_lines[NR][NC],
                          const u_state m[NR][NC]) {
  
  for ( u_xyint i=0; i<NR; i++ ) {
    for ( u_xyint j=0; j<NC; j++ ) {
      
      u_nbcount this_qs[NS]; 
      count_qs(this_qs, i, j, m); 

      // Get line in pre-computed transition probability table
      u_pline line = 0;
      for ( u_state k=0; k<NS; k++ ) {
        // Rcpp::Rcout << "this_qs[" << (int) k << "]: " << (int) this_qs[k] << "\n"; 
        line = line * (1+MAX_NB) + this_qs[k];
      }
      line -= 1;

      // If we have constant number of neighbors, then all_qs only contains the
      // values at each MAX_NB values, so we need to divide by MAX_NB here to fall on the
      // right line in the table of probabilities.
      if ( WRAP ) {
        line = (line+1) / MAX_NB - 1;
      }

      prob_lines[i][j] = line;
      // Rcpp::Rcout << "line: " << (int) line << " pline: " << (int) prob_lines[i][j] << "\n"; 
    }
  }
  
}

// Initialize the count of local densities for each cell in the landscape
void init_local_densities(u_nbcount qs[NR][NC][NS],
                          const u_state m[NR][NC]) {

  for ( uword i=0; i<NR; i++ ) {
    for ( uword j=0; j<NC; j++ ) { 
      count_qs(qs[i][j], i, j, m); 
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
  
  // We loop over the required offsets
  for ( s_xyint o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) { 
    for ( s_xyint o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) { 
        
      if ( WRAP ) { 
        qs[(NR + i + o_r) % NR][(NC + j + o_c) % NC][to] += 
          NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        qs[(NR + i + o_r) % NR][(NC + j + o_c) % NC][from] -= 
          NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        
      } else { 
        // If we don't wrap, then we need to add a bound check to make sure 
        // we are in the matrix
        const s_xyint i_target = i + o_r; 
        const s_xyint j_target = j + o_c; 
        if ( i_target >= 0 && j_target >= 0 && i_target < NR & j_target < NC ) { 
          qs[i_target][j_target][to] += 
            NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
          qs[i_target][j_target][from] -= 
            NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        }
        
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
  // NOTE: we use signed integers here because the pline adjustment can be negative (we 
  // go up in the table). Note that we convert *before* we carry out the difference, 
  // otherwise the integers wrap around 
  s_pline adj = (s_pline) intpow(MAX_NB+1, (NS-1) - to) - 
                  (s_pline) intpow(MAX_NB+1, (NS-1) - from);
  
  // If we have a fixed number of neighbors, then we discard the combinations that 
  // do not sum to this number of neighbors. In the table above for example, this is 
  // not done so you have configurations with, for example, 3 neighbors, that do not 
  // correspond to something that occurs in the grid. This effectively removes every 
  // line but those that fall every 'nb' (the number of neighbors), so we need to divide
  // the adjustment factor accordingly (i.e. instead of going nb-by-nb, we go one by
  // one). 
  adj = WRAP ? adj / ( (s_pline) MAX_NB ) : adj;
  
  // We loop over the required offsets
  for ( s_xyint o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) { 
    for ( s_xyint o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) { 
        
      if ( WRAP ) { 
        pline[(NR + i + o_r) % NR][(NC + j + o_c) % NC] += 
          adj * NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
          
      } else { 
        // If we don't wrap, then we need to add a bound check to make sure 
        // we are in the matrix
        const s_xyint i_changed = i + o_r; 
        const s_xyint j_changed = j + o_c; 
        if ( i_changed > 0 && j_changed > 0 && i_changed < NR & j_changed < NC ) { 
          pline[i_changed][j_changed] += 
            adj * NB_KERNEL[o_r + kernel_semiheight][o_c + kernel_semiwidth]; 
        }
        
      }
    }
  }
  
}

// Compute probability of transitions to other states, and store in tprob_line
void compute_rate(pfloat tprob_line[NS],
                  const u_nbcount qs[NS],
                  const u_pscount ps[NS],
                  const u_nbcount & total_nb,
                  const u_state & from) {
  
  // Set all probs to zero
  memset(tprob_line, 0, sizeof(pfloat)*NS);
  
  double total = 0.0; 
  u_state to = 0; 
  while ( to < NS ) { 
    
    // Check if we will ever transition into the state. This assumes
    // probability is set to zero above, otherwise tprob_line contains undefined or
    // wrong values. This adds a branch but it is usally worth it, especially with 
    // very sparse models. 
    if ( ! transition_matrix[from][to] ) {
      to++; 
      continue;
    }

    // double total = 0.0. 
    // pfloat total = 0.0; 
    s_pline kstart, kend;

    // constant component
    if ( BETA_0_NROW > 0 ) {
      kstart = betas_index[from][to][_beta_0_start];
      kend   = betas_index[from][to][_beta_0_end];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        total += coef_tab_flts[k];
      }
    }

    // q component
    if ( BETA_Q_NROW > 0 && total_nb > 0 ) {
      kstart = betas_index[from][to][_beta_q_start];
      kend = betas_index[from][to][_beta_q_end];
      if ( kstart >= 0 ) {
        // qthis holds the current point at which to extract the value of the coefficient
        // for q. When a transition depends on multiple neighbor expressions, i.e.
        // when P(a -> b) = f(q_1) + f(q_2), then we need to sum both the
        // corresponding coefficient for f(q_1), but also the one for f(q_2).
        // We know they are spaced xpoints apart from each other in the
        // coefficient table, so we just sum values every xpoints
        uword n_coef_to_sum = (kend - kstart + 1) / XPOINTS;
        
        // qpointn is a number used to convert the number of neighbors into 
        // the point at which the value of f(q) is stored in the coefficient tables. 
        // all_qs[ ][NS] holds the total number of neighbors, and is passed as 
        // total_nb. However, here fixed_nb is constexpr so this if() should be 
        // optimized away when we use a fixed number of neighbors as qpointn
        // is then known at compile time. 
        const u_nbcount qpointn = ( XPOINTS - 1 ) / ( FIXED_NB ? MAX_NB : total_nb ); 
          
        for ( uword coef=0; coef<n_coef_to_sum; coef++) {
          uword qthis = qs[ coef_tab_ints[kstart + XPOINTS*coef][_state_1] ] * qpointn;
          total += coef_tab_flts[kstart + qthis + XPOINTS*coef];
        }
      }
    }

    // pp
    if ( BETA_PP_NROW > 0 ) {
      kstart = betas_index[from][to][_beta_pp_start];
      kend   = betas_index[from][to][_beta_pp_end];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat p1 = (pfloat) ps[ coef_tab_ints[k][_state_1] ] / FLOAT_NCELLS;
        const pfloat p2 = (pfloat) ps[ coef_tab_ints[k][_state_2] ] / FLOAT_NCELLS;
        
        total += coef_tab_flts[k] * 
          xy_ab_product(p1, coef_tab_ints[k][_expo_1], 
                        p2, coef_tab_ints[k][_expo_2]);
      }
    }

    // qq
    if ( BETA_QQ_NROW > 0 && total_nb > 0 ) {
      kstart = betas_index[from][to][_beta_qq_start];
      kend   = betas_index[from][to][_beta_qq_end];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat q1 = (pfloat) qs[ coef_tab_ints[k][_state_1] ] / total_nb;
        const pfloat q2 = (pfloat) qs[ coef_tab_ints[k][_state_2] ] / total_nb;
        
        total += coef_tab_flts[k] * 
          xy_ab_product(q1, coef_tab_ints[k][_expo_1], 
                        q2, coef_tab_ints[k][_expo_2]);; 
      }
    }

    // pq
    if ( BETA_PQ_NROW > 0 && total_nb > 0 ) {
      kstart = betas_index[from][to][_beta_pq_start];
      kend   = betas_index[from][to][_beta_pq_end];
      for ( s_pline k=kstart; k<=kend; k++ ) {
        const pfloat p1 = (pfloat) ps[ coef_tab_ints[k][_state_1] ] / FLOAT_NCELLS;
        const pfloat q1 = (pfloat) qs[ coef_tab_ints[k][_state_2] ] / total_nb;
        
        total += coef_tab_flts[k] * 
          xy_ab_product(p1, coef_tab_ints[k][_expo_1], 
                        q1, coef_tab_ints[k][_expo_2]);; 
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
    // while we used to do it separately below. 
    tprob_line[to] = total; 
    
    to++; 
  }

  // Compute cumsum of probabilities. This is already by using a while loop above, 
  // instead of computing the rate of transition to every state, then computing the
  // cumsum separately. 
  // for ( u_state k=1; k<NS; k++ ) {
  //  tprob_line[k] += tprob_line[k-1];
  // }

}


// These memory slots should be contiguous, but maybe the compiler does this already
constexpr u_nbcount empty_vec_qs[NS] = {0}; 
constexpr u_pscount empty_vec_ps[NS] = {0}; 
static const u_nbcount *prev_qs = empty_vec_qs; 
static const u_pscount *prev_ps = empty_vec_ps; 
static u_nbcount prev_total_nb = 0; 
static u_nbcount prev_from = 0; 

// Same as compute_rate, but the last call is memoised, i.e. we will not recompute the 
// probabilities if arguments are the same by updating tprob_line
void compute_rate_memo(pfloat tprob_line[NS],
                       const u_nbcount qs[NS],
                       const u_pscount ps[NS],
                       const u_nbcount & total_nb,
                       const u_state & from) {
  
  // Test if all arguments are the same
  if ( prev_from == from && 
       ( memcmp(qs, prev_qs, sizeof(u_nbcount)*NS) == 0 ) && 
       ( memcmp(ps, prev_ps, sizeof(u_nbcount)*NS) == 0 ) && 
       prev_total_nb == total_nb ) { 
    // Do nothing and return, tprob_line holds the right thing already
    return; 
  }
  
  prev_qs = qs; 
  prev_ps = ps; 
  prev_total_nb = total_nb; 
  prev_from = from; 
  
  compute_rate(tprob_line, qs, ps, total_nb, from);
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
  console_callback(iter, ps_arma, FLOAT_NCELLS);
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
  arma::Col<arma::uword> ps_arma(NS);
  for ( u_state k=0; k<NS; k++ ) {
    ps_arma(k) = (arma::uword) ps[k];
  }
  
  cover_callback(t, ps_arma, FLOAT_NCELLS);
}
