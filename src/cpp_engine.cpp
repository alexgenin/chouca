
#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
    
// Some typedefs for better legibility. ushort typically contains an integer 
// for the state number, which is often way lower than its maximum value (2^16). We 
// cannot use a shorter type here (e.g. unsigned char) because matrices of uchars are not 
// supported by armadillo. We use c-style array in the 'compiled' backend for better 
// performance which do use those, but not here because it adds complexity to the code. 
typedef unsigned short ushort;

// Define some column names for clarity
#define _from 0
#define _to 1
#define _state_1 2
#define _state_2 3
#define _expo_1 4
#define _expo_2 5
#define _qs 3 // only in beta_q, so no overlap with _state_2 above
#define _coef 0

inline void get_local_densities(arma::Col<arma::uword>& qs,
                                const arma::Mat<unsigned short>& m,
                                const arma::uword i,
                                const arma::uword j,
                                const bool wrap,
                                const bool use_8_nb) {
  qs.fill(0);
  arma::uword nc = m.n_cols;
  arma::uword nr = m.n_rows;
  
  // Get neighbors to the left
  if (wrap) {
    
    ushort state_left = m(i, (nc + j - 1) % nc);
    qs(state_left)++; // left
    
    ushort state_right = m(i, (nc + j + 1) % nc);
    qs(state_right)++; // right
    
    ushort state_up = m((nr + i - 1) % nr, j);
    qs(state_up)++; // up
    
    ushort state_down = m((nr + i + 1) % nr, j);
    qs(state_down)++; // down
    
    if (use_8_nb) {
      
      ushort state_upleft = m((nr + i - 1) % nr, (nc + j - 1) % nc);
      qs(state_upleft)++; // upleft
      
      ushort state_upright = m((nr + i - 1) % nr, (nc + j + 1) % nc);
      qs(state_upright)++; // upright
      
      ushort state_downleft = m((nr + i + 1) % nr, (nc + j - 1) % nc);
      qs(state_downleft)++; // downleft
      
      ushort state_downright = m((nr + i + 1) % nr, (nc + j + 1) % nc);
      qs(state_downright)++; // downright
    }
      
  } else {
    
    if (j > 0) {
      ushort state_left = m(i, j - 1);
      qs(state_left)++; // left
    }
    if (j < (nc - 1)) {
      ushort state_right = m(i, j + 1);
      qs(state_right)++; // right
    }
    if (i > 0) {
      ushort state_up = m(i - 1, j);
      qs(state_up)++; // up
    }
    if (i < (nr - 1)) {
      ushort state_down = m(i + 1, j);
      qs(state_down)++; // down
    }
    
    if (use_8_nb) {
      if (i > 0 && j > 0) {
        ushort state_upleft = m(i - 1, j - 1);
        qs(state_upleft)++; // upleft
      }
      if (i > 0 && j < (nc - 1)) {
        ushort state_upright = m(i - 1, j + 1);
        qs(state_upright)++; // upright
      }
      if (i < (nr - 1) && j > 0) {
        ushort state_downleft = m(i + 1, j - 1);
        qs(state_downleft)++; // downleft
      }
      if (i < (nr - 1) && j < (nc - 1)) {
        ushort state_downright = m(i + 1, j + 1);
        qs(state_downright)++; // downright
      }
    }
  }
}

// This is a function that returns the local state counts to R. i and j must be indexed
// the R-way (1-indexing)
// [[Rcpp::export]]
arma::Col<arma::uword> local_dens(const arma::Mat<unsigned short> m,
                                  const arma::uword nstates,
                                  const arma::uword i,
                                  const arma::uword j,
                                  const bool wrap,
                                  const bool use_8_nb) {
  arma::Col<arma::uword> newq(nstates);
  newq.fill(0);
  // m, i and j must be adjusted to take into account the way things are stored in R
  get_local_densities(newq, m, i - 1, j - 1, wrap, use_8_nb);
  
  return (newq);
}

inline void get_local_densities_column(arma::Mat<arma::uword>& qs,
                                       const arma::Mat<unsigned short>& m,
                                       const arma::uword j,
                                       const bool wrap,
                                       const bool use_8_nb) {
                                      
  qs.fill(0);
  arma::uword nc = m.n_cols;
  arma::uword nr = m.n_rows;
  
  if (wrap) {
    for (arma::uword i = 0; i < nr; i++) {
      // left column
      ushort state_left = m(i, (nc + j - 1) % nc);
      qs(i, state_left)++;
      
      // right column
      ushort state_right = m(i, (nc + j + 1) % nc);
      qs(i, state_right)++;
      
      ushort state_up = m((nr + i - 1) % nr, j);
      qs(i, state_up)++; // up
      
      ushort state_down = m((nr + i + 1) % nr, j);
      qs(i, state_down)++; // down
    }
      
    if (use_8_nb) {
      for (arma::uword i = 0; i < nr; i++) {
        ushort state_upleft = m((nr + i - 1) % nr, (nc + j - 1) % nc);
        qs(i, state_upleft)++; // upleft
        
        ushort state_upright = m((nr + i - 1) % nr, (nc + j + 1) % nc);
        qs(i, state_upright)++; // upright
        
        ushort state_downleft = m((nr + i + 1) % nr, (nc + j - 1) % nc);
        qs(i, state_downleft)++; // downleft
        
        ushort state_downright = m((nr + i + 1) % nr, (nc + j + 1) % nc);
        qs(i, state_downright)++; // downright
      }
    }
      
  } else { // no wrap
    
    for (arma::uword i = 0; i < nr; i++) {
      if (j > 0) {
        ushort state_left = m(i, j - 1);
        qs(i, state_left)++; // left
      }
      if (j < (nc - 1)) {
        ushort state_right = m(i, j + 1);
        qs(i, state_right)++; // right
      }
      if (i > 0) {
        ushort state_up = m(i - 1, j);
        qs(i, state_up)++; // up
      }
      if (i < (nr - 1)) {
        ushort state_down = m(i + 1, j);
        qs(i, state_down)++; // down
      }
      
      if (use_8_nb) {
        if (i > 0 && j > 0) {
            ushort state_upleft = m(i - 1, j - 1);
            qs(i, state_upleft)++; // upleft
        }
        if (i > 0 && j < (nc - 1)) {
            ushort state_upright = m(i - 1, j + 1);
            qs(i, state_upright)++; // upright
        }
        if (i < (nr - 1) && j > 0) {
            ushort state_downleft = m(i + 1, j - 1);
            qs(i, state_downleft)++; // downleft
        }
        if (i < (nr - 1) && j < (nc - 1)) {
            ushort state_downright = m(i + 1, j + 1);
            qs(i, state_downright)++; // downright
        }
      }
    }
  }
}

// This is a function that returns the local state counts to R, for the full column of
// a matrix. j must be indexed the R-way (1-indexing)
//[[Rcpp::export]]
arma::Mat<arma::uword> local_dens_col(const arma::Mat<unsigned short> m,
                                      const arma::uword nstates,
                                      const arma::uword j,
                                      const bool wrap,
                                      const bool use_8_nb) {
  arma::Mat<arma::uword> newq(m.n_rows, nstates);
  newq.fill(0);
  // m, i and j must be adjusted to take into account the way things are stored in R
  get_local_densities_column(newq, m, j - 1, wrap, use_8_nb);
  
  return (newq);
}

//
// This file contains the main c++ engine to run CAs
//
// [[Rcpp::export]]
void camodel_cpp_engine(const Rcpp::List ctrl) {
  
  // Unpack control list
  const bool wrap = ctrl["wrap"];
  const arma::Mat<ushort> init = ctrl["init"];
  const ushort ns = ctrl["nstates"];
  const ushort nbs = ctrl["neighbors"]; // Needed to make sure conversion from list is OK
  const bool use_8_nb = nbs == 8 ? true : false;
  const arma::uword substeps = ctrl["substeps"];
  const arma::Col<arma::uword> times = ctrl["times"];
  
  const arma::Mat<ushort> transition_mat = ctrl["transition_mat"];
  
  // Number of samples for qs
  const ushort xpoints = ctrl["xpoints"];
  
  // Extract callbacks which will be used to export data from the simulation 
  const arma::uword console_callback_every = ctrl["console_output_every"];
  const bool console_callback_active = console_callback_every > 0;
  Rcpp::Function console_callback = ctrl["console_callback"];
  
  const arma::uword cover_callback_every = ctrl["save_covers_every"];
  const bool cover_callback_active = cover_callback_every > 0;
  Rcpp::Function cover_callback = ctrl["cover_callback"];
  
  const arma::uword snapshot_callback_every = ctrl["save_snapshots_every"];
  const bool snapshot_callback_active = snapshot_callback_every > 0;
  Rcpp::Function snapshot_callback = ctrl["snapshot_callback"];
  
  const arma::uword custom_callback_every = ctrl["custom_output_every"];
  const bool custom_callback_active = custom_callback_every > 0;
  Rcpp::Function custom_callback = ctrl["custom_callback"];
  
  // Extract things from list. Here we need to split arrays between ints and doubles, as 
  // they cannot be hold by the same array in c++ (but in R it is somewhat possible)
  const arma::Mat<ushort> beta_0_index = ctrl["beta_0_index"];
  const arma::Mat<double> beta_0_vals = ctrl["beta_0_vals"];
  const arma::Mat<ushort> beta_q_index = ctrl["beta_q_index"];
  const arma::Mat<double> beta_q_vals = ctrl["beta_q_vals"];
  const arma::Mat<ushort> beta_pp_index = ctrl["beta_pp_index"];
  const arma::Mat<double> beta_pp_vals = ctrl["beta_pp_vals"];
  const arma::Mat<ushort> beta_pq_index = ctrl["beta_pq_index"];
  const arma::Mat<double> beta_pq_vals = ctrl["beta_pq_vals"];
  const arma::Mat<ushort> beta_qq_index = ctrl["beta_qq_index"];
  const arma::Mat<double> beta_qq_vals = ctrl["beta_qq_vals"];
  
  // Some constants about the simulation
  const arma::uword nr = init.n_rows;
  const arma::uword nc = init.n_cols;
  const double ncells = (double)nr * (double)nc;
  
  // We need to versions of the matrix 
  arma::Mat<ushort> omat = init;
  arma::Mat<ushort> nmat = init;
  
  // Initialize vector with counts of cells in each state in the landscape (used to
  // compute global densities). It is much easier to hold integers that we can 
  // increment/decrement than working on proportions here. 
  arma::Col<int> ps(ns);
  arma::Col<int> delta_ps(ns);
  ps.fill(0);
  delta_ps.fill(0);
  
  // Count the number of cells in each state, and store that in our vector
  for (arma::uword j = 0; j < nc; j++) {
    for (arma::uword i = 0; i < nr; i++) {
      ps(init(i, j))++;
    }
  }
  
  // Allocate some things we will reuse later
  arma::Mat<arma::uword> qs(nr, ns); // number of neighbors in each state for a column
  // arma::Col<double> qs_prop(ns);     // proportion of neighbors in each state
  arma::Col<double> ptrans(ns);      // vector of transition probas
  
  // Allocate essentials things to keep track of time. We use an uword here for time, 
  // which means that the max time step is 2^64 on 64-bits systems (~10e19)
  arma::uword current_t = 0;
  arma::uword last_t = times(times.n_elem - 1);
  arma::uword export_n = 0;
  arma::uword next_export_t = times(export_n);
  
  while (current_t <= last_t) {
    
    // Call callbacks
    if (console_callback_active && current_t % console_callback_every == 0) {
      console_callback(current_t, ps, ncells);
    }
    
    if (next_export_t == current_t) {
        
      if (cover_callback_active && export_n % cover_callback_every == 0) {
        cover_callback(current_t, ps, ncells);
      }
      
      if (snapshot_callback_active && export_n % snapshot_callback_every == 0) {
        snapshot_callback(current_t, omat);
      }
      
      if (custom_callback_active && export_n % custom_callback_every == 0) {
        custom_callback(current_t, omat);
      }
      
      export_n++;
      // Avoids an OOB read when on last iteration that is picked up by ASAN
      if ( export_n < times.n_elem ) {
        next_export_t = times(export_n);
      }
    }
    
    for (arma::uword s = 0; s < substeps; s++) {
        
      for (arma::uword j = 0; j < nc; j++) {
        
        // Get local state counts for the column. Getting it for a column at once is
        // slightly more cache-friendly, and a simple optimization worth going for. 
        get_local_densities_column(qs, omat, j, wrap, use_8_nb);
        
        for (arma::uword i = 0; i < nr; i++) {
          
          // Normalized local densities to proportions
          arma::uword total_nb = accu(qs.row(i));
          
          // Factor to convert the number of neighbors into the point at which the
          // dependency on q is to be found.
          arma::uword qpointn = (xpoints - 1) / total_nb;
          
          // Get current state
          ushort from = omat(i, j);
          
          // Compute probability transitions
          for (ushort to = 0; to < ns; to++) {
              
            // Init probability
            ptrans(to) = 0;
            
            if (transition_mat(from, to) == 0) {
              continue;
            }
            
            // constant component
            for (arma::uword k = 0; k < beta_0_index.n_rows; k++) {
              ptrans(to) += 
                (beta_0_index(k, _from) == from) * 
                (beta_0_index(k, _to) == to) * 
                beta_0_vals(k);
            }
            
            // f(q) component
            for (arma::uword k = 0; k < beta_q_index.n_rows; k++) {
                
                // Lookup which point in the qs function we need to use for the
                // current neighbor situation.
                arma::uword qthis = qs(i, beta_q_index(k, _state_1)) * qpointn;
                
                ptrans(to) += 
                  (beta_q_index(k, _from) == from) * 
                  (beta_q_index(k, _to) == to) *
                  // Given the observed local abundance of this state, which line in
                  // beta_q should be retained ?
                  (beta_q_index(k, _qs) == qthis) * 
                  beta_q_vals(k);
            }
            
            // pp
            for (arma::uword k = 0; k < beta_pp_index.n_rows; k++) {
                const double p1 = ps[beta_pp_index(k, _state_1)] / ncells;
                const double p2 = ps[beta_pp_index(k, _state_2)] / ncells;
                
                ptrans(to) += 
                  (beta_pp_index(k, _from) == from) * // zero when other transition than 
                  (beta_pp_index(k, _to) == to) *     // the one we are dealing with now
                  beta_pp_vals(k, _coef) * 
                  pow(p1, beta_pp_index(k, _expo_1)) * 
                  pow(p2, beta_pp_index(k, _expo_2));
            }
            
            // qq
            for (arma::uword k = 0; k < beta_qq_index.n_rows; k++) {
                const double q1 = (double) qs(i, beta_qq_index(k, _state_1)) / total_nb;
                const double q2 = (double) qs(i, beta_qq_index(k, _state_2)) / total_nb;
                
                ptrans(to) += 
                  (beta_qq_index(k, _from) == from) * // zero when other transition than 
                  (beta_qq_index(k, _to) == to) *     // the one we are dealing with now
                  beta_qq_vals(k, _coef) * 
                  pow(q1, beta_qq_index(k, _expo_1)) * 
                  pow(q2, beta_qq_index(k, _expo_2));
            }
            
            // pq
            for (arma::uword k = 0; k < beta_pq_index.n_rows; k++) {
                const double p1 = ps[beta_pq_index(k, _state_1)] / ncells;
                const double q1 = (double) qs(i, beta_pq_index(k, _state_2)) / total_nb;
                
                ptrans(to) += 
                  (beta_pq_index(k, _from) == from) * // zero when other transition than 
                  (beta_pq_index(k, _to) == to) *     // the one we are dealing with now
                  beta_pq_vals(k, _coef) * 
                  pow(p1, beta_pq_index(k, _expo_1)) * 
                  pow(q1, beta_pq_index(k, _expo_2));
            }
            
            ptrans(to) = ptrans(to) < 0.0 ? 0.0 : ptrans(to);
            
            // NOTE: it seems like a good idea to cap probabilities to 1, but it
            // is really not. The fact that things are above 1 should be handled
            // by increasing the number of substeps, otherwise the rates of
            // probabilities become unbalanced between transitions. For example,
            // if we have P(s -> a) = 2 and P(s -> b) = 0.4, then by setting
            // P(s -> a) to 1.0, we reduce its relative probability compared to
            // P(s -> b) when using substeps >= 2.
            // ptrans(to) = ptrans(to) > 1.0 ? 1.0 : ptrans(to);
          }
          
          // Check if we actually transition. We compute the cumulative sum of 
          // probabilities, then draw a random number (here 'rn'), to see if the 
          // transition occurs or not. 
          // 
          // 0 |-----p0-------(p0+p1)------(p0+p1+p2)------| 1
          //               ^ p0 < rn < (p0+p1) => p1 wins
          //                                           ^ rn > everything => no transition
          //       ^ 0 < rn < p0 => p0 wins
          // Of course the sum of probabilities must be lower than one, otherwise we are
          // making an approximation and may never consider a given transition.
          // 
          // ptrans = cumsum(ptrans); // alternative code, but slower because it needs
          //                          // a copy of the vector
          for (ushort k = 1; k < ns; k++) {
              ptrans(k) += ptrans(k - 1);
          }
          
          ushort new_cell_state = from;
          double rn = Rf_runif(0, 1);
          for (signed short k = (ns - 1); k >= 0; k--) {
            new_cell_state = rn < ptrans(k) ? k : new_cell_state;
          }
          
          if (new_cell_state != from) {
            ushort old_cell_state = nmat(i, j);
            // delta contains the change in ps we need to apply before going to the 
            // next iteration
            delta_ps(new_cell_state)++;
            delta_ps(old_cell_state)--;
            nmat(i, j) = new_cell_state;
          }
        }
      }
      
      // Apply iteration changes in global densities and matrix state
      for (ushort k = 0; k < ns; k++) {
        ps(k) += delta_ps(k);
        delta_ps(k) = 0;
      }
      
      omat = nmat;
    }
    
    current_t++;
    R_CheckUserInterrupt();
  }
}
