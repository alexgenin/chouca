
#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "helpers.h"

// Some typedefs for better legibility. ushort typically contains an integer
// for the state number, which is often way lower than its maximum value (2^16). We
// cannot use a shorter type here (e.g. unsigned char) because matrices of uchars are not
// supported by armadillo. We use c-style array in the 'compiled' backend for better
// performance which do use those, but not here because it adds complexity to the code.
//
// Note that in function signatures we can't use this typedef because they are copied
// by Rcpp in the RcppExports.cpp file, but without typedef and this makes the
// compilation fail on Windows.
typedef unsigned short ushort;

// This is a function that returns the local state counts to R, for the full column of
// a matrix. j must be indexed the R-way (1-indexing)
//[[Rcpp::export]]
arma::Mat<arma::uword> local_dens_col(const arma::Mat<unsigned short>& m,
                                      const arma::uword& nstates,
                                      const arma::uword& j,
                                      const bool& wrap,
                                      const arma::Mat<unsigned short>& kernel) {
  arma::Mat<arma::uword> newq(m.n_rows, nstates);
  newq.fill(0);
  // m, i and j must be adjusted to take into account the way things are stored in R
  get_local_densities_column(newq, m, j - 1, wrap, kernel);

  return(newq);
}

// Get the transition probabilities as a 3D matrix
//[[Rcpp::export]]
arma::cube get_transition_probas_cpp(arma::Mat<ushort>& mat,
                                     Rcpp::List& ctrl) {

  arma::uword nr = mat.n_rows;
  arma::uword nc = mat.n_cols;
  arma::uword ncells = nr * nc;

  // Unpack control list
  const ushort ns = ctrl["nstates"];
  const arma::Mat<ushort> kernel = ctrl["kernel"];
  const bool wrap = ctrl["wrap"];
  // Number of samples for qs
  const ushort xpoints = ctrl["xpoints"];

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

  // Initialize vector with counts of cells in each state in the landscape
  arma::Col<int> ps(ns);
  ps.fill(0);

  // Count the number of cells in each state, and store that in our vector
  for (arma::uword j = 0; j < nc; j++) {
    for (arma::uword i = 0; i < nr; i++) {
      ps( mat(i, j) )++;
    }
  }

  // Define output array
  arma::cube ptrans(mat.n_rows, mat.n_cols, ns);
  ptrans.fill(0);

  for ( arma::uword j=0; j<nc; j++ ) {

    arma::Mat<arma::uword> qs(nr, ns);
    qs.fill(0);

    get_local_densities_column(qs, mat, j, wrap, kernel);

    for ( arma::uword i=0; i<nr; i++ ) {

      ushort from = mat(i, j);

      // Normalized local densities to proportions
      arma::uword total_nb = accu(qs.row(i));
      // Rcpp::Rcout << qs.row(i) << "\n";

      // Factor to convert the number of neighbors into the point at which the
      // dependency on q is to be found.
      // TODO: this will crash when there is no neighbors
      arma::uword qpointn = (xpoints - 1) / total_nb;

      // Compute probability transitions
      for ( ushort to = 0; to < ns; to++ ) {

        // constant component
        for (arma::uword k = 0; k < beta_0_index.n_rows; k++) {
          ptrans(i, j, to) +=
            (beta_0_index(k, _from) == from) *
            (beta_0_index(k, _to) == to) *
            beta_0_vals(k);
        }

        // f(q) component
        for (arma::uword k = 0; k < beta_q_index.n_rows; k++) {

            // Lookup which point in the qs function we need to use for the
            // current neighbor situation.
            arma::uword qthis = qs(i, beta_q_index(k, _state_1)) * qpointn;

            ptrans(i, j, to) +=
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

            ptrans(i, j, to) +=
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

            ptrans(i, j, to) +=
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

            ptrans(i, j, to) +=
              (beta_pq_index(k, _from) == from) * // zero when other transition than
              (beta_pq_index(k, _to) == to) *     // the one we are dealing with now
              beta_pq_vals(k, _coef) *
              pow(p1, beta_pq_index(k, _expo_1)) *
              pow(q1, beta_pq_index(k, _expo_2));
        }

        // ptrans(i, j, to) = ptrans(i, j, to) < 0.0 ? 0.0 : ptrans(i, j, to);
      }

    }
  }

    return ptrans;
}

// This is a transitional function to speed things up, but it will go away
//[[Rcpp::export]]
arma::mat transition_ll(const arma::Mat<ushort>& from_mat,
                     const arma::Mat<ushort>& to_mat,
                     const arma::cube& tprobs) {
  // from/to mat in R integer form

  arma::uword nr = from_mat.n_rows;
  arma::uword nc = from_mat.n_cols;

  arma::mat ll(nr, nc);
  ll.fill(0);
  //   mod <- update_mod(theta, mod)
  //
  // nextp <- next_state_probs(mod, from_mat)
  // ps <- matrix(0, nrow = nrow(from_mat), ncol = ncol(from_mat))
  //
  // for ( i in seq.int(nrow(to_mat)) ) {
  //   for ( j in seq.int(ncol(to_mat)) ) {
  //     ps[i, j] <-
  //       (to_mat[i,j] == from_mat[i,j]) * ( 1 - sum(nextp[i, j, ]) ) + # no change
  //       (to_mat[i,j] != from_mat[i,j]) * ( nextp[i, j, to_mat[i, j]] ) # change
  //   }
  // }
  //
  // sum(log(ps))

  for ( arma::uword i=0; i<nr; i++ ) {
    for ( arma::uword j=0; j<nc; j++ ) {
      ll(i, j) = ( from_mat(i, j) == to_mat(i, j) ) ?
                   1 - accu(tprobs.tube(i, j)) :
                   tprobs(i, j, to_mat(i, j) - 1); // take into account R form
    }
  }

  return(ll);
}
