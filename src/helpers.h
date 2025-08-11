//
// Some helpers used in different functions
//

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

// Define some column names for clarity
#define _from 0
#define _to 1
#define _state_1 2
#define _state_2 3
#define _expo_1 4
#define _expo_2 5
#define _qs 3 // only in beta_q, so no overlap with _state_2 above
#define _coef 0

inline void get_local_densities_column(arma::Mat<arma::uword>& qs,
                                       const arma::Mat<unsigned short>& m,
                                       const arma::uword& j,
                                       const bool& wrap,
                                       const arma::Mat<unsigned short>& kernel) {

  qs.fill(0);
  arma::uword nc = m.n_cols;
  arma::uword nr = m.n_rows;

  const arma::sword kernel_semiheight = ( kernel.n_rows - 1 ) / 2;
  const arma::sword kernel_semiwidth  = ( kernel.n_cols - 1 ) / 2;

  for ( arma::uword i = 0; i < nr; i++ ) {

    // We loop over the required offsets
    for ( arma::sword o_r = - kernel_semiheight; o_r <= kernel_semiheight; o_r++ ) {
      for ( arma::sword o_c = - kernel_semiwidth; o_c <= kernel_semiwidth; o_c++ ) {
        // Rcpp::Rcout << "o_r: " << o_r << " o_c: " << o_c << "\n";

        if ( wrap ) {
          const ushort state = m( (nr + i + o_r) % nr,
                                  (nc + j + o_c) % nc );
          qs(i, state) += kernel(o_r + kernel_semiheight,
                                 o_c + kernel_semiwidth);

        } else {
          // If we don't wrap, then we need to add a bound check to make sure
          // we are in the matrix
          const arma::uword i_target = i + o_r;
          const arma::uword j_target = j + o_c;
          if ( i_target >= 0 && j_target >= 0 && i_target < nr && j_target < nc ) {
            const ushort state = m(i_target, j_target);
            qs(i, state) += kernel(o_r + kernel_semiheight,
                                   o_c + kernel_semiwidth);
          }
        }
      }
    }

  }

}
