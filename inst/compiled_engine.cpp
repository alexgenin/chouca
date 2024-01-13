//
// This is the source code for the 'compiled' chouca engine. Its fields will be
// replaced and the code compiled to run an SCA model.
//
//

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif


// This is required for portability 
#ifdef _OPENMP
# define USE_OMP __USE_OMP__
#else
# define USE_OMP 0
#endif

#if USE_OMP 
# include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// Some typedefs for better legibility
typedef unsigned char uchar;
typedef unsigned short ushort;

// These strings will be replaced by their values
constexpr uword nr = __NR__;
constexpr uword nc = __NC__;
constexpr uchar ns = __NS__;
constexpr bool wrap = __WRAP__;
constexpr bool use_8_nb = __USE_8_NB__;
constexpr uchar n_nb = use_8_nb ? 8 : 4;
constexpr bool fixed_nb = __FIXED_NEIGHBOR_NUMBER__;
constexpr uword substeps = __SUBSTEPS__;
constexpr uword xpoints = __XPOINTS__;
constexpr double ncells = nr * nc;
constexpr uword beta_0_nrow = __BETA_0_NROW__;
constexpr uword beta_q_nrow = __BETA_Q_NROW__;
constexpr uword beta_pp_nrow = __BETA_PP_NROW__;
constexpr uword beta_qq_nrow = __BETA_QQ_NROW__;
constexpr uword beta_pq_nrow = __BETA_PQ_NROW__;
constexpr uword all_qs_nrow = __ALL_QS_NROW__;
constexpr uword cores = __CORES__;

// Transition matrix (true when transition can be done, false when not)
constexpr bool transition_matrix[ns][ns] = __TMATRIX_ARRAY__;
// 5*2 because we have 5 packed tables of coefficients, each of which needing 2
// (where to start, and where to stop for this transition)
constexpr sword betas_index[5 * 2][ns][ns] = __BETAS_INDEX__;

// Whether we want to precompute probabilities or not
#define PRECOMPUTE_TRANS_PROBAS __PRECOMP_PROBA_VALUE__

// The maximum number of neighbors
constexpr arma::uword max_nb = use_8_nb ? 8 : 4;

// Declare rows of coef_tab
constexpr arma::uword coef_tab_nrow = __COEF_TAB_NROW__;

// Coef table
static arma::Mat<ushort> coef_tab_ints(coef_tab_nrow, 6);
static arma::Mat<double> coef_tab_dbls(coef_tab_nrow, 1);

// Include functions and type declarations (this points to 'common.h')
#include "__COMMON_HEADER__"

// Compute transition probabilities between all possible qs states
inline void precompute_transition_probabilites(double tprobs[all_qs_nrow][ns][ns],
                                               const unsigned char all_qs[all_qs_nrow][ns + 1],
                                               const arma::uword ps[ns]) {

  // 
  // This function depends on the matrix all_qs, which contains all possible
  // neighborhood configurations, with an added column at the end with the total
  // number of neighbors (a fixed number under some set of options, e.g. it is always
  // 4 or 8 when we use a toric space as all cells have the same number of neighbors).
  //
  // For example, all_qs looks like this for a model with 3 states and 4 maximum
  // neighbors:
  //       [,1] [,2] [,3] [,4]
  // [1,]    0    0    4    4
  // [2,]    0    1    3    4
  // [3,]    0    2    2    4
  // [4,]    0    3    1    4
  // [5,]    0    4    0    4
  // [6,]    0    4    4    8
  // [7,]    1    0    3    4
  // [8,]    1    1    2    4
  // [9,]    1    2    1    4
  //
  // For all rows of all_qs above, i.e. for all possible neighbor configurations,
  // we need to compute at each iteration the probability of switching from one state to
  // another, which is stored in tprobs. Then once we know the neighbor configuration
  // of a cell, we just need to fetch it from tprobs instead of recomputing it from
  // scratch.
  // 
  // Note that when the neighborhood of a cell changes, we know immediately where to
  // find the new neighborhood configuration all_qs. For example, if we go from the
  // neighborhood (1, 0, 3) to (0, 0, 4), we know we new neighborhood config is 6 rows
  // above in all_qs (and thus in tprobs). This is because in all_qs,
  // the neighborhood configurations have a specific ordering. See adjust_nb_plines()
  // in common.h for more details.
  // 
  for (uword l = 0; l < all_qs_nrow; l++) {

    // Some combinations in all_qs will never be used because the number of
    // neighbors does not sum up to something we will ever encounter. See all_qs above,
    // where some neighborhood configuration sum to 8 (last column), but we use a
    // 4-way neighborhood in that example.
    // Skip those.
    //
    // NOTE: Ideally we would remove those values, and keep only lines in all_qs that
    // sum to valid combinations. This would make all_qs not grow as fast with the number
    // of neighbors and state, thus enabling memoization for more complex models.
    // However, this complexifies a lot the way of keeping
    // track of where to pick the transition probability for each neighbor combination. 
    // 
    // I [Alex] worked out a math formula somewhere but lost it, I just remember the
    // problem simplified into some variant of the bars and stars problem, see
    // https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics). 
    if ( all_qs[l][ns] > n_nb ) {
      continue;
    }

    for ( uchar from = 0; from < ns; from++ ) {

      uword qpointn_factorf; 
      if ( fixed_nb ) { 
        qpointn_factorf = (xpoints - 1) / n_nb;
      } else { 
        qpointn_factorf = (xpoints - 1) / all_qs[l][ns];
      }
      
      // Init probability
      compute_rate(tprobs[l][from],
                   all_qs[l], // qs for this line of all_qs
                   ps, // current ps
                   all_qs[l][ns], // total of neighbors
                   from); // from, to states
    }
  }
}

// [[Rcpp::export]]
void aaa__FPREFIX__camodel_compiled_engine(const arma::Mat<ushort> all_qs_arma,
                                           const Rcpp::List ctrl) {

  // Unpack control list
  const Mat<ushort> init = ctrl["init"]; // this is ushort because init is an arma mat
  const Col<uword> times = ctrl["times"];

  const uword console_callback_every = ctrl["console_output_every"];
  const bool console_callback_active = console_callback_every > 0;
  Rcpp::Function console_callback = ctrl["console_callback"];

  const uword cover_callback_every = ctrl["save_covers_every"];
  const bool cover_callback_active = cover_callback_every > 0;
  Rcpp::Function cover_callback = ctrl["cover_callback"];

  const uword snapshot_callback_every = ctrl["save_snapshots_every"];
  const bool snapshot_callback_active = snapshot_callback_every > 0;
  Rcpp::Function snapshot_callback = ctrl["snapshot_callback"];

  const uword custom_callback_every = ctrl["custom_output_every"];
  const bool custom_callback_active = custom_callback_every > 0;
  Rcpp::Function custom_callback = ctrl["custom_callback"];

  // Extract things from list. These arrays are declared as static above so accessible 
  // from other functions. 
  coef_tab_ints = Rcpp::as<arma::Mat<ushort>>(ctrl["coef_tab_ints"]);
  coef_tab_dbls = Rcpp::as<arma::Mat<double>>(ctrl["coef_tab_dbls"]);

  // Copy some things as c arrays. 
  // Note: we allocate omat/nmat on the heap since they can be big matrices and blow up
  // the size of the C stack beyond acceptable.
  auto old_mat = new uchar[nr][nc];
  auto new_mat = new uchar[nr][nc];
  for (uword i = 0; i < nr; i++) {
    for (uword j = 0; j < nc; j++) {
      old_mat[i][j] = (uchar)init(i, j);
    }
  }
  memcpy(new_mat, old_mat, sizeof(uchar) * nr * nc);

#if PRECOMPUTE_TRANS_PROBAS
  // Convert all_qs to c-style char array
  auto all_qs = new uchar[all_qs_nrow][ns + 1];
  for (uword i = 0; i < all_qs_nrow; i++) {
    for (uword k = 0; k < (ns + 1); k++) {
      all_qs[i][k] = (uchar) all_qs_arma(i, k);
    }
  }
#endif

  // Initialize vectors with counts of cells in each state in the landscape (used to
  // compute global densities)
  uword old_ps[ns];
  uword new_ps[ns];
  memset(old_ps, 0, sizeof(uword) * ns);
  memset(new_ps, 0, sizeof(uword) * ns);

  // Compute global counts of cells
  for (uword i = 0; i < nr; i++) {
    for (uword j = 0; j < nc; j++) {
      old_ps[old_mat[i][j]]++;
    }
  }
  memcpy(new_ps, old_ps, sizeof(uword) * ns);


  // When we do not use memoization, we maintain for each cell of the grid the state
  // of its neighbors (in old/new_qs). Everytime a cell changes state, the values of
  // those matrices are changed for all neighboring cells.
  // 
  // When we use memoization, we do not need to keep track of the local counts of
  // neighbors, we just need to adjust the index (row) in old/new_pline which contains
  // where in tprobs the probability of transition needs to be picked for a cell
  // with that specific neighborhood configuration. There is a relationship between
  // the type of neighborhood change and the number of rows to add or substract to find
  // the new probabilities of transition. This means that when we memoize probabilities we
  // never have to count neighbors, which is a large part of the speed boost.
  // See adjust_nb_plines() in common.h for more information.
  //
#if PRECOMPUTE_TRANS_PROBAS
  // Matrix holding probability line in the precomputed table
  auto old_pline = new uword[nr][nc];
  auto new_pline = new uword[nr][nc];
  initialize_prob_line(old_pline, old_mat);
  memcpy(new_pline, old_pline, sizeof(uword) * nr * nc);

  // Initialize table with precomputed probabilities
  auto tprobs = new double[all_qs_nrow][ns][ns];
  memset(tprobs, 0.0, sizeof(double) * all_qs_nrow * ns * ns);
#else
  // Create tables that hold local densities
  auto old_qs = new uchar[nr][nc][ns];
  auto new_qs = new uchar[nr][nc][ns];
  init_local_densities(old_qs, old_mat);
  memcpy(new_qs, old_qs, sizeof(uchar) * nr * nc * ns);
#endif

  // Initialize random number generator. All values in s must be overwritten with
  // 64 bits. We use doubles, but use memcpy to overwrite the raw bits into s which
  // is defined as integers (uint64_t)
  for (uword i = 0; i < 4; i++) {
    for (uword thread = 0; thread < cores; thread++) {
      // randu returns 64-bit double precision numbers
      double rn = arma::randn<double>();
      memcpy(&s[thread][i], &rn, sizeof(double));
    }
  }

  // Allocate some things we will reuse later. We always need to define this when using
  // multiple cores as otherwise the omp pragma will complain it's undefined when
  // using multiple threads.
#if (! PRECOMPUTE_TRANS_PROBAS) || USE_OMP
  double ptrans[ns];
#endif

  uword current_t = 0.0;
  uword last_t = times(times.n_elem - 1);
  uword export_n = 0;
  uword next_export_t = times(export_n);

  while (current_t <= last_t) {

    // Call callbacks
    if (console_callback_active && current_t % console_callback_every == 0) {
      console_callback_wrap(current_t, old_ps, console_callback);
    }

    if (next_export_t <= current_t) {
      if (cover_callback_active && export_n % cover_callback_every == 0) {
        cover_callback_wrap(current_t, old_ps, cover_callback);
      }

      if (snapshot_callback_active && export_n % snapshot_callback_every == 0) {
        snapshot_callback_wrap(current_t, old_mat, snapshot_callback);
      }

      if (custom_callback_active && export_n % custom_callback_every == 0) {
        custom_callback_wrap(current_t, old_mat, custom_callback);
      }

      export_n++;
      // Avoids an OOB read when on last iteration
      if ( export_n < times.n_elem ) {
        next_export_t = times(export_n);
      }
    }

    for (uword s = 0; s < substeps; s++) {

#if PRECOMPUTE_TRANS_PROBAS
      // Note that we do not need to pass the arrays with model coefficients to the
      // function, because they are known at compile time
      precompute_transition_probabilites(tprobs, all_qs, old_ps);
#endif

#if USE_OMP
      // This parallel block has an implicit barrier at the end, so each thread will
      // process a line and wait for the others to finish before switching to the next
      // line. This ensures that two threads will never write to the same cell in the
      // shared data structures. This slows down things a lot when precomputing
      // probabilities because there is not much to do and threads will spend a lot of 
      // time not trying to step on each others' toes, but when there is a significant 
      // amount of work to do per line (complex model, no memoization of probas), 
      // this works well.
      //
      // Note that this probably assumes that nr >> cores, but it would be a bit of a 
      // pathological case if this was not true.
      // 
      // This has funny indentation but it is the only way to align code in the
      // inside blocks with the non-parallel version
      for (uword c = 0; c < (nr / cores); c++) {
#pragma omp parallel num_threads(cores) default(shared) \
  private(ptrans) reduction(+:new_ps)
      {
      uword i = omp_get_thread_num() * (nr / cores) + c;
      // Make sure we never run the loop if the number of rows is beyond the number
      // of rows in the landscape
      if (i < nr) {
#else
      for (uword i = 0; i < nr; i++) {
#endif
        for (uword j = 0; j < nc; j++) {

          const uchar from = old_mat[i][j];

#if PRECOMPUTE_TRANS_PROBAS
#else
          // Normalized local densities to proportions
          uword qs_total = number_of_neighbors(i, j);

          // Compute probability transitions
          compute_rate(ptrans,
                       old_qs[i][j], // qs
                       old_ps, // ps
                       qs_total, // number of neighbors
                       from); // from (current) state
#endif

#if USE_OMP
          double rn = randunif(omp_get_thread_num());
#else
          double rn = randunif(0);
#endif

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
          uchar new_cell_state = from;
          for (signed char k = (ns - 1); k >= 0; k--) {
#if PRECOMPUTE_TRANS_PROBAS
            uword line = old_pline[i][j];
            new_cell_state = rn < tprobs[line][from][k] ? k : new_cell_state;
#else
            new_cell_state = rn < ptrans[k] ? k : new_cell_state;
#endif
          }

          if ( new_cell_state != from ) {
            new_ps[new_cell_state]++;
            new_ps[from]--;
            new_mat[i][j] = new_cell_state;
#if PRECOMPUTE_TRANS_PROBAS
            adjust_nb_plines(new_pline, i, j, from, new_cell_state);
#else
            adjust_local_density(new_qs, i, j, from, new_cell_state);
#endif
          }

        } // end of loop on j (columns)

#if USE_OMP
      } // closes if() check to make sure lines are within matrix
      } // closes parallel block
#endif
      } // for loop on i

      // Copy old matrix to new, etc.
      memcpy(old_ps, new_ps, sizeof(uword) * ns);
      memcpy(old_mat, new_mat, sizeof(uchar) * nr * nc);
#if PRECOMPUTE_TRANS_PROBAS
      memcpy(old_pline, new_pline, sizeof(uword) * nr * nc);
#else
      memcpy(old_qs, new_qs, sizeof(uchar) * nr * nc * ns);
#endif

    } // end of substep loop

    current_t++;
  }

  // Free up heap-allocated arrays
  delete[] old_mat;
  delete[] new_mat;
#if PRECOMPUTE_TRANS_PROBAS
  delete[] all_qs;
  delete[] tprobs;
  delete[] old_pline;
  delete[] new_pline;
#else
  delete[] old_qs;
  delete[] new_qs;
#endif
}
