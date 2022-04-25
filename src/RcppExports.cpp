// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// local_dens
arma::Col<arma::uword> local_dens(const arma::Mat<ushort> m, const arma::uword nstates, const arma::uword i, const arma::uword j, const bool wrap, const bool use_8_nb);
RcppExport SEXP _chouca_local_dens(SEXP mSEXP, SEXP nstatesSEXP, SEXP iSEXP, SEXP jSEXP, SEXP wrapSEXP, SEXP use_8_nbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<ushort> >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nstates(nstatesSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type j(jSEXP);
    Rcpp::traits::input_parameter< const bool >::type wrap(wrapSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_8_nb(use_8_nbSEXP);
    rcpp_result_gen = Rcpp::wrap(local_dens(m, nstates, i, j, wrap, use_8_nb));
    return rcpp_result_gen;
END_RCPP
}
// local_dens_col
arma::Mat<arma::uword> local_dens_col(const arma::Mat<ushort> m, const arma::uword nstates, const arma::uword j, const bool wrap, const bool use_8_nb);
RcppExport SEXP _chouca_local_dens_col(SEXP mSEXP, SEXP nstatesSEXP, SEXP jSEXP, SEXP wrapSEXP, SEXP use_8_nbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<ushort> >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nstates(nstatesSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type j(jSEXP);
    Rcpp::traits::input_parameter< const bool >::type wrap(wrapSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_8_nb(use_8_nbSEXP);
    rcpp_result_gen = Rcpp::wrap(local_dens_col(m, nstates, j, wrap, use_8_nb));
    return rcpp_result_gen;
END_RCPP
}
// camodel_cpp_engine
void camodel_cpp_engine(const arma::cube trans, const Rcpp::List ctrl, const Rcpp::Function console_callback, const Rcpp::Function cover_callback, const Rcpp::Function snapshot_callback);
RcppExport SEXP _chouca_camodel_cpp_engine(SEXP transSEXP, SEXP ctrlSEXP, SEXP console_callbackSEXP, SEXP cover_callbackSEXP, SEXP snapshot_callbackSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube >::type trans(transSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type ctrl(ctrlSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type console_callback(console_callbackSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type cover_callback(cover_callbackSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type snapshot_callback(snapshot_callbackSEXP);
    camodel_cpp_engine(trans, ctrl, console_callback, cover_callback, snapshot_callback);
    return R_NilValue;
END_RCPP
}
// getline
arma::uword getline(const arma::uvec& qs, const arma::uword& nb, const arma::uword& ns);
RcppExport SEXP _chouca_getline(SEXP qsSEXP, SEXP nbSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(getline(qs, nb, ns));
    return rcpp_result_gen;
END_RCPP
}
// simple_sum
arma::uword simple_sum(const arma::uvec& qs, const arma::uword& nb, const arma::uword& ns);
RcppExport SEXP _chouca_simple_sum(SEXP qsSEXP, SEXP nbSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(simple_sum(qs, nb, ns));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chouca_local_dens", (DL_FUNC) &_chouca_local_dens, 6},
    {"_chouca_local_dens_col", (DL_FUNC) &_chouca_local_dens_col, 5},
    {"_chouca_camodel_cpp_engine", (DL_FUNC) &_chouca_camodel_cpp_engine, 5},
    {"_chouca_getline", (DL_FUNC) &_chouca_getline, 3},
    {"_chouca_simple_sum", (DL_FUNC) &_chouca_simple_sum, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_chouca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
