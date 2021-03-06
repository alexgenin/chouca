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
void camodel_cpp_engine(const arma::Mat<ushort> alpha_index, const arma::Col<double> alpha_vals, const arma::Mat<ushort> pmat_index, const arma::Mat<double> pmat_vals, const arma::Mat<ushort> qmat_index, const arma::Col<double> qmat_vals, const Rcpp::List ctrl, const Rcpp::Function console_callback, const Rcpp::Function cover_callback, const Rcpp::Function snapshot_callback);
RcppExport SEXP _chouca_camodel_cpp_engine(SEXP alpha_indexSEXP, SEXP alpha_valsSEXP, SEXP pmat_indexSEXP, SEXP pmat_valsSEXP, SEXP qmat_indexSEXP, SEXP qmat_valsSEXP, SEXP ctrlSEXP, SEXP console_callbackSEXP, SEXP cover_callbackSEXP, SEXP snapshot_callbackSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<ushort> >::type alpha_index(alpha_indexSEXP);
    Rcpp::traits::input_parameter< const arma::Col<double> >::type alpha_vals(alpha_valsSEXP);
    Rcpp::traits::input_parameter< const arma::Mat<ushort> >::type pmat_index(pmat_indexSEXP);
    Rcpp::traits::input_parameter< const arma::Mat<double> >::type pmat_vals(pmat_valsSEXP);
    Rcpp::traits::input_parameter< const arma::Mat<ushort> >::type qmat_index(qmat_indexSEXP);
    Rcpp::traits::input_parameter< const arma::Col<double> >::type qmat_vals(qmat_valsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type ctrl(ctrlSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type console_callback(console_callbackSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type cover_callback(cover_callbackSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function >::type snapshot_callback(snapshot_callbackSEXP);
    camodel_cpp_engine(alpha_index, alpha_vals, pmat_index, pmat_vals, qmat_index, qmat_vals, ctrl, console_callback, cover_callback, snapshot_callback);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chouca_local_dens", (DL_FUNC) &_chouca_local_dens, 6},
    {"_chouca_local_dens_col", (DL_FUNC) &_chouca_local_dens_col, 5},
    {"_chouca_camodel_cpp_engine", (DL_FUNC) &_chouca_camodel_cpp_engine, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_chouca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
