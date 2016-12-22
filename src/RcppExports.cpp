// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_colors
Rcpp::IntegerVector get_colors(const IntegerVector& pntr, const IntegerVector& idx, const int nvars);
RcppExport SEXP sparseHessianFD_get_colors(SEXP pntrSEXP, SEXP idxSEXP, SEXP nvarsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type pntr(pntrSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const int >::type nvars(nvarsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_colors(pntr, idx, nvars));
    return rcpp_result_gen;
END_RCPP
}
// subst
Rcpp::S4 subst(const Rcpp::NumericMatrix& Y, const Rcpp::IntegerVector& colors, const Rcpp::IntegerVector& jCol, const Rcpp::IntegerVector& ipntr, const double& delta, const int& nvars, const int& nnz);
RcppExport SEXP sparseHessianFD_subst(SEXP YSEXP, SEXP colorsSEXP, SEXP jColSEXP, SEXP ipntrSEXP, SEXP deltaSEXP, SEXP nvarsSEXP, SEXP nnzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type colors(colorsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type jCol(jColSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ipntr(ipntrSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nvars(nvarsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nnz(nnzSEXP);
    rcpp_result_gen = Rcpp::wrap(subst(Y, colors, jCol, ipntr, delta, nvars, nnz));
    return rcpp_result_gen;
END_RCPP
}
