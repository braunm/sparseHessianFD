// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// color
List color(const IntegerVector& pntr, const IntegerVector& idx, const int nvars);
RcppExport SEXP sparseHessianFD_color(SEXP pntrSEXP, SEXP idxSEXP, SEXP nvarsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector& >::type pntr(pntrSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const int >::type nvars(nvarsSEXP);
    __result = Rcpp::wrap(color(pntr, idx, nvars));
    return __result;
END_RCPP
}
// subst_C
S4 subst_C(const NumericMatrix& Y, const IntegerVector& colors, const ListOf<IntegerVector>& W, const ListOf<IntegerVector>& Sp, IntegerVector& colsize_, const double& delta);
RcppExport SEXP sparseHessianFD_subst_C(SEXP YSEXP, SEXP colorsSEXP, SEXP WSEXP, SEXP SpSEXP, SEXP colsize_SEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colors(colorsSEXP);
    Rcpp::traits::input_parameter< const ListOf<IntegerVector>& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const ListOf<IntegerVector>& >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type colsize_(colsize_SEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    __result = Rcpp::wrap(subst_C(Y, colors, W, Sp, colsize_, delta));
    return __result;
END_RCPP
}
// subst2
S4 subst2(const NumericMatrix& Y, const IntegerVector& colors, const ListOf<IntegerVector>& W, const IntegerVector& jCol_, const IntegerVector& ipntr_, const double& delta, const int& nvars, const int& nnz);
RcppExport SEXP sparseHessianFD_subst2(SEXP YSEXP, SEXP colorsSEXP, SEXP WSEXP, SEXP jCol_SEXP, SEXP ipntr_SEXP, SEXP deltaSEXP, SEXP nvarsSEXP, SEXP nnzSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colors(colorsSEXP);
    Rcpp::traits::input_parameter< const ListOf<IntegerVector>& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type jCol_(jCol_SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ipntr_(ipntr_SEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nvars(nvarsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nnz(nnzSEXP);
    __result = Rcpp::wrap(subst2(Y, colors, W, jCol_, ipntr_, delta, nvars, nnz));
    return __result;
END_RCPP
}
