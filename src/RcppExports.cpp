// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getOne
double getOne();
RcppExport SEXP _CompStat_getOne() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getOne());
    return rcpp_result_gen;
END_RCPP
}
// eval_kdensC
NumericVector eval_kdensC(String kcode, NumericVector grid, NumericVector x, double bw);
RcppExport SEXP _CompStat_eval_kdensC(SEXP kcodeSEXP, SEXP gridSEXP, SEXP xSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_kdensC(kcode, grid, x, bw));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CompStat_getOne", (DL_FUNC) &_CompStat_getOne, 0},
    {"_CompStat_eval_kdensC", (DL_FUNC) &_CompStat_eval_kdensC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_CompStat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
