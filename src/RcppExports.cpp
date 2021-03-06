// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// findzerocorr
List findzerocorr(const vec& xell, const vec& x, const double test, const std::string rf_type);
RcppExport SEXP _gcipdr_findzerocorr(SEXP xellSEXP, SEXP xSEXP, SEXP testSEXP, SEXP rf_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type xell(xellSEXP);
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type test(testSEXP);
    Rcpp::traits::input_parameter< const std::string >::type rf_type(rf_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(findzerocorr(xell, x, test, rf_type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gcipdr_findzerocorr", (DL_FUNC) &_gcipdr_findzerocorr, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_gcipdr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
