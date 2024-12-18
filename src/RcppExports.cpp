// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SuSiE_group
Rcpp::List SuSiE_group(arma::vec y, arma::mat X1, arma::mat X2, double L, double sigma2, double delta, int maxIt);
RcppExport SEXP _CPVS_SuSiE_group(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP LSEXP, SEXP sigma2SEXP, SEXP deltaSEXP, SEXP maxItSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type maxIt(maxItSEXP);
    rcpp_result_gen = Rcpp::wrap(SuSiE_group(y, X1, X2, L, sigma2, delta, maxIt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CPVS_SuSiE_group", (DL_FUNC) &_CPVS_SuSiE_group, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_CPVS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
