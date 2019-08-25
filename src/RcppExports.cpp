// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// localRegression
arma::mat localRegression(arma::mat weightmat, arma::mat modelmat, arma::vec xtemp);
RcppExport SEXP _RSEE_localRegression(SEXP weightmatSEXP, SEXP modelmatSEXP, SEXP xtempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type weightmat(weightmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modelmat(modelmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xtemp(xtempSEXP);
    rcpp_result_gen = Rcpp::wrap(localRegression(weightmat, modelmat, xtemp));
    return rcpp_result_gen;
END_RCPP
}
// iterLowess
arma::mat iterLowess(arma::mat weightX, arma::mat weightY, arma::mat modelmat, arma::vec x, int h, double smoothpara, double epsmed);
RcppExport SEXP _RSEE_iterLowess(SEXP weightXSEXP, SEXP weightYSEXP, SEXP modelmatSEXP, SEXP xSEXP, SEXP hSEXP, SEXP smoothparaSEXP, SEXP epsmedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type weightX(weightXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weightY(weightYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modelmat(modelmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type smoothpara(smoothparaSEXP);
    Rcpp::traits::input_parameter< double >::type epsmed(epsmedSEXP);
    rcpp_result_gen = Rcpp::wrap(iterLowess(weightX, weightY, modelmat, x, h, smoothpara, epsmed));
    return rcpp_result_gen;
END_RCPP
}
// lastIterLowess
arma::mat lastIterLowess(arma::mat weightX, arma::mat weightY, arma::mat modelmat, arma::vec x, int h);
RcppExport SEXP _RSEE_lastIterLowess(SEXP weightXSEXP, SEXP weightYSEXP, SEXP modelmatSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type weightX(weightXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weightY(weightYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modelmat(modelmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(lastIterLowess(weightX, weightY, modelmat, x, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RSEE_localRegression", (DL_FUNC) &_RSEE_localRegression, 3},
    {"_RSEE_iterLowess", (DL_FUNC) &_RSEE_iterLowess, 7},
    {"_RSEE_lastIterLowess", (DL_FUNC) &_RSEE_lastIterLowess, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_RSEE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
