// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// load
NumericVector load(NumericVector x);
RcppExport SEXP _odfrid_load(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(load(x));
    return rcpp_result_gen;
END_RCPP
}
// routing_matrix
arma::imat routing_matrix(int s);
RcppExport SEXP _odfrid_routing_matrix(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(routing_matrix(s));
    return rcpp_result_gen;
END_RCPP
}
// test_flat_index
void test_flat_index(int S);
RcppExport SEXP _odfrid_test_flat_index(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type S(SSEXP);
    test_flat_index(S);
    return R_NilValue;
END_RCPP
}
// rod
List rod(NumericVector x);
RcppExport SEXP _odfrid_rod(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rod(x));
    return rcpp_result_gen;
END_RCPP
}
// softmax
NumericVector softmax(const NumericVector& x);
RcppExport SEXP _odfrid_softmax(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(softmax(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_odfrid_load", (DL_FUNC) &_odfrid_load, 1},
    {"_odfrid_routing_matrix", (DL_FUNC) &_odfrid_routing_matrix, 1},
    {"_odfrid_test_flat_index", (DL_FUNC) &_odfrid_test_flat_index, 1},
    {"_odfrid_rod", (DL_FUNC) &_odfrid_rod, 1},
    {"_odfrid_softmax", (DL_FUNC) &_odfrid_softmax, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_odfrid(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
