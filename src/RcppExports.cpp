// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// broadcast_test_1
int broadcast_test_1();
RcppExport SEXP _odfrid_broadcast_test_1() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(broadcast_test_1());
    return rcpp_result_gen;
END_RCPP
}
// broadcast_test_2
int broadcast_test_2();
RcppExport SEXP _odfrid_broadcast_test_2() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(broadcast_test_2());
    return rcpp_result_gen;
END_RCPP
}
// broadcast_test
int broadcast_test();
RcppExport SEXP _odfrid_broadcast_test() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(broadcast_test());
    return rcpp_result_gen;
END_RCPP
}
// routing_matrix
arma::umat routing_matrix(arma::uword s);
RcppExport SEXP _odfrid_routing_matrix(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(routing_matrix(s));
    return rcpp_result_gen;
END_RCPP
}
// model_sample
Rcpp::List model_sample(const arma::umat& ax, const arma::vec& dep_time, const arma::uword sample, const arma::uword warmup, const arma::uword D, const arma::uword print_n);
RcppExport SEXP _odfrid_model_sample(SEXP axSEXP, SEXP dep_timeSEXP, SEXP sampleSEXP, SEXP warmupSEXP, SEXP DSEXP, SEXP print_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type ax(axSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type dep_time(dep_timeSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type print_n(print_nSEXP);
    rcpp_result_gen = Rcpp::wrap(model_sample(ax, dep_time, sample, warmup, D, print_n));
    return rcpp_result_gen;
END_RCPP
}
// load
arma::umat load(arma::umat& x);
RcppExport SEXP _odfrid_load(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(load(x));
    return rcpp_result_gen;
END_RCPP
}
// rod
Rcpp::List rod(const arma::umat& x);
RcppExport SEXP _odfrid_rod(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rod(x));
    return rcpp_result_gen;
END_RCPP
}
// test
void test();
RcppExport SEXP _odfrid_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_odfrid_broadcast_test_1", (DL_FUNC) &_odfrid_broadcast_test_1, 0},
    {"_odfrid_broadcast_test_2", (DL_FUNC) &_odfrid_broadcast_test_2, 0},
    {"_odfrid_broadcast_test", (DL_FUNC) &_odfrid_broadcast_test, 0},
    {"_odfrid_routing_matrix", (DL_FUNC) &_odfrid_routing_matrix, 1},
    {"_odfrid_model_sample", (DL_FUNC) &_odfrid_model_sample, 6},
    {"_odfrid_load", (DL_FUNC) &_odfrid_load, 1},
    {"_odfrid_rod", (DL_FUNC) &_odfrid_rod, 1},
    {"_odfrid_test", (DL_FUNC) &_odfrid_test, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_odfrid(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
