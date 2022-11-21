// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CellCounts
IntegerMatrix CellCounts(List x, List combos);
RcppExport SEXP _COMPASS_CellCounts(SEXP xSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(CellCounts(x, combos));
    return rcpp_result_gen;
END_RCPP
}
// CellCounts_character
IntegerMatrix CellCounts_character(List data, List combinations);
RcppExport SEXP _COMPASS_CellCounts_character(SEXP dataSEXP, SEXP combinationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type combinations(combinationsSEXP);
    rcpp_result_gen = Rcpp::wrap(CellCounts_character(data, combinations));
    return rcpp_result_gen;
END_RCPP
}
// updategammak_cov
SEXP updategammak_cov(SEXP n_s, SEXP n_u, SEXP gammat, SEXP I, SEXP K, SEXP SS, SEXP alphau, SEXP alphas, SEXP alpha, SEXP mk, SEXP Istar, SEXP mKstar, SEXP pp, SEXP pb1, SEXP pb2, SEXP indi, SEXP WK1, SEXP WK0);
RcppExport SEXP _COMPASS_updategammak_cov(SEXP n_sSEXP, SEXP n_uSEXP, SEXP gammatSEXP, SEXP ISEXP, SEXP KSEXP, SEXP SSSEXP, SEXP alphauSEXP, SEXP alphasSEXP, SEXP alphaSEXP, SEXP mkSEXP, SEXP IstarSEXP, SEXP mKstarSEXP, SEXP ppSEXP, SEXP pb1SEXP, SEXP pb2SEXP, SEXP indiSEXP, SEXP WK1SEXP, SEXP WK0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_s(n_sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_u(n_uSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gammat(gammatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type I(ISEXP);
    Rcpp::traits::input_parameter< SEXP >::type K(KSEXP);
    Rcpp::traits::input_parameter< SEXP >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alphau(alphauSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mk(mkSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Istar(IstarSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mKstar(mKstarSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pb1(pb1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pb2(pb2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type indi(indiSEXP);
    Rcpp::traits::input_parameter< SEXP >::type WK1(WK1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type WK0(WK0SEXP);
    rcpp_result_gen = Rcpp::wrap(updategammak_cov(n_s, n_u, gammat, I, K, SS, alphau, alphas, alpha, mk, Istar, mKstar, pp, pb1, pb2, indi, WK1, WK0));
    return rcpp_result_gen;
END_RCPP
}
