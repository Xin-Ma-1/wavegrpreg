// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// L1proj
NumericVector L1proj(NumericVector v, double z, double lambda0);
RcppExport SEXP _wavegrpreg_L1proj(SEXP vSEXP, SEXP zSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    rcpp_result_gen = Rcpp::wrap(L1proj(v, z, lambda0));
    return rcpp_result_gen;
END_RCPP
}
// L12proj
NumericVector L12proj(NumericVector eta, int M, int p, double z, double lambda0);
RcppExport SEXP _wavegrpreg_L12proj(SEXP etaSEXP, SEXP MSEXP, SEXP pSEXP, SEXP zSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    rcpp_result_gen = Rcpp::wrap(L12proj(eta, M, p, z, lambda0));
    return rcpp_result_gen;
END_RCPP
}
// spgaug
arma::vec spgaug(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes, double z, arma::vec theta, arma::vec start, double gamma, double a_min, double a_max, double sig1, double sig2, int Msave, int max1, int max2, double tol, Rcpp::Function L1proj);
RcppExport SEXP _wavegrpreg_spgaug(SEXP WSEXP, SEXP ySEXP, SEXP SigmaSEXP, SEXP sizesSEXP, SEXP zSEXP, SEXP thetaSEXP, SEXP startSEXP, SEXP gammaSEXP, SEXP a_minSEXP, SEXP a_maxSEXP, SEXP sig1SEXP, SEXP sig2SEXP, SEXP MsaveSEXP, SEXP max1SEXP, SEXP max2SEXP, SEXP tolSEXP, SEXP L1projSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_min(a_minSEXP);
    Rcpp::traits::input_parameter< double >::type a_max(a_maxSEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type Msave(MsaveSEXP);
    Rcpp::traits::input_parameter< int >::type max1(max1SEXP);
    Rcpp::traits::input_parameter< int >::type max2(max2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type L1proj(L1projSEXP);
    rcpp_result_gen = Rcpp::wrap(spgaug(W, y, Sigma, sizes, z, theta, start, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj));
    return rcpp_result_gen;
END_RCPP
}
// SPG_gbridge
Rcpp::List SPG_gbridge(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes, double z, double lambda0, arma::vec init, int maxiter, double tol0, double gamma, double a_min, double a_max, double sig1, double sig2, int Msave, int max1, int max2, double tol, Rcpp::Function L1proj);
RcppExport SEXP _wavegrpreg_SPG_gbridge(SEXP WSEXP, SEXP ySEXP, SEXP SigmaSEXP, SEXP sizesSEXP, SEXP zSEXP, SEXP lambda0SEXP, SEXP initSEXP, SEXP maxiterSEXP, SEXP tol0SEXP, SEXP gammaSEXP, SEXP a_minSEXP, SEXP a_maxSEXP, SEXP sig1SEXP, SEXP sig2SEXP, SEXP MsaveSEXP, SEXP max1SEXP, SEXP max2SEXP, SEXP tolSEXP, SEXP L1projSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol0(tol0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_min(a_minSEXP);
    Rcpp::traits::input_parameter< double >::type a_max(a_maxSEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type Msave(MsaveSEXP);
    Rcpp::traits::input_parameter< int >::type max1(max1SEXP);
    Rcpp::traits::input_parameter< int >::type max2(max2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type L1proj(L1projSEXP);
    rcpp_result_gen = Rcpp::wrap(SPG_gbridge(W, y, Sigma, sizes, z, lambda0, init, maxiter, tol0, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj));
    return rcpp_result_gen;
END_RCPP
}
// SPG_glasso
Rcpp::List SPG_glasso(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes, double z, double lambda0, arma::vec init, double gamma, double a_min, double a_max, double sig1, double sig2, int Msave, int max1, int max2, double tol, Rcpp::Function L12proj);
RcppExport SEXP _wavegrpreg_SPG_glasso(SEXP WSEXP, SEXP ySEXP, SEXP SigmaSEXP, SEXP sizesSEXP, SEXP zSEXP, SEXP lambda0SEXP, SEXP initSEXP, SEXP gammaSEXP, SEXP a_minSEXP, SEXP a_maxSEXP, SEXP sig1SEXP, SEXP sig2SEXP, SEXP MsaveSEXP, SEXP max1SEXP, SEXP max2SEXP, SEXP tolSEXP, SEXP L12projSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_min(a_minSEXP);
    Rcpp::traits::input_parameter< double >::type a_max(a_maxSEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type Msave(MsaveSEXP);
    Rcpp::traits::input_parameter< int >::type max1(max1SEXP);
    Rcpp::traits::input_parameter< int >::type max2(max2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type L12proj(L12projSEXP);
    rcpp_result_gen = Rcpp::wrap(SPG_glasso(W, y, Sigma, sizes, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L12proj));
    return rcpp_result_gen;
END_RCPP
}
// SPG_lasso
Rcpp::List SPG_lasso(arma::mat W, arma::vec y, arma::mat Sigma, double z, double lambda0, arma::vec init, double gamma, double a_min, double a_max, double sig1, double sig2, int Msave, int max1, int max2, double tol, Rcpp::Function L1proj);
RcppExport SEXP _wavegrpreg_SPG_lasso(SEXP WSEXP, SEXP ySEXP, SEXP SigmaSEXP, SEXP zSEXP, SEXP lambda0SEXP, SEXP initSEXP, SEXP gammaSEXP, SEXP a_minSEXP, SEXP a_maxSEXP, SEXP sig1SEXP, SEXP sig2SEXP, SEXP MsaveSEXP, SEXP max1SEXP, SEXP max2SEXP, SEXP tolSEXP, SEXP L1projSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_min(a_minSEXP);
    Rcpp::traits::input_parameter< double >::type a_max(a_maxSEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type Msave(MsaveSEXP);
    Rcpp::traits::input_parameter< int >::type max1(max1SEXP);
    Rcpp::traits::input_parameter< int >::type max2(max2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type L1proj(L1projSEXP);
    rcpp_result_gen = Rcpp::wrap(SPG_lasso(W, y, Sigma, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_wavegrpreg_L1proj", (DL_FUNC) &_wavegrpreg_L1proj, 3},
    {"_wavegrpreg_L12proj", (DL_FUNC) &_wavegrpreg_L12proj, 5},
    {"_wavegrpreg_spgaug", (DL_FUNC) &_wavegrpreg_spgaug, 17},
    {"_wavegrpreg_SPG_gbridge", (DL_FUNC) &_wavegrpreg_SPG_gbridge, 19},
    {"_wavegrpreg_SPG_glasso", (DL_FUNC) &_wavegrpreg_SPG_glasso, 17},
    {"_wavegrpreg_SPG_lasso", (DL_FUNC) &_wavegrpreg_SPG_lasso, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_wavegrpreg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
