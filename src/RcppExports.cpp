// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sampleFLC
arma::mat sampleFLC(arma::mat BtY, arma::mat Beta, arma::mat Psi, arma::mat BtB, arma::mat Omega, arma::vec lambda, arma::vec sigmat2);
RcppExport SEXP FDLM_sampleFLC(SEXP BtYSEXP, SEXP BetaSEXP, SEXP PsiSEXP, SEXP BtBSEXP, SEXP OmegaSEXP, SEXP lambdaSEXP, SEXP sigmat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type BtY(BtYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type BtB(BtBSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigmat2(sigmat2SEXP);
    rcpp_result_gen = Rcpp::wrap(sampleFLC(BtY, Beta, Psi, BtB, Omega, lambda, sigmat2));
    return rcpp_result_gen;
END_RCPP
}
// sampleFLC_orthog
arma::mat sampleFLC_orthog(arma::mat BtY, arma::mat Beta, arma::mat Psi, arma::mat Omega, arma::vec lambda, arma::vec sigmat2);
RcppExport SEXP FDLM_sampleFLC_orthog(SEXP BtYSEXP, SEXP BetaSEXP, SEXP PsiSEXP, SEXP OmegaSEXP, SEXP lambdaSEXP, SEXP sigmat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type BtY(BtYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigmat2(sigmat2SEXP);
    rcpp_result_gen = Rcpp::wrap(sampleFLC_orthog(BtY, Beta, Psi, Omega, lambda, sigmat2));
    return rcpp_result_gen;
END_RCPP
}
