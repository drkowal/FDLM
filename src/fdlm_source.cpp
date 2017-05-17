// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Sample the Factor Loading Curves
//'
//' @param BtY (J x T) matrix B.t()*Y for basis matrix B
//' @param Beta (T x K) matrix of factors
//' @param Psi (J x K) matrix of previous factor loading curve coefficients
//' @param BtB (J x J) matrix of B.t()*B
//' @param Omega (J x J) prior precision/penalty matrix
//' @param lamda K-dimensional vector of prior precisions
//' @param sigmat2 T-dimensional vector of time-dependent observation error variances
//' @useDynLib FDLM
//' @import Rcpp
//' @export
// [[Rcpp::export]]
arma::mat sampleFLC(arma::mat BtY, arma::mat Beta, arma::mat Psi, arma::mat BtB,
                           arma::mat Omega, arma::vec lambda, arma::vec sigmat2){
  // Dimensions
  int J = BtY.n_rows;   // Number of basis functions (for each k)
  int K = Beta.n_cols;  // Number of factors/FLCs
  int T = Beta.n_rows;  // Number of time points

  // Terms for constructing sampler for each k
  arma::vec psi_k(J);       // This is the vector that we sample
  arma::mat Quad_sum(J, J); // Quadratic summation term for each k
  arma::vec lin_sum(J);     // Linear summation term for each k
  arma::vec uksamp(J);      // Random sample for psi_k
  arma::mat cholFac(J, J);  // Cholesky factorization of quadratic term

  // Terms for subsetting for each k
  arma::mat beta_j(T, K);               // Matrix of Betas (will subset to T x (K-1) later)
  arma::vec beta_k(T), beta_k_scale(T); // Beta_k and scaled version (i.e., b/variance)
  arma::mat Psi_j(J, K);                // Matrix of FLC coefficients (will subset to J x (K-1) later)

  // Terms for orthogonality constraints
  arma::mat Lcon(J,K-1), BLcon(J,K-1);
  arma::mat psinorm(1,1);

  for(int k = 0; k < K; k++){
    beta_k = Beta.col(k);               // Beta_k
    beta_j = Beta; beta_j.shed_col(k);  // Beta_j for j \ne k
    Psi_j = Psi; Psi_j.shed_col(k);            // d_j for j \ne k

    beta_k_scale = beta_k/sigmat2; // Beta_k/sigmat^2

    // Quadratic term:
    Quad_sum = sum(pow(beta_k, 2.0)/sigmat2)*BtB;
    //Quad_sum = BtB*sum(pow(beta_k, 2.0)/sigmat2);

    // Linear term:
    lin_sum = BtY*beta_k_scale - BtB*Psi_j*beta_j.t()*beta_k_scale;

    // Cholesky factorization of quadratic term:
    cholFac = arma::chol(Quad_sum + lambda(k)*Omega);

    // Sample a normal vector:
    uksamp = rnorm(J);

    // Backsolve/forwardsolve using armadillo:
    psi_k = arma::solve(arma::trimatu(cholFac), arma::solve(arma::trimatl(cholFac.t()), lin_sum) + uksamp);

    // Orthogonality constraints:
    Lcon = BtB*Psi_j;
    BLcon = arma::solve(arma::trimatu(cholFac), arma::solve(arma::trimatl(cholFac.t()), Lcon));

    // psi_k is orthogonal:
    psi_k = psi_k - BLcon*arma::inv_sympd(Lcon.t()*BLcon)*Lcon.t()*psi_k;

    // Normalize:
    psinorm = sqrt(psi_k.t()*BtB*psi_k);

    Psi.col(k) = psi_k/psinorm(0,0);
  }
  return Psi;
}
//' Sample the Factor Loading Curves, For the special case in which BtB = diag(J)
//'
//' @param BtY (J x T) matrix B.t()*Y for basis matrix B
//' @param Beta (T x K) matrix of factors
//' @param Psi (J x K) matrix of previous factor loading curve coefficients
//' @param Omega (J x J) prior precision/penalty matrix
//' @param lamda K-dimensional vector of prior precisions
//' @param sigmat2 T-dim vector of time-dependent observation error variances
//' @export
// [[Rcpp::export]]
arma::mat sampleFLC_orthog(arma::mat BtY, arma::mat Beta, arma::mat Psi,
                    arma::mat Omega, arma::vec lambda, arma::vec sigmat2){
  // Dimensions
  int J = BtY.n_rows;   // Number of basis functions (for each k)
  int K = Beta.n_cols;  // Number of factors/FLCs
  int T = Beta.n_rows;  // Number of time points

  // Terms for constructing sampler for each k
  arma::vec psi_k(J);       // This is the vector that we sample
  arma::mat Quad_sum(J, J); // Quadratic summation term for each k
  arma::vec lin_sum(J);     // Linear summation term for each k
  arma::vec uksamp(J);      // Random sample for psi_k
  arma::mat cholFac(J, J);  // Cholesky factorization of quadratic term

  // Terms for subsetting for each k
  arma::mat beta_j(T, K);               // Matrix of Betas (will subset to T x (K-1) later)
  arma::vec beta_k(T), beta_k_scale(T); // Beta_k and scaled version (i.e., b/variance)
  arma::mat Psi_j(J, K);                // Matrix of FLC coefficients (will subset to J x (K-1) later)

  // Terms for orthogonality constraints
  arma::mat BLcon(J,K-1);
  arma::mat psinorm(1,1);

  for(int k = 0; k < K; k++){
    beta_k = Beta.col(k);               // Beta_k
    beta_j = Beta; beta_j.shed_col(k);  // Beta_j for j \ne k
    Psi_j = Psi; Psi_j.shed_col(k);            // d_j for j \ne k

    beta_k_scale = beta_k/sigmat2; // Beta_k/sigmat^2

    // Quadratic term:
    Quad_sum = sum(pow(beta_k, 2.0)/sigmat2)*Quad_sum.eye();

    // Linear term:
    lin_sum = BtY*beta_k_scale - Psi_j*beta_j.t()*beta_k_scale;

    // Cholesky factorization of quadratic term:
    cholFac = arma::chol(Quad_sum + lambda(k)*Omega);

    // Sample a normal vector:
    uksamp = rnorm(J);

    // Backsolve/forwardsolve using armadillo:
    psi_k = arma::solve(arma::trimatu(cholFac), arma::solve(arma::trimatl(cholFac.t()), lin_sum) + uksamp);

    // Orthogonality constraints:
    BLcon = arma::solve(arma::trimatu(cholFac), arma::solve(arma::trimatl(cholFac.t()), Psi_j));

    // psi_k is orthogonal:
    psi_k = psi_k - BLcon*arma::inv_sympd(Psi_j.t()*BLcon)*Psi_j.t()*psi_k;

    // Normalize:
    psinorm = sqrt(psi_k.t()*psi_k);
    Psi.col(k) = psi_k/psinorm(0,0);
  }
  return Psi;
}
