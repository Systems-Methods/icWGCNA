// [[Rcpp::depends(RcppEigen)]]

// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//' Faster Matrix Multiplication
//'
//' Faster matrix multiplication using C++.
//'
//' @param A,B Matrices.
//' @export
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                     Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C;
  C = A * B;

  return Rcpp::wrap(C);
}
