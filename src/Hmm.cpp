#include <RcppArmadillo.h>
#include <time.h>
#include <stdlib.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double logplus(double a, double b) {
  return std::max(b, a) + log(1 + exp(-std::abs(b - a)));
}

double logplus_approx(double a, double b,
                      double precision=log(pow(10, 20))) {
  // another approximate version of logplus
  // Achieves a relative error approximately e^(-precision).
  if(std::abs(b - a) > precision) {
    return std::max(b, a);
  }else{
    return std::max(b, a) + log(1 + exp(-std::abs(b - a)));
  }
}

// [[Rcpp::export]]
double logsum(arma::vec y) {
  // Returns log of sum(exp(Y))  without overflow or underflow.
  if (y.n_elem == 0) return 1;
  double x = y[0];
  for (int i = 1; i < y.n_elem; i++) {
      x = logplus(x, y[i]);
  }
  return x;
}

// [[Rcpp::export]]
List lAlphaBetaGamma(arma::mat lemission, arma::mat A, arma::colvec Pi) {
  double logll;
  int N_s = lemission.n_rows; // Length of this sequence
  int length = lemission.n_cols; //Number of states
  arma::mat lalpha(N_s, length);
  arma::mat lbeta(N_s, length);
  arma::vec foo(N_s);
  lalpha.col(0) = Pi + lemission.col(0);
  for (int t = 1; t < length ; t++) {
    for (int j = 0; j < N_s; j++) {
      for (int i = 0; i < N_s; i++) {
        foo[i] = lalpha(i, t - 1) + A(i, j);
      }
      lalpha(j, t) = logsum(foo) + lemission(j, t);
    }
  }
  lbeta.col(length - 1) = arma::zeros<arma::vec>(N_s);
  for (int t = length - 2; t > -1; t--) {
    for (int i = 0; i < N_s; i++){
      for (int j = 0; j < N_s; j++){
        foo[j] = A(i, j) + lemission(j, t + 1) + lbeta(j, t + 1);
      }
      lbeta(i, t) = logsum(foo);
    }
  }
  arma::mat lab(N_s, length);
  arma::rowvec L(length);
  arma::mat lGamma(N_s, length);
  arma::mat lXi(N_s, N_s);
  lab = lalpha + lbeta;
  double temp = 0.0;
  for (int t = 0; t < length; t++) {
    L[t] = logsum(lab.col(t));
  }
  for (int i = 0; i < N_s; i++) {
    lGamma.row(i) = lab.row(i) - L;
  }
  for (int i = 0; i < N_s; i++) {
    for (int j = 0; j < N_s; j++) {
      lXi(i, j) = lalpha(i, 0) + A(i, j) + lemission(j, 1) + lbeta(j, 1) - L(1);
    }
  }
    for (int i = 0;  i < N_s; i++) {
      for (int j = 0; j < N_s; j++) {
        for (int t = 1; t < length - 1; t++) {
          lXi(i, j) = logplus(lXi(i, j), lalpha(i, t) + A(i, j) +
                      lemission(j, t + 1) + lbeta(j, t + 1) - L(t + 1));
      }
    }
  }
  A = lXi;
  for (int i = 0; i < N_s; i++){
    temp = logsum(A.row(i).t());
    for (int j = 0; j < N_s; j++){
      A(i, j) = A(i, j) - temp;
    }
  }
  logll = logsum(lalpha.col(length - 1));
  return List::create(Named("lalpha")=lalpha,
                      Named("lbeta")=lbeta,
                      Named("lGamma")=lGamma,
                      Named("lA")=A,
                      Named("logll")=logll);
}

// [[Rcpp::export]]
List AlphaBetaGamma_scaled(arma::mat emission, arma::mat A, arma::colvec Pi){
  int N_s = emission.n_rows; //Number of states
  int length = emission.n_cols; // Length of this sequence
  arma::mat alpha(N_s, length, arma::fill::zeros);
  arma::mat beta(N_s, length, arma::fill::zeros);
  arma::mat Gamma(N_s, length);
  arma::vec Scale(length);
  alpha.col(0) = Pi % emission.col(0);
  Scale[0] = 1 / sum(alpha.col(0));
  alpha.col(0) *= Scale[0];
  for (int t = 1; t < length ; t++){
    for (int j = 0; j < N_s; j++){
      for (int i = 0; i < N_s; i++){
        alpha(j, t) += alpha(i, t - 1) * A(i, j) * emission(j, t);
      }
    }
    Scale[t] = 1 / sum(alpha.col(t));
    alpha.col(t) *= Scale[t];
  }
  beta.col(length - 1) = arma::ones<arma::vec>(N_s);
  beta.col(length - 1) *= Scale[length - 1];
  for (int t = length - 2; t > -1; t--){
    for (int i = 0; i < N_s; i++){
      for (int j = 0; j < N_s; j++){
        beta(i, t) += beta(j, t+1) * A(i, j) * emission(j, t+1);
      }
    }
    beta.col(t) *= Scale[t];
  }
  Gamma = alpha % beta;
  for (int i=0; i<N_s; i++){
    Gamma.row(i) /= Scale.t();
  }
  arma::mat xisum(N_s, N_s, arma::fill::zeros);
  for (int i = 0; i< N_s; i++) {
    for (int j = 0;j < N_s; j++) {
      for (int t = 0;t < length - 1; t++){
        xisum(i, j) += alpha(i, t) * A(i, j) * emission(j, t+1) * beta(j, t+1);
      }
    }
  }
  A = xisum;
  double temp;
  for (int i = 0; i < N_s; i++){
    temp = sum(A.row(i));
    for (int j = 0;j < N_s; j++){
      A(i, j) /= temp;
    }
  }
  return List::create(Named("alpha")=alpha,
                      Named("beta")=beta,
                      Named("Gamma")=Gamma,
                      Named("A")=A,
                      Named("Scale")=Scale);
}

// [[Rcpp::export]]
IntegerVector Viterbi(arma::mat A, arma::colvec Pi, arma::mat emission){
  int N_s = emission.n_rows; //Number of states
  int length = emission.n_cols; // Length of this sequence
  arma::mat T1(N_s, length);
  IntegerMatrix T2(N_s, length);
  arma::colvec maxtemp(N_s);
  for (int j = 0; j < N_s; j++){
    T1(j, 0) = logplus(Pi[j], emission(j, 0));
  }
  for (int i = 1; i < length ; i++){
    for (int j = 0; j < N_s; j++){
      for (int k = 0; k < N_s; k++){
        maxtemp[k] = logplus(T1(k, i - 1), A(k, j));
      }
      T1(j, i) = logplus(maxtemp.max(), emission(j, i));
      T2(j, i) = maxtemp.index_max();
    }
  }
  IntegerVector X(length);
  X[length - 1] = T1.col(length - 1).index_max();
  for (int i = length - 1; i > 0; i--){
    X[i - 1] = T2(X[i], i);
  }
  return X;
}
