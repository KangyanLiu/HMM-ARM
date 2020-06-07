#include <RcppArmadillo.h>
#include <time.h>
#include <stdlib.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
arma::mat UpdateTheta(arma::vec yt, arma::vec logpi, arma::mat logA, arma::vec D, arma::vec CD, arma::vec C,
                double phi, arma::vec mu, arma::mat theta) {
  int n = yt.n_elem;
  arma::vec foo(3);
  arma::vec B1(3);
  arma::vec B2(3);
  arma::mat theta_new(theta.n_rows, theta.n_cols);
  for (int t = 0; t < n; t++) {
    for (int d = 0; d < D[t]; d++) {
      foo = yt[t] / C[t] * pow(phi, (D[t] - 1 - d)) * mu;
      for (int j = 0; j < D[t]; j++) {
        if (d == j) {
          foo = foo - diagvec((pow(phi, (D[t] - 1 - d)) * mu) *
                (pow(phi, (D[t] - 1 - j)) * mu).t()) / C[t] / 2;
        } else {
          foo = foo - (pow(phi, (D[t] - 1 - d)) * mu) *
                (pow(phi, (D[t] - 1 - j)) * mu).t() *
                theta.col(CD[t] + j) / C[t];
        }
      }
      if (d == 0) {
        if (t == 0) {
          B1 = logpi;
        } else {
          B1 = logA * theta.col(CD[t - 1] + D[t - 1] - 1);
        }
      } else {
        B1 = logA * theta.col(CD[t] + d - 1);
      }
      if (d == (D[t] - 1)){
        if (t == (n - 1)){
          B2 = arma::zeros<arma::vec>(3);
        } else {
          B2 = logA.t() * theta.col(CD[t + 1]);
        }
      } else {
        B2 = logA.t() * theta.col(CD[t] + d + 1);
      }
      foo = foo + B1 + B2;
      theta_new.col(CD[t] + d)= exp(foo) / sum(exp(foo));
    }
  }
  return theta_new;
}

// [[Rcpp::export]]
arma::vec UpdateMu(arma::vec yt, arma::vec D, arma::vec C, arma::vec CD,
             double phi, arma::mat theta) {
  arma::mat X = arma::zeros<arma::mat>(3, 3);
  arma::mat Y = arma::zeros<arma::vec>(3);
  for (int t = 0; t < yt.n_elem; t++) {
    for (int i = 0; i < D[t]; i++) {
      for (int j = 0; j < D[t]; j++) {
        if (i == j) {
          X = X + pow(phi, 2 * (D[t] - 1 - i)) /
              C[t] * diagmat(theta.col(CD[t] + i));
        } else {
          X = X + (pow(phi, (D[t] - 1 - i)) *
              theta.col(CD[t] + i)) *
              (pow(phi, (D[t] -1 - j)) * theta.col(CD[t] + j)).t() / C[t];
        }
      }
      Y = Y + pow(phi, (D[t] - 1 - i)) * yt[t] / C[t] * theta.col(CD[t] + i);
    }
  }
  X = X.submat(0, 1, 2, 2);
  return arma::solve(X, Y);
}

// [[Rcpp::export]]
double UpdateS2r(arma::vec yt, arma::vec D, arma::vec CD, arma::vec C,
                 double phi, arma::vec mu, arma::mat theta) {
  double ll = 0.0;
  for (int t = 0; t < yt.n_elem; t++) {
    ll = ll + pow(yt[t], 2) / C[t];
    for (int i = 0; i < D[t]; i++) {
      ll = ll - 2 * yt[t] / C[t] *
           dot(pow(phi, (D[t] - 1 - i)) * mu, theta.col(CD[t] + i));
      for (int j = 0; j < D[t]; j++) {
        if (i == j) {
          ll = ll+ trace((pow(phi, (D[t] - 1 - i)) * mu) *
               (pow(phi, (D[t] - 1 - i)) * mu).t() *
               diagmat(theta.col(CD[t] + i))) / C[t];
        } else {
          ll = ll +  dot((pow(phi, (D[t] - 1 - i)) * mu).t() *
               theta.col(CD[t] + i), (pow(phi, (D[t] - 1 - j)) * mu).t() *
               theta.col(CD[t] + j)) / C[t];
        }
      }
    }
  }
  return ll / yt.n_elem;
  //return ll/sum(sum(theta));
}

// [[Rcpp::export]]
double l1(double phi, arma::vec y, arma::vec x, arma::vec D, arma::vec CD,
          double s2r, arma::vec mu, arma::mat theta) {
  arma::vec C(D.n_elem);
  double cs2r = 0.0;
  double ll = 0.0;
  double yt = 0.0;
  for (int t = 0; t < y.n_elem; t++) {
    yt = y[t] - pow(phi, D[t]) * x[t];
    cs2r = (1 - pow(phi, 2 * D[t])) / (1 - pow(phi, 2)) * s2r;
    ll = ll + pow(yt, 2) / cs2r;
    ll = ll + log(cs2r);
    for (int i = 0; i < D[t]; i++) {
      ll = ll - 2 * yt / cs2r *
           dot(pow(phi, (D[t] - 1 - i)) * mu, theta.col(CD[t] + i));
      for (int j = 0; j < D[t]; j++) {
        if (i == j) {
          ll = ll + trace((pow(phi, (D[t] - 1 - i)) * mu) *
                          (pow(phi, (D[t] - 1 - i)) * mu).t() *
                          diagmat(theta.col(CD[t] + i))) / cs2r;
        } else {
          ll = ll + dot((pow(phi, (D[t] - 1 - i)) * mu).t() *
                              theta.col(CD[t] + i),
                              (pow(phi, (D[t] - 1 - j)) * mu).t() *
                               theta.col(CD[t] + j)) / cs2r;
        }
      }
    }
  }
  return ll;
}

// [[Rcpp::export]]
double l2(arma::vec logpi, arma::mat logA, arma::mat S) {
  double ll = dot(S.col(0), logpi);
  for (int i = 1; i < S.n_cols; i++){
    ll = ll + dot(S.col(i), logA * S.col(i - 1));
  }
  return ll;
}



