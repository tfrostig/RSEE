// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat localRegression(arma::mat weightmat, arma::mat modelmat, arma::vec xtemp) {
  return inv(modelmat.t() * weightmat * modelmat) * modelmat.t() * weightmat * xtemp;
}

// [[Rcpp::export]]
arma::mat iterLowess(arma::mat weightX, arma::mat weightY, arma::mat modelmat, arma::vec x, int h, double smoothpara, double epsmed = 10^(-6)) {
  int n = x.size();
  int p = 2 * h + 1;
  double medres; 
  arma::mat W;
  arma::vec xtemp, res, beta, tempdelta;
  arma::uvec ids;
  W.zeros(p,p); 
  for (int i = h; i < (n - h); i++) {
    xtemp = x(span(i - h, i + h));
    W.diag() = weightX.row(i - h) % weightY.row(i - h);
    beta = localRegression(W, modelmat, xtemp);
    res = abs(xtemp - modelmat * beta) + epsmed;
    medres = median(res.elem(find(W.diag() > 0)));
    tempdelta = pow(1 - pow(res / (medres * smoothpara), 2), 2);
    ids = find(res > smoothpara * medres); // Find indices
    tempdelta.elem(ids).fill(0); 
    weightY.row(i - h) = tempdelta.t(); 
  }
  return weightY;
} 

// [[Rcpp::export]]
arma::mat lastIterLowess(arma::mat weightX, arma::mat weightY, arma::mat modelmat, arma::vec x, int h) {
  int n = x.size();
  int p = 2 * h + 1;
  arma::mat W, betamat; 
  arma::vec xtemp, beta;
  arma::uvec ids;
  W.zeros(p,p); 
  betamat.zeros(n - 2 * h, modelmat.n_cols);
  for (int i = h; i < (n - h); i++) {
    xtemp = x(span(i - h, i + h));
    W.diag() = weightX.row(i - h) % weightY.row(i - h);
    beta = localRegression(W, modelmat, xtemp);
    betamat.row(i - h) = beta.t(); 
  }
  return betamat;
} 

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R


*/
