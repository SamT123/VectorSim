#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat v(arma::colvec a) {
  return a*a.t();}



// [[Rcpp::export]]
arma::cube test_fn(arma::cube M){
  

  
  M.subcube(0,1,1,1,1,1) = M.subcube(0,2,1,1,2,1);  //= M(arma::span::all, 1, 1) ;
  return M;
}

// [[Rcpp::export]]
List test_fn_list(Rcpp::List l){
  for (int i = 0; i<l.size(); i++){
    Rcpp::List inner_l = l[i];
    l[i] = (int)inner_l.size();
  }
  return l;
}

