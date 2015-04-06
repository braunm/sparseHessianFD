#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

using Rcpp::S4;
using Rcpp::List;
using Rcpp::ListOf;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

//[[Rcpp::export]]
NumericMatrix subst2(NumericMatrix Y, IntegerVector colors,
		     ListOf<IntegerVector> W,
		     ListOf<IntegerVector> Sp,
		     double delta) { 
  int nvars = Sp.size();
  NumericMatrix H(nvars, nvars);
  std::fill(H.begin(), H.end(), 0.0);
  
   for (int i=nvars-1; i>=0; i--) {
    int grp = colors(i); // color starts at 1
    for (int j : Sp[i]) {
      if (j <= i) {
	double yi = Y(j, grp);
	double acc = 0;
	for (int k : W[grp]) {	
	  if (k > i) {
	    acc += H(k-1, i);
	  }
	}
	H(i,j) = yi/delta - acc;
	H(j,i) = H(i,j);
      }
    }
  }
	  

  return(H);
}
