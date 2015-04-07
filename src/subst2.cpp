#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <typeinfo>

using Rcpp::S4;
using Rcpp::List;
using Rcpp::ListOf;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::Rcout;
using Rcpp::sapply;
using Eigen::VectorXi;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::SparseMatrix;

template <typename T> std::string type_name();

int get_size(const IntegerVector& x) {
  return(x.size());
}

//[[Rcpp::export]]
S4 subst_C(NumericMatrix Y, IntegerVector colors,
		     ListOf<IntegerVector> W,
		     ListOf<IntegerVector> Sp,
		     double delta) {
  
  const int nvars = Sp.size();
  SparseMatrix<double> M(nvars, nvars);
  const IntegerVector colsize_ = sapply(Sp, get_size);
  const VectorXi colsize = VectorXi::Map(colsize_.begin(), nvars);
  M.reserve(colsize);

  
  for (int i=nvars-1; i>=0; i--) {
    int grp = colors(i)-1; // colors starts at 1
    for (int j : Sp[i]) {
      if (j-1 <= i) {
	double yi = Y(j-1, grp);
	double acc = 0;
	for (int k : W[grp]) {
	  if (k > i) {
	    acc += M.coeff(k-1, i);
	  }
	}
	M.insert(i,j-1) = yi/delta - acc;
      }
    }
  }

  M.makeCompressed();
 
  SparseMatrix<double> H = M.selfadjointView<Lower>();
  return(Rcpp::wrap(H)); 
}
