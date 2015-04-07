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
using Eigen::VectorXd;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::SparseMatrix;
using Eigen::Map;

template <typename T> std::string type_name();

int get_size(const IntegerVector& x) {
  return(x.size());
}

//[[Rcpp::export]]
S4 subst_C(const NumericMatrix& Y,
	   const IntegerVector& colors,
	   const ListOf<IntegerVector>& W,
	   const ListOf<IntegerVector>& Sp,
	   const IntegerVector& colsize_,
	   const double& delta) {
  
  const int nvars = Sp.size();
  SparseMatrix<double> M(nvars, nvars);
  const Map<const VectorXi> colsize = VectorXi::Map(colsize_.begin(), nvars);
  M.reserve(colsize);

  
  for (int i=nvars-1; i>=0; i--) {
    int grp = colors(i)-1; // colors starts at 1
    for (int j : Sp[i]) {
      if (j-1 <= i) {
	double acc = 0;
	for (int k : W[grp]) {
	  if (k > i) {
	    acc += M.coeff(k-1, i);
	  }
	}
	M.insert(i,j-1) = Y(j-1,grp)/delta - acc;
      }
    }
  }

  M.makeCompressed();
 
  SparseMatrix<double> H = M.selfadjointView<Lower>();
  return(Rcpp::wrap(H)); 
}


//[[Rcpp::export]]
S4 subst2(const NumericMatrix& Y,
	  const IntegerVector& colors,
	  const ListOf<IntegerVector>& W,
	  const ListOf<IntegerVector>& Sp,
	  const double& delta,
	  const int& nnz) {

  const int nvars = Sp.size();
  const int ngrp = W.size();

  typedef Eigen::Triplet<double> T;
  std::vector<T> Trips;
  Trips.reserve(nnz*2);
  
  VectorXd B(nvars);
  B.setZero();

  for (int g=0; g<ngrp; g++) {
    Rcout << "group = " << g << "\n";
    B.setZero();
    // for (int i : W[g]) { // make this reverse order?
    //   for (int i=W[g].rbegin(); i != W[g].rend(); ++i) { // must be sorted
    for (int k = W[g].size()-1; k >= 0; --k){
      Rcout << "\tk = " << k << "\n";
      int i = W[g](k)-1;
      Rcout << "\ti = " << i << "\n";
      for (int j : Sp[i]) {
	double z = Y(i, colors(j-1))/delta - B(j-1);
	Rcout << "\t\tj = " << j-1 << "\t" << z << "\n";
	B(j-1) += z;
	Trips.push_back(T(i, j-1, z));
	Trips.push_back(T(j-1, i, z));
      }
    } // move up a row
  } // next group and reset accumulator
  
  SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());
  return(Rcpp::wrap(M)); 
}

