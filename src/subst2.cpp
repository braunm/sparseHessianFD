#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>

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
using Eigen::MatrixXd;
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
	   IntegerVector& colsize_,
	   const double& delta) {
  
  const int nvars = Sp.size();
  SparseMatrix<double> M(nvars, nvars);
  Map<VectorXi> colsize = VectorXi::Map(colsize_.begin(), nvars);
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
	  const IntegerVector& jCol_,
	  const IntegerVector& ipntr_,
	  const double& delta,
	  const int& nvars,
	  const int& nnz) {

  //  auto subst_start = std::chrono::steady_clock::now();
  

  const int ngrp = W.size();

  typedef Eigen::Triplet<double> T;
  std::vector<T> Trips;
  Trips.reserve(nnz*2);
  VectorXd B(nvars);
  B.setZero();

  const Map<const VectorXi> jCol = VectorXi::Map(jCol_.begin(), nvars+1);
  const Map<const VectorXi> ipntr = VectorXi::Map(ipntr_.begin(), nnz);
  
  for (int g=0; g<ngrp; g++) {
    B.setZero();
    for (int k = W[g].size()-1; k >= 0; --k){
      int i = W[g](k)-1;
      for (int m = ipntr(i); m<ipntr(i+1); m++) {
	int j = jCol(m-1)-1;
	double z = Y(j, colors(i)-1)/delta - B(j);
	B(j) += z;
	Trips.push_back(T(i, j, z));
	if (i != j) Trips.push_back(T(j, i, z));
      }
    } // move up a row
  } // next group and reset accumulator
  
  SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());
  return(Rcpp::wrap(M)); 
}

