#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>


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



//[[Rcpp::export]]
Rcpp::S4 subst(const NumericMatrix& Y,
	       const IntegerVector& colors,
	       const ListOf<IntegerVector>& W,
	       const IntegerVector& jCol_,
	       const IntegerVector& ipntr_,
	       const double& delta,
	       const int& nvars,
	       const int& nnz) {
  
  const int ngrp = W.size();

  typedef Eigen::Triplet<double> T;
  std::vector<T> Trips;
  Trips.reserve(nnz*2);
  VectorXd B(nvars);
  B.setZero();

  const double inv_delta = 1/delta;
  const Map<const VectorXi> jCol = VectorXi::Map(jCol_.begin(), nnz);
  const Map<const VectorXi> ipntr = VectorXi::Map(ipntr_.begin(), nvars+1);
  
  for (int g=0; g<ngrp; g++) {
    B.setZero();
    for (int k = W[g].size()-1; k >= 0; --k){
      int i = W[g](k);
      for (int m = ipntr(i); m<ipntr(i+1); m++) {
	int j = jCol(m);
	double z = Y(j, colors(i)) * inv_delta;
	z -= B(j);
	B(j) += z;
	Trips.emplace_back(i, j, z);
	//	if (i != j) Trips.emplace_back(j, i, z);
	
      }
    } // move up a row
  } // next group and reset accumulator
  
  SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());
  M.makeCompressed();
  return(Rcpp::wrap(M)); 
}

