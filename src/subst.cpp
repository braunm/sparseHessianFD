#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>


//[[Rcpp::export]]
Rcpp::S4 subst(const Rcpp::NumericMatrix& Y,
	       const Rcpp::IntegerVector& colors,
	       const Rcpp::ListOf<Rcpp::IntegerVector>& W,
	       const Rcpp::IntegerVector& jCol,
	       const Rcpp::IntegerVector& ipntr,
	       const double& delta,
	       const int& nvars,
	       const int& nnz) {
  
  const int ngrp = W.size();

  typedef Eigen::Triplet<double> T;
  std::vector<T> Trips;
  Trips.reserve(nnz*2);
  std::vector<double> B(nvars);
  const double inv_delta = 1/delta;
  
  for (int g=0; g<ngrp; g++) {
    std::fill(B.begin(), B.end(), 0.0);
    for (int k = W[g].size()-1; k >= 0; --k){
      int i = W[g](k); 
      for (auto j=jCol.begin()+ipntr(i); j != jCol.begin()+ipntr(i+1); j++) {	
	double z = Y(*j, colors(i)) * inv_delta;
	z -= B[*j];
	B[*j] += z;	
	Trips.emplace_back(i, *j, z);
	if (i != *j) Trips.emplace_back(*j, i, z);	
      }
    }// move up a row
  } // next group and reset accumulator
  
  Eigen::SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());
  M.makeCompressed();
  return(Rcpp::wrap(M)); 
}

