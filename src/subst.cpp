#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>




//' @name subst
//' @description Run triangular subsitution algorithm
//' @param Y Matrix of finite differences
//' @param colors Vector of length nvars that identifies color of each variable
//' @param W A list.  Each element represents a color, and contains an integer vector with the indices of the variables with that color.  Indexing is zero-based.
//' @param jCol,ipntr Column indices and row pointers for non-zero elements of lower triangle of Hessian.
//' @param delta Perturbation factor used to compute finite differences.
//' @param nvars Dimension of Hessian (number of variables)
//' @param nnz Number of non-zero elements in the lower triangle of the Hessian.
//' @return A sparse Hessian of class dgCMatrix.
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

