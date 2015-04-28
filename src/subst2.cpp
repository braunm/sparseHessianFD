#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>




//' @title Estimate sparse Hessian
//' @description Estimate Hessian using triangular subsitution algorithm
//' @param Y Matrix of finite differences
//' @param colors Vector of length nvars that identifies color of each variable
//' @param W A list.  Each element represents a color, and contains an integer vector with the indices of the variables with that color.  Indexing is zero-based.
//' @param jCol,ipntr Column indices and row pointers for non-zero elements of lower triangle of Hessian.
//' @param delta Perturbation factor used to compute finite differences.
//' @param nvars Dimension of Hessian (number of variables)
//' @param nnz Number of non-zero elements in the lower triangle of the Hessian.
//' @return A sparse Hessian of class dgCMatrix.
//[[Rcpp::export]]
Rcpp::S4 subst2(const Rcpp::NumericMatrix& Y,
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
  const double inv_delta = 1/delta;
  Rcpp::Rcout << "\n";

 
  // Eigen::MatrixXd B(ngrp, nvars);
  // B.setZero();

  Eigen::MatrixXd H(nvars, nvars);
  H.setZero();
  
  for (int i=nvars-1; i>=0; --i) {
    //   Rcpp::Rcout << "\ni = " << i << "\n";
    int colI = colors(i);
    for (auto j=jCol.begin()+ipntr(i); j != jCol.begin()+ipntr(i+1); j++) {	
      int colJ = colors(*j);
      double yij = Y(i, colJ) * inv_delta;
      double z = yij;
      for (auto k=W[colJ].begin(); k!=W[colJ].end(); ++k) {
      	z -= H(*k, i);
      }
      //   z += B(colJ, i);
      //     Rcpp::Rcout << "\tj = " << *j << "\tz = " << z << "\tyij = " << yij << "\n"; 
      //     B(colJ, *j) += z;
      Trips.emplace_back(i, *j, z);
      H(i, *j) = z;
      //B(colI, *j) += z;
      if (i != *j) {
	Trips.emplace_back(*j, i, z);
      }
    }
  }

  Eigen::SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());
  M.makeCompressed();
 
  
  return(Rcpp::wrap(M)); 
}

