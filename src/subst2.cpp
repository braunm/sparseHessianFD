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
	  const ListOf<IntegerVector>& Sp,
	  const double& delta,
	  const int& nnz) {

  //  auto subst_start = std::chrono::steady_clock::now();
  
  const int nvars = Sp.size();
  const int ngrp = W.size();

  typedef Eigen::Triplet<double> T;
  std::vector<T> Trips;
  Trips.reserve(nnz*2);
  VectorXd B(nvars);
  B.setZero();
  
  //  auto reserve_trips = std::chrono::steady_clock::now();
  
  for (int g=0; g<ngrp; g++) {
    B.setZero();
    // for (int i : W[g]) { // make this reverse order?
    //   for (int i=W[g].rbegin(); i != W[g].rend(); ++i) { // must be sorted
    for (int k = W[g].size()-1; k >= 0; --k){
      int i = W[g](k)-1;
      for (int j : Sp[i]) {
	double z = Y(j-1, colors(i)-1)/delta - B(j-1);
	B(j-1) += z;
	Trips.push_back(T(i, j-1, z));
	if (i != j-1) Trips.push_back(T(j-1, i, z));
      }
    } // move up a row
  } // next group and reset accumulator

  //  auto make_sparse = std::chrono::steady_clock::now();
  
  SparseMatrix<double> M(nvars, nvars);
  M.setFromTriplets(Trips.begin(), Trips.end());

  // auto wrapping = std::chrono::steady_clock::now();
  
  S4 res = Rcpp::wrap(M);
  
  // auto subst_end = std::chrono::steady_clock::now();

  // std::chrono::duration<double> diff1 = reserve_trips - subst_start;
  // std::chrono::duration<double> diff2 = make_sparse - reserve_trips;
  // std::chrono::duration<double> diff3 = wrapping - make_sparse;
  // std::chrono::duration<double> diff4 = subst_end - wrapping;
  // std::chrono::duration<double> diff5 = subst_end - subst_start;

  // Rcout << "setup = " << diff1.count() << " s\n";
  // Rcout << "subst = " << diff2.count() << " s\n";
  // Rcout << "make_sparse = " << diff3.count() << " s\n";
  // Rcout << "wrapping = " << diff4.count() << " s\n";
  // Rcout << "total = " << diff5.count() << " s\n";

  
  return(res); 
}

