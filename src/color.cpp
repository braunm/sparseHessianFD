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
List color(const SparseMatrix<double>& G,
	   const int& nvars,
	   const int& nnz) {

  SparseMatrix<double> G = as<MappedSparseMatrix<double> >(G_); // copy

  int k = 0;

  std::vector<VectorXi> W;
  Eigen::Matrix<bool,Dynamic,1> track_G;
  track_G.resize(nvars);
  track_G.setConstant(true);

  Eigen::Matrix<bool,Dynamic,1> track_S;
  track_S.resize(nvars);
  track_S.setConstant(true);

  VectorXd deg(nvars);
  int r;
  VectorXd Sv(nvars);

  Eigen::Matrix<bool,Dynamic,1> Sr;
  

  while (track_G.any()) {
    k++;
    deg = G.rowwise().sum().array() * track_G.array();
    S = G;
    track_S = track_G;
    W.push_back();
    while (track_S.any()) {
      deg.maxCoeff(&r);
      Sr = S.row(r) * S;
      VectorXb nei = trackS && Sr; // varying length
      for (int i : nei) {
	deg(i) = 0;
	track_S(i) = false;
      }
      W[k-1].push_back(r);
      track_G(r) = false;
      G.col(r).setZero();
    }
  }

  ListOf<IntegerVector> res(k);
  for (int j=0; j<k; j++) {
    res(j).clone(W[k]);
  }
 
  return(Rcpp::wrap(res)); 
}


