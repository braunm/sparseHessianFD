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
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;
using Eigen::Dynamic;

typedef Eigen::Matrix<bool, Dynamic, 1> VectorXb;

//[[Rcpp::export]]
List color(const S4& G_,
	   const int& nvars,
	   const int& nnz) {

  SparseMatrix<double> G = Rcpp::as<MappedSparseMatrix<double> >(G_); // copy

  int k = 0;

  std::vector<std::vector<int> > W;
  Eigen::Matrix<bool,Dynamic,1> track_G;
  track_G.resize(nvars);
  track_G.setConstant(true);

  Eigen::Matrix<bool,Dynamic,1> track_S;
  track_S.resize(nvars);
  track_S.setConstant(true);

  VectorXd deg(nvars);
  int r;


  Eigen::Matrix<bool,1,Dynamic> Sr;
  

  while (track_G.any()) {
    k++;
    for (int i=0; i<nvars; i++) {
      deg(i) = G.row(i).sum() * track_G(i);
    }
    auto S = G;
    track_S = track_G;
    W.push_back(std::vector<int>());
    while (track_S.any()) {
      deg.maxCoeff(&r);
      Sr = S.row();
      //    Sr = S.row(r) * S; 
      for (int i=0; i<nvars; i++) {
	bool srs = any(Sr.transpose() && S.col(i));
	if (track_S(i) && srs) {
	  deg(i) = 0;
	  track_S(i) = false;
	}
      }   
      W[k-1].push_back(r);
      track_G(r) = false;
      for (SparseMatrix<double>::InnerIterator jt(G, r); jt; ++jt) {
	jt.valueRef() = 0;
      }
      //      G.col(r).setZero();
    }
  }

  List res(k);
  
  for (int j=0; j<k; j++) {
    res[j] = Rcpp::wrap(W[k]);
  }
 
  return(Rcpp::wrap(res)); 
}


