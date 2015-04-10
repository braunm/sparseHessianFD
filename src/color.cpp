#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <common_R.hpp>

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
typedef Eigen::SparseMatrix<bool, Eigen::RowMajor> SparseMatrixXbRow;
typedef Eigen::SparseMatrix<bool, Eigen::ColMajor> SparseMatrixXbCol;

//[[Rcpp::export]]
List color(const IntegerVector& rows,
	   const IntegerVector& cols,
	   const int& nvars) {

  //  Rcout << "Starting\n";
  int nnz = rows.size();

  SparseMatrixXbRow GR(nvars, nvars);
  SparseMatrixXbCol GC(nvars, nvars);

  typedef Eigen::Triplet<bool> T;
  std::vector<T> tripletList;
  tripletList.reserve(nnz);
  for(int i=0; i<nnz; i++) {
    tripletList.push_back(T(rows(i), cols(i), true));
  }
  
  GR.setFromTriplets(tripletList.begin(), tripletList.end());
  GR.makeCompressed();
  GC.setFromTriplets(tripletList.begin(), tripletList.end());
  GC.makeCompressed();
  
  
  //  Rcout << "GR:\n" << GR << "\n";
  // Rcout << "GC:\n" << GC << "\n";

  //  Rcout << "Sparse matrix is filled\n";
  
  int k = 0;

  std::vector<std::set<int> > W;

  std::set<int> track_G(nvars);
  std::iota(track_G.begin(), track_G.end(), 0);
  std::set<int> track_S(nvars);
  std::iota(track_S.begin(), track_S.end(), 0);

  VectorXi deg(nvars);
  std::set<int> common(nvars);
  int r;

  //  const std::vector<int> nm(nvars);
  //std::iota(nm.begin(), nm.end(), 0);
  //  VectorXb nei(nvars);

  Eigen::Matrix<bool,1,Dynamic> Srow(nvars);
  
  // Rcout << "Starting algorithm\n";
  while (!track_G.empty()) {
    k++;
    check_interrupt();
    //  Rcout << "k = " << k << "\n";
    //  Rcout << "track_G:\n" << track_G << "\n";
  
    for (int i=0; i<nvars; i++) {
      VectorXb Gi = GR.row(i);
      std::set<int> Gi; // set of elements in row i (from pointer)
      std::set_intersection(track_G.begin(), track_G.end(),
			    Gi.begin(), Gi.end(),
			    std::back_inserter(common));
      deg(i) = common.size();
    }

    std::vector<std::set> SC = GC;
    std::vector<std::set> SR = GR;
    track_S = track_G;
    W.push_back(std::set<int>());
    while (!track_S.empty()) {
      check_interrupt();
      deg.maxCoeff(&r);
      nei = any_common(SR[r], SC) %intersect% track_S;
      for (int i : nei) {
	  deg(i) = 0;
	  track_S.erase(i); // remove neighbors of r from S
	}	
      }
      W[k-1].insert(r+1);
      track_G.erase(r); // remove r from uncolored set
  
      for (int i=0; i<nvars; i++) {
      	GR[i].erase(r); // zero out column r in GR
      	GC[i].erase(r); // zero out row r in GC 
      }
    }
  }

  // Rcout << "Finished algorithm\n";
  
  List res(k);
  // Rcout << "Wrapping\n";
  for (int j=0; j<k; j++) {
    res[j] = W[j];
  }
  
  return(Rcpp::wrap(res)); 
}


