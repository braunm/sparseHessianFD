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

typedef std::set<int> S;


//[[Rcpp::export]]
List color(const IntegerVector& ipntr,
	   const IntegerVector& jCols,
	   const IntegerVector& jpntr,
	   const IntegerVector& iRows,
	   const int nvars) {

  //  Rcout << "Starting\n";
  int nnz = jCols.size()-1;
  assert(nnz == iRows.size()-1);
  
  std::vector<std::set<int> > W; // colorings

  S uncolored(nvars);
  std::iota(uncolored.begin(), uncolored.end(), 0);
  S A(nvars);
  std::iota(A.begin(), A.end(), 0);

  VectorXi deg(nvars);
  std::set<int> common(nvars);
  int r;


  //  Eigen::Matrix<bool,1,Dynamic> Srow(nvars);
  
  // Rcout << "Starting algorithm\n";
  int k = 0;
  while (!uncolored.empty()) {
    k++;
    //  Rcout << "k = " << k << "\n";
    W.push_back(std::set<int>());
    A = uncolored; 
    for (auto g = uncolored.begin(); g != uncolored.end(); ) {
      int gstart = ipntr(g);
      int gend = ipntr(g+1);
      S gcols = S(jCols.begin()+gstart; jCols.begin()+gend()); // cols in row g
      S A = uncolored;
      while (!A.empty()) {
	deg.maxCoeff(&r);
	S rcols = S(jCols.begin()+ipntr(r); jCols.begin()+ipntr(r+1));
	for (auto j=A.begin(); j != A.end(); ) {
	  S jrows = S(iRows.begin()+jpntr(*j), iRows.begin()+jpntr(*j +1));
	  S nei(A.size()); // set of neighborhood indices
	  set_intersect(rcols.begin(), rcols.end(),
			jrows.begin(), jrows.end(),
			std::back_inserter(nei));
	  


	}


      }
      

      
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
 
    while (!track_S.empty()) {
      for (auto Si=track_S.begin(); != track_S.end() ;) { // advance interator in loop
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

 SparseMatrixXbRow GR(nvars, nvars);
  SparseMatrixXbCol GC(nvars, nvars);


  GR.setFromTriplets(tripletList.begin(), tripletList.end());
  GR.makeCompressed();
  GC.setFromTriplets(tripletList.begin(), tripletList.end());
  GC.makeCompressed();
  
  typedef Eigen::Triplet<bool> T;
  std::vector<T> tripletList;
  tripletList.reserve(nnz);
  for(int i=0; i<nnz; i++) {
    tripletList.push_back(T(rows(i), cols(i), true));
  }
  

  // Rcout << "Finished algorithm\n";
  
  List res(k);
  // Rcout << "Wrapping\n";
  for (int j=0; j<k; j++) {
    res[j] = W[j];
  }
  
  return(Rcpp::wrap(res)); 
}


