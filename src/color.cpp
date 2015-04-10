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

template<typename T>
void print_set(const std::set<T>& X) {

  int z = 0;
  for (auto i : X) {
    Rcout << i << "\t";
    z++;
    if (fmod(z,5) == 0) Rcout << "\n";
  }
  Rcout << "\n";

}

//[[Rcpp::export]]
List color(const IntegerVector& pntr, //row/col pointer
	   const IntegerVector& idx, // col/row index
	   const int nvars) {

  // All pointers and indices are ZERO-BASED
  // and for the FUUL MATRICES (not just LT)
  //  Rcout << "Starting\n";
  int nnz = idx.size();

  
  std::vector<std::set<int> > W; // colorings

  std::vector<int> ix(nvars);
  std::iota(ix.begin(), ix.end(), 0);
  
  S uncolored(ix.begin(), ix.end());
  VectorXi deg(nvars);
  int r;


  //  Eigen::Matrix<bool,1,Dynamic> Srow(nvars);
  
  

  Rcout << "Computing initial degrees\n";
  std::vector<std::set<int> > P(nvars);  
  for (int m=0; m < nvars; m++) {
    P[m] = S(idx.begin()+pntr(m), idx.begin()+pntr(m+1));
    deg(m) = P[m].size() - 1;
  }

  Rcout << "Starting algorithm\n";
  int k = 0;
  while (!uncolored.empty() && (k < 6)) {
    Rcout << "k = " << k << ".  Uncolored:\n";
    print_set(uncolored);
    W.push_back(std::set<int>()); // initialize bin for color k    
    S A(uncolored); // all uncolored are candidates
    while (!A.empty()) { // while there are still candidates
      Rcout << "Candidates-start:\n";
      print_set(A);
      deg.maxCoeff(&r); // chose r,  candidate with highest degree
      W[k].insert(r+1); // put r in color k (returning 1-based indices)
      S rcols = S(idx.begin()+pntr(r), idx.begin()+pntr(r+1)); // cols in row r
      Rcout << "\tr = " << r << "\n";
      S nei;
      for (auto j : A) { // fill the neighborhood
	Rcout << "\tj = " << j << "\n";
	S jrows = S(idx.begin()+pntr(j), idx.begin()+pntr(j+1)); // rows in col j
	S in_nei;
	in_nei.clear();
	set_intersection(rcols.begin(), rcols.end(),
			 jrows.begin(), jrows.end(),
			 std::inserter(in_nei, in_nei.begin()));
	Rcout << "in_nei:\n";
	print_set(in_nei);
	//	Rcout << "\tIs " << j << " in neighborhood? " << in_nei.empty() << "\n";
	if (!in_nei.empty()) nei.insert(j);
      }
      for (auto j : nei) {
	  deg(j) = 0;
	  A.erase(j); // remove from candidate set and move to next candidate
      } 
    }
    // remove r as neighbor for everyone else, and recompute degrees
    deg.setZero();
    for (auto i : uncolored) {
      for (auto rr : W[k]) {
	P[i].erase(rr);
      }
    }
    for (int i : uncolored) {
      deg[i] = P[i].size(); // self has already been removed
    }
    Rcout << "Candidates-end:\n";
    print_set(A);
    uncolored.erase(r);     
    k++; // advance to next color
  } // end loop on uncolored
  Rcout << "Finished algorithm\n";

  List res(k);
  for (int j=0; j<k; j++) {
    res[j] = W[j];
  }
  Rcout << "Returning\n";
  return(Rcpp::wrap(res));  
}





      
//       VectorXb Gi = GR.row(i);
//       std::set<int> Gi; // set of elements in row i (from pointer)
//       std::set_intersection(track_G.begin(), track_G.end(),
// 			    Gi.begin(), Gi.end(),
// 			    std::back_inserter(common));
//       deg(i) = common.size();
//     }

//     std::vector<std::set> SC = GC;
//     std::vector<std::set> SR = GR;
//     track_S = track_G;
 
//     while (!track_S.empty()) {
//       for (auto Si=track_S.begin(); != track_S.end() ;) { // advance interator in loop
// 	deg.maxCoeff(&r);
// 	nei = any_common(SR[r], SC) %intersect% track_S;
// 	for (int i : nei) {
// 	  deg(i) = 0;
// 	  track_S.erase(i); // remove neighbors of r from S
// 	}	
//       }
//       W[k-1].insert(r+1);
//       track_G.erase(r); // remove r from uncolored set
      
//       for (int i=0; i<nvars; i++) {
//       	GR[i].erase(r); // zero out column r in GR
//       	GC[i].erase(r); // zero out row r in GC 
//       }
//     }
//   }

//  SparseMatrixXbRow GR(nvars, nvars);
//   SparseMatrixXbCol GC(nvars, nvars);


//   GR.setFromTriplets(tripletList.begin(), tripletList.end());
//   GR.makeCompressed();
//   GC.setFromTriplets(tripletList.begin(), tripletList.end());
//   GC.makeCompressed();
  
//   typedef Eigen::Triplet<bool> T;
//   std::vector<T> tripletList;
//   tripletList.reserve(nnz);
//   for(int i=0; i<nnz; i++) {
//     tripletList.push_back(T(rows(i), cols(i), true));
//   }
  

 
// }


