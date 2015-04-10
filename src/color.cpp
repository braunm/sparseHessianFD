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
  Eigen::Matrix<bool,Dynamic,1> track_G;
  track_G.resize(nvars);
  track_G.setConstant(true);

  Eigen::Matrix<bool,Dynamic,1> track_S;
  track_S.resize(nvars);
  track_S.setConstant(true);

  VectorXi deg(nvars);
  int r;

  VectorXb nei(nvars);

  Eigen::Matrix<bool,1,Dynamic> Srow(nvars);
  
  // Rcout << "Starting algorithm\n";
  while (track_G.any()) {
    k++;
    //    Rcout << "k = " << k << "\n";
    //  Rcout << "track_G:\n" << track_G << "\n";
    //  Rcout << "GR:\n" << GR << "\n";
    for (int i=0; i<nvars; i++) {
      VectorXb Gi = GR.row(i);
      //   Rcout << "Gi :\n" << Gi << "\n"; 
      deg(i) = Gi.count() * track_G(i);
      //      Rcout << "\ti = " << i << "\tdeg = " << deg(i) << "\n";
    }
    SparseMatrixXbRow SR = GR;
    SparseMatrixXbCol SC = GC;
    track_S = track_G;
    W.push_back(std::set<int>());
    //   Rcout << "track_S:\n";
    while (track_S.any()) {
      check_interrupt();
      //      Rcout << "track_S start:\n" << track_S << "\n\n";
      deg.maxCoeff(&r);
      Srow = SR.row(r);      
      nei = ((Srow * SC).transpose().array() *  track_S.array()).matrix();
      //     Rcout << "\tr = " << r << "\n";
      //   Rcout << "\tDeg\tSrow\ttrack_S\tNei\n";
      //   for (int jj=0; jj<nvars; jj++) {
      //	Rcout << "\t" << deg(jj) << "\t" <<  Srow(jj) << "\t" << track_S(jj) << "\t" << nei(jj) << "\n";
      //  }
      for (int i=0; i<nvars; i++) {
	if (nei(i)) {
	  deg(i) = 0;
	  track_S(i) = false;	    
	}	
      }
      W[k-1].insert(r+1);
      track_G(r) = false;
  
      for (SparseMatrixXbRow::InnerIterator jt(GR, r); jt; ++jt) {
      	jt.valueRef() = false;
      }
      for (SparseMatrixXbCol::InnerIterator it(GC, r); it; ++it) {
      	it.valueRef() = false;
      }
      for (int i=0; i<nvars; i++) {
	GR.coeffRef(i,r) = false;
	GC.coeffRef(r,i) = false;
      }
      //    Rcout << "\ttrack_S:\n" << track_S << "\n\n";
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


