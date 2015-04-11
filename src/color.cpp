#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <common_R.hpp>
#include <chrono>

using Rcpp::List;
using Rcpp::IntegerVector;
using Rcpp::Rcout;
using Eigen::VectorXi;


typedef std::set<int> S;

template<typename T>
void print_container(const T& X) {
  int z = 0;
  for (auto i : X) {
    Rcout << i << "   ";
    z++;
    if (fmod(z,7) == 0) Rcout << "\n";
  }
  Rcout << "\n";
}

template<typename T>
void print_vec(const Eigen::MatrixBase<T>& X) {
  int z = 0;
  int n = X.size();
  for (int i=0; i<n; i++) {
    Rcout << X(i) << "   ";
    z++;
    if (fmod(z,7) == 0) Rcout << "\n";
  }
  Rcout << "\n";
}

//[[Rcpp::export]]
List color(const IntegerVector& pntr, //row/col pointer
	   const IntegerVector& idx, // col/row index
	   const int nvars) {

  using namespace std::chrono;
  typedef std::chrono::duration<double> Dur;
  
  //  auto start = high_resolution_clock::now();
  
  // All pointers and indices are ZERO-BASED
  // and for the FULL MATRICES (not just LT)
  //  Rcout << "Starting\n";
  std::vector<std::set<int> > W; // colorings
  
  S uncolored;
  VectorXi deg(nvars);
  int r;

  std::vector<std::set<int> > P(nvars);

  for (int m=0; m < nvars; m++) {
    uncolored.insert(m);
    P[m] = S(idx.begin()+pntr(m), idx.begin()+pntr(m+1));
    deg(m) = P[m].size();
  }

  //  auto init = high_resolution_clock::now();

  S A, rcols, jrows, in_nei, Wk;
  int k = 0;
  while (!uncolored.empty()) {
    Wk.clear();
    A = S(uncolored); // all uncolored are candidates
    while (!A.empty()) { // while there are still candidates
      deg.maxCoeff(&r); // chose r,  candidate with highest degree
      Wk.insert(r); // put r in color k
      uncolored.erase(r); // remove r from uncolored      
      rcols = S(P[r].begin(), P[r].end());
      for (auto j = A.begin(); j != A.end(); ) {
	jrows = S(idx.begin()+pntr(*j), idx.begin()+pntr(*j+1)); // rows in col j
	in_nei.clear();
	set_intersection(rcols.begin(), rcols.end(),
			 jrows.begin(), jrows.end(),
			 std::inserter(in_nei, in_nei.begin()));
	if (!in_nei.empty()) {
	  deg(*j) = 0;
	  j = A.erase(j);
	} else {
	  j++;
	}
      }
      P[r].clear();
    } // until no more candidates
    // remove r as neighbor for everyone else, and recompute degrees
    deg.setZero();
    for (auto i : uncolored) {
      for (auto rr : Wk) P[i].erase(rr);       
      deg[i] = P[i].size();      
    }
    W.push_back(Wk);
    k++; // advance to next color
  } // end loop on uncolored

  //  auto alg = high_resolution_clock::now();
  
  List res(k);
  for (int j=0; j<k; j++) {
    res[j] = W[j]; // returning ZERO-BASED INDEXES
  }

  // auto wrapper = high_resolution_clock::now();

  // Dur diff1 = init - start;
  // Dur diff2 = alg - init;
  // Dur diff3 = wrapper - alg;

  // Rcout << "Diff1 = " << diff1.count() << "s \n";
  // Rcout << "Diff2 = " << diff2.count() << "s \n";
  // Rcout << "Diff3 = " << diff3.count() << "s \n";

  return(Rcpp::wrap(res));  
}



