#include <Rcpp.h>


using Rcpp::List;
using Rcpp::IntegerVector;
using Rcpp::Rcout;


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

/*
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
*/

//' @title Color graph of sparse Hessian
//' @description Generate valid cyclic coloring of variables that is consistent with estimating sparse Hessian with a triangle substitution method.
//' @param pntr,idx row pointers and column indices (CSC or CSR format; same since Hessian matrix is symmetric). Must use zero-based indexing.
//' @param nvars Dimension of Hessian (number of variables)
//' @param LT true (default) if ordering with degree in lower triangle, and false if using full Hessian.  Included only for testing, to see if it matters (I don't think it does).
//' @return A list.  Each element of the list represents a color, and contains an integer vector with the indices of the variables with that color.  Indices are zero-based.
//[[Rcpp::export]]
List color_graph(const IntegerVector& pntr, //row/col pointer
		 const IntegerVector& idx, // col/row index
		 const int nvars,
		 const bool LT = true) {
  
  // All pointers and indices are ZERO-BASED
  // and for the FULL MATRICES (not just LT)

  std::vector<std::set<int> > W; // colorings
  S uncolored;
  std::vector<int> deg(nvars);
  
  std::vector<std::set<int> > P(nvars), jrows(nvars);

  for (int m=0; m < nvars; m++) {
    uncolored.insert(m);
    P[m] = S(idx.begin()+pntr(m), idx.begin()+pntr(m+1)); // rows
    jrows[m] = P[m];
    if (LT) {
      deg[m] = std::distance(P[m].begin(), P[m].upper_bound(m));
    } else {
      deg[m] = P[m].size();
    }
  }

  int k = 0;
  while (!uncolored.empty()) {
    int r;
    S Wk;
    S A = S(uncolored); // all uncolored are candidates

    while (!A.empty()) { // while there are still candidates
      r = std::distance(deg.begin(), max_element(deg.begin(), deg.end()));
      Wk.insert(r); // put r in color k
      uncolored.erase(r); // remove r from uncolored      
      for (auto j = A.begin(); j != A.end(); ) {
	S in_nei;
	set_intersection(P[r].begin(), P[r].end(),
			 jrows[*j].begin(), jrows[*j].end(),
			 std::inserter(in_nei, in_nei.begin()));	
	if (!in_nei.empty()) {
	  deg[*j] = 0;
	  j = A.erase(j);
	} else {
	  ++j;
	}
      }
      P[r].clear();      
    } // until no more candidates
    
    // remove r as neighbor for everyone else, and recompute degrees
    fill(deg.begin(), deg.end(), 0);
    for (auto i : uncolored) {
      for (auto rr : Wk) {
	P[i].erase(rr);
	if (LT) {
	  deg[i] = std::distance(P[i].begin(), P[i].upper_bound(i));
	} else {
	  deg[i] = P[i].size();
	}
      }
    }
    W.push_back(Wk);
    k++; // advance to next color
  } // end loop on uncolored
  
  List res(k);
  for (int j=0; j<k; j++) {
    res[j] = W[j]; // returning ZERO-BASED INDEXES
  }
  
  return(Rcpp::wrap(res));  
}



