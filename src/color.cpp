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
    if (fmod(z,10) == 0) Rcout << "\n";
  }
  Rcout << "\n";
}

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

void checkInterrupt() {
   if( ! R_ToplevelExec(chkIntFn, NULL) )
     Rcpp::stop( "user interuption" );
}


//' @title Vertex coloring of sparse symmetric matrix
//' @description Generate proper coloring of a sparse symmetric matrix.
//' @param pntr,idx row pointers and column indices (CSC or CSR format; same since  matrix is symmetric). Must use zero-based indexing.
//' @param nvars Dimension of matrix (number of variables)
//' @return A list.  Each element of the list represents a color, and contains an integer vector with the indices of the variables with that color.  Indices are zero-based.
//[[Rcpp::export]]
Rcpp::IntegerVector get_colors(const IntegerVector& pntr, //row/col pointer
			 const IntegerVector& idx, // col/row index
			 const int nvars) {
  
  // All pointers and indices are ZERO-BASED
  // and for the FULL MATRICES (not just LT)
  
  std::vector<std::set<int> > P(nvars);
  std::vector<std::set<int> > forb(nvars);
  Rcpp::IntegerVector colors(nvars);
  std::set<int> used;
  std::set<int> valid;

  for (int m=0; m < nvars; m++) {
    P[m] = S(idx.begin()+pntr(m), idx.begin()+pntr(m+1)); // rows
  }

  
  int max_color = 0;
  used.insert(0);
  for (int i=0; i<nvars; i++) { 
    if (forb[i].empty()) {
      colors[i] = 0;
    } else {
      valid.clear();
      set_difference(used.begin(), used.end(),
		     forb[i].begin(), forb[i].end(),
		     std::inserter(valid,valid.begin()));
      if (valid.empty()) { // add new color
	max_color++;
	used.insert(max_color);
	colors[i] = max_color;
      } else {
	colors[i] = *valid.begin();
      }
    }

    for (auto j : P[i]) {
      forb[j].insert(colors[i]);
    }     
  }

  return(Rcpp::wrap(colors));
}



