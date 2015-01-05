// Rfunc.cpp.  This file is part of sparseHessianFD a contributed package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun.  See LICENSE file for details.


#ifndef __SPARSEHESSIANFD_FUNC__
#define __SPARSEHESSIANFD_FUNC__

#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <except.h>
#include <assert.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;

extern "C"  
{
  void dssm_(int *, int *, int *, int *, int *, int *, int *, int *, int *,
	     int *, int *, int *, int *, int *);
}

extern "C"  
{
  void fdhs_(int *, int *, int *, int *, int *, int *, int *,
	     int *, int *, double *, double *, double *, int *);
}


class Rfunc {

  typedef std::vector<int> ivec;
  typedef std::vector<double> dvec;
  typedef Eigen::Triplet<double> TT;


public:
  
  Rfunc(const int, const Rcpp::Function, const Rcpp::Function);
  ~Rfunc();

  double get_f(const NumericVector&);
  
  NumericVector get_df(const NumericVector&);
  
  Rcpp::List get_fdf(const NumericVector&);
  
  
  Rcpp::S4 get_hessian(const NumericVector&);
  
  
  void hessian_init(const IntegerVector&,
		    const IntegerVector&,
		    int, double);
  
  int get_nnz(); 
  
  
  int nvars; 
  const Rcpp::Function fn;
  const Rcpp::Function gr;
  
private:

  void get_hessian_CSC(const NumericVector&,
		       ivec&,
		       ivec&,
		       dvec&);

  
  
  ivec iRow; // row indices of nonzero elements
  ivec jCol; // col indices of nonzero elements
  ivec listp; // permutation used for DSSM/FDHS
  ivec ngrp; // number of groups for DSSM/FDHS
  ivec ipntr; // for each row, pointer to first element
  ivec jpntr; // for each col, pointer to first element
  dvec fhes; // values for nonzero elements
  dvec fd_eps_vec; // eps used for finite differencing 
  NumericMatrix pert; // perturbation for finite differencing
  NumericMatrix fd; // the finite differences
  int mingrp, maxgrp;
  
  int dssm_info;

  int fd_method;

  int nnz;
  double eps;

  int DSSM_wrap();  // process Hessian structure
  void FDHS_wrap();

  void sort_CSC_cols(ivec&,
		     ivec&,
		     dvec&
		     );
  

  void sort_CSC_cols(ivec&,
		     ivec&
		     ); // for sorting only the indices
  
  void compute_hessian_fd(const NumericVector&);

  NumericVector tmp1;
  NumericVector tmp2;

  ivec irnTmp;
  ivec jclTmp;
  dvec valTmp;
  Eigen::SparseMatrix<double> BkTmp;

 
};

Rfunc::Rfunc(const int nvars_,
	     const Rcpp::Function fn_,
	     const Rcpp::Function gr_) :
  nvars(nvars_), fn(fn_), gr(gr_), nnz(0)
{
   Rcout << "Constructor\n";
  
}

Rfunc::~Rfunc() {
  Rcout << "Destructing...";
  Rcout << "   Complete\n ";
}


int Rfunc::get_nnz() {
  return(nnz);
}


double Rfunc::get_f(const NumericVector& P) {
  
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n",
					 __FILE__, __LINE__);
  NumericVector res = fn(P);
  return(res(0));

}

NumericVector Rfunc::get_df(const NumericVector& P) {
  
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n",
					 __FILE__, __LINE__);
   NumericVector grad = gr(P);
  return(grad);  
}

Rcpp::List Rfunc::get_fdf(const NumericVector& P)
{

  double f = get_f(P);
  NumericVector df = get_df(P);
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("val") = f,
				      Rcpp::Named("grad") = Rcpp::wrap(df)
				      );
    
  return (res);   
}

Rcpp::S4 Rfunc::get_hessian(const NumericVector& P) {


  // Get Hessian using sparse finite differencing from a hessObj
  
  if (fd_method<0) throw MyException("Error:  Hessian is not initialized",
				     __FILE__, __LINE__);

  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n",
					 __FILE__, __LINE__);

 
  get_hessian_CSC(P, irnTmp, jclTmp, valTmp);

  Rcout << "Copying Hessian\n";

  std::vector<TT> Trips;
  Trips.resize(nnz);

  Eigen::SparseMatrix<double> out;
  out.resize(nvars, nvars);
  
  
  // copy hessian to sparse structure elements
  int ind, nels;
  for (int j=0; j<nvars; j++) {
    ind = jclTmp[j];
    nels = jclTmp[j+1] - ind;
    for (int i=0; i<nels; i++) {
      //      BkTmp.coeffRef(irnTmp[ind+i], j) = valTmp[ind+i];
      Trips.push_back(TT(irnTmp[ind+i], j, valTmp[ind+1]));
    }
  }

  out.setFromTriplets(Trips.begin(), Trips.end());
  Rcout << out <<"\n\n";
  
  //  Eigen::SparseMatrix<double> out = BkTmp.selfadjointView<Eigen::Lower>();


  
  Rcout << "Hessian copied\n";
  return(Rcpp::wrap(out));

  
}

/*
Below this point, functions to compute sparse hessian using FD
 */


void Rfunc::hessian_init(const IntegerVector& hess_iRow,
			 const IntegerVector& hess_jCol,
			 int fd_method_, double eps_)
{
			
// copy indices.  iRow and jCol are destroyed during DSSM
    fd_method = fd_method_;
   eps = eps_;
   using std::endl;
   if (fd_method >= 0) {  // use fd_method = -1 for no Hessian
     listp.resize(nvars);
     ngrp.resize(nvars);
     ipntr.resize(nvars+1);
     jpntr.resize(nvars+1);
     
     tmp1 = NumericVector(nvars, 0.0);
     tmp2 = NumericVector(nvars, 0.0);    
     fd_eps_vec.resize(nvars);
     std::fill(fd_eps_vec.begin(), fd_eps_vec.end(), 0.0);
     nnz = hess_iRow.size();
     fhes.resize(nnz);
     std::fill(fhes.begin(), fhes.end(), 0.0);
     iRow.resize(nnz);
     std::copy(hess_iRow.begin(), hess_iRow.end(), iRow.begin());
     jCol.resize(nnz);
     std::copy(hess_jCol.begin(), hess_jCol.end(), jCol.begin());
       dssm_info = DSSM_wrap(); // convert structure information
       if (dssm_info < 0) {
       Rcout << "Problem with hessian structure.  Check column ";
       Rcout  << -dssm_info << "." << endl;
       throw MyException ("Exception thrown. ", __FILE__, __LINE__);
     }
     if (dssm_info == 0) {
       throw MyException ("DSSM_info = 0 (internal problem).",
			  __FILE__, __LINE__);
     }
      pert = NumericMatrix(nvars, maxgrp);
     std::fill(pert.begin(), pert.end(), 0.0);
     fd = NumericMatrix(nvars, maxgrp); // maxgrp is set by DSSM_wrap
     std::fill(fd.begin(), fd.end(), 0.0);
      for (int i=0; i<nvars; i++) {
       pert(i, ngrp[i]-1) = 1.;  // construct perturbation matrix from DSSM results
     }
   }
    
   irnTmp.resize(nnz);
   jclTmp.resize(nvars+1);
   valTmp.resize(nnz);
    
   // BkTmp.resize(nvars, nvars);
   // BkTmp.reserve(nnz);
}
  

void Rfunc::compute_hessian_fd(const NumericVector& P) {
  
  /*
    fd.col = f(x + dx) - f(x)
    pert identifies which groups should be perturbed.  1 for yes and 0 for no.
    each column is a color group.
    For dense, one-column-at-a-time estimation, all elements of pert are zero, except one.
    
    fd is the output matrix, and each row represents the row of the output hessian.
    difference is NOT divided by eps
  */

 if (fd_method<0) throw MyException("Error:  Hessian is not initialized",
				    __FILE__, __LINE__);
  
  tmp1 = get_df(P);  // returns current gradient to tmp1
  std::fill(fd_eps_vec.begin(), fd_eps_vec.end(), 1e-5);
  
  for (int i=0; i<nvars; i++) {
    tmp2(i) = P(i) + fd_eps_vec[i];
    fd_eps_vec[i] = tmp2(i) - P(i);
  }

    
  // It will be worthwhile to create a parallel version of this

  for (int j=0; j<maxgrp; j++) {
    for (int i=0; i<nvars; i++) {
      tmp2(i) = P(i) + fd_eps_vec[i] * pert(i,j);
    }
    fd(Rcpp::_, j) = get_df(tmp2);
    fd(Rcpp::_, j) = fd(Rcpp::_, j) - tmp1;
  }

  FDHS_wrap();  // call FDHS routine
  
}


  void Rfunc::get_hessian_CSC(const NumericVector& P,
			      ivec& irn,
			      ivec& jcl,
			      dvec& vals
			      )
{
  
  compute_hessian_fd(P); // gets CSC format, but unsorted within columns
  
  
  // copy output.  will then be sorted.
  irn = iRow;
  vals = fhes;
  jcl = jpntr;
  
  sort_CSC_cols(irn, jcl, vals);

  // convert to 0-based indexing

  for (int i=0; i<irn.size(); i++) {
    irn[i] = irn[i] - 1;
    jcl[i] = jcl[i] - 1;
  }
  
}


void Rfunc::sort_CSC_cols(ivec& irn,
			  ivec& jcl
			  )
{

   int p0, p1, nels;


  for (int col=0; col<nvars; col++){
    p0 = jcl[col]-1;  // index of first element in column 
    p1 = jcl[col+1]-1; // index of first element in next column
    nels = p1 - p0;
    std::sort(irn.begin()+p0, irn.begin()+p1-1);
  }
}


void Rfunc::sort_CSC_cols(ivec& irn,
			  ivec& jcl,
			  dvec& vals
			  )
{

  
  int row, p0, p1, nels;

  /* Rcout << "irn:\n"; */
  /* for (auto i : irn) { */
  /*   Rcout << i << "\n"; */
  /* } */
  
  
  /* Rcout << "\njcl:\n"; */
  /* for (auto j : jcl) { */
  /*   Rcout << j << "\n"; */
  /* } */

  for (int col=0; col<nvars; col++){
    p0 = jcl[col] - 1;  // index of first element in column 
    p1 = jcl[col+1] - 1; // index of first element in next column
    nels = p1 - p0;
    for (int z=0; z<nels; z++) {
      auto M = std::min_element(irn.begin()+p0+z, irn.begin()+p0+nels);
      row = std::distance(irn.begin()+p0+z, M);
      std::swap(irn[p0+z], irn[p0+z+row]);
      std::swap(vals[p0+z], vals[p0+z+row]);
    }
  }
}


void Rfunc::FDHS_wrap() {

  std::vector<int> iwa(nvars);
  
  int numgrp;
  
  for (numgrp = 1; numgrp <= maxgrp; numgrp++) {
    
    //    double * fhesd_ptr = fd.col(numgrp-1).data();  //assumes fd is col major, and each grp is a column

    double * fhesd_ptr = &fd(0, numgrp-1);
    
    /* fdhs_(&nvars, iRow.data(), jpntr.data(), jCol.data(), ipntr.data(), */
    /* 	  listp.data(), ngrp.data(), &maxgrp, &numgrp,  */
    /* 	  fd_eps_vec.data(), fhesd_ptr, fhes.data(), iwa.data()); */

    fdhs_(&nvars, &iRow[0], &jpntr[0], &jCol[0], &ipntr[0],
	  &listp[0], &ngrp[0], &maxgrp, &numgrp, 
	  &fd_eps_vec[0], fhesd_ptr, &fhes[0], &iwa[0]);
    
  }
  
  return;
  
}

int Rfunc::DSSM_wrap() {
  
  // converts structure information into format needed for FDHS

  int liwa = 6 * nvars; 
  
  std::vector<int> iwa(liwa);
  
  int info=0;
  
  /* dssm_(&nvars, &nnz, iRow.data(), jCol.data(), &fd_method, */
  /* 	listp.data(), ngrp.data(), &maxgrp, &mingrp, */
  /* 	&info, ipntr.data(), jpntr.data(), iwa.data(), &liwa); */

    dssm_(&nvars, &nnz, &iRow[0], &jCol[0], &fd_method,
	&listp[0], &ngrp[0], &maxgrp, &mingrp,
	&info, &ipntr[0], &jpntr[0], &iwa[0], &liwa);
  
  return info;
}


#endif









