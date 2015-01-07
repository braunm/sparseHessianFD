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
using Rcpp::as;
using Rcpp::wrap;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::MatrixBase;

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

  template<typename TP, typename TG>
    void get_df_(const MatrixBase<TP>&, const MatrixBase<TG>&);
  
private:
  
  
  VectorXi iRow; // row indices of nonzero elements
  VectorXi jCol; // col indices of nonzero elements
  VectorXi listp; // permutation used for DSSM/FDHS
  VectorXi ngrp; // number of groups for DSSM/FDHS
  VectorXi ipntr; // for each row, pointer to first element
  VectorXi jpntr; // for each col, pointer to first element
  VectorXd fhes; // values for nonzero elements
  VectorXd fd_eps_vec; // eps used for finite differencing 
  MatrixXd pert; // perturbation for finite differencing
  MatrixXd fd; // the finite differences
  int mingrp, maxgrp, dssm_info, fd_method;

  int nnz;
  double eps;

  int DSSM_wrap();  // process Hessian structure
  void FDHS_wrap();
  
  void compute_hessian_fd(const NumericVector&);

  VectorXi hess_iRow;
  VectorXi hess_jCol;

  NumericVector tmp1, tmp2; 
};

Rfunc::Rfunc(const int nvars_,
	     const Rcpp::Function fn_,
	     const Rcpp::Function gr_) :
  nvars(nvars_), fn(fn_), gr(gr_), nnz(0)
{}

Rfunc::~Rfunc() {
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

template<typename TP, typename TG>
  void Rfunc::get_df_(const MatrixBase<TP>& P_,
		     const MatrixBase<TG>& out_) {

  NumericVector P(nvars);
  std::copy(P_.derived().data(), P_.derived().data() + nvars, P.begin());
  
  MatrixBase<TG>& out = const_cast<MatrixBase<TG>&>(out_);

  
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n",
					 __FILE__, __LINE__);
  NumericVector grad = gr(P);

  out = as<VectorXd>(grad);
  
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


void Rfunc::hessian_init(const IntegerVector& hess_iRow_,
			 const IntegerVector& hess_jCol_,
			 int fd_method_, double eps_) 
{


  hess_iRow = Rcpp::as<VectorXi>(hess_iRow_);
  hess_jCol = Rcpp::as<VectorXi>(hess_jCol_);
			
  // copy indices.  iRow and jCol are destroyed during DSSM
  fd_method = fd_method_;
  eps = eps_;

  
  if (fd_method >= 0) {  // use fd_method = -1 for no Hessian
    listp.resize(nvars);
    ngrp.resize(nvars);
    ipntr.resize(nvars+1);
    jpntr.resize(nvars+1);
    
    tmp1 = NumericVector(nvars, 0.0);
    tmp2 = NumericVector(nvars, 0.0);    
    fd_eps_vec.resize(nvars);
    fd_eps_vec.setConstant(eps);
    
    nnz = hess_iRow.size();
    fhes.setZero(nnz);

    iRow = hess_iRow;
    jCol = hess_jCol;

    dssm_info = DSSM_wrap(); // convert structure information
    if (dssm_info < 0) {
      Rcout << "Problem with hessian structure.  Check column ";
      Rcout  << -dssm_info << ".\n";
      throw MyException ("Exception thrown. ", __FILE__, __LINE__);
    }
    if (dssm_info == 0) {
      throw MyException ("DSSM_info = 0 (internal problem).",
			 __FILE__, __LINE__);
    }
    Rcout << "dssm_info = " << dssm_info << "\n";
    
    pert.setZero(nvars, maxgrp);
    fd.setZero(nvars, maxgrp);

    for (int i=0; i<nvars; i++) {
      pert(i, ngrp(i)-1) = 1.;  // construct perturbation matrix from DSSM results
    }

    /* int ind, nels; */
    /* for (int j=0; j<nvars; j++) { */
    /*   ind = jpntr[j]-1; */
    /*   nels = jpntr[j+1] -1 - ind; */
    /*   for (int i=0; i<nels; i++) { */
    /* 	jCol[ind+i] = j; */
    /*   } */
    /* } */
  }
}


Rcpp::S4 Rfunc::get_hessian(const NumericVector& P) {


  if (fd_method<0) throw MyException("Error:  Hessian is not initialized",
				     __FILE__, __LINE__);

  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n",
					 __FILE__, __LINE__);

 Rcout << "Computing Hessian\n";
  
 compute_hessian_fd(P);

 Rcout << "iRow, jCol, fhes\n";
 for (int i=0; i<nnz; i++) {
   Rcout << iRow(i) << "\t" << jCol(i) << "\t" << fhes(i) << "\n";
 }

 
 std::vector<TT> Trips;
 Trips.resize(nnz);
 Eigen::SparseMatrix<double> sp;
 sp.resize(nvars, nvars);
 
 Rcout << "Copying Hessian\n";
 
 // copy hessian to sparse structure elements
 int ind, nels;
 for (int j=0; j<nvars; j++) { 
   ind = jpntr(j)-1; 
   nels = jpntr(j+1) -1 - ind;
   Rcout << "j = " << j << "  ind = " << ind << "  nels = " << nels << "\n";
   for (int i=0; i<nels; i++) {
     Rcout <<  "\ti = " << i << "  iRow = " << iRow(ind+i) << "  jCol = " << jCol(ind+i) << "  fhes = " << fhes(ind+i) << "\n";
     Trips.push_back(TT(iRow(ind+i)-1, j, fhes(ind+i))); 
   } 
 } 
 

 /* for (int k=0; k<nnz; k++) { */
 /*   Trips.push_back(TT(iRow(k)-1, jCol(k)-1, fhes(k))); */
 /* } */
 
 sp.setFromTriplets(Trips.begin(), Trips.end());

 Rcout << sp << "\n";
 
 Eigen::SparseMatrix<double> out = sp.selfadjointView<Eigen::Lower>();
 
 Rcout << "Returning Hessian\n";
 
 return(Rcpp::wrap(out));
}

/*
  Below this point, functions to compute sparse hessian using FD
 */





void Rfunc::compute_hessian_fd(const NumericVector& P) {
  
  /*
    fd.col = f(x + dx) - f(x)
    pert identifies which groups should be perturbed.  1 for yes and 0 for no.
    each column is a color group.
    For dense, one-column-at-a-time estimation, all elements of pert are zero, except one.
    
    fd is the output matrix, and each row represents the row of the output hessian.
    difference is NOT divided by eps
  */

  VectorXd vars = as<VectorXd>(P);
  
  Rcout << "eps = " << eps << "\n";
  
  if (fd_method<0) throw MyException("Error:  Hessian is not initialized",
				     __FILE__, __LINE__);
  
  VectorXd tmp1(nvars);
  get_df_(vars, tmp1);  // returns current gradient to tmp1
  fd_eps_vec.setConstant(eps);


  VectorXd tmp2 = vars + fd_eps_vec;
  fd_eps_vec = tmp2 - vars;
  
  /* for (int i=0; i<nvars; i++) { */
  /*   tmp2(i) = P(i) + fd_eps_vec(i); */
  /*   fd_eps_vec(i) = tmp2(i) - P(i); */
  /* } */

 

  Eigen::PermutationMatrix<Eigen::Dynamic> tmpP(listp - VectorXi::Constant(nvars, 1));
  Eigen::PermutationMatrix<Eigen::Dynamic> tmpPinv = tmpP.inverse();
  MatrixXd PP = tmpP.toDenseMatrix().cast<double>();
  MatrixXd PPinv = tmpPinv.toDenseMatrix().cast<double>();
  Rcout << "\nPP:\n" << tmpP.indices() << "\n\n";
  Rcout << "\nPPinv:\n" << tmpPinv.indices() << "\n\n";

  Rcout << "\npert:\n" << pert << "\n\n";
  
  
    
  // It will be worthwhile to create a parallel version of this

  for (int j=0; j<maxgrp; j++) {
    /* for (int i=0; i<nvars; i++) { */
    /*   tmp2(i) = P(i) + fd_eps_vec[i] * pert(i,j); */
    /* } */
    //   tmp2 = P + fd_eps_vec * pert(Rcpp::_, j);
    VectorXd tmp3(nvars);
    tmp2 = vars + eps * pert.col(j);
    get_df_(tmp2, tmp3);
    fd.col(j) = tmp3 - tmp1;
    
    /* NumericVector v2 = wrap(tmp2); */
    /* VectorXd tmp3 = as<VectorXd>(get_df(v2)); */
    
    /* fd.col(j) = Rcpp::as<VectorXd>(get_df(wrap(tmp2))); */
    /* fd.col(j) -= Rcpp::as<VectorXd>(tmp1); */
  }

  //  fd = PPinv * fd; // Permutation

  FDHS_wrap();  // call FDHS routine
  
  
  Rcout << "listp, ngrp:\n";
  for (int i=0; i<nvars; i++) {
    Rcout << listp(i) << "\t" << ngrp(i) << "\n";
  }


  MatrixXd fdtmp = (fd.array() / eps).matrix();
  Rcout << "fdtmp:\n" << fdtmp << "\n\n";


  /* MatrixXd pat(nvars, nvars);; */
  /* pat.setZero(); */
  /* for (int i=0; i<nnz; i++){ */
  /*   pat(hess_iRow(i)-1, hess_jCol(i)-1) = 1.; */
  /* } */

  /* Rcout << "Pattern:\n" << pat << "\n\n"; */
  

  
  MatrixXd HH = ((pert * fdtmp.transpose()).array()).matrix();
  HH = HH;

  Rcout << "HH:\n" << HH << "\n\n";


  
}

void Rfunc::FDHS_wrap() {

  VectorXi iwa(nvars);  
  int numgrp;
  
  for (numgrp = 1; numgrp <= maxgrp; numgrp++) {

    VectorXd htmp(nvars);
    htmp = fd.col(numgrp-1);
    
    //   double * fhesd_ptr = &fd(0, numgrp-1);    
    fdhs_(&nvars, iRow.data(), jpntr.data(), jCol.data(), ipntr.data(),
	  listp.data(), ngrp.data(), &maxgrp, &numgrp, 
	  fd_eps_vec.data(), htmp.data(), fhes.data(), iwa.data());    
  }
  
  return;
  
}

int Rfunc::DSSM_wrap() {
  
  // converts structure information into format needed for FDHS

  int liwa = 6 * nvars;   
  VectorXi iwa(liwa);  
  int info = 0;

  dssm_(&nvars, &nnz, iRow.data(), jCol.data(), &fd_method,
	listp.data(), ngrp.data(), &maxgrp, &mingrp,
	&info, ipntr.data(), jpntr.data(), iwa.data(), &liwa);


  /* Rcout << "After DSSM (iRow, jCol)\n"; */
  /* for (int i=0; i<nnz; i++) { */
  /*   Rcout << iRow(i) << "\t" << jCol(i) << "\n"; */
  /* } */
  /* Rcout << "\n(jpntr, ipntr):\n"; */
  /* for (int i=0; i<nvars+1; i++) { */
  /*   Rcout << jpntr(i) << "\t" << ipntr(i) << "\n\n"; */
  /* } */

  
  return info;
}

#endif









