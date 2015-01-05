// Rfunc.cpp.  This file is part of sparseHessianFD a contributed package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun.  See LICENSE file for details.


#ifndef __SPARSEHESSIANFD_RFUNC__
#define __SPARSEHESSIANFD_RFUNC__

#include <RcppEigen.h>
#include <common_R.hpp>

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Lower;

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


public:
  
  Rfunc(const int, const Rcpp::Function, const Rcpp::Function);
  
  ~Rfunc();
  
  template <typename Tpars>
  void get_f(const Eigen::MatrixBase<Tpars>&, const double&);
  
  template <typename Tpars, typename Tgrad>
  void get_df(const Eigen::MatrixBase<Tpars>&, const Eigen::MatrixBase<Tgrad>&);
  
  template <typename Tpars, typename Tgrad>
  void get_fdf(const Eigen::MatrixBase<Tpars>&, const double&, const Eigen::MatrixBase<Tgrad>&);
  
  
  template<typename Tpars, typename Tout>
  void get_hessian(const Eigen::MatrixBase<Tpars>&, const Eigen::SparseMatrixBase<Tout>&);
  
  
  template<typename TP, typename TX, typename Tout>
  void get_hessian_CSC(const MatrixBase<TP>&,
		       const MatrixBase<TX>&,
		       const MatrixBase<TX>&,
		       const MatrixBase<Tout>&);
  
  
  template<typename Tin>
  void hessian_init(const MatrixBase<Tin>&,
		    const MatrixBase<Tin>&,
		    int, double);
  
  int get_nnz(); 
  
  
  int nvars; 
  const Rcpp::Function fn;
  const Rcpp::Function gr;
  
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
  int mingrp, maxgrp;
  
  int dssm_info, fdhs_info;

  int fd_method;

  int nnz;
  double eps;

  int DSSM_wrap();  // process Hessian structure
  void FDHS_wrap();

  template<typename TX, typename Tout>
    void sort_CSC_cols(const MatrixBase<TX>&,
		       const MatrixBase<TX>&,
		       const MatrixBase<Tout>&
		       );

  template<typename TX>
    void sort_CSC_cols(const MatrixBase<TX>&,
		       const MatrixBase<TX>&
		       ); // for sorting only the indices

  void compute_hessian_fd(const VectorXd&);

  VectorXd tmp1;
  VectorXd tmp2;

  VectorXi irnTmp;
  VectorXi jclTmp;
  VectorXd valTmp;
  SparseMatrix<double> BkTmp;

 
};

Rfunc::Rfunc(const int nvars_,
	     const Rcpp::Function fn_,
	     const Rcpp::Function gr_) :
  nvars(nvars_), fn(fn_), gr(gr_), nnz(0)
{
}


int Rfunc::get_nnz() {
  return(nnz);
}


template<typename Tpars>
void Rfunc::get_f(const MatrixBase<Tpars>& P_, const double& f_) {
  
  Eigen::MatrixBase<Tpars>& P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n", __FILE__, __LINE__);
 
  double &f = const_cast<double&>(f_);

  Rcpp::NumericVector pars(P.derived().data(), P.derived().data() + P.derived().size());
  
  double res = Rcpp::as<double>(fn(pars));
  f = res;
  return;

}


template<typename Tpars, typename Tgrad>
void Rfunc::get_df(const MatrixBase<Tpars>& P_, const MatrixBase<Tgrad>& df_) {
  

  using Rcpp::NumericVector;

  Eigen::MatrixBase<Tpars>& P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
  Eigen::MatrixBase<Tgrad> & df = const_cast<Eigen::MatrixBase<Tgrad>& >(df_);
  
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n", __FILE__, __LINE__);
  if (df.size()!=nvars) throw MyException("Incorrect gradient length\n", __FILE__, __LINE__);
  
  NumericVector pars(P.derived().data(), P.derived().data() + P.size());
  
  NumericVector grad_  = gr(pars);
  
  VectorXd grad = VectorXd::Map(grad_.begin(), nvars);
  
  df = grad;
  
  return;  
}

template<typename Tpars, typename Tgrad>
void Rfunc::get_fdf(const Eigen::MatrixBase<Tpars>& P_, const double& f_,
		    const Eigen::MatrixBase<Tgrad>& df_)
{

  get_f(P_, f_);
  get_df(P_, df_);
  return;    
  
}

template<typename Tpars, typename Tout>
void Rfunc::get_hessian(const Eigen::MatrixBase<Tpars>& P_,
			const Eigen::SparseMatrixBase<Tout>& out_) {

  // Get Hessian using sparse finite differencing from a hessObj
  
  if (fd_method<0) throw MyException("Error:  Hessian is not initialized", __FILE__, __LINE__);

  Eigen::SparseMatrixBase<Tout> & out = const_cast<Eigen::SparseMatrixBase<Tout>& >(out_);
  
  VectorXd P = P_;

  get_hessian_CSC(P, irnTmp, jclTmp, valTmp);

  // copy hessian to sparse structure elements
  int ind, nels;
  for (int j=0; j<nvars; j++) {
    ind = jclTmp(j);
    nels = jclTmp(j+1) - ind;
    for (int i=0; i<nels; i++) {
      BkTmp.coeffRef(irnTmp(ind+i),j) = valTmp(ind+i);
    }
  }

  out = BkTmp.selfadjointView<Lower>();

  
}

/*
Below this point, functions to compute sparse hessian using FD
 */


template<typename Tin>
void Rfunc::hessian_init(const MatrixBase<Tin>& hess_iRow,
			 const MatrixBase<Tin>& hess_jCol,
			 int fd_method_, double eps_)
{
			
// copy indices.  iRow and jCol are destroyed during DSSM


  fd_method = fd_method_;
  eps = eps_;

  using std::endl;
  if (fd_method>=0) {  // use fd_method = -1 for no Hessian
    
    listp.setZero(nvars);
    ngrp.setZero(nvars);
    ipntr.setZero(nvars+1);
    jpntr.setZero(nvars+1);
    tmp1.setZero(nvars);
    tmp2.setZero(nvars);

    fd_eps_vec.resize(nvars);

    nnz = hess_iRow.size();

    fhes.setZero(nnz);

    iRow = hess_iRow;
    jCol = hess_jCol;
    
    dssm_info = DSSM_wrap(); // convert structure information
    if (dssm_info < 0) {
      TRUST_COUT << "Problem with hessian structure.  Check column " << -dssm_info << "." << endl;
      //    throw MyException ("Exception thrown. ", __FILE__, __LINE__);
    }
    if (dssm_info == 0) {
      throw MyException ("DSSM_info = 0 (internal problem).", __FILE__, __LINE__);
    }
        
    pert.setZero(nvars, maxgrp);
    fd.setZero(nvars, maxgrp); // maxgrp is set by DSSM_wrap
    
    for (int i=0; i<nvars; i++) {
      pert(i,ngrp(i)-1) = 1.;  // construct perturbation matrix from DSSM results
    }
  }

  irnTmp.setZero(nnz);
  jclTmp.setZero(nvars+1);
  valTmp.setZero(nnz);

  BkTmp.resize(nvars, nvars);
  BkTmp.reserve(nnz);
}
  

void Rfunc::compute_hessian_fd(const VectorXd& P) {
  
  /*
    fd.col = f(x + dx) - f(x)
    pert identifies which groups should be perturbed.  1 for yes and 0 for no.
    each column is a color group.
    For dense, one-column-at-a-time estimation, all elements of pert are zero, except one.
    
    fd is the output matrix, and each row represents the row of the output hessian.
    difference is NOT divided by eps
  */

 if (fd_method<0) throw MyException("Error:  Hessian is not initialized", __FILE__, __LINE__);
  
  get_df(P, tmp1);  // returns current gradient to tmp1
  fd_eps_vec.setConstant(eps);

  tmp2 = P + fd_eps_vec;
  fd_eps_vec = tmp2.col(0) - P;
  
  // It will be worthwhile to create a parallel version of this

  for (int i=0; i<maxgrp; i++) {
    
    tmp2.array() = P.array() + fd_eps_vec.array()*pert.col(i).array();
    get_df(tmp2, fd.col(i));
    fd.col(i) -= tmp1;
  }
  
  FDHS_wrap();  // call FDHS routine
  
}

template<typename TP, typename TX, typename Tout>
  void Rfunc::get_hessian_CSC(const MatrixBase<TP>& P_,
					    const MatrixBase<TX>& irn_,
					    const MatrixBase<TX>& jcl_,
					    const MatrixBase<Tout>& vals_
					    )
{
  
  VectorXd P = P_;  // copy P for perturbations

  MatrixBase<TX>& irn = const_cast<MatrixBase<TX>&>(irn_);
  MatrixBase<TX>& jcl = const_cast<MatrixBase<TX>&>(jcl_);
  MatrixBase<Tout>& vals = const_cast<MatrixBase<Tout>&>(vals_);
  
  compute_hessian_fd(P); // gets CSC format, but unsorted within columns
  
  
  // copy output.  will then be sorted.
  irn = iRow;
  vals = fhes;
  jcl = jpntr;

  sort_CSC_cols(irn, jcl, vals);

  // convert to 0-based indexing

  irn.array() -= 1;
  jcl.array() -= 1;
}


template<typename TX>
void Rfunc::sort_CSC_cols(const MatrixBase<TX>& irn_,
					const MatrixBase<TX>& jcl_
					)
{
  
  MatrixBase<TX>& irn = const_cast<MatrixBase<TX>&>(irn_);
  MatrixBase<TX>& jcl = const_cast<MatrixBase<TX>&>(jcl_);
  
  
  int row, p0, p1, mc, nels;


  for (int col=0; col<nvars; col++){
    p0 = jcl(col)-1;  // index of first element in column 
    p1 = jcl(col+1)-1; // index of first element in next column
    nels = p1 - p0;
    for (int z=0; z<nels; z++) {
      irn.segment(p0+z,nels-z).minCoeff(&row, &mc);
      std::swap(irn(p0+z),irn(p0+z+row));
    }
  }
}

template<typename TX, typename Tout>
  void Rfunc::sort_CSC_cols(const MatrixBase<TX>& irn_,
					  const MatrixBase<TX>& jcl_,
					  const MatrixBase<Tout>& vals_
					  )
{
  
  MatrixBase<TX>& irn = const_cast<MatrixBase<TX>&>(irn_);
  MatrixBase<TX>& jcl = const_cast<MatrixBase<TX>&>(jcl_);
  MatrixBase<Tout>& vals = const_cast<MatrixBase<Tout>&>(vals_);
  
  int row, p0, p1, mc, nels;
  
  for (int col=0; col<nvars; col++){
    p0 = jcl(col)-1;  // index of first element in column 
    p1 = jcl(col+1)-1; // index of first element in next column
    nels = p1 - p0;
    for (int z=0; z<nels; z++) {
      irn.segment(p0+z,nels-z).minCoeff(&row, &mc);
      std::swap(irn(p0+z),irn(p0+z+row));
      std::swap(vals(p0+z),vals(p0+z+row));  
    }
  }
}

void Rfunc::FDHS_wrap() {
  
  
  VectorXi iwa(nvars);
  
  int numgrp;
  
  for (numgrp=1; numgrp<=maxgrp; numgrp++) {
    
    double * fhesd_ptr = fd.col(numgrp-1).data();  //assumes fd is col major, and each grp is a column
    
    fdhs_(&nvars, iRow.data(), jpntr.data(), jCol.data(), ipntr.data(),
	  listp.data(), ngrp.data(), &maxgrp, &numgrp, 
	  fd_eps_vec.data(), fhesd_ptr, fhes.data(), iwa.data());
  }
  
  return;
  
}

int Rfunc::DSSM_wrap() {
  
  // converts structure information into format needed for FDHS
  
  int liwa = 6*nvars; 
  
  VectorXi iwa(liwa);
  
  int info;
  
  dssm_(&nvars, &nnz, iRow.data(), jCol.data(), &fd_method,
	listp.data(), ngrp.data(), &maxgrp, &mingrp,
	&info, ipntr.data(), jpntr.data(), iwa.data(), &liwa);
  
  return info;
}

Rfunc::~Rfunc(){
}

#endif




