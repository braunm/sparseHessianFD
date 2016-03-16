// Rinterface.cpp.  Part of the sparseHessianFD package for the R programming language.
//
// Copyright (C) 2012 Michael Braun.  See LICENSE file for details.


#ifndef __SPARSE_HESSIAN_FD
#define __SPARSE_HESSIAN_FD

#include <RcppEigen.h>
#include <common_R.hpp>

inline bool my_ret_bool(bool x) {return(x);}

#define my_assert(x) do { \
    if(!my_ret_bool(x) ) throw MyException(EIGEN_MAKESTRING(x), __FILE__, __LINE__); \
} \
while (false);

#ifdef eigen_assert
#undef eigen_assert
#endif
#define eigen_assert(x) my_assert(x)

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Rfunc.cpp>


using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;
using Rcpp::Function;
using Rcpp::as;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::SparseMatrix;


RcppExport SEXP sparseHessianObj__new(SEXP nvars_, SEXP fn_, SEXP gr_, 
				      SEXP hs_, SEXP fd_method_, SEXP eps_) {
  
  int nvars = as<int>(nvars_);
  int fd_method = as<int>(fd_method_);
  double fd_eps = as<double>(eps_);

  Function fn(fn_);
  Function gr(gr_);

  List hs = as<List>(hs_);

  IntegerVector iRow_((SEXP)hs["iRow"]);
  IntegerVector jCol_((SEXP)hs["jCol"]);
 
  int nnz = iRow_.size();

  Map<VectorXi> iRow(iRow_.begin(),nnz);
  Map<VectorXi> jCol(jCol_.begin(),nnz);

  Rfunc* func = new Rfunc(nvars, fn, gr);

  func->hessian_init(iRow, jCol, fd_method, fd_eps);

  Rcpp::XPtr<Rfunc> ptr(func, true);

  return ptr;

}



RcppExport SEXP sparseHessianObj__fn(SEXP ptr_, SEXP x_)
{
  BEGIN_R_INTERFACE

  Rcpp::XPtr<Rfunc> func(ptr_);
  NumericVector x2(x_);
  int nvars = x2.size();

  Map<VectorXd> x(x2.begin(),nvars);

  double val;

  func->get_f(x, val);

  return(Rcpp::wrap(val));

  END_R_INTERFACE
}


RcppExport SEXP sparseHessianObj__gr(SEXP ptr_, SEXP x_)
{
  BEGIN_R_INTERFACE

  Rcpp::XPtr<Rfunc> func(ptr_);
  NumericVector x2(x_);
  int nvars = x2.size();

  Map<VectorXd> x(x2.begin(),nvars);

  double val;
  VectorXd grad(nvars);

  func->get_fdf(x, val, grad);

  return(Rcpp::wrap(grad));

  END_R_INTERFACE
}


RcppExport SEXP sparseHessianObj__hessian(SEXP ptr_, SEXP x_)
{
  BEGIN_R_INTERFACE

  Rcpp::XPtr<Rfunc> func(ptr_);
  NumericVector x2(x_);
  int nvars = x2.size();
  int nnz = func->get_nnz();

  Map<VectorXd> x(x2.begin(),nvars);

 
  SparseMatrix<double> hess(nvars, nvars);
  hess.reserve(nnz);

  func->get_hessian(x, hess);

  return(Rcpp::wrap(hess));

  END_R_INTERFACE
}



RcppExport SEXP sparseHessianObj__fngr(SEXP ptr_, SEXP x_)
{
  BEGIN_R_INTERFACE

  Rcpp::XPtr<Rfunc> func(ptr_);
  NumericVector x2(x_);
  int nvars = x2.size();

  Map<VectorXd> x(x2.begin(),nvars);

  double val;
  VectorXd grad(nvars);

  func->get_fdf(x, val, grad);

  List res = List::create(Rcpp::Named("fn") = val,
			  Rcpp::Named("gr") = Rcpp::wrap(grad)
			  );

  return(res);

  END_R_INTERFACE
}




RcppExport SEXP sparseHessianObj__all(SEXP ptr_, SEXP x_)
{
  BEGIN_R_INTERFACE

  Rcpp::XPtr<Rfunc> func(ptr_);
  NumericVector x2(x_);
  int nvars = x2.size();
  int nnz = func->get_nnz();

  Map<VectorXd> x(x2.begin(),nvars);

  double val;
  VectorXd grad(nvars);
  SparseMatrix<double> hess(nvars, nvars);
  hess.reserve(nnz);

  func->get_fdf(x, val, grad);
  func->get_hessian(x, hess);
  
  List res = List::create(Rcpp::Named("fn") = val,
			  Rcpp::Named("gr") = Rcpp::wrap(grad),
			  Rcpp::Named("hs") = Rcpp::wrap(hess)
			  );

  return(res);

  END_R_INTERFACE
}



#endif




