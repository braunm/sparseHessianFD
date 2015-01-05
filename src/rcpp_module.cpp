
#include <func.h>

/* Not exposing functions related to sparse Hessian, because the Hessian is dense */

RCPP_MODULE(sparseHessianFD){  
  Rcpp::class_< Rfunc >("sparseHessianFD")
    
    .constructor<const int, const Rcpp::Function, const Rcpp::Function>()
    
    .method( "fn", & Rfunc::get_f)
    .method( "gr", & Rfunc::get_df)
    .method( "fngr", & Rfunc::get_fdf)
    .method( "hessian", & Rfunc::get_hessian)
    .method( "hessian.init", & Rfunc::hessian_init)
    .method( "get.nnz", & Rfunc::get_nnz)
    ;
}
