
#include <func.h>



RCPP_MODULE(sparseHessianFD){  
  Rcpp::class_< Rfunc >("sparseHessianFD")
    
    .constructor<const int, const Rcpp::Function, const Rcpp::Function>()
    
    .method( "fn", & Rfunc::get_f, "function value")
    .method( "gr", & Rfunc::get_df, "gradient")
    .method( "fngr", & Rfunc::get_fdf, "list containing function value and gradient")
    .method( "hessian", & Rfunc::get_hessian, "sparse Hessian as dgCMatrix object")
    .method( "hessian.init", & Rfunc::hessian_init, "used internally to initialize object with sparsity pattern")
    .method( "nnz", & Rfunc::get_nnz, "number of non-zero elements in lower triangle")
    .method( "nvars", & Rfunc::get_nvars, "length of parameter vector")
    ;
}
