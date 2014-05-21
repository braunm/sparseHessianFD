
// common_R.hpp -- this file is part of sparseHessianFD, a contributed package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun.  See LICENSE file for details.


#ifndef COMMON_R
#define COMMON_R

#ifndef TRUST_COUT
#define TRUST_COUT Rcpp::Rcout
#endif

#define ERROR_HANDLER R_Interface_Error_Handler

#include <Rcpp.h>
#include <R_ext/Utils.h>
#include "exceptions.hpp"

#define BEGIN_R_INTERFACE try {

  //

#define END_R_INTERFACE  } catch (  const MyException& ex) { \
       ::Rf_error(ex.what());			\
  } catch( const std::exception& __ex__ ) {		\
          forward_exception_to_r( __ex__ );		\
       } catch(...) {				\
    TRUST_COUT << "Unknown error\n";				\
    ::Rf_error( "c++ exception (unknown reason)" );	\
  }  \
  return R_NilValue;

template<typename T>
void R_Interface_Error_Handler(const T & ex) {

  // takes exception object and does R-friendly things to it
  ex.print_message();
  Rf_error("R error\n");

}


static inline void check_interrupt_impl(void* /*dummy*/) {
 R_CheckUserInterrupt();
}

inline bool check_interrupt() {
  return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

#endif
