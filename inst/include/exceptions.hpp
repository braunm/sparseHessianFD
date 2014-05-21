// exceptions.hpp.  Part of the sparseHessianFD package for the R programming language.
//
// Copyright (C) 2013 Michael Braun.  See LICENSE file for details.

#ifndef MY_EXCEPTION
#define MY_EXCEPTION


#include <iostream>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <cstring>

#include <common_R.hpp>

class MyException : public std::exception {

public:

  const std::string reason;
  const std::string file;
  const int line;
  
  std::string message;
  

 MyException(const std::string reason_,
	     const std::string file_,
	      const int line_) :
    reason(reason_), file(file_), line(line_)
  {
    
    std::ostringstream oss;
	
	oss << "\nException thrown from File " << file << "  at Line " << line <<".\n";
	oss << "Reason : " << reason << "\n";
	
	message = oss.str();

  }

  virtual ~MyException() throw() {};

  virtual const char* what() const throw()
  {
    return message.c_str();
 
  }

  void print_message(){

    TRUST_COUT << message << std::endl;
  }


};




#endif
