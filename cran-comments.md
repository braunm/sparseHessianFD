##Submission notes for sparseHessianFD, version 0.2.0


### Resubmission notes

-  modified the deprecated function new.sparse.hessian.obj to be
   backwards compatible with the bayesGDS package (including
   exporting it).  bayesGDS now passes R CMD check --as-cran.
   Admittedly, the current version of bayesGDS is a mess, and is due
   for a complete overhaul.  The plan was/is to get its dependencies,
   including this package, working correctly first.

-  I understand that a package should be listed in only one of the
   Depends or LinkingTo fields (Section 1.1.3 of Writing R
   Extensions).  But I cannot figure out how to get the package to
   work without having Rcpp in both of them.  Having Rcpp in both
   fields works, and I am not aware of any harm this causes.  Of
   course, I am open to suggestions on how to fix this, if it is even necessary.


### Test environments

-  local OS X 10.10.2 install, R 3.1.2, CRAN compiled binary
-  win_builder, both R-release and R-devel

### R CMD check results
There were no ERRORs or WARNINGs

The NOTE is related to the LICENSE file.  See below.
Possible misspelled words are either proper nouns, acronyms, or names
of other R packages.

### LICENSE
The LICENSE file has not changed from previous versions.  ACM holds the
copyright on the original Fortran code, and allows redistribution
under their standard license.
