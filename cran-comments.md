##Submission notes for sparseHessianFD, version 0.3.0


### Resubmission notes

- This version is a complete re-write of the package.  Specifically,
  all ACM-copyrighted code has been removed, and replaced with an
  original implementation of the underlying algorithms.  The package
  is now licensed under MPL.  The LICENSE file has been removed.

### Test environments

-  local OS X 10.11.3 install, R 3.2.4, CRAN compiled binary
-  win_builder, both R-release and R-devel

### R CMD check results

There were no ERRORs or WARNINGs.  The only note is related to the
change to the license.

### Compatibility with bayesGDS package

This package is a reverse suggests with bayesGDS.  The tests and
vignettes in bayesGDS 0.6.1 (but not the package functions themselves)
use some sparseHessianFD functions that are deprecated. Thus, the
tests in bayesGDS will generate a warning.  A major update to
bayesGDS is in the works, but it is not yet ready, and will need
sparseHessianFD 0.3.0 anyway.

Once sparseHessianFD 0.3.0 is on CRAN, I can upload a patch to bayesGDS,
but the patched version will not work with the old version of sparseHessianFD.
