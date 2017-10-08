##Submission notes for sparseHessianFD, version 0.3.3

### Changes from last version

-  Fix registration of native routines, related to changes in Rcpp.

### Test environments

-  local macOS 10.12.6 install
-  R 3.4.2 (CRAN compiled binary) and R 3.4.2 patched (binary from r.research.att.com).
-  win_builder, both R-release and R-devel

### R CMD check results

There were no ERRORs or WARNINGs. The only NOTE flagged the word
"Hessians" in the DESCRIPTION file as a possible mis-spelling. It is not.

