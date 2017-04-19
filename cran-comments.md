##Submission notes for sparseHessianFD, version 0.3.3

- This is a minor update that corresponds to the final version
submitted to the Journal of Statistical Software.

- Explicit registration of native routines, as required by R 3.4.0.

### Test environments

-  local macOS 10.12.4 install
-  R 3.3.3 (CRAN compiled binary) and R 3.4.0 RC  (r72542, binary from r.research.att.com).
-  win_builder, both R-release and R-devel

### R CMD check results

There were no ERRORs or WARNINGs. The only NOTE flagged the word
"Hessians" in the DESCRIPTION file as a possible mis-spelling. It is not.

