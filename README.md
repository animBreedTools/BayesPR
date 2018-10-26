# BayesPR
- For a data of 1Kx3K, and for 10,000 iterations, it takes:
  - ~3 minutes for single-trait 
  - ~9 minutes for two-trait analysis
- I think it is quite fast now, but there is room for further improvement. 
  - Possibly by pre-allocating full RHS and LHS, and updating it
- When I use a matrix for two traits, BLAS.axpy!() masses up
  - i.e., BLAS.axpy!(a,X,Y[:,t]) in a for loop, where t is the trait, does not work as expected (if I dont mass it up somewhere else in the code) This is why I hardcoded it for two traits
  - i.e., BLAS.axpy!(a,X[:,locus],ycorr1) and BLAS.axpy!(a,X[:,locus],ycorr2), where ycorr is the corrected phenotypes
