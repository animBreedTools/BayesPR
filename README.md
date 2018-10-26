# BayesPR
- For a data of 1Kx3K, and for 10,000 iterations, it takes:
  - ~3 minutes for single-trait 
  - ~9 minutes for two-trait analysis
- I think it is quite fast now, but there is room for further improvement. 
  - Possibly by pre-allocating full RHS and LHS, and updating it
