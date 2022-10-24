# DTVGM_PML
Functions for hydrological simulation using the DTVGM-PML

# DTVGM_PML.R
R function for DTVGM_PML model

# DTVGM_PML.cpp
C++ function exporting to R using Rcpp package  
using the following code to load DTVGM_PML function  
```
  library(Rcpp)
  sourceCpp("./DTVGM_PML.cpp")
  DTVGM_PML(...)
```

# DTVGM_0.1.0.tar.gz
R package for DTVGM_PML model  
using the following code to install "DTVGM" package  
```
  install.packages("./DTVGM_0.1.0.tar.gz", repo=NULL, type=”source”)
  library(DTVGM)
  DTVGM_PML(...)
```
