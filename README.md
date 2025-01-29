# RetroFun-RVS

## General Description
RetroFun-RVS is a family-based burden test permitting the integration of functional annotations. One of the critical feature of the method is to consider only affected members across families. The package can be used whenever you have standard pedigree structures and data. 

## Package installation

Before installing RetroFunRVS, ACAT and RVS packages used as dependencies should be installed using:
``
BiocManager::install("RVS")
devtools::install_github("yaowuliu/ACAT")
``

RetroFunRVS can be installed using 
``
devtools::install_github("abureau/RetroFun-RVS", force=T)
``

## Example 

An use case to illustrate core functions has been provided, please refer to the wiki or the package vignette.  

## Contact 
For any questions, comments or suggestions, please contact Alexandre Bureau at alexandre.bureau@fmed.ulaval.ca
## Citation 
Mangnier L. et al. (2025) RetroFun-RVS: a retrospective family-based framework for rare variant analysis incorporating functional annotations.
 Genetic Epidemiology 49:e70001, http://dx.doi.org/10.1002/gepi.70001

