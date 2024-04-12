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
devtools::install_github("lmangnier/RetroFun-RVS", force=T)
``

## Example 

An use case to illustrate core functions has been provided, please refer to the wiki or the package vignette.  

## Contact 
For any questions, comments or suggestions, please contact Loic Mangnier at loic.mangnier@gmail.com
## Citation 
RetroFun-RVS: a retrospective family-based framework for rare variant analysis incorporating functional annotations
Loic Mangnier, Alexandre Bureau, bioRxiv 2022.06.21.497085; doi: https://doi.org/10.1101/2022.06.21.497085
