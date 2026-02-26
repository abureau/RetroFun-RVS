#' Compute the p-values associated with each functional annotation
#'
#'This function computes both the ACAT-combined and Fisher's method p-values
#'
#'@param null.value.by.fam is a dataframe with four colums (FamID, Expected, Variance and Covariance) returned by the compute.null function
#'@param null.covar.by.indiv is a list of variances and covariances for each pair of affected subjects in every family returned by the compute.null.by.indiv function
#'@return The expected genotype value, variance and covariance for each pedigree within a data.frame
#'@export

replace.covar.by.fam = function(null.value.by.fam,null.covar.by.indiv)
{
  # Calcule la somme des covariances individuelles
  covar.by.fam = sapply(null.covar.by.indiv,function(covar) sum(covar$`Covar(ij,i'j')`))
  names(covar.by.fam) = sapply(null.covar.by.indiv,function(l) l$FamId)
  null.value.by.fam[,"CoVar"] = covar.by.fam[null.value.by.fam$FamID]
  null.value.by.fam
}