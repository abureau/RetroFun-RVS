#' Compute the p-values associated with each functional annotation
#'
#'This function computes both the ACAT-combined and Fisher's method p-values
#'
#'@param null.value.by.fam is a dataframe with four colums (FamID, Expected, Variance and Covariance) returned by the compute.null function
#'@param null.covar.by.indiv is a list of variances and covariances for each pair of affected subjects in every family returned by the compute.null.by.indiv function
#'@return The expected genotype value, variance and covariance for each pedigree within a data.frame
#'@export

update.covar.by.fam = function(null.value.by.fam,null.covar.by.indiv)
{
  # Calcule 2 x la somme des covariances individuelles
  #covar.by.fam = sapply(null.covar.by.indiv,function(covar) 2*sum(covar$`Covar(ij,i'j')`[lower.tri(covar$`Covar(ij,i'j')`)])+sum(diag(covar$`Covar(ij,i'j')`)))
  # Cette version prend 2 fois la somme de la matrice et soustrait la somme de la diagonale qu'il faut compter une seule fois.
  # Ça évite d'extraire la diagonale inférieure.
  covar.by.fam = sapply(null.covar.by.indiv,function(covar) 2*sum(covar$`Covar(ij,i'j')`)-sum(diag(covar$`Covar(ij,i'j')`)))
  names(covar.by.fam) = sapply(null.covar.by.indiv,function(l) l$FamId)
  null.value.by.fam[,"CoVar"] = covar.by.fam[null.value.by.fam$FamID]
  null.value.by.fam
}