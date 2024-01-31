#' Function computing bootstrap scores
#'
#' This function computes the burden score on non-parametric bootstrap samples. Genotypes are resampled based on
#' sharing probabilities
#' @param agg.genos.by.fam a list return by agg.genos.by.fam
#' @param n.unique.config.by.fam a list with the corresponding genotype configurations in each family
#' @param Z_annot a p*q matrix of functional annotations. The first column should be composed only with ones
#' @param W a vector of weights. Should be a unitary vector when unweighted version is needed.
#' @param prob.sharing.by.fam a list with the corresponding sharing probabilities for each genotype configuration
#' @param nullval a dataframe with four colums (FamID, Expected, Variance and Covariance) return by compute.null function
#' @param nrep the number of bootstrap samples to generate
#' @return a matrix of scores where each column is an annotation and each row is a bootstrap sample
#'
#'@export
bootRetroFunRVS = function(agg.genos.by.fam, n.unique.config.by.fam, Z_annot=NULL, W=NULL, prob.sharing.by.fam,nullval,nrep=1000)
{
  if(is.null(Z_annot)) Z_annot = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
  if(is.null(W)) W = rep(1, nrow(Z_annot)) 

  score = matrix(NA,nrep,ncol(Z_annot))
  # Loop over the bootstrap replicates
  for (r in 1:nrep)
  {
    ech = resample.genos.by.fam(agg.genos.by.fam, n.unique.config.by.fam=n.unique.config.by.fam, prob.sharing.by.fam=prob.sharing.by.fam)
    score[r,] = compute.Burden.by.Annot(nullval, ech,  Z_annot = Z_annot, W = W)$B
  }
  score
}
