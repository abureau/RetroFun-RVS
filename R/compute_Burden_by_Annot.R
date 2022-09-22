#' Compute the Burden Statistic for q functional annotations
#'
#'This function computes the functional annotation burden statistic. It assumes that the first parameter is
#'the data.frame with null values for each family with the corresponding variance and covariance. The second parameter should be
#'a list with aggregated genotypes by family with related variant index. Statistics are computed based on Z_annot which is a matrix where first column is only ones.
#'Weights are allowed through the W parameter
#'
#'@param null.value.by.fam is a dataframe with four colums (FamID, Expected, Variance and Covariance) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z_annot is a p*q matrix of functional annotations. The first column should be composed only with ones
#'@param W is a vector of weights. Should be a unitary vector when unweighted version is needed.
#'
#'@return a list with each score by annotation
#'@export
#'
compute.Burden.by.Annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W){

  #Replace NA by zero in Z
  Z_annot[is.na(Z_annot)] = 0

  convert.neg.to.pos = function(col){
    if(min(col)<0) col = col - min(col)
    else col
  }

  if(ncol(Z_annot) > 1){
    #Make negative scores positive
    Z_annot = apply(Z_annot, 2, convert.neg.to.pos)
  } else{
    if(min(Z_annot)< 0) Z_annot = Z_annot-min(Z_annot)
    else Z_annot = Z_annot
  }

  ped_agg = aggregate.geno.by.fam$ped_agg
  
  W_mat = diag(W, nrow=length(W), ncol=length(W))
  W_sub = W_mat[aggregate.geno.by.fam$index_variants,aggregate.geno.by.fam$index_variants]
  
  Expected = null.value.by.fam[,c("FamID", "Expected")]

  split_G_agg_by_fam = split(ped_agg, 1:nrow(ped_agg))

  diff_obs_expected_by_fam = lapply(split_G_agg_by_fam, function(x) {
    pedigree = x$pedigree

    x_tmp = x[,-1]
    x_tmp[x_tmp==0] = NA
    doe = x_tmp - Expected[Expected$FamID==pedigree,"Expected"]
  })

  diff_obs_expected_all_fam = do.call("rbind", diff_obs_expected_by_fam)
  diff_obs_expected_all_fam[is.na(diff_obs_expected_all_fam)] = 0

  S_by_var = colSums(diff_obs_expected_all_fam)
  Wz_sub = W_sub%*%Z_annot[aggregate.geno.by.fam$index_variants,]


  S_Wz = S_by_var%*%Wz_sub
  Burden_by_annot = S_Wz^2

  return(list("B"=Burden_by_annot))

}
