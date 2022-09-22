#' Compute the Variance for q functional annotations
#'This function computes the functional annotation variance statistic. It assumes that the first parameter is
#'the data.frame with null values for each family with the corresponding variance and covariance. The second parameter should be
#'a list with aggregated genotypes by family with related variant index. Statistics are computed based on Z_annot which is a matrix where first column is only ones.
#'Weights are allowed through the W parameter. Haplotype structure can be assumed with the independence parameter.
#'
#'@param null.value.by.fam is a dataframe with four colums (FamID, Expected, Variance and Covariance) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z_annot is a p*q matrix of functional annotations. The first column should be composed only with ones
#'@param W is a vector of weights. Should be a unitary vector when unweighted version is needed.
#'@param independence is a boolean. Default value is FALSE. If variant independence can be assumed, use independence=TRUE.
#'
#'@return a list with each variance by annotation
#'@export
#'
compute.Var.by.Annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W, independence=F){
  list_var_Annot = list()

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
  W_mat = diag(W, nrow=length(W), ncol=length(W))
  Wz = W_mat%*%Z_annot
  Wz_sub = Wz[aggregate.geno.by.fam$index_variants,]

  for(col_A in 1:ncol(Wz_sub)){
    Score_label = paste0("Score", col_A)
    Wz_sub_annot = Wz_sub[,col_A]

    cov_annot = c()

    for(r in 1:nrow(aggregate.geno.by.fam$ped_agg)){
      ped = aggregate.geno.by.fam$ped_agg[r,"pedigree"]

      x = aggregate.geno.by.fam$ped_agg[r,-1]
      diff_x = unique(x[x>0])

      i = which(x>0)
      Wz_sub_annot_fam = Wz_sub_annot[i]

      cov_tmp = NA
      var_tmp = sum(Wz_sub_annot_fam^2 * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"])

      if(independence ==T){
        cov_annot = c(cov_annot,var_tmp)
      }
      else{
        if(length(i)==1){
          cov_tmp = var_tmp } else{
            cw = combn(Wz_sub_annot_fam,2)
            if(length(diff_x) ==1){
              cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw), function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"]
            } else{
              cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw),function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"CoVar"]
            }
          }

        cov_annot = c(cov_annot, cov_tmp)
      }
    }
    list_var_Annot[[Score_label]] = sum(cov_annot)
  }
  return(list_var_Annot)
}
