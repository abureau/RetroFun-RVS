#' Function resampling genotypes by fam
#'
#' This functions proceeds to non-parametric bootstrap for each family. Genotypes can be resampled based on
#' sharing probabilities
#' @param agg.genos.by.fam a list return by agg.genos.by.fam
#' @param n.unique.config.by.fam a list with the corresponding genotype configurations in each family
#' @param prob.sharing.by.fam a list with the corresponding sharing probabilities for each genotype configuration
#'
#' @return resample data.frame
#'
#'@export

resample.genos.by.fam = function(agg.genos.by.fam, n.unique.config.by.fam, prob.sharing.by.fam=NULL){

  index_non_null = apply(agg.genos.by.fam$ped_agg[,-1,drop=FALSE],1,function(x) which(x>0))
#  n_non_null = apply(agg.genos.by.fam$ped_agg[,-1],1,function(x) length(which(x>0)))
  n_non_null = sapply(index_non_null,length)
    
  agg_tmp = agg.genos.by.fam
  agg_tmp_ped_agg = agg.genos.by.fam$ped_agg[,-1,drop=FALSE]


  for(x in 1:length(agg_tmp$ped_agg$pedigree)){
    famid = as.character(agg_tmp$ped_agg$pedigree[x])

    if(is.null(prob.sharing.by.fam[[famid]])){
      sample_geno = sample(n.unique.config.by.fam[[famid]],n_non_null[x], replace=T)
    }

    else {
      sample_geno = sample(n.unique.config.by.fam[[famid]],n_non_null[x], replace=T, prob = prob.sharing.by.fam[[famid]])

    }

    agg_tmp_ped_agg[x,index_non_null[[x]]] = sample_geno
  }
  agg_tmp_ped_agg = data.frame("pedigree"=agg_tmp$ped_agg[,1],agg_tmp_ped_agg)
  agg_tmp$ped_agg = agg_tmp_ped_agg

  return(agg_tmp)
}
