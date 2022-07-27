#' Function resampling genotypes by fam
#'
#' This functions proceeds to non-parametric bootstrap for each family. Genotypes can be resampled based on
#' sharing probabilities
#' @param agg.genos.by.fam a list return by agg.genos.by.fam
#' @param prob.sharing.by.famid a list with the corresponding sharing probabilities for each genotype configurations in families
#'
#' @return resample data.frame
#'
#'@export

resample.genos.by.fam = function(agg.genos.by.fam, prob.sharing.by.famid, probs=NULL){


  index_non_null = apply(agg.genos.by.fam$ped_agg[,-1],1,function(x) which(x>0))
  n_non_null = apply(agg.genos.by.fam$ped_agg[,-1],1,function(x) length(which(x>0)))

  agg_tmp = agg.genos.by.fam
  agg_tmp_ped_agg = agg.genos.by.fam$ped_agg[,-1]


  for(x in 1:length(agg_tmp$ped_agg$pedigree)){
    famid = agg_tmp$ped_agg$pedigree[x]

    if(is.null(probs)){
      sample_geno = sample(1:length(prob.sharing.by.famid[[famid]]),n_non_null[x], replace=T)
    }

    else {
      sample_geno = sample(1:length(prob.sharing.by.famid[[famid]]),n_non_null[x], replace=T, prob = probs)

    }

    agg_tmp_ped_agg[x,index_non_null[[x]]] = sample_geno
  }
  agg_tmp_ped_agg = data.frame("pedigree"=agg_tmp$ped_agg[,1],agg_tmp_ped_agg)
  agg_tmp$ped_agg = agg_tmp_ped_agg

  return(agg_tmp)
}
