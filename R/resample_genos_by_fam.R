#' Function resampling genotypes by fam
#'
#' This functions proceeds to non-parametric bootstrap for each family. Genotypes can be resampled based on
#' sharing probabilities
#' @param aggregate.geno.by.fam a list return by aggregate.geno.by.fam
#' @param prob_sharing_by_famid a list with the corresponding sharing probabilities for each genotype configurations in families
#'
#' @return resample data.frame
#'
#'@export

resample.genos.by.fam = function(aggregate.geno.by.fam, prob_sharing_by_famid=NULL){


  index_non_null = apply(aggregate.geno.by.fam$ped_agg[,-1],1,function(x) which(x>0))
  n_non_null = apply(aggregate.geno.by.fam$ped_agg[,-1],1,function(x) length(which(x>0)))

  agg_tmp = aggregate.geno.by.fam
  agg_tmp_ped_agg = aggregate.geno.by.fam$ped_agg[,-1]


  for(x in 1:length(agg_tmp$ped_agg$pedigree)){
    famid = agg_tmp$ped_agg$pedigree[x]

    if(is.null(resample.genos.by.fam)){
      sample_geno = sample(1:length(prob_sharing_by_famid[[famid]]),n_non_null[x], replace=T)
    }

    else{
      sample_geno = sample(1:length(prob_sharing_by_famid[[famid]]),n_non_null[x], replace=T, prob = prob_sharing_by_famid)

    }

    agg_tmp_ped_agg[x,sample(1:ncol(agg_tmp_ped_agg),n_non_null[x])] = sample_geno
  }
  agg_tmp_ped_agg = data.frame("pedigree"=agg_tmp$ped_agg[,1],agg_tmp_ped_agg)
  agg_tmp$ped_agg = agg_tmp_ped_agg

  return(agg_tmp)
}
