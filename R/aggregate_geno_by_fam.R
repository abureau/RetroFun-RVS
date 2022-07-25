usethis::use_package("dplyr")
usethis::use_package("plyr")
#'Aggregate genotypes among affected individuals for each family
#'
#'This function aggregates variants across family members. It takes a pedfile as first parameter.
#'The six first columns must be family id, individual id, father and mother ids, sex and phenotype
#'Homozyguous variants are replaced by their heterozyguous configurations.
#'
#'@param pedfile is a genotype file in ped format: A .ped file
#'@return A list with the ped file corrected and aggregated by family and index each variants observed in families
#'@export

aggregate.geno.by.fam = function(pedfile){
  p = read.table(pedfile, header = F)
  fam = p[,1:6]
  affected = which(fam$V6==2)
  genos = p[,7:ncol(p)]
  
  genos[genos==1] = 0
  genos[genos==2] = 1
  
  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))
  
  df.genos.affected = df.genos[affected,]
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  df.genos.affected[df.genos.affected==2] = 1
  df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected)==0)]
  
  df.genos.affected$pedigree = fam[affected,"V1"]
  
  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum)
  
  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,]
  
  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1])))
  
  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))
  
}
