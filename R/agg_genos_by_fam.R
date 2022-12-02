# usethis::use_package("dplyr")
# usethis::use_package("plyr")

#'Aggregate genotypes among affected individuals for each family
#'
#'This function aggregates variants across family members. It takes a pedfile as first parameter.
#'The six first columns must be family id, individual id, father and mother ids, sex and phenotype
#'Homozyguous variants are replaced by their heterozyguous configurations.
#'
#'@param pedfile is a genotype file in ped format: A path to a .ped file
#'@param correction a string corresponding to the applied correction, none when no correction is applied
#'replace when homozyguous configurations are replaced by their corresponding heterozyguous configurations, remove when variants with
#'at least one homozyguous configurations
#'@return A list with the ped file corrected and aggregated by family and index each variants observed in families
#'@export

agg.genos.by.fam = function(pedfile, correction=c("none","replace","remove")){
  p = read.table(pedfile, header = FALSE)
  fam = p[,1:6]
  affected = which(fam$V6==2)
  genos = p[,7:ncol(p)]

  genos[genos==1] = 0
  genos[genos==2] = 1

  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))

  colnames(df.genos) = paste0("X", 1:ncol(df.genos))

  #Keeping only affected individuals
  df.genos.affected = df.genos[affected,,drop=FALSE]

  #Keeping only one variants per haplotype
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))


  #Option 1: No correction applied
  if(correction == "none"){
    df.genos.affected = df.genos.affected
  }

  #Option 2: Homozyguous --> Heterozyguous
  else if(correction == "replace"){
    df.genos.affected[df.genos.affected==2] = 1
  }

  #Option 3: Homozyguous are remove
  else if(correction == "remove"){
    df.genos.affected = df.genos.affected %>% dplyr::select(-which(colSums(df.genos.affected==2, na.rm = TRUE)>0))
  }

  else{
    stop("Please provide a correct value for the correction parameter...")
  }


  df.genos.affected = df.genos.affected %>% dplyr::select(-which(colSums(df.genos.affected, na.rm = TRUE)==0))
  df.genos.affected = df.genos.affected %>% dplyr::select(-which(colSums(is.na(df.genos.affected))==nrow(df.genos.affected)))

  df.genos.affected$pedigree = fam[affected,"V1"]

  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum, na.rm=TRUE ,na.action = NULL)

  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1,drop=FALSE])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,,drop=FALSE]

  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1,drop=FALSE])))

  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))

}
