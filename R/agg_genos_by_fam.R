# usethis::use_package("dplyr")
# usethis::use_package("plyr")

#'Aggregate genotypes among affected individuals for each family
#'
#'This function aggregates variants across family members. It takes a pedfile as first parameter.
#'The six first columns must be family id, individual id, father and mother ids, sex and phenotype
#'Homozyguous variants are replaced by their heterozyguous configurations.
#'
#'@param pedfile.path is a genotype file in ped format: A path to a .ped file
#'@param pedfile a data.frame corresponding to genotype by family
#'@param Z_annot is a p*q matrix of functional annotations. The first column should be composed only with ones
#'@param exclude_annot vector of indices of the columns of the Z_annot matrix to exclude when computing 
#'the maximal annotation for each variant (default = 1)
#'@param correction a string corresponding to the applied correction, none when no correction is applied
#'replace when homozyguous configurations are replaced by their corresponding heterozyguous configurations, remove when variants with
#'at least one homozyguous configurations
#'@return A list with the ped file corrected and aggregated by family and index each variants observed in families
#'@export

agg.genos.by.fam = function(pedfile.path=NULL, pedfile=NULL, Z_annot=NULL, exclude_annot=1, correction=c("none","replace","remove")){
  if(is.null(pedfile) & !is.null(pedfile.path)){
    p = read.table(pedfile.path, header = FALSE)
  }
  
  else if(!is.null(pedfile) & is.null(pedfile.path)){
    p=pedfile
  }
  
  
  fam = p[,1:6]
  affected = which(fam[,6]==2)
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
  if (!is.null(Z_annot))
  {
  if (ncol(Z_annot)>1)
  {
    dupvar = which(duplicated(t(df.genos.affected)))
    treated = logical(ncol(df.genos.affected))
    
    #Replace NA by zero in Z
    Z_annot[is.na(Z_annot)] = 0
    
    var_to_keep = numeric(0)
    for (i in dupvar)
    {
      # Si le variant n'a pas déjà été traité
      if (!treated[i])
      {
        # Déterminer les copies de ce variant dupliqué
        copies = which(apply(df.genos.affected==df.genos.affected[,i],2,all))
        # Calculer l'annotation maximale par variant (on exclut certaines annotations comme l'ordonnée à l'origine)
        Zmax.copies = apply(Z_annot[copies,-exclude_annot, drop=FALSE],1,max)
        # On garde le variant avec l'annotation maximale (ou le 1er variant avec l'annotation maximale s'il y en a plusieurs)
        var_to_keep = c(var_to_keep,copies[which.max(Zmax.copies)[1]])
        treated[copies] = TRUE
      }
    }
    logical_to_keep = (1:ncol(df.genos.affected))%in%var_to_keep
    # Conserver les variants retenus et ceux non impliqués dans une duplication
    df.genos.affected = df.genos.affected[,logical_to_keep|(!treated)]
  }
    # Sinon il y a une seule annotation qui est une colonne de 1, donc pas de sélection basée sur les annotations
    else df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  }
  # Sinon Z_annot est nulle, donc pas de sélection basée sur les annotations
  else df.genos.affected = data.frame(t(unique(t(df.genos.affected))))


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

  attributes(df.genos.agg.by.fam)$correction = correction
  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))

}


