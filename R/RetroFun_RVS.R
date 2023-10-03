#' Compute the p-values associated with each functional annotation
#'
#'This function computes both the ACAT-combined and Fisher's method p-values
#'
#'@param null.value.by.fam is a dataframe with four colums (FamID, Expected, Variance and Covariance) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z_annot is a p*q matrix of functional annotations. The first column should be composed only with ones
#'@param W is a vector of weights. Should be a unit vector when an unweighted approach is needed.
#'@param independence is a boolean. Default value is FALSE. If variant independence can be assumed, use independence=TRUE.
#'@param contrib.by.fam a boolean. If TRUE the function returns the weighted statistic by family
#'@param contrib.by.var a boolean. If TRUE the function returns the weighted statistic by variant
#'@return a data.frame of p-values for each score by annotation. ACAT and Fisher combined p-values are computed.
#'@export

RetroFun.RVS = function(null.value.by.fam, aggregate.geno.by.fam, Z_annot, W, independence=FALSE, contrib.by.fam=FALSE, contrib.by.var=FALSE){

  # if(.check.parameters.agg.null(null.value.by.fam, aggregate.geno.by.fam)%in%c(1,2)){
  #   warning("Correction parameters are not consistent with options provided for cryptic relatedness or consanguinity, please check ...")
  # }

  run.Burden = compute.Burden.by.Annot(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W, contrib.by.fam = contrib.by.fam, contrib.by.var = contrib.by.var)
  Burden.Stat = run.Burden$B
  Var.Stat = unlist(compute.Var.by.Annot(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W,independence = independence))

  p = pchisq(Burden.Stat/Var.Stat,1,lower.tail = FALSE)

  df.p=data.frame(p)
  colnames(df.p) = paste0("Score_",colnames(Z_annot))

  #colnames(df.p) = paste0("Score",1:length(p))
  
  if(length(p)==1){
    colnames(df.p) = "Burden"
    df.p$ACAT = ACAT::ACAT(df.p$Burden)

    df.p$Fisher = pchisq(-2*sum(log(df.p$Burden)),2* length(df.p$Burden), lower.tail = F)
    }
  
  else{
    colnames(df.p) = paste0("Score_",colnames(Z_annot))

  df.p$ACAT = apply(df.p,1,function(x) ACAT::ACAT(x[!is.nan(x)]))

  df.p$Fisher = apply(df.p[,-ncol(df.p)],1, function(x){
    pchisq(-2*sum(log(x[!is.nan(x)])),2* length(x[!is.nan(x)]), lower.tail = F)
  })

  if(contrib.by.fam==TRUE & contrib.by.var==TRUE){
    return(list("Resuts"=df.p, "Stat.by.fam" = run.Burden$Stat.by.fam , "Stat.by.var" = run.Burden$Stat.by.var))
  }

  else if(contrib.by.fam==TRUE & contrib.by.var==FALSE){
    return(list("Resuts"=df.p, "Stat.by.fam" = run.Burden$Stat.by.fam ))
  }

  else if(contrib.by.fam==FALSE & contrib.by.var==TRUE){
    return(list("Resuts"=df.p, "Stat.by.var" = run.Burden$Stat.by.var ))
  }
  else{
    return(df.p)
  }
    }
  return(df.p)
}
