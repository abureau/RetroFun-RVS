#'Function computing the theoretical variance and covariance for each family
#'
#' @param configs.with.probs a list with configurations and corresponding probabilities by family
#' @param distinguishHomo a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise
#' #' @param cryptic.relatedness a boolean: TRUE when cryptic relatedness is expected, FALSE otherwise
#' @return a dataframe with the genotype variance for each family
#'
compute.var.by.fam = function(configs.with.probs, distinguishHomo=FALSE, cryptic.relatedness=FALSE){
  
  
  if(distinguishHomo==FALSE & cryptic.relatedness==FALSE){
    
    
    var.all.fam = sapply(configs.with.probs, function(x){
      config.fam = x$configs
      proba.fam = x$probs
      
      var = sum(config.fam^2*proba.fam) -
        sum(config.fam * proba.fam)^2
      
    })
    
    
    covar.all.fam = sapply(configs.with.probs, function(x){
      
      config.fam = x$configs
      proba.fam = x$probs
      
      joint = outer(proba.fam,proba.fam,pmin)/sum(outer(proba.fam,proba.fam,pmin))
      marginal = apply(joint,1,sum)
      moy = sum(config.fam*marginal)
      sum(outer(config.fam,config.fam,"*")*joint) - moy^2
    })
    
    covar.all.fam=pmin(covar.all.fam, var.all.fam)
    
    df_var = data.frame("FamID"=names(configs.with.probs), "Var"=var.all.fam)
    df_covar = data.frame("FamID"=names(configs.with.probs), "CoVar"=covar.all.fam)
    
    return(merge(df_var, df_covar, by="FamID"))
  }
  
  else if(distinguishHomo==TRUE | cryptic.relatedness==TRUE){
    
    
    var.all.fam = sapply(configs.with.probs, function(x) {
      config = unlist(sapply(x, function(y) y$configs))
      proba = unlist(sapply(x, function(y) y$probs))
      sum(config ^ 2 * proba) - sum(config * proba)^2
    })
    
    covar.all.fam = sapply(configs.with.probs, function(x) {
      
      config = unlist(sapply(x, function(y) y$configs))
      proba = unlist(sapply(x, function(y) y$probs))
      # En présence de probabilités ex aequo, le calcul de la somme de la matrice fonctionne
      # seulement avec la méthode "average" et le calcul de la somme de chaque ligne
      # (distribution marginale) seulement avec la méthode "min". Ce n'est pas clair pourquoi.
      rang = rank(proba,ties.method="average")
      rang.min = rank(proba,ties.method="min")
      lproba = length(proba)
      # Somme de la matrice du minimum par paire d'éléments de proba (outer(proba,proba,pmin))
      sm = sum((2*(lproba-rang)+1)*proba)
      # La prochaine instruction n'est pas conservée. Est-ce que des NaN peuvent arriver?
      #      joint[is.nan(joint)] = 0
      # Somme des ligne de la matrice du minimum par paire d'éléments de proba
      marginal = numeric(lproba)
      for (l in 1:lproba)
      {
        marginal[l] = ((lproba-rang.min[l]+1)*proba[l] + sum(proba[rang.min<rang.min[l]]))/sm
      }
      moy = sum(config*marginal)
      # Somme du produit "extérieur" de config multiplié élément-par-élément avec la matrice joint (sum(outer(config,config,"*")*joint))
      covar = 0
      for (l in 1:lproba)
        covar = covar + sum(config[l]*config * ifelse(proba[l]<proba,proba[l],proba)/sm)
      covar = covar - moy^2
      covar
    })
    
    covar.all.fam = pmin(var.all.fam, covar.all.fam)
    df_var = data.frame("FamID"=names(configs.with.probs), "Var"=var.all.fam)
    df_covar = data.frame("FamID"=names(configs.with.probs), "CoVar"=covar.all.fam)
    
    return(merge(df_var, df_covar, by="FamID"))
    
  }
  
}
