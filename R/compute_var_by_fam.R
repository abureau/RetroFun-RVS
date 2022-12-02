#'Function computing the theoretical variance and covariance for each family
#'
#' @param configs.with.probs a list with configurations and corresponding probabilities by family
#' @param inbreeding a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise
#' @return a dataframe with the genotype variance for each family
#'
compute.var.by.fam = function(configs.with.probs, inbreeding=FALSE){


  if(inbreeding==FALSE){


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

  else if(inbreeding==TRUE){


    var.all.fam = sapply(configs.with.probs, function(x) {
      config = unlist(sapply(x, function(y) y$configs))
      proba = unlist(sapply(x, function(y) y$probs))
      sum(config ^ 2 * proba) - sum(config * proba)^2
    })

    covar.all.fam = sapply(configs.with.probs, function(x) {

      config = unlist(sapply(x, function(y) y$configs))
      proba = unlist(sapply(x, function(y) y$probs))
      joint = outer(proba,proba,pmin)/sum(outer(proba,proba,pmin))
      joint[is.nan(joint)] = 0
      marginal = apply(joint,1,sum)
      moy = sum(config*marginal)
      covar = sum(outer(config,config,"*")*joint) - moy^2
      covar
    })

    covar.all.fam = pmin(var.all.fam, covar.all.fam)
    df_var = data.frame("FamID"=names(configs.with.probs), "Var"=var.all.fam)
    df_covar = data.frame("FamID"=names(configs.with.probs), "CoVar"=covar.all.fam)

    return(merge(df_var, df_covar, by="FamID"))

  }

  else{
    print("Please provide TRUE or FALSE for the inbreeding parameter...")
  }

}
