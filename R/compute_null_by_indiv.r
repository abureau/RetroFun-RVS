#' Compute the null value for pedigrees
#'
#' This function computes the genotype null value for one or more families. It assumes that the
#' the first parameter is a list with all the possible family configurations while the second
#' parameter is a list with all the sharing probabilities for the all family configurations
#'
#' @param pedigree is a object of pedigree format
#' @param eta multiplicative factor for the joint carrier probability for computing the covariance between two variants for a pairs of subjects (between 0 and 1)
#' @param etap multiplicative factor for the marginal carrier probability for computing the covariance between two variants in the same subject (between 0 and 1)
#' @return List over pedigrees of lists of the individual variance-covariance matrix for one variant and the individual covariance matrix for two variants 
#' @export
compute.null.by.indiv <- function(pedigree,eta=0.95,etap=0.95) {
  
  ### Configurations possibles et probabilités de partage
  configs.with.probs <- lapply(names(pedigree), function(fam){ 
    
    # Identification des sujets atteints 
    carriers <- pedigree[[fam]]$id[pedigree[[fam]]$affected == 1]
    
    # Génération des configurations possibles
    carrier.sets <- list()
    for(i in 1:length(carriers)){
      carrier.sets <- c(carrier.sets, combn(carriers, i, simplify = FALSE))
    }
    
    # Matrice des configurations
    configs <- t(sapply(carrier.sets, function(set) as.integer(carriers %in% set)))
    colnames(configs) <- carriers
    
    # Probabilité de partage
    probs <- sapply(carrier.sets, function(vec)
      RVS::RVsharing(pedigree[[fam]], carriers = vec, useAffected = TRUE))
    
    # Données retrournées
    list(configs = configs, probs = probs)
  })
  
  names(configs.with.probs) <- names(pedigree)
  
  ## Appel de la fonction compute.var.cov.by.indiv pour calculer la variance et les covariances
  
  var.covar <- compute.var.cov.by.indiv(configs.with.probs,eta=eta,etap=etap)
  
  return(var.covar)
}

