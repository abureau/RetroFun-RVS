#' Compute the null value for pedigrees using parallel programming
#'
#' This function computes the genotype null value for one or more families. It assumes that the
#' the first parameter is a list with all the possible family configurations while the second
#' parameter is a list with all the sharing probabilities for the all family configurations
#'
#' @param pedigree is a object of pedigree format
#' @param distinguishHomo a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise
#' @param cryptic.relatedness a boolean: TRUE when cryptic relatedness is expected, FALSE otherwise
#' @param kinshipCoeff mean kinship coefficient among the pedigree founders
#' @param fam.ids a vector of family ids where inbreeding or cryptic relatedness may be computed
#' @param out.prob.dist a boolean: TRUE when the probability distributions need to be a part of the output
#' @return The expected genotype value, variance and covariance for each pedigree within a data.frame
#' @export

compute.null.parallel = function(pedigree, distinguishHomo = FALSE, cryptic.relatedness=FALSE, kinshipCoeff=NULL, fam.ids=NULL, out.prob.dist = FALSE){

  if(!distinguishHomo%in%c(TRUE,FALSE) | !cryptic.relatedness%in%c(TRUE,FALSE) ){
    stop("distinguishHomo or cryptic.relatedness parameters were not well set, please check ...")
  }

  #Check for consanguinity and correction parameters provided by the user

  if(class(pedigree) == "pedigreeList" | class(pedigree) =="list"){
    if(length(.check.consanguinity(pedigree)) > 0 & distinguishHomo == FALSE){
      warning(paste0("Consanguinity within pedigrees ", .check.consanguinity(pedigree), " distinguishHomo = TRUE should be provided ..."))

    }
    else if(length(.check.consanguinity(pedigree)) == 0 & distinguishHomo == TRUE & cryptic.relatedness == FALSE){
      warning("No consanguinity and no cryptic relatedness within pedigrees, please set distinguishHomo = FALSE ...")
    }
  }

  else if(class(pedigree) == "pedigree"){
    if(.check.consanguinity(pedigree) == TRUE & distinguishHomo == FALSE){
      warning(paste0("Consanguinity within pedigrees ", .check.consanguinity(pedigree), " distinguishHomo = TRUE should be provided ..."))

    }
    else if(.check.consanguinity(pedigree) ==FALSE & distinguishHomo == TRUE & cryptic.relatedness == FALSE){
      warning("No consanguinity and no cryptic relatedness within pedigrees, please set distinguishHomo = FALSE ...")
    }
  }


  if(distinguishHomo == FALSE & cryptic.relatedness == TRUE){
    warning("In presence of cryptic relatedness, the parameter choice distinguishHomo = FALSE is not optimal, please set distinguishHomo = TRUE ... ")
  }
  
  #Function to combine parallel results when distinguishHomo==TRUE
  fun_to_comb <- function(x) {
    len <- length(x); list <- vector(mode = "list", length = (len/2))
    for (i in 1:(len/2)){ list[[i]] <- list("configs" = x[[(2*i)-1]], "probs" = x[[2*i]]) }
    return(list)
  }
  
  if(class(pedigree) == "list"){
    famid = names(pedigree)
    
    l = lapply(famid, function(fam){
      
      nf = length(RVS::processPedigree(pedigree[[fam]])$founders)
      
      #Subset on affected status
      carriers = pedigree[[fam]]$id[pedigree[[fam]]$affected==1]
      carrier.sets = list()
      
      #Compute all the possible carrier configurations
      for(i in 1:length(carriers)){
        carrier.sets = c(carrier.sets, combn(carriers,i,simplify=FALSE))
      }
      
      #Scenario 1: No consanguinity
      if(distinguishHomo == FALSE){
        
        configs = sapply(carrier.sets, length)
        
        #Scenario 1-1: No consanguinity + No cryptic relatedness
        if(cryptic.relatedness == FALSE){
          
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE)
          }
          
          #probs = sapply(carrier.sets, fun_to_apply)
          probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
        }
        
        #Scenario 1-2: No consanguinity + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){
          
          if(!is.null(kinshipCoeff)){
            fun_to_apply <- function(vec) {
              number.of.carriers = length(vec)
              RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE, kinshipCoeff=kinshipCoeff, kinshipOrder=nf%/%2+1, distinguishHomo = TRUE)
            }
            
            #probs = sapply(carrier.sets, fun_to_apply)
            probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          }
          else{
            stop("Please provide a correct value for the kinship coeff ...")
          }
          
        }
        
        config_probs = list("configs" = configs, "probs" = probs)
        
      }
      
      #Scenario 2: Consanguinity
      else if (distinguishHomo == TRUE){
        
        #Scenario 2-1: Consanguinity + No cryptic relatedness
        if(cryptic.relatedness==FALSE){
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            
            prob.distinguish.homo = RVsharing(pedigree[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE)
            names.config = names(prob.distinguish.homo)
            
            config.homo = sapply(names.config,
                                 function(n) {sum(sapply(1:nchar(n), function(x){as.numeric(substr(n,x,x))})) + number.of.carriers})
            return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
          }
          
          #config_probs = lapply(carrier.sets, fun_to_apply)
          config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          config_probs <- fun_to_comb(x = config_probs)
          return(config_probs)
        }
        #Scenario 2-2: Consanguinity + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){
          
          if(!is.null(kinshipCoeff)){
            fun_to_apply <- function(vec) {
              number.of.carriers = length(vec)
              
              prob.distinguish.homo = RVsharing(myPedTrimmed[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE, kinshipCoeff = kinshipCoeff, kinshipOrder = nf%/%2+1)
              names.config = names(prob.distinguish.homo)
              
              config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
                as.numeric(substr(n,x,x))
              })) + number.of.carriers})
              return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
            }
            
            #config_probs = lapply(carrier.sets, fun_to_apply)
            config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
            config_probs <- fun_to_comb(x = config_probs)
          }
          else{
            stop("Please provide a correct value for the kinship coefficient ...")
          }
          return(config_probs)
        }
      }})
    
    names(l) = famid #names(pedigree)
    
  } else if(class(pedigree)=="pedigreeList"){
    famid = unique(pedigree$famid)
    
    l = lapply(famid, function(fam){
      
      nf = length(RVS::processPedigree(pedigree[fam])$founders)
      
      #Subset on affected status
      carriers = pedigree[fam]$id[pedigree[fam]$affected==1]
      carrier.sets = list()
      
      #Compute all the possible carrier configurations
      for(i in 1:length(carriers)){
        carrier.sets = c(carrier.sets, combn(carriers,i,simplify=FALSE))
      }
      
      #Scenario 1: No consanguinity
      if(distinguishHomo == FALSE){
        
        configs = sapply(carrier.sets, length)
        
        #Scenario 1-1: No consanguinity + No cryptic relatedness
        if(cryptic.relatedness == FALSE){
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            RVsharing(pedigree[fam], carriers=vec,useAffected = TRUE)
          }
          
          #probs = sapply(carrier.sets, fun_to_apply)
          probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          
        }
        
        #Scenario 1-2: No consanguinity + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){
          
          if(!is.null(kinshipCoeff)){
            fun_to_apply <- function(vec) {
              number.of.carriers = length(vec)
              RVsharing(pedigree[fam], carriers=vec,useAffected = TRUE, kinshipCoeff=kinshipCoeff, kinshipOrder=nf%/%2+1, distinguishHomo = TRUE)
            }
            
            #probs = sapply(carrier.sets, fun_to_apply)
            probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          }
          else{
            stop("Please provide a correct value for the kinship coeff ...")
          }
          
        }
        
        config_probs = list("configs" = configs, "probs" = probs)
      }
      
      #Scenario 2: Consanguinity
      else if (distinguishHomo == TRUE){
        
        #Scenario 2-1: Consanguinity + No cryptic relatedness
        if(cryptic.relatedness==FALSE){
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            
            prob.distinguish.homo = RVsharing(pedigree[fam], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE)
            names.config = names(prob.distinguish.homo)
            
            config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
              as.numeric(substr(n,x,x))
            })) + number.of.carriers})
            return(list("configs"=config.homo,"probs"= prob.distinguish.homo))
          }
          
          #config_probs = lapply(carrier.sets, fun_to_apply)
          config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          config_probs <- fun_to_comb(x = config_probs)
          return(config_probs)
        }
        #Scenario 2-2: Consanguinity + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){
          
          if(!is.null(kinshipCoeff)){
            fun_to_apply <- function(vec) {
              number.of.carriers = length(vec)
              
              prob.distinguish.homo = RVsharing(pedigree[fam], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE, kinshipCoeff = kinshipCoeff, kinshipOrder = nf%/%2+1)
              names.config = names(prob.distinguish.homo)
              
              config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
                as.numeric(substr(n,x,x))
              })) + number.of.carriers})
              return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
            }
            
            #config_probs = lapply(carrier.sets, fun_to_apply)
            config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
            config_probs <- fun_to_comb(x = config_probs)
            
          }
          else{
            stop("Please provide a correct value for the kinship coefficient ...")
          }
        }
      }})
    
    names(l) = famid #names(pedigree)
    
  } else if(class(pedigree)=="pedigree"){
    
    famid = unique(pedigree$famid)
    l = list()
    
    nf = length(RVS::processPedigree(pedigree)$founders)
    
    #Subset on affected status
    carriers = pedigree$id[pedigree$affected==1]
    carrier.sets = list()
    
    #Compute all the possible carrier configurations
    for(i in 1:length(carriers)){
      carrier.sets = c(carrier.sets, combn(carriers,i,simplify=FALSE))
    }
    
    #Scenario 1: No consanguinity
    if(distinguishHomo == FALSE){
      
      configs = sapply(carrier.sets, length)
      
      #Scenario 1-1: No consanguinity + No cryptic relatedness
      if(cryptic.relatedness == FALSE){
        fun_to_apply <- function(vec) {
          number.of.carriers = length(vec)
          RVsharing(pedigree, carriers=vec,useAffected = TRUE)
        }
        
        #probs = sapply(carrier.sets, fun_to_apply)
        probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
      }
      
      #Scenario 1-2: No consanguinity + Cryptic relatedness
      else if(cryptic.relatedness==TRUE){
        
        if(!is.null(kinshipCoeff)){
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            RVsharing(pedigree, carriers=vec,useAffected = TRUE, kinshipCoeff=kinshipCoeff, kinshipOrder=nf%/%2+1, distinguishHomo = TRUE)
          }
          
          #probs = sapply(carrier.sets, fun_to_apply)
          probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
        }
        else{
          stop("Please provide a correct value for the kinship coeff ...")
        }
        
      }
      
      config_probs = list("configs" = configs, "probs" = probs)
      
    }
    
    #Scenario 2: Consanguinity
    else if (distinguishHomo == TRUE){
      
      #Scenario 2-1: Consanguinity + No cryptic relatedness
      if(cryptic.relatedness==FALSE){
        fun_to_apply <- function(vec) {
          number.of.carriers = length(vec)
          
          prob.distinguish.homo = RVsharing(pedigree, carriers=vec, distinguishHomo = TRUE, useAffected = TRUE)
          names.config = names(prob.distinguish.homo)
          
          config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
            as.numeric(substr(n,x,x))
          })) + number.of.carriers})
          return(list("configs"=config.homo,"probs"= prob.distinguish.homo))
        }
        
        #config_probs = lapply(carrier.sets, fun_to_apply)
        config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
        config_probs <- fun_to_comb(x = config_probs)
      }
      #Scenario 2-2: Consanguinity + Cryptic relatedness
      else if(cryptic.relatedness==TRUE){
        
        if(!is.null(kinshipCoeff)){
          fun_to_apply <- function(vec) {
            number.of.carriers = length(vec)
            
            prob.distinguish.homo = RVsharing(pedigree, carriers=vec, distinguishHomo = TRUE,useAffected = TRUE, kinshipCoeff = kinshipCoeff, kinshipOrder = nf%/%2+1)
            names.config = names(prob.distinguish.homo)
            
            config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
              as.numeric(substr(n,x,x))
            })) + number.of.carriers})
            return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
          }
          
          #config_probs = lapply(carrier.sets, fun_to_apply)
          config_probs <- foreach(set = carrier.sets, .combine = 'c') %dopar% {fun_to_apply(set)}
          config_probs <- fun_to_comb(x = config_probs)
        }
        else{
          stop("Please provide a correct value for the kinship coefficient ...")
        }
      }
    }
    
    l[[famid]] = config_probs
    #names(l) = famid #names(pedigree)
  }
  
  
  if(is.null(fam.ids)){
    
    l = l
    
  } else if(!is.null(fam.ids)){
    
    l = l[fam.ids]
  }
  
  
  expected = compute.expected.by.fam(l, distinguishHomo=distinguishHomo, cryptic.relatedness=cryptic.relatedness)
  var.covar = compute.var.by.fam.parallel(l, distinguishHomo=distinguishHomo, cryptic.relatedness=cryptic.relatedness)
  
  df.expected.var.covar = merge(expected, var.covar, by="FamID")
  attributes(df.expected.var.covar)$distinguishHomo = distinguishHomo
  
  if (out.prob.dist){
    output <- list("distributions" = l, "expected.values" = df.expected.var.covar)
  } else {
    output <- df.expected.var.covar 
  }
  return(output)
}