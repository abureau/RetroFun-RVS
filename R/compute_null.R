# usethis::use_package("RVS")

#' Compute the null value for pedigrees
#'
#' This function computes the genotype null value for one or more families. It assumes that the
#' the first parameter is a list with all the possible family configurations while the second
#' parameter is a list with all the sharing probabilities for the all family configurations
#'
#' @param pedigree is a object of pedigree format
#' @param inbreeding a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise
#' @param cryptic.relatedness a boolean: TRUE when cryptic relatedness is expected, FALSE otherwise
#' @param cryptic.coefficient an integer corresponding to the cryptic degree among founders
#' @param fam a vector of family ids where inbreeding or cryptic relatedness may be computed
#' @return The expected genotype value, variance and covariance for each pedigree within a data.frame
#' @export

compute.null = function(pedigree, inbreeding = FALSE, cryptic.relatedness=FALSE, kinshipCoeff=NULL, fam.ids=NULL){

  if(is.null(fam.ids)){

    l = lapply(names(pedigree), function(fam){

      nf = length(RVS::processPedigree(pedigree[[fam]])$founders)

      #Subset on affected status
      carriers = pedigree[[fam]]$id[pedigree[[fam]]$affected==1]
      carrier.sets = list()

      #Compute all the possible carrier configurations
      for(i in 1:length(carriers)){
        carrier.sets = c(carrier.sets, combn(carriers,i,simplify=FALSE))
      }
      #Scenario 1: No inbreeding
      if(inbreeding==FALSE){

        configs = sapply(carrier.sets, length)

        #Scenario 1-1: No inbreeding + No cryptic relatedness
        if(cryptic.relatedness==FALSE){

          probs = sapply(carrier.sets, function(vec) {
            number.of.carriers = length(vec)
            RVS::RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE)
          })
        }

        #Scenario 1-2: No inbreeding + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){

          if(!is.null(kinshipCoeff)){
            probs = sapply(carrier.sets, function(vec) {
              number.of.carriers = length(vec)
              RVS::RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE, kinshipCoeff=kinshipCoeff,kinshipOrder=nf%/%2+1)
            })
          }
          else{
            print("Please provide a correct value for the kinship coeff ...")
          }

        }

        config_probs = list("configs" = configs, "probs" = probs)

      }

      #Scenario 2: Inbreeding
      else if (inbreeding == TRUE){

        #Scenario 2-1: Inbreeding + No cryptic relatedness
        if(cryptic.relatedness==FALSE){
          config_probs = lapply(carrier.sets, function(vec) {
            number.of.carriers = length(vec)

            prob.distinguish.homo = RVS::RVsharing(pedigree[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE)
            names.config = names(prob.distinguish.homo)

            config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
              as.numeric(substr(n,x,x))
            })) + number.of.carriers})
            return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
          })
          return(config_probs)
        }
        #Scenario 2-2: Inbreeding + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){

          if(!is.null(kinshipCoeff)){
            config_probs = lapply(carrier.sets, function(vec) {
              number.of.carriers = length(vec)

              prob.distinguish.homo = RVS::RVsharing(myPedTrimmed[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE, kinshipCoeff = kinshipCoeff, kinshipOrder = nf%/%2+1)
              names.config = names(prob.distinguish.homo)

              config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
                as.numeric(substr(n,x,x))
              })) + number.of.carriers})
              return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
            })
          }
          else{
            print("Please provide a correct value for the kinship coeff ...")
          }
          return(config_probs)
        }
      }})
    names(l) = names(pedigree)
  }


  else if(!is.null(fam.ids)){

    l = lapply(fam.ids, function(fam){

      nf = length(RVS::processPedigree(pedigree[[fam]])$founders)

      #Subset on affected status
      carriers = pedigree[[fam]]$id[pedigree[[fam]]$affected==1]
      carrier.sets = list()

      #Compute all the possible carrier configurations
      for(i in 1:length(carriers)){
        carrier.sets = c(carrier.sets, combn(carriers,i,simplify=FALSE))
      }
      #Scenario 1: No inbreeding
      if(inbreeding==FALSE){

        configs = sapply(carrier.sets, length)

        #Scenario 1-1: No inbreeding + No cryptic relatedness
        if(cryptic.relatedness==FALSE){

          probs = sapply(carrier.sets, function(vec) {
            number.of.carriers = length(vec)
            RVS::RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE)
          })
        }

        #Scenario 1-2: No inbreeding + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){

          if(!is.null(kinshipCoeff)){
            probs = sapply(carrier.sets, function(vec) {
              number.of.carriers = length(vec)
              RVS::RVsharing(pedigree[[fam]], carriers=vec,useAffected = TRUE, kinshipCoeff=kinshipCoeff,kinshipOrder=nf%/%2+1)
            })
          }
          else{
            print("Please provide a correct value for the kinship coeff ...")
          }

        }

        config_probs = list("configs" = configs, "probs" = probs)

      }

      #Scenario 2: Inbreeding
      else if (inbreeding == TRUE){

        #Scenario 2-1: Inbreeding + No cryptic relatedness
        if(cryptic.relatedness==FALSE){
          config_probs = lapply(carrier.sets, function(vec) {
            number.of.carriers = length(vec)

            prob.distinguish.homo = RVS::RVsharing(pedigree[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE)
            names.config = names(prob.distinguish.homo)

            config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
              as.numeric(substr(n,x,x))
            })) + number.of.carriers})
            return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
          })
          return(config_probs)
        }
        #Scenario 2-2: Inbreeding + Cryptic relatedness
        else if(cryptic.relatedness==TRUE){

          if(!is.null(kinshipCoeff)){
            config_probs = lapply(carrier.sets, function(vec) {
              number.of.carriers = length(vec)

              prob.distinguish.homo = RVS::RVsharing(myPedTrimmed[[fam]], carriers=vec, distinguishHomo = TRUE,useAffected = TRUE, kinshipCoeff = kinshipCoeff, kinshipOrder = nf%/%2+1)
              names.config = names(prob.distinguish.homo)

              config.homo = sapply(names.config, function(n) {sum(sapply(1:nchar(n), function(x){
                as.numeric(substr(n,x,x))
              })) + number.of.carriers})
              return(list("configs"=config.homo,"probs"=prob.distinguish.homo))
            })
          }
          else{
            print("Please provide a correct value for the kinship coeff ...")
          }
          return(config_probs)
        }
      }})
    names(l) = names(pedigree)
  }

  if(inbreeding==TRUE){

    expected = compute.expected.by.fam(l, inbreeding=TRUE)
    var.covar = compute.var.by.fam(l, inbreeding=TRUE)
  }

  else if (inbreeding==FALSE){
    expected = compute.expected.by.fam(l, inbreeding=FALSE)
    var.covar = compute.var.by.fam(l, inbreeding=FALSE)
  }

  return(merge(expected, var.covar, by="FamID"))
}
