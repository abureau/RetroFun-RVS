#'Function computing the genotype expected value for each family
#'
#' @param configs.with.probs a list with configurations and corresponding probabilities by family
#' @param inbreeding a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise 
#' @return a dataframe with the genotype expected value for each family
#' 
#' 


compute.expected.by.fam = function(configs.with.probs, inbreeding=FALSE){
  
  if(inbreeding==FALSE){
    expected = sapply(configs.with.probs, function(x) sum(x$configs* x$probs))
    
    return(data.frame("FamID" = names(configs.with.probs), "Expected" = expected))
  }
  
  else if(inbreeding==TRUE){
    expected = sapply(configs.with.probs, function(x){
      sum(sapply(x, function(y){
        sum(y$configs*y$probs)
      }))
    })
    
    return(data.frame("FamID" = names(configs.with.probs), "Expected" = expected))
  }
  else{
    print("Please provide TRUE or FALSE for the inbreeding parameter...")
  }
  
  
}
