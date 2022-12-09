#'Function computing the genotype expected value for each family
#'
#' @param configs.with.probs a list with configurations and corresponding probabilities by family
#' @param distinguishHomo a boolean: TRUE when inbreeding among family members is expected, FALSE otherwise
#' #' @param cryptic.relatedness a boolean: TRUE when cryptic relatedness is expected, FALSE otherwise
#' @return a dataframe with the genotype expected value for each family
#'


compute.expected.by.fam = function(configs.with.probs, distinguishHomo=FALSE, cryptic.relatedness=FALSE){

  if(distinguishHomo==FALSE & cryptic.relatedness==FALSE ){
    expected = sapply(configs.with.probs, function(x) sum(x$configs* x$probs))

    return(data.frame("FamID" = names(configs.with.probs), "Expected" = expected))
  }

  else if(distinguishHomo==TRUE | cryptic.relatedness == TRUE){
    expected = sapply(configs.with.probs, function(x){
      sum(sapply(x, function(y){
        sum(y$configs*y$probs)
      }))
    })

    return(data.frame("FamID" = names(configs.with.probs), "Expected" = expected))
  }

}
