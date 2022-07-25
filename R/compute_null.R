#' Compute the null value for pedigrees
#'
#' This function computes the genotype null value for one or more families. It assumes that the
#' the first parameter is a list with all the possible family configurations while the second
#' parameter is a list with all the sharing probabilities for the all family configurations
#'
#' @param pedigree.configurations is all possible configurations for people having the risk allele: A list
#' @param pedigree.probas is all probabilities associated with all pedigree configurations: A list
#' @return The expected genotype value, variance and covariance for each pedigree within a data.frame
#' @export

compute.null = function(pedigree.configurations, pedigree.probas){
  #pedigree.configurations and pedigree.probas must have the same length
  if(length(pedigree.configurations) != length(pedigree.probas)) print("Please use objects with same length")
  else{
    l = length(pedigree.configurations)
    expected_values = sapply(1:l, function(x) sum(pedigree.configurations[[x]] * pedigree.probas[[x]] ))

    df_expected_by_fam = data.frame("FamID" = names(pedigree.configurations), "Expected" = expected_values)
    df_var_covar_by_fam = compute.var.by.fam(pedigree.configurations, pedigree.probas)

    df_output = merge(df_expected_by_fam, df_var_covar_by_fam, by="FamID")
    return(df_output)
  }

}
