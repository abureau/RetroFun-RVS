.check.parameters.agg.null = function(null.value.by.fam, aggregate.geno.by.fam) {
  
  
  if(null.value.by.fam$distinguishHomo == TRUE & (attributes(aggregate.geno.by.fam$ped_agg)$correction == "replace" | attributes(aggregate.geno.by.fam$ped_agg)$correction == "remove")){
    return(2)
  }
  
  else if(null.value.by.fam$distinguishHomo == FALSE & (attributes(aggregate.geno.by.fam$ped_agg)$correction == "none")){
    return(1)
  }
  
  else{
    return(0)
  }
}
