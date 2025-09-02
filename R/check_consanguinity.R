#' usethis::use_package("kinship2")

#' Check whether there is consanguinity among pedigrees
#'
#' @param pedigrees is a pedigree of a pedigreeList format or list of pedigrees
#'

.check.consanguinity = function(pedigrees){


  if(class(pedigrees) == "list"){
    ped.consanguinity = c()

    for(ped in 1:length(pedigrees)){
      df.ped = kinship2::as.data.frame.pedigree(pedigrees[[ped]])

      couple = unique(df.ped[df.ped$dadid!=0,c("dadid","momid")])
      kinship.mat = kinship2::kinship(pedigrees[[ped]])

      consanguinity.coeff = c()
      for(i in 1:nrow(couple)){
        dadid = as.character(couple$dadid[i])
        momid = as.character(couple$momid[i])

        consanguinity.coeff = c(consanguinity.coeff,kinship.mat[dadid,momid ])
      }

      ped.consanguinity = c(ped.consanguinity, any(consanguinity.coeff>0))
    }

    return(which(ped.consanguinity))
  }

  else if (class(pedigree) == "pedigreeList"){
  else if (class(pedigrees) == "pedigreeList"){

    ped.consanguinity = c()
    famid = unique(pedigrees$famid)

    for(ped in famid){
      df.ped = kinship2::as.data.frame.pedigree(pedigrees[as.character(ped)])

      couple = unique(df.ped[df.ped$dadid!=0,c("dadid","momid")])
      kinship.mat = kinship2::kinship(pedigrees[as.character(ped)])

      consanguinity.coeff = c()
      for(i in 1:nrow(couple)){
        dadid = as.character(couple$dadid[i])
        momid = as.character(couple$momid[i])

        consanguinity.coeff = c(consanguinity.coeff,kinship.mat[dadid,momid])
      }

      ped.consanguinity = c(ped.consanguinity, any(consanguinity.coeff>0))
    }

    return(which(ped.consanguinity))
  }

  else if(class(pedigrees) == "pedigree"){

    df.ped = kinship2::as.data.frame.pedigree(pedigrees)

    couple = unique(df.ped[df.ped$dadid!=0,c("dadid","momid")])
    kinship.mat = kinship2::kinship(pedigrees)

    consanguinity.coeff = c()
    for(i in 1:nrow(couple)){
      dadid = as.character(couple$dadid[i])
      momid = as.character(couple$momid[i])

      consanguinity.coeff = c(consanguinity.coeff,kinship.mat[dadid,momid])
    }

    return(any(consanguinity.coeff>0))
  }


}


