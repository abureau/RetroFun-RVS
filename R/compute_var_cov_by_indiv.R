# Fonction de calcul de la variance et des covariances 
# Le paramètre configs.with.probs est une liste qui contient :
#   - les configurations possibles des individus porteurs d’un variant
#   - ainsi que leurs probabilités de partage
# Le paramètre eta correspond à la probabilité que deux individus distincts 
# d’une même famille portent chacun un variant différent.

compute.var.cov.by.indiv <- function(configs.with.probs, eta,etap) {
  
  resultat <- lapply(names(configs.with.probs), function(fam) {
    
    ## Extraction des configurations possibles et de leurs probabilités pour chaque famille
    
    configs <- configs.with.probs[[fam]]$configs
    probs   <- configs.with.probs[[fam]]$probs
    nb.indiv <- ncol(configs)

    ## Calcul de la variance individuelle et de la covariance qu’un même individu porte simultanément deux variants différents.
    p<-numeric(nb.indiv)
   var.ind<-numeric(nb.indiv)
   covar2<-numeric(nb.indiv)
   for (i in 1:nb.indiv) {
     p[i]<-sum(probs[configs[, i] == 1], na.rm = TRUE)
     
     #Calcul de la variance individuelle
     var.ind[i]<-p[i]-p[i]^2
     
     #covariance qu’un même individu porte simultanément deux variants différents.
     covar2[i]<-etap*p[i]-p[i]^2
   }
    
    ## Calcul des covariances covar1 et covar3
    covar1 <- matrix(0, nb.indiv, nb.indiv)  ## Calcul de la covariance entre deux individus (i et ip) concernant le même variant j
    
    covar3 <- matrix(0, nb.indiv, nb.indiv)  ## Calcul de la covariance entre deux individus (i et ip) concernant deux variants distincts j et j'
    
     for (i in 1:(nb.indiv-1)) {
       for (ip in (i+1):nb.indiv) {
         
         #### Calcul des probabilités que deux individus partagent simultanément un même variant.
        
         
         p11j <- sum(probs[configs[, i] == 1 & configs[, ip] == 1], na.rm = TRUE)
         

         ## Covariance entre deux individus (i et ip) concernant le même variant j.
         covar1[i,ip] <- p11j - p[i]^2
         covar1[ip, i] <- covar1[i, ip] 
         
         # Covariance entre deux individus (i et ip) concernant deux variants distincts j et j'.
         covar3[i, ip] <- eta*p11j - p[i]^2
         covar3[ip, i] <- covar3[i, ip]  
       }
     }
   
   diag(covar1)<-var.ind
   
   diag(covar3)<-covar2
  
    ##Mise sous forme triangulaire des matrices de covariance
    
    covar1[upper.tri(covar1)] <- 0
    covar3[upper.tri(covar3)] <- 0
    
    ## Resultats par famille 
    return(list("FamId" =fam,"Covar(ij,i'j)"=covar1,"Covar(ij,i'j')"=covar3))
  })
  
  return(resultat)
}

