# Fonction de calcul de la variance et des covariances 
# Le paramètre configs.with.probs est une liste qui contient :
#   - les configurations possibles des individus porteurs d’un variant
#   - ainsi que leurs probabilités de partage
# Le paramètre eta correspond à la probabilité que deux individus distincts 
# d’une même famille portent chacun un variant différent.

compute.var.cov.by.indiv <- function(configs.with.probs,distinguishHomo=FALSE, cryptic.relatedness=FALSE, eta,etap) {
  
  resultat <- lapply(names(configs.with.probs), function(fam) {
    
    ## Extraction des configurations possibles et de leurs probabilités pour chaque famille
    
    configs <- configs.with.probs[[fam]]$configs
    probs   <- configs.with.probs[[fam]]$probs
    nb.indiv <- ncol(configs)

   p<-numeric(nb.indiv)
   var.ind<-numeric(nb.indiv)
   covar2<-numeric(nb.indiv)
   covar1 <- matrix(0, nb.indiv, nb.indiv)  ## Calcul de la covariance entre deux individus (i et ip) concernant le même variant j
   covar3 <- matrix(0, nb.indiv, nb.indiv)  ## Calcul de la covariance entre deux individus (i et ip) concernant deux variants distincts j et j'

    if (distinguishHomo=FALSE) {
     ## Calcul de la variance individuelle et de la covariance qu’un même individu porte simultanément deux variants différents.
     for (i in 1:nb.indiv) {
     p[i]<-sum(probs[configs[, i] == 1], na.rm = TRUE)
     
     #Calcul de la variance individuelle
     var.ind[i]<-p[i]-p[i]^2
     
     #covariance qu’un même individu porte simultanément deux variants différents.
     covar2[i]<-etap*p[i]-p[i]^2
   }
    
    ## Calcul des covariances covar1 et covar3
    
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
    }
    else if(distinguishHomo==TRUE) {
      
      # Calcul de la matrice de probabilités de 1 ou 2 copies du variant pour chaque paire d'individu
      mp = array(0,c(nb.indiv-1,nb.indiv,2,2))
      vp = matrix(0,nb.indiv,2)
      for (i in 1:(nb.indiv-1)) {
        for (ip in (i+1):nb.indiv) {
          vi = which(configs[, i] == 1 & configs[, ip] == 1)
          for (ci in vi)
          {
              #pattern = as.numeric(strsplit(R.utils::intToBin(patternIndex),"")[[1]]) + 1
              #pattern = c(rep(1,nb.indiv-length(pattern)),pattern)
              #if(pattern[i]==1&pattern[ip]==1) 
              mp[i,ip,1,1] = mp[i,ip,1,1] + sum(probs[[ci]][sapply(names(probs[[ci]]),function(bin) substr(bin,i,i)=="0" & substr(bin,ip,ip)=="0")])
              mp[i,ip,1,2] = mp[i,ip,1,2] + sum(probs[[ci]][sapply(names(probs[[ci]]),function(bin) substr(bin,i,i)=="0" & substr(bin,ip,ip)=="1")])
              mp[i,ip,2,1] = mp[i,ip,2,1] + sum(probs[[ci]][sapply(names(probs[[ci]]),function(bin) substr(bin,i,i)=="1" & substr(bin,ip,ip)=="0")])
              mp[i,ip,2,2] = mp[i,ip,2,2] + sum(probs[[ci]][sapply(names(probs[[ci]]),function(bin) substr(bin,i,i)=="1" & substr(bin,ip,ip)=="1")])
          }
          # On a besoin de calculer les probabilités marginales par sujet une seule fois
          if (i==1)
          {
          # Probabilité que le sujet 1 porte 1 ou 2 copies
          if (ip==1) vp[i,] = mp[i,ip,,1] + mp[i,ip,,2]
          # Probabilité que le sujet 2 porte 1 ou 2 copies
          vp[ip,] = mp[i,ip,1,] + mp[i,ip,2,]
          }
          
          ## Covariance entre deux individus (i et ip) concernant le même variant j.
          covar1[i,ip] <- mp[i,ip,1,1] + 2*mp[i,ip,1,2] + 2*mp[i,ip,2,1] + 4*mp[i,ip,2,2]  - (vp[i,1] + 2*vp[i,2])*(vp[ip,1] + 2*vp[ip,2])
          covar1[ip, i] <- covar1[i, ip] 
          
          # Covariance entre deux individus (i et ip) concernant deux variants distincts j et j'.
          covar3[i, ip] <- eta*mp[i,ip,1,1] + 2*mp[i,ip,1,2] + 2*mp[i,ip,2,1] + 4*mp[i,ip,2,2]  - (vp[i,1] + 2*vp[i,2])*(vp[ip,1] + 2*vp[ip,2])
          covar3[ip, i] <- covar3[i, ip]
          
        }
      ## Calcul de la variance individuelle et de la covariance qu’un même individu porte simultanément deux variants différents.
      var.ind[i] = vp[i,1] + 4*vp[i,2] - (vp[i,1] + 2*vp[i,2])^2
      covar[i] = etap*(vp[i,1] + 4*vp[i,2]) - (vp[i,1] + 2*vp[i,2])^2
      }
    }
    ## Resultats par famille 
    return(list("FamId" =fam,"Covar(ij,i'j)"=covar1,"Covar(ij,i'j')"=covar3))
  })
  
  return(resultat)
}

