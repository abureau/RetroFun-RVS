## ---- include = FALSE---------------------------------------------------------
knitr::opts_knit$set(root.dir = "C:\\Users\\loicm\\Documents\\recherche\\Github\\RetroFun-RVS\\data")


## ----setup--------------------------------------------------------------------
library(RetroFunRVS)

set.seed(1234)
pedfile.clean = agg.genos.by.fam("sample_ped.ped")
pedfile.clean

## -----------------------------------------------------------------------------
load("SPAPsimpleprob.RData")
null = compute.null(forsim.N.list, forsim.pattern.prob.list)
head(null)

## -----------------------------------------------------------------------------

#Fist annotation with equal weights
first.annot = sample(c(0,1),510, replace = T)
#Second annotation with unequal weights
second.annot = sample(c(0,1), 510, replace = T, prob=c(0.8, 0.2))

#Annotation matrix, the first column should be composed only with ones
Z = matrix(c(rep(1,510),first.annot, second.annot),ncol=3,nrow=510)

#Equal weights, assuming correlation between variants
RetroFun.RVS(null, pedfile.clean, Z, diag(1, nrow=510,ncol=510), independence = F)

#Here we can assume weights depending on MAF
beta.weights = diag(dbeta(runif(510, 0.00001, 0.01), 1,25), ncol=510, nrow=510)

RetroFun.RVS(null, pedfile.clean, Z, beta.weights, independence = F) 

#No annotation is predictive with the trait
third.annot = sample(c(0,1),510, replace = T, prob = c(0.9,0.1))
fourth.annot = sample(c(0,1),510, replace = T, prob = c(0.95,0.05))

Z.nonpred = matrix(c(rep(1,510),third.annot, fourth.annot),ncol=3,nrow=510)

RetroFun.RVS(null, pedfile.clean, Z.nonpred, diag(1, nrow=510,ncol=510), independence = F)

## -----------------------------------------------------------------------------
prob.sharing.by.famid = lapply(1:length(forsim.pattern.prob.list), function(x) tapply(forsim.pattern.prob.list[[x]], forsim.N.list[[x]], sum))

names(prob.sharing.by.famid) = names(forsim.pattern.prob.list)

resample.genos = resample.genos.by.fam(pedfile.clean, prob.sharing.by.famid)

## -----------------------------------------------------------------------------

#Bootstrap burden statistic
compute.Burden.by.Annot(null,resample.genos, Z, diag(1, nrow=510,ncol=510))

#Observed burden statistic
compute.Burden.by.Annot(null,pedfile.clean, Z, diag(1, nrow=510,ncol=510))

