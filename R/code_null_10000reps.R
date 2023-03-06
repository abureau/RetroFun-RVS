library(RetroFunRVS)

pedfiles_null = list.files("C:\\\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Last_sim_null_10000reps\\TAD", full.names = T)


load("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\SPAPsimpleprob.RData")
null.extended.ped = compute.null(forsim.N.list, forsim.pattern.prob.list)

names(null.extended.ped)

agg.geno.by.fam.for.sim = function(pedfile, FamID=NULL){
  p = read.table(pedfile, header = F)
  fam = p[,1:6]
  affected = which(fam$V6==2)
  genos = p[,7:ncol(p)]

  genos[genos==1] = 0
  genos[genos==2] = 1

  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))

  df.genos.affected = df.genos[affected,]
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  df.genos.affected[df.genos.affected==2] = 1
  df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected)==0)]

  df.genos.affected$pedigree = fam[affected,"V1"]

  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum)

  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,]

  df.genos.agg.by.fam$pedigree = ifelse(df.genos.agg.by.fam$pedigree %in% FamID, df.genos.agg.by.fam$pedigree,sub('^[A-Z]', '', df.genos.agg.by.fam$pedigree))

  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1])))

  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))

}


pedfiles_agg_null = lapply(1:10000, function(rep) agg.geno.by.fam.for.sim(pedfiles_null[[rep]],null.extended.ped$FamID))



pvalues_null = lapply(1:10000, function(rep) RetroFun.RVS(null.extended.ped, pedfiles_agg_null[[rep]], Z,W))
pvalues_null_without_4 = lapply(1:10000, function(rep) RetroFun.RVS(null.extended.ped, pedfiles_agg_null[[rep]], Z_without_4,W))

Z_without_4 = Z[,-4]

hist(sapply(pvalues_null, function(x) x$ACAT))
gaston::qqplot.pvalues(sapply(pvalues_null, function(x) x$Score4))

which.min(sapply(pvalues_null, function(x) x$ACAT))

RetroFun.RVS(null.extended.ped, pedfiles_agg_null[[6879]], Z_without_4,W)

compute.Burden.by.Annot(null.extended.ped, pedfiles_agg_null[[6879]], Z,W)
compute.Var.by.Annot(null.extended.ped, pedfiles_agg_null[[3]], Z,W)

load("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\code\\pedsimple.RData")

mean(as.numeric(table(pedsimple$famid)[names(table(pedsimple$famid))%in%pedfiles_agg_null[[6879]]$ped_agg[,1]]))

max(sapply(1:10000, function(rep) mean(as.numeric(table(pedsimple$famid)[names(table(pedsimple$famid))%in%pedfiles_agg_null[[rep]]$ped_agg[,1]]))))

mean(rowSums(pedfiles_agg_null[[6879]]$ped_agg[,-1]))
nrow(pedfiles_agg_null[[6879]]$ped_agg)
Z[pedfiles_agg_null[[6879]]$index_variants,]
sapply(1:10000, function(x) nrow(pedfiles_agg_null[[x]]$ped_agg))
boxplot(sapply(1:10000, function(rep) mean(rowSums(pedfiles_agg_null[[rep]]$ped_agg[,-1]))))

which(Z[,4]==1)
pedfiles_agg_null[[6879]]$ped_agg$X225
pedfiles_agg_null[[6879]]$ped_agg$pedigree
null.extended.ped

((5-1.337513)*19.92426)^2
W[225,225]
null.extended.ped


compute.Burden.by.Annot(null.extended.ped, pedfiles_agg_null[[6879]], Z,W)
compute.Var.by.Annot(null.extended.ped, pedfiles_agg_null[[6879]], Z,W)


gaston::qqplot.pvalues(sapply(pvalues_null_without_4, function(x) x$ACAT))
colSums(Z)
