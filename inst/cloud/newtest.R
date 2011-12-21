

date()

f5 = function() {
 library(parallel)
 nc = detectCores()
 options(mc.cores=nc)
 library(GGtools)

nsFilter2 = function(sms, var.cutoff=.5) {
 alliq = apply(exprs(sms),1,IQR)
 qs = quantile(alliq,var.cutoff)
 sms[ which(alliq > qs), ]
}
s20 = getSS("GGdata", "20")
fndind = which(s20$mothid==0 & s20$fathid==0)

 
run8 = get.cis.eQTLs( "GGdata", ~male, geneApply=mclapply, chromApply=mpi.parLapply,
  exTransform=function(x) MAFfilter(nsFilter2(clipPCs(x, 1:10)[,fndind], var.cutoff=.5),
  lower=0.05),
  folderstem="run8" )
run8
save(run8, file="run8.rda")

}

library(Rmpi)

mpi.spawn.Rslaves(nsl=7)

options(verbose=TRUE)

#set up the workers

mpi.parLapply(1:7, function(x) 
 {library(parallel); suppressPackageStartupMessages(library(GGtools))})

# execute the big run

f5()

date()
 
mpi.close.Rslaves()
mpi.quit()
 
