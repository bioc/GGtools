
detScreen = function(gge=c20GGceu, psn="206918_s_at", 
  chrmeta=chr20meta, chr="chr20", gran=50, gene="CPNE1", ...) {
opar = par()
cpn = regseq(gge, psn, 
  seq(1,ncol(gge@phenoData@pData),gran), chrmeta, ...)
par(mfrow=c(1,2))
plot(cpn$locs, -log10(cpn$pva), main=paste(psn,chr), xlab="position",
ylab="-log10 p Ho:B=0")
bot = which.min(cpn$pva)
ggrplot(gge, gene, names(bot))
par(opar)
invisible(list(bot=bot, cpn=cpn))
}

