masterSnps = function(mgw, n=50, auto=TRUE, orgdb="org.Hs.eg.db",
   minl10=5, gstart=0, gend=3e9, genomesize=3e9, pcex=1, 
   pal=rainbow(20), numxax=FALSE, ...) {
#
# visualize a multiGwSnpScreenResult (filtered)
#
 palette(pal)
 gn = geneIds(mgw@geneset)
 if (auto) mgw@.Data = lapply(mgw, function(x)x[1:22])
 allp = lapply(mgw, lapply, function(x) -log10(sort(p.value(x))[1:n]))
 names(allp) = gn
# deal with gene locations
 require(orgdb, character.only=TRUE)
 orgtag = gsub("\\.db", "", orgdb)
 smap = get(paste(orgtag, "SYMBOL", sep=""))
 lmap = get(paste(orgtag, "CHRLOC", sep=""))
 lens = get(paste(orgtag, "CHRLENGTHS", sep=""))
 nlens = names(lens)
 lenn = gsub("chr", "", nlens)
 names(lens) = lenn
 drc = union(grep("hap", lenn), union( grep("X", lenn), 
      union(grep("M", lenn), grep("Y", lenn))))
 lens = lens[-drc]
 nl = names(lens)
 nln = as.numeric(nl)
 lens = lens[order(nln)]
 rsmap = revmap(smap)
 egid = na.omit(unlist(mget(gn, rsmap, ifnotfound=NA)))
 locs = abs(unlist(lapply(mget(egid, lmap), "[", 1))) # deal with multiloc
 nl = names(locs)
 chr = gsub(".*\\.", "", nl)
 hap = grep("hap", chr)
 ngn = mget(egid, smap)  # should be clean
 names(locs) = ngn
 dr = is.na(locs)
 dr[hap] = TRUE
 if (auto) { 
           drc = which(chr %in% c("X", "Y"))
           if (length(drc)>0) dr[drc] = TRUE
           }
 if (sum(dr)>0) {
  allp = allp[-which(dr)]
  locs = locs[-which(dr)]
  chr = chr[-which(dr)]
  }
 off = cumsum(c(0,lens))
 nchr = as.numeric(chr)
 gloc = locs + off[nchr]
#
 oksnp = lapply(allp, function(x) {uu = unlist(x); uu[uu>minl10]})
 names(oksnp) = names(allp)
 empt = sapply(oksnp, function(x)length(x)==0)
 OKS = oksnp
 if (sum(empt)>0) OKS = oksnp[-which(empt)]
 nn = sapply(OKS,names)
 us = unique(unlist(nn))
 ULOC = snpLocs.Hs(rsid(as.character(us)))
 if (any(bc <- is.na(ULOC[1,]))) {
    badcol = which(bc)
    ULOC = ULOC[,-badcol]
 }
 rsids = paste("rs", ULOC[1,], sep="")
 gwlocs = ULOC[2,]
 names(gwlocs) = rsids
# snplocs = lapply(oksnp, function(x) {
#     cat(1); if (length(x)>0) snpLocs.Hs(rsid(names(x)))})
 snplocs = NULL
oksn = lapply(OKS, names)
oksnloc = lapply(oksn, function(x) gwlocs[x])
allglocs = gloc
snwlocg = gwlocs
okglocs = allglocs[names(oksnloc)]
okgx = list()
for (i in 1:length(okglocs)) okgx[[i]] = rep(okglocs[i], length(oksnloc[[i]]))
names(okgx) = names(oksnloc)
ans = list(allp=allp, locs=locs, chr=chr, lens=lens, gloc=okglocs, off=off,
    oksnp=OKS, snplocs=snplocs, gwlocs=gwlocs, okgx=okgx, okgy=oksnloc)
#
# now do visualization -- first layout everything
#
layout(matrix(1:2,nr=1), widths=c(.15,.85))
#
# now a fake plot to hold legend
#
par(mar=c(0,0,0,0))
plot(0,0,pch=" ", xlab=" ", ylab=" ", axes=FALSE,
 xlim=c(0,10), ylim=c(0,100))
#
# legend preferably without a box, telling us where the
#  genes are and their names
#
pchs = 1:length(names(ans$gloc))
pchs = pchs %% 19
pchs[pchs==0] = 9
legend(0,100, legend=names(ans$gloc), 
    pch=pchs,
    col=1:length(ans$gloc), box.col=par("bg"))
#
# now the data plot
#
par(mar=c(4,4,4,4))
plot(ans$okgx[[1]], ans$okgy[[1]], xlim=c(gstart,gend), ylim=c(0,genomesize), 
   pch=1, cex=pcex, axes=FALSE, 
   xlab= ifelse(numxax, "genomewide location", "chromosome harboring gene"), 
   ylab="chromosome harboring SNP", ...)
#
# put in a guide to help see cis-assoc
#
abline(0,1,col="gray")
#
# show where the genes are, with colors
# 
#for (i in 1:length(ans$gloc)) {
#  segments(ans$gloc[i]+.8e7*(i-1), genomesize, 
#    ans$gloc[i], (.9)*genomesize, col=i, lwd=3)
#  }
#
# axes -- an ugly business
#
midpts = (off[-1] + off[-length(off)])/2
if (!numxax) axis(1, at=off[-1], lab=1:22, cex.lab=.8)
   else axis(1, cex.lab=.8)
axis(2, at=off[-1], lab=1:22, cex.lab=.8)
#
# following from 1 to overplot
#
for (i in 1:length(ans$okgx)) points(ans$okgx[[i]], ans$okgy[[i]], 
    pch=ifelse(i %% 19, i %% 19, 9), cex=pcex, col=i)
abline(v=0, col="gray")
invisible(ans)
}

