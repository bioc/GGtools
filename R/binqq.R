
qqhex = function(sco, p1, p2, p3, fdr, nxbins=20,
  thrs = c(0, .001, .005, .01, .05)) {
indset = lapply(thrs, function(x) which(fdr<=x))
scothrs = lapply(indset, function(x) min(sco[x]))
qqp = qqplot( sco, c(p1, p2, p3) , plot.it=FALSE)
list(hb=hexbin( qqp$x, qqp$y, xbins=nxbins ), thrs=thrs, scothrs=unlist(scothrs) )
}

binqq = function(  qqob, ylim=c(0, 76), xlim=c(0, 30), end45=5,  ... ) {
# par(mar=c(5,4,4,4))
 plot( qqob$hb@ycm, qqob$hb@xcm , pch= " ", xlab="Permuted", ylab="Observed",
     ylim=ylim, xlim=xlim )
 #points( qqob$hb@ycm, qqob$hb@xcm , pch= 19, col="gray")
 ok = which(qqob$hb@xcm < qqob$scothrs[1] )
 text( qqob$hb@ycm[ok], qqob$hb@xcm[ok] , labels=qqob$hb@count[ok], pch=.4 )
 text( qqob$hb@ycm[length(ok)+1], qqob$hb@xcm[length(ok)+1] , labels=sum(qqob$hb@count[-ok]), pch=.4 )
 segments( 0, 0, end45, end45, lwd = 3 )
 abline(h=qqob$scothrs, lty=1:length(qqob$scothrs) )
 axis(4, at=qqob$scothrs, labels=qqob$thrs, cex.axis=.75, las=2 )
}

