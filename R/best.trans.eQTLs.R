best.trans.eQTLs = function(smpack, rhs, genechrnum, snpchrnum,
	K = 20, targdirpref="tsco", batchsize=200, radius=2e6,
        genequeryprefix="", snploadprefix="chr", snplocprefix="chr", geneannopk,
	snpannopk, exFilter=function(x)x, smFilter=function(x)x,
	geneApply=lapply) {
   transScores( smpack=smpack, snpchr=paste( snploadprefix, snpchrnum, sep=""), 
		rhs=rhs, K=K, 
		targdirpref=targdirpref, chrnames=genechrnum,
		gchrpref=genequeryprefix, schrpref=snplocprefix,
		radius=radius, shortfac=10, wrapperEndo=smFilter,
		geneannopk=geneannopk, snpannopk=snpannopk )
}
		
