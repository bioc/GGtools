#! /Users/stvjc/Desktop/ExternalSoft/R-devel/bin/R
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#

tag(html)
tag(head)
tag('body style="color: rgb(0, 0, 0);" alink="#358cff" link="#358cff" 
vlink="#990099"')
grepsymbols<-function(string){
	library(hgfocus)
	hgsymbols <- as.list(hgfocusSYMBOL)
	symbols=names(hgsymbols)
	pathprobe <- as.list(hgfocusPATH2PROBE)
	chrs <- as.list(hgfocusCHR)
	nameschrs=names(chrs)
	
	inds=grep(string,hgsymbols,value=T,ignore.case =T)
	#named vector of symbols
 
	namesinds=names(inds)
	#probeids

	syms=matrix(ncol=2,nrow=length(inds))
	
	for(i in 1:length(inds)){
		syms[i,1]=inds[i];
		index=which(nameschrs==namesinds[i])
		syms[i,2]=as.character(chrs[index])}

syms=syms[order(syms[,1]),]		
	tag('body style="color: rgb(0, 0, 0)" alink="#358cff" 
link="#358cff" 
	vlink="#990099"')
	tag(center)
	cat("\n")

	tag("table border=0")
	cat("\n")
	tag(tr)
	tag(td)
	cat("SYMBOL")
	untag("td")
	tag("td")
	cat("Chromosome")
	untag("td")
	tag("td")
	cat("Select")
	untag("td")
	untag("tr")
	cat("\n")
	
	for(i in 1:length(inds)){
		tag("tr")
		for(j in 1:2){
			tag("td")
			cat(syms[i,j])
			untag("td")
		}
	tag("td")
	b=paste('<a 

href="../R.cgi/nsnps.R?gene=',syms[i,1],'&chromosone=',syms[i,2],'">select</a>',sep="")
	cat(b)

	untag("td")
	untag("tr")
	cat("\n")
	}
	untag(table)
	cat("\n")

}
		

grepsymbols(formData$string)
