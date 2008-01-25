#! @RPATH@
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#

tag(html)
tag(head)
tag('body style="color: rgb(0, 0, 0);" alink="#358cff" link="#358cff" 
vlink="#990099"')
poptable<-function(string){
	
a=read.delim("bboutf.txt",header=F,sep=":",colClasses=c("character","character"))
	library(hgfocus)
	hgsymbols <- as.list(hgfocusSYMBOL)
	symbols=names(hgsymbols)
	pathprobe <- as.list(hgfocusPATH2PROBE)
	chrs <- as.list(hgfocusCHR)
	nameschrs=names(chrs)
	ids=pathprobe[grep(a[grep(string,a[,2],ignore.case = 
TRUE),1],names(pathprobe),ignore.case = TRUE)]
	ids2=unlist(ids)
	syms=matrix(data="",ncol=3,nrow=length(ids2))
	colnames(syms)=c("symbol","chromosome","symbolunchanged")
	for(i in 1:length(ids2)){ # should be vectorized!
		ind=which(symbols==ids2[i])
		syms[i,1]=hgsymbols[[ind]]
		syms[i,3]=hgsymbols[[ind]] 
		}
	for(i in 1:length(ids2)){
		ind=which(nameschrs==ids2[i])[1]
		syms[i,2]=chrs[[ind]][1] # watch for dupd genes
                }
	toupper(syms[,1])
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
	
	for(i in 1:length(ids2)){
		tag("tr")
		for(j in 1:2){
			tag("td")
			cat(syms[i,j])
			untag("td")
		}
	tag("td")
	b=paste('<a href="../R.cgi/nsnps.R?gene=',syms[i,3],'&chromosome=',syms[i,2],'">select</a>',sep="")
	cat(b)

	untag("td")
	untag("tr")
	cat("\n")
	}
	untag(table)
	cat("\n")
}
		

poptable(formData$string)
