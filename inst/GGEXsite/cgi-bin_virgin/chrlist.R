#! @RPATH@
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#
tag(html)
tag(head)
tag('body style="color: rgb(0, 0, 0);" alink="#358cff" link="#358cff" 
vlink="#990099"')

getsymbolsbychr<-function(chr){
library(hgfocus.db)

	hgsymbols <- as.list(hgfocusSYMBOL)
	symbolnames=names(hgsymbols);
	chrs <- as.list(hgfocusCHR)
	symbolstemp=names(subset(chrs,chrs==as.character(chr)))
	symbols=matrix(ncol=2,nrow=length(symbolstemp));
	for(i in 1:length(symbolstemp)){
		index=which(symbolnames==symbolstemp[i]);
		symbols[i,1]=as.character(hgsymbols[index]);
		symbols[i,2]=chr;
	}
symbols=symbols[order(symbols[,1]),]
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
	
	for(i in 1:length(symbolstemp)){
		tag("tr")
		for(j in 1:2){
			tag("td")
			cat(symbols[i,j])
			untag("td")
		}
	tag("td")
			b=paste('<a href="../R.cgi/nsnps.R?gene=',symbols[i,1],'&chromosome=',symbols[i,2],'">select</a>',sep="")
	cat(b)

	untag("td")
	untag("tr")
	cat("\n")
	}
	untag(table)
	cat("\n")

}

getsymbolsbychr(formData$string)

