#! /Users/stvjc/Desktop/ExternalSoft/R-devel/bin/R
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#

tag(html)
tag(head)
tag('body style="color: rgb(0, 0, 0);" alink="#358cff" link="#358cff" vlink="#990099"')
popkegg2<-function(string){
	
a=read.delim("bbout.txt",header=F,sep=":",colClasses=c("character","character"))
a[,2]=tolower(a[,2])
a=a[order(a[,2]),]
	a=a[grep(string,a[,2],ignore.case = TRUE),2];
	tag('body style="color: rgb(0, 0, 0)" alink="#358cff" 
link="#358cff" 
	vlink="#990099"')
	tag(center)
	cat("\n")

	tag("table border=0")
	cat("\n")
	tag(tr)
	tag(td)
	cat("Path")
	untag("td")
	tag("td")
	cat("Select")
	untag("td")
	untag("tr")
	cat("\n")
	
	
	for(i in 1:length(a)){	
		d=gsub('[[:punct:]]',"",a[i])
		p='[[:space:]]';	
		d=gsub(p,"",d)
			tag("tr")
			tag("td")
			cat(a[i])
			untag("td")
		
	tag("td")
	b=paste('<a href="../R.cgi/keggsearch2.R?string=',d,'">select</a>',sep="");
	cat(b)

	untag("td")
	untag("tr")
	cat("\n")
	}
	untag(table)
	cat("\n")
}

popkegg2(formData$string)
