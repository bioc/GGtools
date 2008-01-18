#! /Users/stvjc/Desktop/ExternalSoft/R-devel/bin/R
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#

tag(html)
tag(head)
tag('body style="color: rgb(0, 0, 0);" alink="#358cff" link="#358cff"
vlink="#990099"')

tag('form action="../R.cgi/ggtools3.R" method ="head"')
a=paste('<input name="chromosome" type="hidden" value="',formData$chromosome,'">',sep="")
cat(a)
untag(input)
a=paste('<input name="gene" type="hidden" value="',formData$gene,'">',sep="")
cat(a)
untag(input)
cat("Number of SNPs to report:")
tag('input type="text" name="nsnps"')
untag(input)
tag('input type="submit" value="Send"')
untag(input)
