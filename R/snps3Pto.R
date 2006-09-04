snps3PrimeTo = function(gn, rad=50000) {
 chr = as.character(geneLocs[gn,"chr"])
 metaDat = paste("chr", chr, "meta", sep="")
 md = get(metaDat)
 snpdf =  get("meta", md@meta)
 pos = snpdf[,"pos"]
 ini = geneLocs[gn,"end"]
 las = ini+50000
 snpID(rownames(snpdf)[which(pos > ini & pos < las)])
}
 
 
 
snps5PrimeTo = function(gn, rad=50000) {
 chr = as.character(geneLocs[gn,"chr"])
 metaDat = paste("chr", chr, "meta", sep="")
 md = get(metaDat)
 snpdf =  get("meta", md@meta)
 pos = snpdf[,"pos"]
 ini = geneLocs[gn,"beg"]
 las = ini-50000
 snpID(rownames(snpdf)[which(pos > las & pos < ini)])
}

snpsNear = function(gn, rad=50000) {
 chr = as.character(geneLocs[gn,"chr"])
 metaDat = paste("chr", chr, "meta", sep="")
 md = get(metaDat)
 snpdf =  get("meta", md@meta)
 pos = snpdf[,"pos"]
 ini = geneLocs[gn,"beg"]
 las = ini-50000
 tail = geneLocs[gn,"end"]
 las2 = tail+50000
 snpID(rownames(snpdf)[which(pos > las & pos < las2)])
}
