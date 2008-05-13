parse2SM = function(pop="CEU", rel="r23a", bld="b36") {
 chrtok = as.character(c(1:22, "X", "Y"))
 tmpl = function(chr, pop=pop, rel=rel, bld=bld)
   paste("file://genotypes_chr", chr, "_", pop, "_", rel,
     "_nr.", bld, ".txt.gz", sep="")
 allfi = sapply(chrtok, tmpl, pop, rel, bld)
 lapply(allfi, read.HapMap.data)
}
