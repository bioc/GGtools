#
#makeLocNCDF = function(smSupp, fn, report=TRUE, doRSID=FALSE) {
# # take a snpMatrix hapmap 'support' data frame
# # and convert relevant information to NetCDF
## major dimension
# require(org.Hs.eg.db)
# ns = nrow(smSupp)
## compute cumlocs using CHRLENGTHS in org.Hs.eg.db
# cat("computing cumulative locations...")
# allp = smSupp$Posi
# allc  = as.character(smSupp$Chrom)
# allc[allc=="chrY"] = "24"
# allc[allc=="chrX"] = "23"
# allc = gsub("chr", "", allc)
# allc = as.numeric(allc)
# sp = split(allp, allc)
## cmax = sapply(sp,max)
## cmp = cumsum(as.numeric(cmax))
# cmp = cumsum(as.numeric(org.Hs.egCHRLENGTHS[1:24]))
# for (i in 2:length(sp)) sp[[i]] = sp[[i]] + cmp[i-1]
# cumlocs = unlist(sp)
# cat("done.\n")
## do NCDF work
# dimSNP = dim.def.ncdf( "allsnps", "snp", 1:ns )
# varCumloc = var.def.ncdf( "cumloc", "bases", dimSNP, -1,
#    longname="Cumulative base count", prec="double")
# varChr = var.def.ncdf( "chr", "token", dimSNP, -1,
#    longname="chr num (x=23, y=24)", prec="short")
# varl = list( varCumloc, varChr )
# if (doRSID) { # also see below for put
#  # convert rs numbers to numeric
#   allrs = rownames(smSupp)
#   rsc = gsub("rs", "", allrs)
#   rsint = as.numeric(rsc)
#   varRSN = var.def.ncdf( "rsnum", "token", dimSNP, -1,
#      longname="dbsnp ID dropping rs prefix", prec="integer")
#   varl[[3]] = varRSN
#   }
# cat("writing to disk...")
# ncnew = create.ncdf( fn, varl )
# put.var.ncdf( ncnew, varCumloc, cumlocs, start=1, count=ns )
# put.var.ncdf( ncnew, varChr, as.integer(allc), start=1, count=ns )
# if (doRSID) put.var.ncdf( ncnew, varRSN, as.integer(rsint), start=1, count=ns )
# ans = close.ncdf(ncnew)
# cat("done.\n")
# if (report) {
#   tmp = open.ncdf(fn)
#   print(tmp)
#   close.ncdf(tmp)
#   }
# invisible(ans)
#}
# 
