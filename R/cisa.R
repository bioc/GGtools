sampsInVCF = function(tf) {
#
# probe into VCF file to determine sample names
#  perhaps you don't need the chr?  can determine from header?
#
 samples(scanVcfHeader(tf))
}

cisAssoc = function( summex, vcf.tf, rhs=~1, nperm=3, cisradius=1000, 
    stx=force, vtx=force, snfilt=function(x) gsub("chr", "", x),
    genome="hg19", assayind=1, lbmaf=1e-6 ) {
 #
 # take all features from SummarizedExperiment
 # harmonize samples between summex and vcf.tf (TabixFile)
 # obtain genotypes of variants cis to features in summex -- using only SNVs!!
 # compute associations between stx(features) and vtx(genotypes)
 # store in a GRanges ordered by variant address with relevant metadata
 #
 # obtain sample IDs and harmonize genotypes and molec phenotype assay data
 #
 usn = unique(seqnames(summex))
 if(length(usn)>1) stop("current implementation insists that length(unique(seqnames(summex)))==1 as VCF assumed chr-specific")
 sampidsInSumm = colnames(summex)
 sampidsInVCF = sampsInVCF(vcf.tf, chrToTry=snfilt(usn))
 oksamp = intersect(sampidsInSumm, sampidsInVCF)
 stopifnot(length(oksamp)>0)
 summex = summex[, oksamp]
 #
 # harmonize annotation for seqnames -- could use style methods here
 #
 sn = snfilt(as.character(seqnames(summex)))
 stopifnot(length(ctouse <- unique(sn))==1)
 #
 # generate cis search space for assay probes
 #
 cisr = rowData(summex)+cisradius
 seqlevels(cisr) = snfilt(seqlevels(cisr)) # must use VCF-oriented seqlevels
 #
 # first pass at genotype data retrieval
 #
 vp = ScanVcfParam(fixed="ALT", info=NA, geno="GT", 
     samples=oksamp, which=cisr)  # which will sort variants into groups defined by probes
 vdata = readVcf(vcf.tf, genome=genome, param=vp) # compressed
 #
 # retain only SNVs with MAF > lbmaf
 #
 rdd = rowData(vdata)
 ok = nchar(rdd$REF)== 1
 ok = ok & (sapply(lapply(rdd$ALT,force),nchar) == 1)
 vdata = vdata[which(ok),]
 gtdata = genotypeToSnpMatrix(vdata)
 csumm = col.summary(gtdata[[1]])
 inmafs = csumm[,"MAF"]
 bad = which(inmafs<lbmaf)
 if (length(bad)>0) {
   vdata = vdata[-bad,]
   gtdata = genotypeToSnpMatrix(vdata)
   }
 varrd = rowData(vdata)  # would like to use this as the backbone of test result report
 #
 # use a list mapping probes to SNVs in cis to organize the testing
 #
 uo = unique(varrd$paramRangeID)
 vdl = split(varrd, varrd$paramRangeID)[uo] # will reorder without uo
 snpbyprobe = sapply(vdl, names)
 probes2test = names(snpbyprobe)
 numdata = assays(summex)[[assayind]][probes2test,]
 #
 # force the formula to have form ex~[rhs]
 #
 infmla = as.formula(paste("ex", paste(as.character(rhs), collapse="")))
 tsts = vector("list", length(probes2test) )
 summs = col.summary(gtdata$genotypes)
 mafs = summs[,"MAF"]
 names(mafs) = rownames(summs)
 #
 # loop over cis map to collect tests
 #
 suppressWarnings({
 for (i in 1:length(probes2test)) {
   ex = numdata[ probes2test[i], ]
   tsts[[i]] = snp.rhs.tests( formula=infmla, 
           snp.data=gtdata$genotypes[, snpbyprobe[[ probes2test[i] ]] ], family="gaussian",
           data=data.frame(ex=ex, as(colData(summex), "data.frame") ), uncertain=TRUE )
   }
  })
 #
 # test under permutation for plug-in FDR
 #
 perms = vector("list", nperm)
 permit = function(x) x[sample(1:length(x), size=length(x), replace=FALSE)]
 for (j in 1:nperm) {
   perms[[j]] = vector("list", length(probes2test))
   suppressWarnings({
   for (i in 1:length(probes2test)) {
     ex = permit(as.numeric(numdata[ probes2test[i], ]))
     perms[[j]][[i]] = snp.rhs.tests( formula=infmla, 
             snp.data=gtdata$genotypes[, snpbyprobe[[ probes2test[i] ]] ], family="gaussian",
             data=data.frame(ex=ex, as(colData(summex), "data.frame") ), uncertain=TRUE )
     }
    })
   }
 names(tsts) = probes2test
 #
 # bind test results to varrd GRanges instance -- note that ALT may be either
 # CharacterList or DNAStringSetList and this may need attention on collection
 #
 chisqs = unlist(lapply(tsts, chi.squared))
 varrd$chisq = chisqs
 pnames = paste0("permScore_", 1:nperm)
 for (i in 1:nperm) {
   mcols(varrd)[ ,pnames[i] ] = unlist(lapply(perms[[i]], chi.squared))
   }
 varrd$snp = names(varrd)
 varrd$MAF = as.numeric(mafs[varrd$snp])
 varrd$probeid = varrd$paramRangeID
 names(varrd) = NULL
 varrd
}