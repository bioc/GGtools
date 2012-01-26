setGeneric("bindScores", function(obj, locs, isSNPlocs)standardGeneric("bindScores"))

.bindScores.cwt = function(cwt, gr, isSNPlocs=TRUE) {
# cwt is gwSnpTests result
# gr is GRanges for locations
# isSNPlocs will trigger modification of RefSNP_id metadata to include "rs" prefix
 scored = cwt[[1]]@snp.names
 scores = -log10(p.value(cwt[[1]]))
 bad = which(is.na(scores))
 if (any(bad)) {
   scores = scores[-bad]
   scored = scored[-bad]
   }
 names(scores) = scored
 if (isSNPlocs) {
    locd = elementMetadata(gr)$RefSNP_id
    locd = paste("rs", locd, sep="")
    } else locd = names(gr)
 oksnp = intersect(scored, locd)
 oksco = scores[oksnp]
 okloc = gr[ na.omit(match(names(oksco), locd)) ]
 elementMetadata(okloc)$score = oksco
 okloc
}

setMethod("bindScores", c("cwSnpScreenResult", "GRanges", "logical"),
   function(obj, locs, isSNPlocs) .bindScores.cwt(obj, locs, isSNPlocs))
