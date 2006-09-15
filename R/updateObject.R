# by Martin Morgan (14 Sep 2006)

setMethod("updateObject",
          signature(object="racExSet"),
          function(object, ..., verbose=FALSE) {
              if (verbose)
                message("updateObject(object = racExSet)")
              if (!isS4(object)) {
                  new("racExSet",
                      phenoData=updateObject(phenoData(object), verbose=verbose),
                      experimentData=updateObject(experimentData(object), verbose=verbose),
                      annotation=updateObject(annotation(object), verbose=verbose),
                      exprs=exprs(object),
                      racs=snps(object),
                      SNPalleles=SNPalleles(object)
                      )
              } else object
          })

