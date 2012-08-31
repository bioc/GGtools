
heavyTest = function() {

# verifies that transScores agrees with snp.rhs.tests to obtain best
# trans scores

tenOn2021 = 
c( "GI_4557248-S",  "GI_15451784-S", "GI_4557678-S", "GI_9951914-S", "GI_21327679-S",
 "GI_7669476-I", "GI_7669478-A", "GI_4557290-I", "GI_4557294-A", "GI_41406053-S") 

suppressPackageStartupMessages(library(GGtools))

mt1 = meta.transScores(c("GGdata", "GGdata"), rhs=list(~1, ~1), snpchr="22", chrnames=as.character(c(20,21)),
#
# all tests are trans, to verify buffering approach
#
        radius = 2e+06,  K=4, targdir="uiu2",
    probesToKeep = tenOn2021, batchsize = 200, 
    geneannopk = "illuminaHumanv1.db", 
    snpannopk = "SNPlocs.Hsapiens.dbSNP.20111119", gchrpref = "", SMFilterList = list( function(x) x[probeId(tenOn2021),],
         function(x) x[probeId(tenOn2021), ]),
    schrpref = "ch", exFilter = list( function(x) x, function(x) x)) 

transTab(mt1)
}
