# verifies that transScores agrees with snp.rhs.tests to obtain best
# trans scores

tenOn2021 = 
c( "GI_4557248-S",  "GI_15451784-S", "GI_4557678-S", "GI_9951914-S", "GI_21327679-S",
 "GI_7669476-I", "GI_7669478-A", "GI_4557290-I", "GI_4557294-A", "GI_41406053-S") 

tconf = new("TransConfig")
radius(tconf) = 2000000L
smpack(tconf) = "GGdata"
rhs(tconf) = ~1
snpchr(tconf) = "22"
chrnames(tconf) = c("20", "21")
gbufsize(tconf) = 4
smFilter(tconf) = function(x) x[probeId(tenOn2021),]
snpannopk(tconf) = snplocsDefault()
schrpref(tconf) = "ch"
exFilter(tconf) = function(x)x
gchrpref = ""
batchsize(tconf) = 200L

suppressPackageStartupMessages(library(GGtools))

t1 = transScores(tconf)

#"GGdata", rhs=~1, snpchr="22", chrnames=as.character(c(20,21)),
#
# all tests are trans, to verify buffering approach
#
#        radius = 2e+06,  K=4,
#    probesToKeep = tenOn2021, batchsize = 200, 
#    geneannopk = "illuminaHumanv1.db", 
#    snpannopk = "SNPlocs.Hsapiens.dbSNP.20111119", gchrpref = "", 
#    schrpref = "ch", exFilter = function(x) x) 

tt1 = transTab(t1)

c22 = getSS("GGdata", "22")

exl = lapply(tenOn2021, function(x) exprs(c22)[x,])
rhst = lapply(1:length(exl), function(g) {
     ex = exl[[g]]
     snp.rhs.tests(ex~1, snp.data=smList(c22)[[1]], fam="gaussian", uncertain=TRUE) })
csnp1 = sapply(rhst, function(x)chi.squared(x)[1])
csnp50 = sapply(rhst, function(x)chi.squared(x)[50])


all(abs(floor(sort(csnp1, decreasing=TRUE)[1:4]*10)/10 - tt1[1:4,2]) < .01)
all(abs(floor(sort(csnp50, decreasing=TRUE)[1:4]*10)/10 - tt1[197:200,2]) < 0.01)



