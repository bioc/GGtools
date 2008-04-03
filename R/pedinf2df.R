pedinf2df = function(fn, ...) {
# create a data frame based on a hapmap pedigree info file
 V1 = V2 = V3 = V4 = V5 = V6 = V7 = 0
 x = read.table(fn, h=FALSE)
 attach(x)
 on.exit(detach(x))
 famid = V1
 persid = V2
 fathid = V3
 mothid = V4
 sampid = gsub("..*Sample:", "", as.character(V7))
 sampid = gsub(":.*$", "", sampid)
 isFounder = (V3 == 0 & V4 == 0)
 male = (V5 == 1)
 isAmom = persid %in% mothid
 isAdad = persid %in% fathid
 data.frame(famid=famid, persid=persid, mothid=mothid, fathid=fathid,
    sampid=sampid, isFounder=isFounder, male=male, isAmom = isAmom,
    isAdad = isAdad)
}
