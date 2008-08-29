
dropMono = function(x) {
 # output of parsePhPairs
 allh = unlist(x[[1]])
 saveh = t(sapply(allh, function(x) strsplit(x, "")[[1]]))
 uu = unique(allh)
 uus = t(sapply(uu, function(x) strsplit(x, "")[[1]]))
 mono = apply(uus,2,function(x)length(unique(x))==1)
 saveh = saveh[,-which(mono),drop=FALSE]
 rownames(uus) = NULL
 list(polyu=uus[,-which(mono),drop=FALSE], polyfull=saveh, monoinds=which(mono))
}

getTags = function(polyu) {
 df = data.frame(polyu)
 require(rpart, quietly=TRUE)
 Y = factor(1:nrow(df))
 tr1 = rpart(Y~., minsplit=1, data=df)
 tags = as.character(unique(tr1$frame$var))
 tags[tags!="<leaf>"]
}
 
personalTags = function(x) {
 m1 = dropMono(x)
 t1 = getTags(m1$polyu)
 inds = as.numeric(gsub("X", "", t1))
 m1$polyfull[, inds, drop=FALSE]
}

personalHap = function(x) {
 tt = personalTags(x)
 fixn = rownames(tt)
 fixn = fixn[seq(1,length(fixn),2)]
 fixn = gsub("1$", "", fixn)
 ttt = apply(tt,1,paste,collapse="")
 tttt = matrix(ttt,nr=2)
 ans = apply(tttt,2,function(z)paste(sort(z),collapse=":"))
 names(ans) = fixn
 ans
}

