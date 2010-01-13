
dropMono = function(x) {
 # output of parsePh.out
 allh = unlist(x[[1]])
 saveh = t(sapply(allh, function(x) strsplit(x, "")[[1]]))
 uu = unique(allh)
 uus = t(sapply(uu, function(x) strsplit(x, "")[[1]]))
 mono = apply(uus,2,function(x)length(unique(x))==1)
 rownames(uus) = NULL
 if (any(mono)) {
     saveh = saveh[,-which(mono),drop=FALSE]
     polyu = uus[,-which(mono),drop=FALSE]
     }
 else polyu = uus
 list(polyu=polyu, polyfull=saveh, monoinds=which(mono), rsid=x$rsid)
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
 if (length(m1$monoinds)>0) keeprs = m1$rsid[-m1$monoinds]
 else keeprs = m1$rsid
 list(tagdata=m1$polyfull[, inds, drop=FALSE], rsid=keeprs[inds])
}

personalHap = function(x) {
 tmp = personalTags(x)
 tt = tmp$tagdata
 fixn = rownames(tt)
 fixn = fixn[seq(1,length(fixn),2)]
 fixn = gsub("1$", "", fixn)
 ttt = apply(tt,1,paste,collapse="")
 if (ncol(tt) == 1) dim(ttt) = dim(tt) # in case we have a single nuc tag
 tttt = matrix(ttt,nr=2)
 ans = apply(tttt,2,function(z)paste(sort(z),collapse=":"))
 names(ans) = fixn
 list(hapdata=ans, rsid=tmp$rsid)
}

