countRare = function(x) {
 unph = strsplit(x,"")
 nunph = lapply(unph, function(x) if(any(x == "N")) return(c(NA,NA)) else return(x))
 ac = table(unlist(nunph))
 comm = names(ac)[which.min(ac)]
 sapply(nunph, function(x) sum(x == comm))
}
getRare = function(x) {
 unph = strsplit(x,"")
 nunph = lapply(unph, function(x) if(any(x == "N")) return(c(NA,NA)) else return(x))
 ac = table(unlist(nunph))
 rare = names(ac)[which.min(ac)]
 rare
}
