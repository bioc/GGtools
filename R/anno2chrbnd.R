anno2chrbnd = function(tag) {
 x = get(paste(tag,"CHRLENGTHS",sep=""))
 cumsum(x)
}
