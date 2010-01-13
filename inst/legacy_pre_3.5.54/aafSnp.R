# we overload some of the annaffy methods to
# render information about SNP

setClass("aafSNP", "character", prototype = character(0))

aafSNP = function(x) {
 result = vector("list", length(x))
 attrs = list(class="aafSNP")
 for (i in 1:length(x)) {
  ann = x[[i]]
  attributes(ann) = attrs
  result[[i]]  = ann
  }
 class(result) = "aafList"
 asS4(result)
}

setMethod("getURL", "aafSNP", function(object) {

    url = "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs="
    urlsuffix <- ""

    if( !length(object) )
        return(character(0))
    return(paste(url, object, urlsuffix, sep = ""))
})


setMethod("getHTML", "aafSNP", function(object) {

    if( !length(object) )
        return("")
    if( length(url <- getURL(object)) )
        return(paste(paste("<a href=\"", url, "\">", object, "</a>", sep = ""), collapse = " "))
    else
        return(text)
})

