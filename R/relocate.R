
relocate = function(old, new, obj, ffind=1) {
# do the fffile component
 inst = obj@fffile
 fref = attr(attributes(inst)[["physical"]], "filename")
 ans = gsub(old, new, fref)
 attr(attributes(inst)[["physical"]], "filename") = ans
 obj@fffile = inst
# do the summary component (MAF, RAF) -- still a list but could simplify
 inst = obj@summaryList[[ffind]]
 fref = attr(attributes(inst)[["physical"]], "filename")
 ans = gsub(old, new, fref)
 attr(attributes(inst)[["physical"]], "filename") = ans
 obj@summaryList[[ffind]] = inst
# return
 obj
}
