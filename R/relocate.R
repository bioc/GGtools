
relocate = function(old, new, obj, ffind=1) {
# do the fflist component
 inst = obj@fflist[[ffind]]
 fref = attr(attributes(inst)[["physical"]], "filename")
 ans = gsub(old, new, fref)
 attr(attributes(inst)[["physical"]], "filename") = ans
 obj@fflist[[ffind]] = inst
# do the summary component (MAF, RAF)
 inst = obj@summaryList[[ffind]]
 fref = attr(attributes(inst)[["physical"]], "filename")
 ans = gsub(old, new, fref)
 attr(attributes(inst)[["physical"]], "filename") = ans
 obj@summaryList[[ffind]] = inst
# return
 obj
}
