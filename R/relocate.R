
relocate = function(old, new, obj, ffind=1) {
 inst = obj@fflist[[ffind]]
 fref = attr(attributes(inst)[["physical"]], "filename")
 ans = gsub(old, new, fref)
 attr(attributes(inst)[["physical"]], "filename") = ans
 obj@fflist[[ffind]] = inst
 obj
}
