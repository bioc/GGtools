

interleave2cis = function( cisp, permcisp ) {
 if (length(cisp) != length(permcisp)) stop("inputs must have same length")
 newl = list()
 for (i in seq(2, 2*length(cisp), 2)) {
 newl[[i-1]] = cisp[[i/2]]
 newl[[i]] = permcisp[[i/2]]
 }
 names(newl) = rep(names(cisp), each=2)
 names(newl)[ seq(2,2*length(cisp),2) ] =
   paste( "p_", names(newl)[ seq(2,2*length(cisp),2) ], sep="" )
 newl
}

