gffprocess = function(basename="fullyri100k", n_in=44, headpatt="_1A", tmpForSort="/freshdata/tmp") {
#
# given a collection of gff3 files representing chunked cis-eQTL search results, order them
# for concatenation and tabix indexing -- requires external classic unix grep and sort utilities
# along with tabix and bgzip
#
# we assume that -T argument to sort is followed by a folder path to be used for temporary files
# more modern sort utilities allow --temporary-directory= 
#
 allgff = dir(patt="gff3$")
 stopifnot(length(allgff)==n_in)
 topind = grep(headpatt, allgff)
 stopifnot(length(topind)==1)
 tocat1 = allgff[-topind]
 tf1 = tempfile()
 tf2 = tempfile()
 run1 = system(paste0("cat ", paste(tocat1, collapse=" "), " | grep -v '^#' > ", tf1))
 stopifnot(run1 == 0)
 run2 = system(paste0("cat ", allgff[topind], " ", tf1, " > ", tf2 ))
 stopifnot(run2 == 0)
#!/bin/bash
 buildgrep = function(infile, outfile) {
#
# sorts and zips also
#
   templ = "(grep ^'#' %INF%; grep -v ^'#' %INF% | sort -T %TDIR% -k1,1 -k4,4n) | bgzip > %OUTF%"
   templ = gsub("%INF%", infile, templ)
   templ = gsub("%OUTF%", outfile, templ)
   templ = gsub("%TDIR%", tmpForSort, templ)
   templ
}
 outname = paste0(basename, ".gff3.gz")
 gcomm = buildgrep(tf2, outname) 
 run3 = system(gcomm) # now outname is bgzipped and sorted
 stopifnot(run3 == 0)
 run4 = system(paste0("tabix -p gff ", outname))
 stopifnot(run4 == 0)
 NULL
}
 

 
