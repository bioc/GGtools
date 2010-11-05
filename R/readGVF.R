
readGVF = function(gzfname, ncheck=100, ndata=NULL,
   els=c("ID", "Variant_seq", "Reference_seq", "Total_reads", "Genotype")) {
  con = gzfile(gzfname, "r")
  on.exit(close(con))
  grab = readLines(con, n=ncheck)
  numprag = length(grep("^##", grab))
  numcomm = length(grep("^#[^#]", grab))
  numbl = length(grep("^\ |^$", grab))
  nstuff = list(numprag = numprag, numcomm = numcomm, numbl = numbl)
  toskip = sum(unlist(nstuff))
  close(con)
  con = gzfile(gzfname, "r")
  for (i in 1:toskip)
    tmp = readLines(con, n=1)
  #
  if (is.null(ndata)) {
   nrec = system(paste("zcat ", gzfname, "|wc"), intern=TRUE)
   nrec = gsub("^ *", "", nrec)
   nrec = as.numeric(strsplit(nrec, " ")[[1]][1])
   ndata = nrec - toskip
   }
  start = integer(ndata)
  end = integer(ndata)
  feattype = character(ndata)
  strand = character(ndata)
  source = character(ndata)
  seqname = character(ndata)
  #valdata = matrix("", nc=ndata, nr=length(els))
  #rownames(valdata) = els
  valdl = list()
  for (i in 1:length(els)) valdl[[els[i]]] = character(ndata)
  # first 8 columns are GFF3, 9th is a specific key=value collection, semicolon-separated
  for (i in 1:ndata) {
     if (i %% 1000 == 0) cat(i)
     tmp = strsplit(readLines(con, n=1), "\t")[[1]]
     seqname[i] = tmp[1]
     source[i] = tmp[2]
     feattype[i] = tmp[3]
     start[i] = tmp[4]
     end[i] = tmp[5]
     strand[i] = tmp[7]
     meat = lapply(strsplit(tmp[9], ";")[[1]], function(x)strsplit(x, "=")[[1]])
     vals = sapply(meat, "[", 2)
     names(vals) = sapply(meat, "[", 1)
     #valdata[, i ] = as.character(vals[rownames(valdata)])   
     for (j in 1:length(els)) valdl[[els[j]]][i] = vals[els[j]]
     }
  ans = GRanges( ranges=IRanges(start=as.integer(start), end=as.integer(end)), seqnames=seqname, 
       strand = Rle(strand),
       Variant_seq=Rle(valdl[["Variant_seq"]]), Reference_seq = Rle(valdl[["Reference_seq"]]),
       Total_reads = as.integer(valdl[["Total_reads"]]), Genotype=Rle(valdl[["Genotype"]])  )
  #list(start=start, valdata=valdata, ans=ans)
  ans
}

