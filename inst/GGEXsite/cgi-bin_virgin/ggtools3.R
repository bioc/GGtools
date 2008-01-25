#! @RPATH@
#
# YOU NEED TO BE SURE ABOVE IS A PATH TO AN R including CGIwithR and GGtools
# and GGdata
#

# THESE NEED TO BE WRITEABLE AND CONSISTENT WITH APACHE SETUP

        globGraphDir <- "@GLOB_GRAPH_DIR@"
        graphDirURLroot <- "@GRAPH_DIR_URL_ROOT@"
        graphURLroot <- "@GRAPH_URL_ROOT@"

# NO MORE CHANGES


tag(html)
tag(head)

library(GGtools)
library(GGdata)
tag('body style="color: rgb(0, 0, 0)" alink="#358cff" link="#358cff" 
vlink="#990099"')
tag(center)

uniqueNames = function(n) {
 tags = sample(letters, n, replace=TRUE)
 un = tempfile(tags)
 gsub("\\/", "", un)
}

meta=paste("chr",formData$chromosome,"meta",sep="");
genes=paste("chr",formData$chromosome, "GGceuRMA",sep="");
gene=formData$gene;
data(geneLocs_hsa)
 
chrString = paste("chr", formData$chromosome, sep="")

xx <- unlist(as.list(hgfocusSYMBOL))
 if(any(gene %in% xx)){
  if(as.numeric(formData$nsnps)<=50){

    num=as.numeric(formData$nsnps)+2
    uNames = uniqueNames(num)
file1= paste(globGraphDir, "/", uNames[2], ".txt", sep="")
    nums=floor(rnorm(num, mean=1000000, sd=158))
    #mainjpg=paste("mygraph",nums[1],".jpg",sep="");
    
    data(list=meta)
    data(list=genes)
    scr1 = snpScreen(get(genes), get(meta), genesym(gene), ~., 
               fastAGMfitter)
    webPNG( paste(uNames[1],"png", sep="."), graphDir = globGraphDir )
    plot_mlp(scr1, get(meta), geneLocDF=geneLocs_hsa)
    
    pvals=scr1[[4]];
    pvals=sort(pvals);
    snpnames=names(pvals);

tag(big)
cat("Parameters:")
untag(big)
tag(br)
tag(br)
cat("Chromosome: ")
cat(formData$chromosome)
tag(br)
cat("Gene: ")
g=paste('<a 
href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',formData$gene,'">',formData$gene,'</a>',sep="")
cat(g)
tag(br)
cat("Number of SNPs: ")
cat(formData$nsnps)
tag(br)
tag(br)

tag(big)
cat("Model Fitting Info:")
untag(big)
tag(br)
tag(br)
print(scr1)
tag(br)
tag(br)
tag(big)
cat("Results:")
untag(big)
a=paste('<a href="', graphDirURLroot, '',uNames[2],'.txt">TEXT 
</a>',sep="");
cat(a) 
a=paste('<a href="', graphDirURLroot, '',uNames[2],'.rdf">RDF </a>',sep="");
cat(a) 
tag("table border=0")
tag(tr)
tag(td)
imgname=paste(uNames[1],".png",sep="")
#print(imgname)
img(imgname)
untag(td)
tag(td)
tag("TABLE border=1")
tag(tr)
tag(td)
cat("SNP ID")
untag(td)
tag(td)
cat("P Value")
untag(td)
untag(tr)

# SNP link table

for(i in 1:as.numeric(formData$nsnps)){
	tag(tr)
	tag(td)
	a=snpnames[i]
	b=strsplit(a,"rs")[[1]][2]
	b=paste('<a 
href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',b,'">',a,'</a>',sep="")
	cat(b)
	untag(td)
	tag(td)
	a=pvals[i]
	cat(a)
	untag(td)
	untag(tr)
	}
untag(table)
untag(td)
untag(tr)
untag(table)

# SNP pictures

nrows=ceiling(as.numeric(formData$nsnps)/2);
tag("TABLE border=0")
tag(br)
tag(br)
k = 2
for(i in 1:nrows){
	tag("TR")
	for(j in 1:2){
		if(((i-1)*2+j)<=as.numeric(formData$nsnps)){
		tag(td)
		graphnum=(i-1)*2+j+2
		fname=paste("mygraph",i,j,graphnum,sep="")		
		file1= globGraphDir;
                k = k+1
		webPNG(curfn <- paste(uNames[k], ".png", sep=""), graphDir = file1)

		 plot_EvG(get(genes), genesym(gene), snpnames[(i-1)*2+j], "hgfocus")



		dir=graphDirURLroot;
		img( curfn , width="375", height="300")
		untag(td)
		}
	untag(tr)	
	}
}

# text report generation

file1= paste(globGraphDir, "/", uNames[2], ".txt", sep="")
zz <- file(file1, "w")
cat("GG Tools Report\n***************\n\nParameters:\n\nGene: ", file = 
zz) 
cat(formData$gene, file = zz) 
cat("\nChromosome: ", file = zz) 
cat(formData$chromosome , file = zz) 
cat("\n", file = zz) 
cat("\n****\nData structure:\n", file=zz)
thold = textConnection("sap", "w")
sink(thold)
show(get(genes))
sink()
writeLines(sap, con=zz)
cat("\n*****\n", file = zz) 
cat("Number of SNPS: ", file = zz)
cat(formData$nsnps, file = zz)
cat("\n",file=zz)
cat("Call:\n ",file=zz)
f=scr1@call
f=as.character(f)
f2=paste('snpScreen(racExSet=',f[2],', snpMeta=',f[3],', 
gene=',f[4],',"\n", formTemplate=',f[5],', fitter = fastAGM, gran = 
1)',sep="")
cat(f2,file=zz)
cat("\n",file=zz)
cat("Number of fits attempted: ",file=zz)
f=length(scr1[[1]]);
cat(f,file=zz)
cat("\n",file=zz)
cat("Number of Fits Successful: ",file=zz)
f=sum(!is.na(scr1[[1]]))
cat(f,file=zz);
cat("\n\nResults:\n\n",file=zz)
cat("SNP ID : P-value",file=zz)
cat("\n",file=zz)
for(i in 1:as.numeric(formData$nsnps)){
	a=snpnames[i]
	cat(a,file=zz)
	a=pvals[i]
	cat(" : ",file=zz)
	cat(a,file=zz)
	cat("\n",file=zz)
	}
close(zz)

######################
##RDF GENERATION
######################
file2= paste(globGraphDir, "/", uNames[2], ".rdf", sep="")
zzz <- file(file2, "w")
cat('<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
             xmlns:ggexp="http://www.ggexplorer.org/terms/#">

<rdf:Description rdf:about="http://www.ggexplorer.org/terms#config">
        <ggexp:gene>',file=zzz)
cat(formData$gene,file=zzz)
cat('</ggexp:gene>
	<ggexp:chromosome>',file=zzz)
cat(formData$chromosome,file=zzz)
cat('</ggexp:chromosome>
	<ggexp:nsnps>',file=zzz)
cat(formData$nsnps,file=zzz)
cat('</ggexp:nsnps>
</rdf:Description>

<rdf:Description rdf:about="http://www.ggexplorer.org/terms#experiment">
	<ggexp:experimenter>Cheung VG</ggexp:experimenter>
	<ggexp:lab>Department of Pediatrics, University of Pennsylvania, 
Philadelphia, Pennsylvania 19104, USA. 
vcheung@mail.med.upenn.edu</ggexp:lab>
	<ggexp:snpfeatures>117417</ggexp:snpfeatures>
	<ggexp:snpsamples>58</ggexp:snpsamples>
	<ggexp:genefeatures>8793</ggexp:genefeatures>
	<ggexp:genesamples>58</ggexp:genesamples>
	<ggexp:exptitle>Title: Mapping determinants of human gene 
expression by regional and genome-wide association.</ggexp:exptitle> 
</rdf:Description>


<rdf:Description rdf:about="http://www.ggexplorer.org/terms#modelinfo">
	<ggexp:call>snpScreen(racExSet=racExSet, snpMeta=snpMeta, 
gene=gene,"", formTemplate=formTemplate, fitter = fastAGM, gran = 
1)</ggexp:call>',file=zzz)
cat('<ggexp:fitattempts>',file=zzz)
f=length(scr1[[1]]);
cat(f,file=zzz)
cat('</ggexp:fitattempts>
	<ggexp:successfulfits>',file=zzz)
f=sum(!is.na(scr1[[1]]))
cat(f,file=zzz)
cat('</ggexp:successfulfits>
</rdf:Description>

<rdf:Description rdf:about="http://www.ggexplorer.org/terms#data">
	<ggexp:topsnps>
	<rdf:Seq>',file=zzz)
for(i in 1:as.numeric(formData$nsnps)){
	a=snpnames[i];
	b=paste('<rdf:li rdf:resource="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',a,'"/>',"\n",sep="");
	cat(b,file=zzz)}
cat('</rdf:Seq>
	</ggexp:topsnps>
</rdf:Description>',file=zzz)
cat("\n",file=zzz)
for(i in 1:as.numeric(formData$nsnps)){
	a=snpnames[i];
	b=paste('<rdf:Description rdf:about="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',a,'">',"\n",sep="");
	cat(b,file=zzz)
	cat('<ggexp:pvalue>',file=zzz)
	a=pvals[i]
	cat(a,file=zzz)
	cat('</ggexp:pvalue>',file=zzz)
	cat("\n",file=zzz)
cat('</rdf:Description>',file=zzz)
cat("\n",file=zzz)}
cat('</rdf:RDF>',file=zzz)
close(zzz)
#########################################
#END RDF
#########################################

untag(table)
}else{cat("Please choose less than 50 snps")}
}else{cat("Please re-enter the symbol")}
untag(body)
untag(html)





