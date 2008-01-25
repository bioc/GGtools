GGEsite -- README.txt

This is the infrastructure for ggexplorer.org

This assumes that a recent (ca. 2007) version of apache is in use.

Makevars in GGEsite must be modified to reflect the locations
of various R and CGI-related resources

0) R_PATH gives the absolute path to the executable for an R
that has CGIwithR and Bioconductor GGtools and GGdata packages 
installed 

1) HTDOCS_ABSOLUTE dir holds html and png files used for the
presentation layer

2) CGI_ABSOLUTE holds R code for the cgi-bin layer

3) RGRAPHS_RELATIVE is the path to a folder where graphs
can be written by R, relative to HTDOCS_ABSOLUTE

4) We assume you are capable of
running everything in the trivial.* sources supplied with
CGIwithR.

To prepare and install, 
a) modify Makevars to be accurate for your R and Apache 
b) run make
c) with sufficient privileges, run make install
d) point browser to http://localhost/ggexpl.html
