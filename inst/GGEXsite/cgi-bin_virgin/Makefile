# DO NOT CHANGE THIS FILE UNLESS YOU ABSOLUTELY KNOW WHAT YOU ARE DOING
# MAKE ANY CHANGES TO ../Makevars

include ../Makevars

RPATH=$(R_PATH)
SUBS1="s.@RPATH@.$(RPATH)." 
SUBS2="s.@GLOB_GRAPH_DIR@.$(HTDOCS_ABSOLUTE)$(RGRAPHS_RELATIVE)."
SUBS3="s.@GRAPH_DIR_URL_ROOT@.$(RGRAPHS_RELATIVE)."
SUBS4="s.@GRAPH_URL_ROOT@.$(RGRAPHS_RELATIVE)."

RFILES=chrlist.R       ggtools3.R      keggsearch2.R   symbolsearch.R \
	ggtools.R       keggsearch1.R   nsnps.R         trivial.R

.sedfi:
	echo $(SUBS1) > .sedfi
	echo $(SUBS2) >> .sedfi
	echo $(SUBS3) >> .sedfi
	echo $(SUBS4) >> .sedfi

substitute: .sedfi
	for X in $(RFILES) ; do \
		echo $$X ; \
		sed -i bak -f .sedfi $$X ; \
		done ;
	-rm -rf .sedfi
