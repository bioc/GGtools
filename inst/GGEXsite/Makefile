include Makevars

subst:
	-rm -rf cgi-bin
	cp -r cgi-bin_virgin cgi-bin
	cd cgi-bin; make substitute

install:
	cp cgi-bin/* cgi-bin/.Rprofile $(CGI_ABSOLUTE)
	cp htdocs/* $(HTDOCS_ABSOLUTE)
