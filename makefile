
PACKAGE = spectral
VER     = 0.8.1

check:	;
	R CMD check $(PACKAGE)

build:	;
	R CMD INSTALL $(PACKAGE)
	R CMD BUILD $(PACKAGE) 

clear:	;
	rm $(PACKAGE)/src/*.o $(PACKAGE)/src/*.so $(PACKAGE)/src/symbols.rds $(PACKAGE)/src/.DS_Store

