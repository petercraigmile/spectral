
PACKAGE = spectral
VER     = 0.8

check:	;
	R CMD check $(PACKAGE)

build:	;
	R CMD INSTALL $(PACKAGE) 

clear:	;
	rm $(PACKAGE)/src/*.o $(PACKAGE)/src/*.so $(PACKAGE)/src/symbols.rds $(PACKAGE)/src/.DS_Store

