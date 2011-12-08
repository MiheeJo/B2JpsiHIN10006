ROOTCFLAGS	=	$(shell root-config --cflags)
ROOTGLIBS		=	$(shell root-config --glibs)

CPP					=	g++
CPPFLAGS		=	-g -fPIC -Wno-deprecated -O2 -ansi
LD					=	g++
LDFLAGS			=	-g
SOFLAGS			=	-shared

CPPFLAGS		+= $(ROOTCFLAGS)
NGLIBS			=	$(ROOTGLIBS)
NGLIBS			+= -L/afs/cern.ch/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_4/external/slc5_amd64_gcc434/lib -lMathMore -lMinuit -lRooFit -lRooFitCore -lFoam 
GLIBS				= $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR	=	./
CPP					+= -I$(INCLUDEDIR)
CPP					+= -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include
OUTLIB			= ../../../lib

.SUFFIXES:	.cc,.C,.hh,.h
.PREFIXES:	../../../lib

RooHistPdfConv.o: $(INCLUDEDIR)/RooHistPdfConv.cpp
	$(CPP) $(CPPFLAGS) -c -o $(OUTLIB)/libRooHistPdfConv.o $(NGLIBS) $<

Tree2Datasets:	$(INCLUDEDIR)tree2Datasets.cpp
	$(CPP) $(CPPFLAGS) -o Tree2Datasets $(GLIBS) $ $<

Fit2DDataPbPb:	$(INCLUDEDIR)fit2DData_pbpb.cpp
	$(CPP) $(CPPFLAGS) -o Fit2DDataPbPb $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DDataPeePbPb:  $(INCLUDEDIR)fit2DData_pee_pbpb.cpp
	$(CPP) $(CPPFLAGS) -o Fit2DDataPeePbPb $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DJpsiPEE:    $(INCLUDEDIR)fit2DJpsi_PEE.cxx   
		$(CPP) $(CPPFLAGS) -o Fit2DJpsiPEE $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DDatapp:	$(INCLUDEDIR)fit2DData_pp.cpp
	$(CPP) $(CPPFLAGS) -o Fit2DDatapp $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DData:	$(INCLUDEDIR)fit2DJpsi_PEE.cxx
	$(CPP) $(CPPFLAGS) -o Fit2DData $(OUTLIB)/*.o $(GLIBS) $ $<

#Toy2DData:	$(INCLUDEDIR)toy2DData.cpp
#	$(CPP) $(CPPFLAGS) -o Toy2DData $(OUTLIB)/*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIB)*.so
