ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lFoam -lTMVA
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lFoam -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES -I..
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD		= g++
LDFLAGS		= -g -O -Wall -fPIC
SOFLAGS		= -shared

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= $(ROOTLIBS) -lEG -I.. -L. -L../TopTreeProducer/src 
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES         = $(wildcard Tools/src/*.cc StatProcedure/src/*.cc BkgEstimationMethods/src/*.cc  Selection/src/*.cc Reconstruction/src/*.cc MCInformation/src/*.cc tinyxml/*.cc KinFitter/src/*.cc JESMeasurement/src/*.cc WHelicities/src/*.cc)
HEADERS         = $(wildcard Tools/interface/*.h StatProcedure/interface/*.h BkgEstimationMethods/interface/*.h  Selection/interface/*.h Reconstruction/interface/*.h MCInformation/interface/*.h tinyxml/*.h Kinfitter/interface/*.h JESMeasurement/interface/*.h WHelicities/interface/*.h)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))
#for libTopTreeAnaContent.so
SOURCESDIC	= $(wildcard Reconstruction/src/Observables.cc Reconstruction/src/MEzCalculator.cc Content/src/*.cc ../TopTreeProducer/src/TRoot*.cc JESMeasurement/src/Monster.cc JESMeasurement/src/LightMonster.cc WHelicities/src/WTree.cc)
HEADERSDIC	= $(wildcard Content/interface/*.h ../TopTreeProducer/interface/TRoot*.h JESMeasurement/interface/Monster.h JESMeasurement/interface/LightMonster.h WHelicities/interface/WTree.h)
OBJECTSDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

# headers and sources for btag eff analysis lib
SOURCESBTAGDIC	= $(wildcard BtagEffAnalysis/src/TRoot*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc)
HEADERSBTAGDIC	= $(wildcard BtagEffAnalysis/interface/TRoot*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h)
OBJECTSBTAGDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESBTAG         = $(wildcard BtagEffAnalysis/src/*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc Tools/src/PlottingTools.cc)
HEADERSBTAG         = $(wildcard BtagEffAnalysis/interface/*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h Tools/interface/PlottingTools.h)
OBJECTSBTAG		= $(SOURCESBTAG:.$(SrcSuf)=.$(ObjSuf))


all:  libTopTreeAna42.so libTopTreeAnaContent42.so ;  cp libTopTreeAna42.so ~/lib/ ; cp libTopTreeAnaContent42.so ~/lib/

btag: libBtagAnalysis42.so; cp libBtagAnalysis42.so ~/lib/

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(OBJECTSDIC) $(OBJECTSBTAG) $(OBJECTBTAGDIC) $(DEPENDS) *Dict.* core 

.SUFFIXES: .$(SrcSuf) .C .o .so

###

Dict.$(SrcSuf): $(HEADERSDIC) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERSDIC) ./LinkDef.h

libTopTreeAna42.so: $(OBJECTS) 
	@echo "Building libTopTreeAna42..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

libTopTreeAnaContent42.so: $(OBJECTSDIC)  Dict.o  
	@echo "Building libTopTreeAnaContent42..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

# specific stuff for btag eff analysis ONLY

BtagDict.$(SrcSuf): $(HEADERSBTAGDIC) ./BtagLinkDef.h
	@echo "Generating dictionary BtagDict..."
	@rootcint -f BtagDict.$(SrcSuf) -c $(DEFINES) $(HEADERSBTAGDIC) ./BtagLinkDef.h

libBtagAnalysis42.so: $(OBJECTSBTAG) BtagDict.o
	@echo "Building libBtagAnalysis..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@


macros/Demo.exe: macros/Demo_binning.cc config/Datasets.cc $(HEADERS) libTopTreeAna42.so libTopTreeAnaContent42.so
	$(LD) -lTopTreeAna -lTopTreeAnaContent $(LIBS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

macros/Full.exe: macros/FullChainEstimation.cc $(HEADERS) libTopTreeAna42.so libTopTreeAnaContent42.so
	$(LD) -lTopTreeAna -lTopTreeAnaContent $(LIBS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

macros/Cross.exe: macros/CrossSectionMeasurement.cc $(HEADERS) libTopTreeAna42.so libTopTreeAnaContent42.so
	$(LD) -lTopTreeAna -lTopTreeAnaContent $(LIBS) -I`root-config --incdir` $< $(LDFLAGS) -o $@
