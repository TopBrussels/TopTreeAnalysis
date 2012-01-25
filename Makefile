ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
UNAME = $(shell uname -s)
ifeq ($(UNAME), Darwin)
DllSuf        = dylib
endif
ifeq ($(UNAME), Linux)
DllSuf        = so
endif

OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lFoam -lTMVA
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lFoam -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES -I.. -I./TMVA/include
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD		= g++
LDFLAGS		= -g -O -Wall -fPIC
ifeq ($(UNAME), Darwin)
SOFLAGS         = -dynamiclib
endif
ifeq ($(UNAME), Linux)
SOFLAGS         = -shared
endif

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= -I./TMVA/include -L./TMVA/lib $(ROOTLIBS) -lEG -I.. -L. -L../TopTreeProducer/src 
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES         = $(wildcard Tools/src/*.cc StatProcedure/src/*.cc BkgEstimationMethods/src/*.cc  Selection/src/*.cc Reconstruction/src/*.cc MCInformation/src/*.cc tinyxml/*.cc KinFitter/src/*.cc JESMeasurement/src/*.cc WHelicities/src/*.cc TopFCNC/src/*.cc)
HEADERS         = $(wildcard Tools/interface/*.h StatProcedure/interface/*.h BkgEstimationMethods/interface/*.h  Selection/interface/*.h Reconstruction/interface/*.h MCInformation/interface/*.h tinyxml/*.h Kinfitter/interface/*.h JESMeasurement/interface/*.h WHelicities/interface/*.h TopFCNC/interface/*.h)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))
#for libTopTreeAnaContent.so
SOURCESDIC	= $(wildcard Reconstruction/src/Observables.cc Reconstruction/src/MEzCalculator.cc Content/src/*.cc ../TopTreeProducer/src/TRoot*.cc JESMeasurement/src/Monster.cc JESMeasurement/src/LightMonster.cc WHelicities/src/WTree.cc TopFCNC/src/*.cc)
HEADERSDIC	= $(wildcard Content/interface/*.h ../TopTreeProducer/interface/TRoot*.h JESMeasurement/interface/Monster.h JESMeasurement/interface/LightMonster.h WHelicities/interface/WTree.h TopFCNC/interface/*.h)
OBJECTSDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

# headers and sources for btag eff analysis lib
SOURCESBTAGDIC	= $(wildcard BtagEffAnalysis/src/TRoot*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc)
HEADERSBTAGDIC	= $(wildcard BtagEffAnalysis/interface/TRoot*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h)
OBJECTSBTAGDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESBTAG         = $(wildcard BtagEffAnalysis/src/*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc Tools/src/PlottingTools.cc)
HEADERSBTAG         = $(wildcard BtagEffAnalysis/interface/*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h Tools/interface/PlottingTools.h)
OBJECTSBTAG		= $(SOURCESBTAG:.$(SrcSuf)=.$(ObjSuf))


all:  libTopTreeAna42.$(DllSuf) libTopTreeAnaContent42.$(DllSuf)
	cp libTopTreeAna42.$(DllSuf) ~/lib/ ; cp libTopTreeAnaContent42.$(DllSuf) ~/lib/

btag: libBtagAnalysis42.$(DllSuf)
	cp libBtagAnalysis42.$(DllSuf) ~/lib/

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(OBJECTSDIC) $(OBJECTSBTAG) $(OBJECTBTAGDIC) $(DEPENDS) macros/*.exe *Dict.* *.$(DllSuf) core 

.SUFFIXES: .$(SrcSuf) .C .o .$(DllSuf)

###

Dict.$(SrcSuf): $(HEADERSDIC) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERSDIC) ./LinkDef.h

libTopTreeAna42.$(DllSuf): $(OBJECTS) libTopTreeAnaContent42.$(DllSuf)
	@echo "Building libTopTreeAna42..."
	$(LD) $(LIBS) -lTopTreeAnaContent42 $(SOFLAGS) $(LDFLAGS) $+ -o $@

libTopTreeAnaContent42.$(DllSuf): $(OBJECTSDIC)  Dict.o  
	@echo "Building libTopTreeAnaContent42..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

# specific stuff for btag eff analysis ONLY

BtagDict.$(SrcSuf): $(HEADERSBTAGDIC) ./BtagLinkDef.h
	@echo "Generating dictionary BtagDict..."
	@rootcint -f BtagDict.$(SrcSuf) -c $(DEFINES) $(HEADERSBTAGDIC) ./BtagLinkDef.h

libBtagAnalysis42.$(DllSuf): $(OBJECTSBTAG) BtagDict.o
	@echo "Building libBtagAnalysis..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO

macros/Demo.exe: macros/Demo_binning.cc config/Datasets.cc $(HEADERS) libTopTreeAna42.$(DllSuf) libTopTreeAnaContent42.$(DllSuf)
	$(LD) -lTopTreeAna -lTopTreeAnaContent $(LIBS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

macros/%.exe: macros/%.cc $(HEADERS) libTopTreeAna42.$(DllSuf) libTopTreeAnaContent42.$(DllSuf)
	$(LD) -lTopTreeAna42 -lBtagAnalysis42 -lTopTreeAnaContent42 $(LIBS) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

SOURCES_MACROS = $(wildcard macros/*.cc)

macros: $(SOURCES_MACROS:.cc=.exe)

