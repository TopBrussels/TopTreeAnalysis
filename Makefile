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
ROOTLIBS      = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES -I..
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)  -I./TMVA/include
ifeq ($(UNAME), Darwin)
CXXFLAGS        += -I/opt/local/include
endif
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
ifeq ($(UNAME), Darwin)
LIBS            += -I/opt/local/include
endif
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES         = $(wildcard Tools/src/*.cc StatProcedure/src/*.cc BkgEstimationMethods/src/*.cc  Selection/src/*.cc Reconstruction/src/*.cc MCInformation/src/*.cc tinyxml/*.cc KinFitter/src/*.cc JESMeasurement/src/*.cc WHelicities/src/*.cc TopFCNC/src/*.cc InclFourthGenSearch/src/*.cc StopSearchesBG/src/*.cc)
HEADERS         = $(wildcard Tools/interface/*.h StatProcedure/interface/*.h BkgEstimationMethods/interface/*.h  Selection/interface/*.h Reconstruction/interface/*.h MCInformation/interface/*.h tinyxml/*.h Kinfitter/interface/*.h JESMeasurement/interface/*.h WHelicities/interface/*.h TopFCNC/interface/*.h InclFourthGenSearch/interface/*.h StopSearchesBG/interface/*.h)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))
#for libTopTreeAnaContent.so
SOURCESDIC	= $(wildcard Reconstruction/src/Observables.cc Reconstruction/src/MEzCalculator.cc Content/src/*.cc ../TopTreeProducer/src/TRoot*.cc JESMeasurement/src/Monster.cc JESMeasurement/src/LightMonster.cc WHelicities/src/WTree.cc TopFCNC/src/TopFCNC_Evt.cc InclFourthGenSearch/src/InclFourthGenTree.cc StopSearchesBG/src/*.cc BkgEstimationMethods/src/VJetEstimation.cc)
HEADERSDIC	= $(wildcard Content/interface/*.h ../TopTreeProducer/interface/TRoot*.h JESMeasurement/interface/Monster.h JESMeasurement/interface/LightMonster.h WHelicities/interface/WTree.h TopFCNC/interface/TopFCNC_Evt.h InclFourthGenSearch/interface/InclFourthGenTree.h StopSearchesBG/interface/*.h BkgEstimationMethods/interface/VJetEstimation.h)
OBJECTSDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

# headers and sources for btag eff analysis lib
SOURCESBTAGDIC	= $(wildcard BtagEffAnalysis/src/TRoot*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc)
HEADERSBTAGDIC	= $(wildcard BtagEffAnalysis/interface/TRoot*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h)
OBJECTSBTAGDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESBTAG         = $(wildcard BtagEffAnalysis/src/*.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc Tools/src/PlottingTools.cc)
HEADERSBTAG         = $(wildcard BtagEffAnalysis/interface/*.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h Tools/interface/PlottingTools.h)
OBJECTSBTAG		= $(SOURCESBTAG:.$(SrcSuf)=.$(ObjSuf))

# headers and sources for mtop analysis lib
SOURCESMTOPDIC	= $(wildcard JESMeasurement/src/LightMonster.cc Content/src/Dataset.cc Content/src/AnalysisEnvironment.cc)
HEADERSMTOPDIC	= $(wildcard JESMeasurement/interface/LightMonster.h Content/interface/Dataset.h Content/interface/AnalysisEnvironment.h)
OBJECTSMTOPDIC	= $(SOURCESMTOPDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESMTOP   = $(wildcard JESMeasurement/src/LightMonster.cc Tools/src/MultiSamplePlot.cc Content/src/Dataset.cc Tools/src/PlottingTools.cc Content/src/AnalysisEnvironment.cc MCInformation/src/*ReWeighting.cc MCInformation/src/ResolutionFit.cc)
HEADERSMTOP   = $(wildcard JESMeasurement/interface/LightMonster.h Tools/interface/MultiSamplePlot.h Content/interface/Dataset.h Tools/interface/PlottingTools.h Content/interface/AnalysisEnvironment.h MCInformation/interface/*ReWeighting.h MCInformation/interface/ResolutionFit.h)
OBJECTSMTOP		= $(SOURCESMTOP:.$(SrcSuf)=.$(ObjSuf))


all:  libTopTreeAna42.$(DllSuf) libTopTreeAnaContent42.$(DllSuf)
	cp libTopTreeAna42.$(DllSuf) ~/lib/ ; cp libTopTreeAnaContent42.$(DllSuf) ~/lib/

btag: libBtagAnalysis42.$(DllSuf)
	cp libBtagAnalysis42.$(DllSuf) ~/lib/

mtop: libMtopAnalysis42.$(DllSuf)
	cp libMtopAnalysis42.$(DllSuf) ~/lib/

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

# specific stuff for mtop analysis ONLY

MtopDict.$(SrcSuf): $(HEADERSMTOPDIC) ./MtopLinkDef.h
	@echo "Generating dictionary MtopDict..."
	@rootcint -f MtopDict.$(SrcSuf) -c $(DEFINES) $(HEADERSMTOPDIC) ./MtopLinkDef.h

libMtopAnalysis42.$(DllSuf): $(OBJECTSMTOP) MtopDict.o
	@echo "Building libMtopAnalysis..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@


ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO

macros/%.exe: macros/%.cc $(HEADERS) libTopTreeAna42.$(DllSuf) libTopTreeAnaContent42.$(DllSuf)
	$(LD) -lTopTreeAna42 -lTopTreeAnaContent42 $(LIBS) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

SOURCES_MACROS = $(wildcard macros/*.cc)

macros: $(SOURCES_MACROS:.cc=.exe)

