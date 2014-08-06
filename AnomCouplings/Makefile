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

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
ROOTLIBS_NoTMVA = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam
ROOTGLIBS       = $(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA

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

CXXFLAGS	     += $(ROOTCFLAGS)
LIBS		        = -I./TMVA/include -L./TMVA/lib $(ROOTLIBS) -lEG -I.. -L. -L../TopTreeProducer/src 
LIBS_NoTMVA     = $(ROOTLIBS_NoTMVA) -lEG -I.. -L. -L../TopTreeProducer/src 
ifeq ($(UNAME), Darwin)
LIBS           += -I/opt/local/include
endif
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
#SOURCES         = $(wildcard Tools/src/*.cc StatProcedure/src/*.cc BkgEstimationMethods/src/*.cc  Selection/src/*.cc Reconstruction/src/*.cc MCInformation/src/*.cc tinyxml/*.cc KinFitter/src/*.cc JESMeasurement/src/*.cc WHelicities/src/*.cc InclFourthGenSearch/src/*.cc)
#HEADERS         = $(wildcard Tools/interface/*.h StatProcedure/interface/*.h BkgEstimationMethods/interface/*.h  Selection/interface/*.h Reconstruction/interface/*.h MCInformation/interface/*.h tinyxml/*.h Kinfitter/interface/*.h JESMeasurement/interface/*.h WHelicities/interface/*.h InclFourthGenSearch/interface/*.h)
#OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
#DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
#SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

# headers and sources for ANOMCPL analysis lib
#SOURCESANOMCPLDIC	= $(wildcard AnomCouplings/src/AnomCouplTree.cc)
#HEADERSANOMCPLDIC	= $(wildcard AnomCouplings/interface/AnomCouplTree.h)
#OBJECTSANOMCPLDIC	= $(SOURCESANOMCPLDIC:.$(SrcSuf)=.$(ObjSuf))
SOURCESANOMCPLDIC	= $(wildcard AnomCouplings/src/*.cc)
HEADERSANOMCPLDIC	= $(wildcard AnomCouplings/interface/*.h)
OBJECTSANOMCPLDIC	= $(SOURCESANOMCPLDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESANOMCPL   = $(wildcard AnomCouplings/src/*.cc)
HEADERSANOMCPL   = $(wildcard AnomCouplings/interface/*.h)
OBJECTSANOMCPL	 = $(SOURCESANOMCPL:.$(SrcSuf)=.$(ObjSuf))

anomcpl: libANOMCPLAnalysis53.$(DllSuf)
	cp libANOMCPLAnalysis53.$(DllSuf) ~/lib/

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(OBJECTSDIC) $(OBJECTSBTAG) $(OBJECTBTAGDIC) $(OBJECTSMTOP) $(OBJECTSMTOPDIC) $(OBJECTSTOPFCNC) $(OBJECTSTOPFCNCDIC) $(DEPENDS) macros/*.exe *Dict.* *.$(DllSuf) core 

.SUFFIXES: .$(SrcSuf) .C .o .$(DllSuf)

###

# specific stuff for AnomCouplings analysis ONLY

ANOMCPLDict.$(SrcSuf): $(HEADERSANOMCPLDIC) ./AnomCouplLinkDef.h
	@echo "Generating dictionary ANOMCPLDict..."
	@rootcint -f ANOMCPLDict.$(SrcSuf) -c $(DEFINES) $(HEADERSANOMCPLDIC) ./AnomCouplLinkDef.h

libANOMCPLAnalysis53.$(DllSuf): $(OBJECTSANOMCPL) ANOMCPLDict.o
	@echo "Building libVLQAnalysis..."
	$(LD) $(LIBS_NoTMVA) $(SOFLAGS) $(LDFLAGS) $+ -o $@

ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO

macros/%.exe: macros/%.cc $(HEADERS) libTopTreeAna53.$(DllSuf) libTopTreeAnaContent53.$(DllSuf)
	$(LD) -lTopTreeAna53 -lTopTreeAnaContent53 $(LIBS) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

SOURCES_MACROS = $(wildcard macros/*.cc)

macros: $(SOURCES_MACROS:.cc=.exe)

