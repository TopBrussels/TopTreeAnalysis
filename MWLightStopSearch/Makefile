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
#SOURCES         = $(wildcard WHelicities/src/*.cc)
#HEADERS         = $(wildcard WHelicities/interface/*.h)
#OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
#DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
#SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

# headers and sources for madweight_stop analysis lib
SOURCESSTOPDIC	= $(wildcard WHelicities/src/*.cc)
HEADERSSTOPDIC	= $(wildcard WHelicities/interface/*.h)
OBJECTSSTOPDIC	= $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))

SOURCESSTOP     = $(wildcard WHelicities/src/*.cc)
HEADERSSTOP     = $(wildcard WHelicities/interface/*.h)
OBJECTSSTOP	= $(SOURCESSTOP:.$(SrcSuf)=.$(ObjSuf))


all:  libStopAnalysis53.$(DllSuf)
	cp libStopAnalysis53.$(DllSuf) ~/lib/

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTSSTOP) $(OBJECTSTOPDIC)  $(DEPENDS) macros/*.exe *Dict.* *.$(DllSuf) core 

.SUFFIXES: .$(SrcSuf) .C .o .$(DllSuf)

# specific stuff for madweight_stop analysis ONLY

StopDict.$(SrcSuf): $(HEADERSSTOPDIC) ./WHelicityLinkDef.h
	@echo "Generating dictionary StopDict..."
	@rootcint -f StopDict.$(SrcSuf) -c $(DEFINES) $(HEADERSSTOPDIC) ./WHelicityLinkDef.h

libStopAnalysis53.$(DllSuf): $(OBJECTSSTOP) StopDict.o
	@echo "Building libStopAnalysis..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

#ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO
#
#macros/%.exe: macros/%.cc $(HEADERS) libTopTreeAna53.$(DllSuf) libTopTreeAnaContent53.$(DllSuf)
#	$(LD) -lTopTreeAna53 -lTopTreeAnaContent53 $(LIBS) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@
#
#SOURCES_MACROS = $(wildcard macros/*.cc)
#
#macros: $(SOURCES_MACROS:.cc=.exe)

