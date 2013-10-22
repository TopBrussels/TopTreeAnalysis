#when an argument (mymacro.cc) is given ($1) do the compiling
if [ "${1}" != "" ]; then 
	#export LD_LIBRAY_PATH=$LD_LIBRARY_PATH:.
	export LD_LIBRAY_PATH=$LD_LIBRARY_PATH:~/lib
	LDFLAGS= `root-config --ldflags`
	ROOTINCDIR= `root-config --incdir`
	#ROOTLIBS= "`root-config --libs` -lTreePlayer"
	ROOTLIBSS= `root-config --libs`
	
	MACRO=$1
	
	EXECUTABLE=`echo ${1} | cut -d'.' -f1`
	
	echo "Will compile macro ${MACRO}"
	
	#compiling the macro: 
	#g++ ${LDFLAGS} -I../../../ -I../../ -I.. -L${HOME}/lib	-lTopFcncAnalysis53 -I ${ROOTINCDIR} ${ROOTLIBS} ${MACRO} -o ${EXECUTABLE}
	
	g++ -g -L ~/lib -I../../../ -I../../ -I.. -L${HOME}/lib	-lTopTreeAnaContent53 -lTopTreeAna53 -lMLP -lTreePlayer -XMLIO -I ${ROOTINCDIR} ${ROOTLIBSS} -I../../ -L. -L../TopTreeProducer/src  ${MACRO} -o	${EXECUTABLE}
	
	if [ -e ${EXECUTABLE} $; then
		echo "${EXECUTABLE} was succesfully compiled"
	else
		echo "Compilation failed"
	fi


#when no macro is given, ask for it 
else
	echo "Please specify a macro (eg: mymacro.cc)"
fi
