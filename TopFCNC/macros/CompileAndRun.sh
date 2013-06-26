#!/bin/bash   
if [ "${1}" != "" ]; then

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
  #:${HOME}/lib

  #echo ${LD_LIBRARY_PATH}
  LDFLAGS=`root-config --ldflags`
  ROOTINCDIR=`root-config --incdir`
  ROOTLIBS="`root-config --libs` -lTreePlayer"
  
  MACRO=$1

  EXECUTABLE=`echo ${1} | cut -d'.' -f1`
    
  echo "Will compile macro ${MACRO}"

  #g++ ${LDFLAGS} -I../../../ -I../../ -I.. -L${HOME}/lib -lTopTreeAnaContent53 -lTopTreeAna53 -I `root-config --incdir` `root-config --libs` ${MACRO} -o ${EXECUTABLE}
  #g++ ${LDFLAGS} -I../../../ -I../../ -I.. -L${HOME}/lib -lTopFcncAnalysis53 -I `root-config --incdir` `root-config --libs` -lTreePlayer ${MACRO} -o ${EXECUTABLE}
  g++ ${LDFLAGS} -I../../../ -I../../ -I.. -L${HOME}/lib -lTopFcncAnalysis53 -I ${ROOTINCDIR} ${ROOTLIBS} ${MACRO} -o ${EXECUTABLE}

  if [ -e ${EXECUTABLE} ]; then
  #  echo "Done. Will now run the associated executable ${EXECUTABLE}"
    echo "${EXECUTABLE} was sucessfully compiled"
  #./TopFCNC_EventSelection | tee $2
  else
    echo "Compilation failed!"
  fi
else
  echo "Please, specify a macro (ex : foo.cc)"
fi
