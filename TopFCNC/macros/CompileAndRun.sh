#!/bin/bash   
if [ "${1}" != "" ]; then

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/lib

  #echo ${LD_LIBRARY_PATH}

  MACRO=$1

  EXECUTABLE=`echo ${1} | cut -d'.' -f1`
    
  echo "Will compile macro ${MACRO}"

  g++ -m32 -I../../../ -I../../ -I.. -L${HOME}/lib -lTopTreeAnaContent42 -lTopTreeAna42 -I `root-config --incdir` `root-config --libs` ${MACRO} -o ${EXECUTABLE}

  echo "Done. Will now run the associated executable ${EXECUTABLE}"

  #./TopFCNC_EventSelection | tee $2
else
  echo "Please, specify a macro (ex : foo.cc)"
fi
