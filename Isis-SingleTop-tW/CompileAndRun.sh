export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib

g++ -g -L ~/lib -lTopTreeAnaContent53 -lTopTreeAna53 -lMLP  -lTreePlayer -lXMLIO -I `root-config --incdir` `root-config --libs` -I../../ -L. -L../TopTreeProducer/src $1.cc -o $1

if [ $? == "0" ];then

    ./$1 $(echo $* | awk -F "$1 " {print'$2'})

fi

