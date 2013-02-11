export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib

stamp="$1_$(date +%s)"

g++ -g -L ~/lib -lTopTreeAnaContent53 -lTopTreeAna53 -lMLP  -lTreePlayer -lXMLIO -I `root-config --incdir` `root-config --libs` -I../../ -L. -L../TopTreeProducer/src $1.cc -o $stamp

#echo $* | awk -F "$1 " {print'$2'} 

cp -v $stamp $1

if [ $? == "0" ];then

    echo "Running ./$stamp $(echo $* | awk -F "$1 " {print'$2'})"
    sleep 1
    #./$* | awk -F "$1 " {print'$2'}
    ./$stamp $(echo $* | awk -F "$1 " {print'$2'})

fi

echo "removing temp executable $stamp"

rm $stamp -v
