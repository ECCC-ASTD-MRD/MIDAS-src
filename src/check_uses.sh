cd modules
grep category *.f*90 |grep MODULE|grep prefix |grep -iv obsspacedata > ../tempfile.txt
cd ../

while IFS= read -r ll
do
    modName=`echo $ll | awk '{print $1}' | sed 's/://g' | sed 's/.f90//g' | sed 's/.ftn90//g'`
    prfName=`echo $ll | awk '{print $5}' | sed 's/(prefix=//g' | sed "s/'//g" | sed 's/"//g'`
    echo "Checking for unneeded 'uses' of **${modName}** with prefix **${prfName}**"

    for i in `grep -li "^ *use * $modName" modules/*.f*90 programs/*.f90`
    do
	searchPrefix=`grep -li "${prfName}_" $i`
	searchStruct=`grep -li "struct_${prfName}" $i`
	if [[ -z ${searchPrefix}${searchStruct} ]]
	then
	    echo "        ***************" $i "****************"
	fi
    done
done < tempfile.txt

rm tempfile.txt
