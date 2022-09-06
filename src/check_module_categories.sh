
thisCategory=""
allPrevious=""
firstTime=true

#for modNameFile in `cat module_categories.txt`;
while IFS= read -r modNameFile
do

    if [[ ${modNameFile} == *"==="* ]]; then
       continue
    fi
    if [[ ${modNameFile} == " *" ]]; then
       continue
    fi
    if [[ -z ${modNameFile} ]]; then
       continue
    fi

    #echo ...${modNameFile}...

    if [[ ${modNameFile} == [0-9]* ]]; then
	if [ "${firstTime}" = false ]; then
	    echo ""
	    echo "====START OF A NEW CATEGORY==============="
	fi
	allPrevious="${allPrevious} ${thisCategory}"
	allPrevious=`echo ${allPrevious} | tr " " "\n" | sort -u | tr "\n" " "`
	thisCategory=""
	firstTime=false
    else
	thisCategory="${thisCategory} ${modNameFile}"
	useList=`grep -i '^ *use .*_mod' modules/${modNameFile} | awk '{print $2}' | sed 's/,//g' | sort -u`
	for useName in `echo ${useList}`
	do
	    useFileName=`grep -il "^ *module *${useName}" modules/*.f*90 | sed 's/modules\///g'`
	    if [ -z "${useFileName}" ]; then
		echo "SKIP checking module: ${useName} used in ${modNameFile}"
	    else
		useInPreviousGrep=`echo ${allPrevious} | grep ${useFileName}`
		if [[ ! -z "${useInPreviousGrep}" ]]; then
		    echo "OUTOFORDER modNameFile = ${modNameFile}, useFileName = ${useFileName}"
		fi
	    fi
	done	       
    fi

done < module_categories.txt
