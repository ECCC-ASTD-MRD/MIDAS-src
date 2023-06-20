
echo ""
echo "Purpose: Scan modules in the 'module_categories.txt' file and identify"
echo "         any module that is used by another module that appears at a"
echo "         lower level in the dependency hierarchy. The number of such"
echo "         modules should be small, however some occurences are inevitable."
echo "         This scripts identifies such occurences with:"
echo "           'OUTOFORDER modNameFile = __lower-level_module_including_the_use__, ..."
echo "                       useFileName = __higher-level_used_module__'"

thisCategory=""
allPrevious=""
categoryNumber=0

while IFS= read -r modNameFile
do

    # skip unimportant lines in the file
    if [[ ${modNameFile} == *"========"* ]]; then
       continue
    fi
    if [[ ${modNameFile} == " *" ]]; then
       continue
    fi
    if [[ -z ${modNameFile} ]]; then
       continue
    fi

    if [[ ${modNameFile} == [0-9]* ]]; then
	# this line is the header for a new category

	let categoryNumber=categoryNumber+1
	echo ""
	echo "===============START OF CATEGORY ${categoryNumber}==============="

	allPrevious="${allPrevious} ${thisCategory}"
	allPrevious=`echo ${allPrevious} | tr " " "\n" | sort -u | tr "\n" " "`
	thisCategory=""
    else
	# this line is a module within the current category

	thisCategory="${thisCategory} ${modNameFile}"
	useList=`grep -i '^ *use .*_mod' ${modNameFile} | awk '{print $2}' | sed 's/,//g' | sort -u`
	for useName in `echo ${useList}`
	do
	    useFileName=`grep -il "^ *module *${useName}" *.f*90`
	    if [ -z "${useFileName}" ]; then
		# no file corresponding to this module because in external library
		echo "SKIPMODULE modNameFile = ${modNameFile}, useName = ${useName}"
	    else
		# check if module file appears in any previous (higher-level) category
		useInPreviousGrep=`echo ${allPrevious} | grep ${useFileName}`
		if [[ ! -z "${useInPreviousGrep}" ]]; then
		    echo "OUTOFORDER modNameFile = ${modNameFile}, useFileName = ${useFileName}"
		fi
	    fi
	done	       
    fi

done < module_categories.txt
