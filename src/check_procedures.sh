#
# check_procedures.sh:
#
#     example: check_procedures.sh "^ *implicit *none"
#
searchString="$1"

echo
echo "********************************************"
echo "*** search string = "\"${searchString}\" "***"
echo "********************************************"

# loop over all files
for file in */*.f90; 
do
    echo
    echo " ----------- "
    echo $file
    echo " ----------- "

    # loop over all subroutines
    for proc in `grep -hi "^ *subroutine" $file |sed 's/ //g' |sed 's/subroutine/ /gI'| cut -d'(' -f1| awk '{ print tolower($0) }'`
    do
	# check for some things in this subroutine
	echo "subroutine $proc"
	cat $file | awk '{ print tolower($0) }' | awk "/^ *subroutine *${proc} *\(|^ *subroutine *${proc} *$/,/^ *end *subroutine *${proc} *$/" |grep -i "${searchString}"
    done

    # loop over all functions
    for proc in `grep -hi "^ *function" $file |sed 's/ //g' |sed 's/function/ /gI'| cut -d'(' -f1| awk '{ print tolower($0) }'`
    do
	# check for some things in this function
	echo "function $proc"
	cat $file | awk '{ print tolower($0) }' | awk "/^ *function *${proc} *\(|^ *function *${proc} *$/,/^ *end *function *${proc} *$/" |grep -i "${searchString}"
    done
done
