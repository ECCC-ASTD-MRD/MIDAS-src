#!/bin/bash

#set -e

# ARGUMENTS SPECIFY DIRECTORIES OF MIDAS SOURCE CODE AND DESTINATION

codedir=${1:-../../src}
docdir=${2:-${PWD}}

echo "Generating namelist information for the source directory ${codedir}"
echo "Output will appear in the directory ${docdir}"

cat <<EOF > $docdir/namelists_html
  <HTML style="font-family:arial">
  <HEAD>
     <TITLE> MIDAS fortran software namelists </TITLE>
  </HEAD>
  <H1> Namelists used in MIDAS fortran software </H1>
EOF

cd $codedir
directories=""
for filename in ./* ; do 
  if [[ -d $filename ]]; then
    directories="$directories $filename"
  fi
done
for filename in ./*/* ; do 
  if [[ -d $filename ]]; then
    directories="$directories $filename"
  fi
done
echo "Directories found: $directories"

namelistnames=''
for directory in $directories ; do
firsttime="yes"
for filename in $directory/*.f* ; do
  for namelistname in `grep -hi 'namelist */' $filename |grep -iv [a-zA-Z0-9]namelist |grep -iv namelist[a-zA-Z0-9] |grep -v '^ *!'|cut -d"/" -f2`; do
    gotitalready=`echo "$namelistnames" | grep -i $namelistname`
    if [[ -z "${gotitalready// }" ]]; then
      echo " "
      echo "starting work on namelistname = $namelistname"
      if [[ "$firsttime" == "yes" ]]; then
        echo "<H2> Code in directory: $directory </H2>" >> $docdir/namelists_html
        firsttime="no"
      fi
      namelistname2=`echo $namelistname | tr '[:lower:]' '[:upper:]'`
      filename2=`basename $filename`
      echo "<H3>NAMELIST: $namelistname2 (in file $filename2)</H3>" >> $docdir/namelists_html

      namelistnames="$namelistnames $namelistname"
      namelistvars=`grep -hi "namelist */ *${namelistname} */" $filename |grep -iv [a-zA-Z0-9]namelist |grep -iv namelist[a-zA-Z0-9] |grep -v '^ *!'|cut -d"/" -f3`
      if [[ ! -z `echo $namelistvars | rev |grep ' *\&'` ]]; then
        echo "!---------------------------------------------------------------------!"
        echo "!!!line continuation in namelist definition, trying to deal with it !!!"
        echo "!---------------------------------------------------------------------!"
        echo "namelistvars = $namelistvars "
        namelistvars_first=$namelistvars
        for index in `seq 1 100`; do
          namelistvars_next=`grep -A${index} -hi "${namelistvars_first}" $filename |tail -1`
          echo "namelistvars_next = $namelistvars_next "
          namelistvars="$namelistvars $namelistvars_next"
          if [[ -z `echo $namelistvars_next | rev |grep ' *\&'` ]]; then
            break
          fi
        done
      fi
      namelistvars=`echo ${namelistvars//'\n'/ }`
      namelistvars=`echo ${namelistvars//,/ }`
      namelistvars=`echo ${namelistvars//&/}`
      namelistvars=`echo $namelistvars | tr ' ' '\n' |sort -u |tr '\n' ' '`
      namelistvars="$namelistvars "
      namelistvars2=''
      for index in `seq 1 100`; do
        namelistvar=`echo "$namelistvars" | cut -d' ' -f$index`
        namelistvar=`echo ${namelistvar// /}`
        if [[ -z "${namelistvar// }" ]]; then 
          break
        fi
        namelistvars2="$namelistvars2 $namelistvar"
      done

      echo "namelistvars2= $namelistvars2"
      echo '<p style="padding:6px; color: black; background-color: white; border: black 2px solid">' >> $docdir/namelists_html
      echo '<table style="width:100%">' >> $docdir/namelists_html
      for namelistvar in $namelistvars2 ; do
        echo "greping for type definition of variable $namelistvar in $filename"
        numMatchesWithDescrip=`grep -iEw "$namelistvar" "$filename" | grep -iE "^ *(integer|real|logical|character|type).*${namelistvar}.*!" |wc -l`
	if (( "$numMatchesWithDescrip" > 0 )); then
	    echo "matches with description found"
	    greppedline=`grep -iEw "$namelistvar" "$filename" | grep -iE "^ *(integer|real|logical|character|type).*${namelistvar}.*!" | head -1 | sed -e 's/!/<b>!/'`
	    greppedline="${greppedline}</b>"
	else
	    echo "no matches with description found"
	    greppedline=`grep -iEw "$namelistvar" "$filename" | grep -iE "^ *(integer|real|logical|character|type).*${namelistvar}.*" | head -1`
	fi
        greppedline="$(echo -e "${greppedline}" | sed -e 's/^[[:space:]]*//')"
        echo "<tr>" >> $docdir/namelists_html
        echo "<td>$namelistvar</td>" >> $docdir/namelists_html
        echo "<td>$greppedline</td>" >> $docdir/namelists_html
        echo "</tr>" >> $docdir/namelists_html
      done
      echo '</table>' >> $docdir/namelists_html
      echo '</p>' >> $docdir/namelists_html

    fi
  done
done
done

cat <<EOF >> $docdir/namelists_html
</HTML>
EOF

mv $docdir/namelists_html $docdir/namelists.html
