#!/bin/bash
codedir=$1
docdir=${2:-${PWD}/doc}

if [ -z $codedir ]; then
  echo " "
  echo "   create_f90doc.sh: must specify the main code directory!" 
  echo " "
  exit
fi

cd $codedir

if [[ ! -d programs ]]; then
  echo ""
  echo "The selected directory does not seem to be the MIDAS 'src' directory!"
  echo "ABORTING..."
  echo ""
  exit
fi

codedir=$PWD
cd ../

cd scripts
scriptdir=$PWD
cd ../

[ ! -d "${docdir}" ] && mkdir ${docdir}
cd ${doc}

echo " "
echo "Generating documentation for source code: $codedir"
echo "Documentation will be located in:         $docdir"
echo "Using scripts in:                         $scriptdir"
echo " "

read -p "Press ENTER to continue"

f90doc=$scriptdir/f90doc/f90doc-0.4.0/f90doc

cd $docdir
rm -f *.html */*.html

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
      namelistvars="$namelistvars "
      namelistvars2=''
      for index in `seq 1 100`; do
        namelistvar=`echo "$namelistvars" | cut -d' ' -f$index`
        namelistvar=`echo ${namelistvar// /}`
        #echo "starting work on namelistvar=$namelistvar"
        if [[ -z "${namelistvar// }" ]]; then 
          break
        fi
        gotitalready=`echo "$namelistvars2" | grep -i $namelistvar`
        if [[ -z "${gotitalready// }" ]]; then
          #echo "namelistvar = $namelistvar"
          namelistvars2="$namelistvars2 $namelistvar"
        fi
      done

      echo "namelistvars2= $namelistvars2"
      echo '<p style="padding:6px; color: black; background-color: white; border: black 2px solid">' >> $docdir/namelists_html
      echo '<table style="width:100%">' >> $docdir/namelists_html
      for namelistvar in $namelistvars2 ; do
        echo "greping for type definition of variable $namelistvar in $filename"
        greppedline=`grep -wi $namelistvar $filename | grep -i -E -w 'integer|real|logical|character' |grep -v '^ *!' |head -1`
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

cd $docdir

for i in $codedir/programs/*.f*90 ; do 
  if [ -f $i ]; then
    echo "Processing file $i"
    $f90doc $i
    echo "return code = $?"
  fi
  echo " "
done

for i in $codedir/modules/*.f*90 ; do 
  echo "Processing file $i"
  $f90doc $i
  echo "return code = $?"
  echo " "
done

for i in $codedir/modules/* ; do 
  if [ -d $i ]; then
    ii=`basename $i`
    echo "Entering subdirectory $ii"
    mkdir -p $ii
    found_files=0
    for i2 in $codedir/modules/$ii/*.f*90 ; do
      if [ -f $i2 ]; then
        echo "Processing file $i2"
        found_files=1
        $f90doc -o temp_ $i2
        echo "return code = $?"
        ib=`echo $i2 | cut -f1 -d.`
        mv temp_* $ii/
        cd $ii
        for iii in temp_* ; do
          iiii=`echo $iii|cut -c6-`
          mv $iii $iiii
        done
        cd ..        
        echo " " 
      fi
    done
    if [ $found_files == 0 ]; then
      rmdir $ii
    fi
  fi
done

cd $codedir
revision=$(git describe --always --dirty=_M 2>/dev/null)
cd $docdir

cat <<EOF > index_html
  <HTML style="font-family:arial">
  <HEAD>
     <TITLE> MIDAS fortran software </TITLE>
  </HEAD>
  <H1> Documentation of MIDAS fortran software </H1>
  <b>Documentation created on:</b> `TZ='America/New_York' date`<br><br>
  <b>Code Version:</b> $revision<br><br>
  <b><A style="text-decoration:none" HREF="namelists.html">Page with NAMELIST documentation</A><br><br>
  <H2> Programs directory </H2>
EOF

numcol=3
echo '<p style="padding:6px; color: black; background-color: white; border: black 2px solid">' >> index_html
echo '<table style="width:100%">' >> index_html
FILES=`grep -li program *.html`
count=0
echo "<tr>" >> index_html
for i in $FILES ; do 
  if [ -f $i ]; then
    count=$((count+1))
    ib=`basename $i .html`
    echo "<td><A style=\"text-decoration:none\" HREF=""${i}"">${ib}</A></td>" >> index_html
    if [ $((count % numcol)) == 0 ]; then 
      echo "</tr><tr>" >> index_html
    fi
  fi
done
echo "</tr>" >> index_html
echo '</table>' >> index_html
echo '</p>' >> index_html

echo "<H2> Modules directory </H2>" >> index_html
echo '<p style="padding:6px; color: black; background-color: white; border: black 2px solid">' >> index_html
echo '<table style="width:100%">' >> index_html
FILES=`grep -Li program *.html`
count=0
echo "<tr>" >> index_html
for i in $FILES ; do 
  if [ -f $i ]; then
    count=$((count+1))
    ib=`basename $i .html`
    echo "<td><A style=\"text-decoration:none\" HREF=""${i}"">${ib}</A></td>" >> index_html
    if [ $((count % numcol)) == 0 ]; then 
      echo "</tr><tr>" >> index_html
    fi
  fi
done
echo "</tr>" >> index_html
echo '</table>' >> index_html
echo '</p>' >> index_html

for i in * ; do 
  if [ -d $i ]; then
    dd=`basename $i`    
    echo "<H2> Modules subdirectory: $dd </H2>" >> index_html
    echo "entering directory: $dd"
    cd $dd
    echo '<p style="padding:6px; color: black; background-color: white; border: black 2px solid">' >> ../index_html
    FILES=`grep -li program *.html`
    for i2 in $FILES ; do
      if [ -f $i2 ]; then
        ib=`basename $i2 .html`
        echo "<A style=\"text-decoration:none\" HREF=""${i2}"">${ib}</A> (Program)<br>" >> ../index_html
      fi
    done
    echo '<table style="width:100%">' >> ../index_html
    FILES=`grep -Li program *.html`
    count=0
    echo "<tr>" >> ../index_html
    for i2 in $FILES ; do
      if [ -f $i2 ]; then
        count=$((count+1))
        ib=`basename $i2 .html`
        echo "<td><A style=\"text-decoration:none\" HREF=""${i2}"">${ib}</A></td>" >> ../index_html
        if [ $((count % numcol)) == 0 ]; then 
          echo "</tr><tr>" >> ../index_html
        fi
      fi
    done
    echo '</table>' >> ../index_html
    echo '</p>' >> ../index_html
    cd ../
  fi
done

cat <<EOF >> index_html
</HTML>
EOF

mv $docdir/index_html $docdir/index.html
mv $docdir/namelists_html $docdir/namelists.html

# Move all html files back into main directory (no subdirectories)
for i in * ; do 
  if [ -d $i ]; then
    dd=`basename $i`    
    cd $dd
    for i2 in *.html ; do
      if [ -f $i2 ]; then
        mv $i2 ../
      fi
    done
    cd ../
    rmdir $i
  fi
done

echo "The documentation has been generated in ${docdir}."
