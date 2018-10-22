#!/bin/bash

filein=$1

# obtain list of tables in the file
echo ".tables" > gettables.sh
tablelist=`sqlite3 $filein < gettables.sh`
rm gettables.sh

# loop over all tables
for table in $tablelist; do
  echo "select * from $table ;" > dumptable.sh
  sqlite3 -header -column $filein < dumptable.sh > ${filein}_${table}.dat
done
rm dumptable.sh
