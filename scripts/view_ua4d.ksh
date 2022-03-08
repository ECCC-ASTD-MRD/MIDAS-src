#!/bin/bash
# this script modifies the structure of headers for UA 4D files to treat them with a common sqlr_readSqlite fortran routine 
# Attention ! Make sure that the sql table resume containing date, time and run exists in the output files
sqlfile=$1
FAM='ua'

KEY=" concat(tsonde,hlat,hlon,ifnull(vcoord,-999.99),id_obs) "
ELEMS="12004,12203,10051,10004,11011,11012,12001,11001,11002,12192,10194"

cd $TMPDIR
rdbgen -f uadsplit.rdb -type ua
cat  <<EOF > command
attach 'uadsplit.rdb' as db1;

create temporary table  data2 as select  $KEY  idobs ,id_data,vcoord,varno,vcoord_type,obsvalue,flag,omp,oma,OBS_ERROR,FG_ERROR,VERT_SOUNDING_SIGN,fso from data natural join  header    
where  varno in( $ELEMS ) ;
create temporary table ids as select distinct idobs from data2;


create temporary table  Header2 as select distinct   $KEY   idobs,id_stn,strftime('%Y%m%d', tsonde) as date,strftime('%H%M%S', tsonde) as time,status,id_stn,hlat lat,hlon lon,codtyp,elev  from header natural join data where $KEY  in (select idobs from ids) and id_data in (select id_data from data2)   ;

create temporary table  Header3 as select rowid idobsint,idobs, id_stn,lat,lon,date,time,status,codtyp,elev from  header2;
create temporary table Data3 as select idobsint,id_data, varno,vcoord,vcoord_type,obsvalue,flag from  Header3 natural join Data2;

insert into db1.header  (id_obs,id_stn,lat,lon,date,time,status,codtyp,elev)   select idobsint, id_stn,lat,lon,date,time,status,codtyp,elev from header3;
insert into db1.data (id_obs,id_data, varno,vcoord,vcoord_type,obsvalue,flag ) select idobsint,id_data, varno,vcoord,vcoord_type,obsvalue,flag  from  data3;
create table db1.rdb4_schema( schema  varchar(9) );
insert into  db1.rdb4_schema values( '${FAM}' );
EOF
d.sqlite ${sqlfile} < command
mv uadsplit.rdb ${sqlfile}_ua4d
