#!/bin/ksh

touch revisionnumber.cdk
revnum=`svnversion -n`
echo "      character(len=10) :: crevision='${revnum}'" > revisionnumber.cdk
