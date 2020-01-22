#! /bin/sh

###########################################################
##
##  USER CONFIGURATION
##
###########################################################

BACKEND=daley
FRONTEND=eccc-ppp4
JOBNAME=midasCompilation
DIR_BUILD=../compiledir
NCORES=8
VERBOSE=2
CLEAN=false
DIRECT_FRONTEND_COMPILE=false

###########################################################
##  compilation and SSM needed for compilation
##
## -- should not change that
DOT_CONFIG=./programs/commons/compile_setup.sh
