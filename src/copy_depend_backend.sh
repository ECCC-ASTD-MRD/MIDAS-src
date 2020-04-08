#! /bin/sh
##  sourcing user configuration
source ./config.dot.sh

##  sourcing compilation configuration and SSM packages
source ${DOT_CONFIG}

##  sourcing utilitary functions
source ./func.dot.sh

copy_depend ${BACKEND}
