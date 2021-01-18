#! /bin/sh


## **Temporary** retrocompatibility
if [ ! -z ${COMPILE_MIDAS_ADD_DEBUG_OPTIONS+x} ] || \
   [ ! -z ${COMPILEDIR_MIDAS_MAIN+x} ] || \
   [ ! -z ${COMPILEDIR_MIDAS_ID+x} ] ;
then

    cat << EndOfMessage

╔══╡ WARNING ╞═════════════════════════════════════════════════════════════════╗
║                                                                              ║
║  Please modify your profile to account for this change:                      ║
║                                                                              ║
║  All  COMPILE*_MIDAS_*  variable have been renamed  MIDAS_COMPILE_*          ║
║                                                                              ║
║  *  MIDAS_COMPILE_DIR_MAIN  <-  COMPILEDIR_MIDAS_MAIN                        ║
║  *  MIDAS_COMPILE_ADD_DEBUG_OPTIONS  <-  COMPILE_MIDAS_ADD_DEBUG_OPTIONS     ║
║  *  MIDAS_COMPILE_DIR_ID  <-  COMPILEDIR_MIDAS_ID                            ║
║                                                                              ║
║                                                                              ║
║  Retro-compatibility is maintained **until** next release!                   ║
║                                                                              ║
║  Thank you for your understanding.                                           ║
║           - your friendly CMDA developper                                    ║
║                                                                              ║
║  PRESS ANY KEY TO CONTINUE                                                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

EndOfMessage
    read
    set -x
    source ${__toplevel}/src/programs/commons/retroComp_export.sh
    set +x
fi
