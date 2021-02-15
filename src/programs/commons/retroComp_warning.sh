#! /bin/sh

typeset __do_retro=false

## **Temporary** retrocompatibility
if [ ! -z ${COMPILE_MIDAS_ADD_DEBUG_OPTIONS+x} ]
then
    ## check if new variable also present and then if equal
    if [ ! -z ${MIDAS_COMPILE_ADD_DEBUG_OPTIONS+x} ]; then
        echo "Both MIDAS_COMPILE_ADD_DEBUG_OPTIONS and \
COMPILE_MIDAS_ADD_DEBUG_OPTIONS defined"
        if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS}" != "${COMPILE_MIDAS_ADD_DEBUG_OPTIONS}" ]; then
            echo "... but have conflicting values! exiting."
            exit 1
        fi
    else
        __do_retro=true
    fi
fi


if [ ! -z ${COMPILEDIR_MIDAS_MAIN+x} ]
then
    ## check if new variable also present and then if equal
    if [ ! -z ${MIDAS_COMPILE_DIR_MAIN+x} ]; then
        echo "Both MIDAS_COMPILE_DIR_MAIN and COMPILEDIR_MIDAS_MAIN defined"
        if [ "${MIDAS_COMPILE_DIR_MAIN}" != "${COMPILEDIR_MIDAS_MAIN}" ]; then
            echo "... but have conflicting values! exiting."
            exit 1
        fi
    else
        __do_retro=true
    fi
fi


if [ ! -z ${COMPILEDIR_MIDAS_ID+x} ]
then
    ## check if new variable also present and then if equal
    if [ ! -z ${MIDAS_COMPILE_DIR_ID+x} ]; then
        echo "Both MIDAS_COMPILE_DIR_ID and COMPILEDIR_MIDAS_ID defined"
        if [ "${MIDAS_COMPILE_DIR_ID}" != "${COMPILEDIR_MIDAS_ID}" ]; then
            echo "... but have conflicting values! exiting."
            exit 1
        fi
    else
        __do_retro=true
    fi
fi


if ${__do_retro}
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
║  If both are present, they must be equals to prevent this message from       ║
║  showing.  New variables have precedence.                                    ║
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
