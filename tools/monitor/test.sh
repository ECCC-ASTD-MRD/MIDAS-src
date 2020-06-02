#!/bin/bash

## Ce script sert a tester le programme 'midas.monitor.Abs'.

set -e

MONITOR_PGM=${1:-./midas.monitor.Abs}

echo "Testing program '${MONITOR_PGM}'"
cat > run_after_rebm_done.sh <<EOF
#!/bin/ksh

echo "On a trouve le fichier 'VAR3D_STATUS.dot' qui contenait 'VAR3D_STATUS=REBM_DONE'"
echo ORDENV_PLAT=\${ORDENV_PLAT}
EOF
chmod +x run_after_rebm_done.sh

[ -f VAR3D_STATUS.dot ] && rm VAR3D_STATUS.dot

function launch_monitor {
    set -e

    ${MONITOR_PGM} VAR3D_STATUS.dot ./run_after_rebm_done.sh
    sleep 10

    cat > VAR3D_STATUS.dot <<EOF
VAR3D_STATUS=VAR3D_BEG
EOF
    sleep 10

    cat > VAR3D_STATUS.dot <<EOF
EOF
    sleep 10

    cat > VAR3D_STATUS.dot <<EOF
VAR3D_STAT
EOF
    sleep 10

    cat > VAR3D_STATUS.dot <<EOF
VAR3D_STATUS=REBM_DON
EOF
    sleep 10

    cat > VAR3D_STATUS.dot <<EOF
VAR3D_STATUS=REBM_DONE
EOF

} ## End of function 'launch_monitor'

launch_monitor 2>&1 > monitor_stdout.txt | tee monitor_stderr.txt

diff_status=0
cat <<EOF | diff - monitor_stderr.txt || diff_status=1
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
EOF

cat <<EOF | diff - monitor_stdout.txt || diff_status=1
On a trouve le fichier 'VAR3D_STATUS.dot' qui contenait 'VAR3D_STATUS=REBM_DONE'
ORDENV_PLAT=${ORDENV_PLAT}
file=VAR3D_STATUS.dot,cmd=run_after_rebm_done.sh
Executing:run_after_rebm_done.sh VAR3D_STATUS.dot
EOF

if [ "${diff_status}" -eq 0 ]; then
    echo "The 'midas.monitor.Abs' tests are successull!!!"
else
    echo "The 'midas.monitor.Abs' tests did not pass!!!"
    exit 1
fi
