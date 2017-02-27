#!/bin/ksh

## Ce script sert a tester le programme 'oavar_monitor_rebm'.

set -ex

cat > run_after_rebm_done.sh <<EOF
#!/bin/ksh

echo "On a trouve le fichier 'VAR3D_STATUS.dot' qui contenait 'VAR3D_STATUS=REBM_DONE'"
which editfst
echo ORDENV_PLAT=\${ORDENV_PLAT}
EOF
chmod +x run_after_rebm_done.sh

[ -f VAR3D_STATUS.dot ] && rm VAR3D_STATUS.dot
./oavar_monitor_rebm VAR3D_STATUS.dot run_after_rebm_done.sh

sleep 10

cat > VAR3D_STATUS.dot <<EOF
VAR3D_STATUS=VAR3D_BEG
EOF

sleep 10

cat > VAR3D_STATUS.dot <<EOF
VAR3D_STATUS=REBM_DONE
EOF
