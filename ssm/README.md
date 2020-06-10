## Create your own SSM domain

First, to build the SSM packages, one must compile all the programs using:
```bash
cd src/programs
./compile_all_plat.sh
cd ../../tools/monitor
make install PGM=${PWD}/../../compiledir/midas_abs/midas.monitor_$(../../midas.version.sh).Abs
cd ../..
```

Then, SSM publication is done in two steps:
 * SSM packaging
 * publish packages in SSM domain

### SSM packaging

Then, SSM packaging is done with:
```bash
VERSION=$(./midas.version.sh | cut -c3-)
## to set variables 'MACHINE_PPP' and 'MACHINE_SUPER'
. maestro/suites/midas_system_tests/set_machine_list.dot
cd ssm
SSM_PACKAGES=${SSM_PACKAGES:-${HOME}/data_maestro/ords/SSM/midas/packages}
./package --midas-abs ${PWD}/../compiledir/midas_abs --packages ${SSM_PACKAGES}/${VERSION} --frontend ${MACHINE_PPP} --backend ${MACHINE_SUPER}
cd ..
```

You can specify `${SSM_PACKAGES}` to a directory where you want the packages to be copied. They will be used in the next step.

### Publish packages in a SSM domain

Then, SSM publish is done with:
```bash
cd ssm
DOMAIN_BASE=${DOMAIN_BASE:-${HOME}/data_maestro/ords/SSM/midas}
./publish --packages ${SSM_PACKAGES}/${VERSION} --post-install ${PWD}/post-install --workdir ${PWD}/publish-workdir-$$ --domainpath ${DOMAIN_BASE}/${VERSION}
```
Then you can use the SSM domain published with:
```bash
. ssmuse-sh -d ${DOMAIN_BASE}/${VERSION}
```
