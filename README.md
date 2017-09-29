# OAVAR Fortran coding standards (under construction):

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/Assimilation/oavar_Coding_Standards)

# OAVAR Fortran code documentation:

* [Master branch](http://hpfx.science.gc.ca/~mab001/f90doc/master/)

# OAVAR test suite

You can install a maestro suite with a serie of tests to evaluate the
changes made to the code.

On the `science.gc.ca` network, you can install the suite with the commande
```bash
maestro/suites/oavar_system_tests/install_suite.sh
```

Once the `xflow` appears, just launch the node `/Tests`.

You can change to the binary just just compiled by changing the
variable `ENVAR_abs` in the file
`maestro/suites/oavar_system_tests/experiment.cfg`.  It could have the
following form
```bash
ENVAR_abs=${HOME}/data_maestro/ords/oavar_abs/oavar_${ORDENV_PLAT}-${OAVAR_version}.Abs
ENVAR_ominusf_abs=${ABS_DIR}/ominusf_${ORDENV_PLAT}-${OAVAR_version}.Abs
```
if you compiled the programs with the commands
```bash
ssh eccc-ppp2
cd ${WHERE YOU CODE IS}
cd src
./compile_oavar.sh
cd ominusf
./compile_ominusf.sh
```
and
```bash
ssh brooks
cd ${WHERE YOU CODE IS}
cd src
./compile_oavar.sh
```

# Updating the results

When the changes introduced modify the result, then one must update the
results to that the tests pass.  After checking carefully the listing
and the results obtained, you can update them with the command:
```bash
maestro/suites/oavar_system_tests/update_results.sh
```

You can give as argument a single test but by default it collects the
results for all the tests.  It also saves the listing.

This will save the new references results locally on your account and
update the system test suite to point to these results.

