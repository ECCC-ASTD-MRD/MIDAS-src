# OAVAR Fortran coding standards (under construction):

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/Assimilation/oavar_Coding_Standards)

# OAVAR Fortran code documentation:

* [Master branch](http://iweb.cmc.ec.gc.ca/~armabue/f90doc/master/index.html)

# OAVAR test suite

You can install a maestro suite with a serie of tests to evaluate the changes made to the code.

On the `science.gc.ca` network, you can install the suite with the commande
```bash
oavar_system_tests_science/install_suite.sh
```

Once the `xflow` appears, just launch the node `/Tests`.

You can change to the binary just just compiled by changing the variable `ENVAR_abs` in the file `oavar_system_tests_science/experiment.cfg`.  It could have the following form
```bash
ENVAR_abs=${HOME}/data_maestro/ords/oavar_abs/oavar_${ORDENV_PLAT}-${OAVAR_version}.Abs
```
