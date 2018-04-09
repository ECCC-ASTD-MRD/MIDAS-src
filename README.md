# MIDAS Fortran coding standards (under construction):

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards)

# MIDAS Fortran code documentation:

* [Master branch](http://hpfx.science.gc.ca/~erv000/midas-doc/latest)

# Compiling the programs

To compile the programs used in this code, use the commands
```bash
ssh eccc-ppp2  ## or eccc-ppp1
cd ${WHERE YOUR CODE IS}
cd src/programs
yes '' | ./compile_all.sh
```
and
```bash
ssh brooks  ## or hare
cd ${WHERE YOUR CODE IS}
cd src/programs
yes '' | ./compile_all.sh
```

# MIDAS test suite

You can install a maestro suite with a serie of tests to evaluate the
changes made to the code.

On the `science.gc.ca` network, you can install the suite with the command
```bash
maestro/suites/midas_system_tests/install_suite.sh
```

Once the `xflow` appears, just launch the node `/Tests`.

The suite is configured to use by default the programs you just
compiled.

# Updating the results

*This script is not working for now!*

When the changes introduced modify the result, then one must update the
results to that the tests pass.  After checking carefully the listing
and the results obtained, you can update them with the command:
```bash
maestro/suites/midas_system_tests/update_results.sh
```

You can give as argument a single test but by default it collects the
results for all the tests.  It also saves the listing.

This will save the new references results locally on your account and
update the system test suite to point to these results.

