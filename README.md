# MIDAS Fortran coding standards:

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards)

# MIDAS Fortran code documentation:

The documentation for officially supported branches is available:
* [`master` branch](http://hpfx.science.gc.ca/~sanl000/midas-doc/latest)
* [`master` branch (new sphinx prototype)](http://hpfx.science.gc.ca/~sanl000/midas-sphinx-doc/latest-master)
* [`v_3.3` branch (new sphinx prototype)](http://hpfx.science.gc.ca/~sanl000/midas-sphinx-doc/latest-v_3.3)

# Contributing

We strongly suggest anyone considering to contribute to the MIDAS
 project, to follow the workflow documented in the [contributing
 guide](CONTRIBUTING.md).

# Getting a local copy of the code

To simply get a local copy of the code from an existing branch
associated with an issue, we suggest the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/1.2
clone_projet --no-central -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
or if one is interested in the latest version of the master branch
```bash
clone_projet --no-central -c master git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-master
```

## Getting code related to operational system
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/1.2
clone_projet --no-central -c v_3.2 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.2
```

## Getting code related to the latest final cycle version or installed in parallel run
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/1.2
clone_projet --no-central -c v_3.3 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.3
```

# Compiling a single program

To compile a program for a given platform, one has to do:
```bash
ssh ${host}  ## '${host}' can be 'hare', 'brooks', 'eccc-ppp1', eccc-ppp2, 'gpsc*'
cd ${WHERE YOUR CODE IS}
cd src/programs
./compile_program.sh ${program}
```
where `program` may be one the file with extention `.f90` in the
sub-directory `src/programs`.

The listing of the compilation will let you know where the program
binary is.

Note, if the compilation fails due to a change in the dependencies between
the main program and modules or among the modules, the `src_files` files
can be automatically updated by running the script `make_src_files.sh` that is 
located in the same directory as the main programs.

# Compiling all the programs

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

## Compiling all programs on both platforms

A script, `compile_all_plat.sh`, has been written to compile all
programs of this project on supported platforms:
`ubuntu-14.04-amd64-64` and `sles-11-broadwell-64-xc40`.  The
compiling is done in parallel.  You can call it with:
```bash
cd src/programs
./compile_all_plat.sh
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

## Updating the results

The results can be updated by running the task `UnitTest/update` for
the wanted test.  You have to specify a path to store the new results
with the variable `UnitTest_reference_update` in the test
configuration file.

The listings will be collected at the same time.

Once the new results have been collected, you can update the variable
`UnitTest_reference` with the path to the new results.

You can then commit the changes to the configuration.

When you will ask a merge-request, the new results will be copied in a
safe directory with all other reference results by one of the
maintainer.

# SSM

The [CI](CI.md) has been configured to produce a SSM domain under
```
/fs/ssm/eccc/mrd/rpn/anl/midas
```
automatically when a tag is pushed.

You can produce your own SSM domain using these commands:
```bash
cd ssm
export DOMAIN_PATH=${PATH TO THE SSM DOMAIN THAT WILL BE PRODUCED}
export MIDAS_REVISION=${name of a branch or a tag}
./publish
```
Then, you can use this domain with:
```bash
. ssmuse-sh -d ${DOMAIN_PATH}
```

# Automatic Testing using GitLab-CI

An automatic system of tests has been developed.  For each push in the
`master` branch the system tests are launched to guarantee that the
all the tests pass for the `master` branch.  The [instructions for
automatic testing using GitLab-CI are available in a separate
file](CI.md).
