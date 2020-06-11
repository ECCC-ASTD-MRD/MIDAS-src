# MIDAS Fortran coding standards:

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards)
* [Automatic documentation standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Documentation_Standards)

# MIDAS Fortran code documentation:

The documentation for officially supported branches is available:
* [`master` branch](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-master)
* [`v_3.3` branch](http://hpfx.science.gc.ca/~sanl000/midas-sphinx-doc/latest-v_3.3)
* [`v_3.4` branch](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.4)

# Contributing

We strongly suggest anyone considering to contribute to the MIDAS
 project, to follow the workflow documented in the [contributing
 guide](CONTRIBUTING.md).

# Getting a local copy of the code

To simply get a local copy of the code from an existing branch
associated with an issue, we suggest the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.1
clone_projet --no-central -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
or if one is interested in the latest version of the master branch
```bash
clone_projet --no-central -c master git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-master
```

## Getting code related to operational system
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.1
clone_projet --no-central -c v_3.4 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.4
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
ssh eccc-ppp3  ## or eccc-ppp4
cd ${WHERE YOUR CODE IS}
cd src/programs
yes '' | ./compile_all.sh
```
and
```bash
ssh banting    ## or daley
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

## Hosts used to run the test suite

When running `install_suite.sh`, links are created under `hub` and
`listings` just like any `maestro` suite.  If you want to control the
hosts used, you can put the list of hosts in the environment variable
```bash
MIDAS_MAKE_LINKS_MACHINE_LIST
```

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
```bash
/fs/ssm/eccc/mrd/rpn/anl/midas
```
automatically when a tag is pushed.

## Create your own SSM domain

You can create your own SSM domain using the script `ssm/domaingen`
which takes two optional arguments:
 1. `DOMAIN_BASE`: a directory where the SSM domain will be published
   * default: `${HOME}/data_maestro/ords/SSM/midas`
 2. `SSM_PACKAGES`: a directory where packages will be copied before published in the SSM domain
   * default: `${DOMAIN_BASE}/packages`

# Tools

Several tools related to MIDAS are included in the codebase.  Those
tools have a code separated from the main code in MIDAS.

## `midas_scripts`

Those are the helper scripts which launch MIDAS programs.

Refer to the [`midas_scripts/README.md`](tools/midas_scripts/README.md) for more details.

## `midas.monitor`

This program monitors a file to react to its content.

See [`monitor/README.md`](tools/monitor/README.md) for more details.

## `midas.findTrials`

This scripts finds the trial name extensions in an assimilation window.

See [`findTrials/README.md`](tools/findTrials/README.md) for more details.

# Automatic Testing using GitLab-CI

An automatic system of tests has been developed.  For each push in the
`master` branch the system tests are launched to guarantee that the
all the tests pass for the `master` branch.  The [instructions for
automatic testing using GitLab-CI are available in a separate
file](CI.md).
