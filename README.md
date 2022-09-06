# Breaking change!

If you compile MIDAS code after tag `v_3.7.2` on HPCR-U2, then you
must update your login profile to version `1.19.0`:
```bash
ln -svi /fs/ssm/eccc/mrd/ordenv/profile/1.19.0 ~/.profile_1.19.0
rm -v ~/.profile && ln -svi .profile_1.19.0 ~/.profile
```

This change is backward compatible for your suites but you absolutely
need to update your profile to compile any MIDAS code after version
`v_3.7.2`.

To know if your code is after `v_3.7.2`, you can execute:
```bash
git describe
```

# MIDAS Fortran coding standards:

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards)
* [Automatic documentation standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Documentation_Standards)

# MIDAS releases general and code documentations:

The documentation for officially supported branches is available:
* `main` branch
  * [General documentation (`README.md`) - this page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/main/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-main)
* `v_3.7` branch
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.7/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/v_3.7.2)

# Contributing

We strongly suggest anyone considering to contribute to the MIDAS
 project, to follow the workflow documented in the [contributing
 guide](CONTRIBUTING.md).

# Getting a local copy of the code

To simply get a local copy of the code from an existing branch
associated with an issue, we suggest the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
clone_projet --no-central -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
or if one is interested in the latest version of the `main` branch
```bash
clone_projet --no-central -c main git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-main
```

## Getting code related to IC-3 system on HPCR-U2

```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
clone_projet --no-central -c v_3.7 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.7
```

If you created a new branch with the GitLab web UI, then the branch
has been created using the default branch which is `main`.  One must
reset it to the release branch.  One can simply do:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
clone_projet --no-central -c v_3.7 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
cd midas-${ISSUE_NUMBER}
git checkout -b ${ISSUE_NUMBER}-complete-the-name-of-the-branch-as-on-GitLab
git push origin ${ISSUE_NUMBER}-complete-the-name-of-the-branch-as-on-GitLab --force
```

# Compiling MIDAS

[`src/midas_build`](./src/README.md)
is now the official compilation tool to build MIDAS.
To proceed to compilation tasks, you should be in the `src/` directory.
`midas_build` compiles by default on both platforms.

## Compiling a single program
To compile a single program on both platforms, do the following from the 
frontnode:
```bash
cd ${where_your_code_is}
cd src
./midas_build ${program_basename}.Abs
```
where `program_basename` is the **basename** of one the files with extention `.f90` in
the sub-directory `src/programs` or [`splitobs`](./src/README.md#splitobs-an-external-program).
If you installed the [auto-completion feature](./src/README.md#auto-completion)
you can browse all install targets by pressing `<TAB>` following `./midas_build`.

By default the binary will be installed in 
`${HOME}/data_maestro/ords/midas-bld/midas_abs/` 
(this can be [configured by environment variables](./src/README.md#configuring-the-compilation-and-linking-process)).

## Compiling all programs
To compile all programs (`src/programs/*.f90` as well as
[`splitobs`](./src/README.md#splitobs-an-external-program)),
simply do:
```bash
cd ${where_your_code_is}
cd src
./midas_build
```

## Complete documentation on using `midas_build` and `make`
If you are [contributing a new program, changing external dependencies](./src/README.md#adding-a-new-program-or-changing-external-dependencies),
recompiling a lot or debugging the code,
you should take the time to read the detailed instructions found in
[`src/README.md`](./src/README.md).

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

## Interactive debugging

You can debug interactively a MIDAS program by launching a job and
login into it.  For that, you can set
```bash
UnitTest_stop_for_interactive_work=yes
```
in `maestro/suites/midas_system_tests/experiment.cfg`.

to prepare the working directory for interactive debugging.  The first
one is the suggested method.

After setting this, launch the `/${pathToTest}/UnitTest/run` task
itself and look at the log messages.  You will see something like
this:
```
We stop here to let you continue to work interactively.
You can use:
    ${SEQ_EXP_HOME}/modules/UnitTest/scripts/launch_interactive.sh -exp ${SEQ_EXP_HOME} -node /${pathToTest}/UnitTest/run -date ${SEQ_DATE}
to launch an interactive job.
```

For example, if `pathToTest=Tests/thinning/IASI` and you are working
on issue 493, then you will have something:
```
We stop here to let you continue to work interactively.
You can use:
    ${HOME}/.suites/midas-493/modules/UnitTest/scripts/launch_interactive.sh -exp ${HOME}/.suites/midas-493 -node /Tests/thinning/IASI/UnitTest/run -date 20211117000000
to launch an interactive job.
```

The script `launch_interactive.sh` is launching an interactive job and
it will give you the instructions on how to login to the job and start
working.

#### Debugging tools

This setup is building an interactive environment to launch again and
again the program you are debugging.  This is the script
`launch_program.sh`.  You launch the program with:
```bash
./launch_program.sh ${path_to_program}
```

Sometimes, it is very helpful to use some debugging tools.  So the
script `launch_program.sh` is supporting two tools for debugging:
 * [`gdb`](https://www.gnu.org/software/gdb) and
 * [`DDT`](https://portal.science.gc.ca/confluence/display/SCIDOCS/DDT)

which can be activated by using respectively the options `--gdb` and
`--ddt` when calling `launch_program.sh` like this:
```bash
./launch_program.sh ${path_to_program} --ddt
```
For DDT, the GUI is automatically started for you.

When you use those tools, it is suggested to compile with debugging
options enabled by using:
```bash
export MIDAS_COMPILE_ADD_DEBUG_OPTIONS=yes
```

#### Preparing interactive mode without changing any configuration (advanced)

You can prepare the interactive mode by using this command:
```bash
maestro -n /${pathToTest}/UnitTest/run -d ${SEQ_DATE} -s submit -o -args 'UnitTest_stop_for_interactive_work=yes'
```
instead of modifying `experiment.cfg` as mentioned in the previous
section.

You could also modify the file
```
maestro/suites/midas_system_tests/abs.dot
```
since it is ignored by Git.

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

# Using `midas.splitobs.Abs`

To gain access to the program `midas.splitobs.Abs`, one can load one
of the MIDAS SSM domain with:
```bash
. ssmuse-sh -d eccc/mrd/rpn/anl/midas/3.6.8
```
or point directly to the program in the SSM domain:
```
/fs/ssm/eccc/mrd/rpn/anl/midas/3.6.8/${ORDENV_PLAT}/bin/midas.splitobs.Abs
```

An English help message is printed to screen by calling
`midas.splitobs.Abs` one these ways:
```bash
midas.splitobs.Abs -h
midas.splitobs.Abs -help
midas.splitobs.Abs --help
midas.splitobs.Abs ## without any argument
```
and a French documentation is printed when the program is called with:
```bash
midas.splitobs.Abs -aide
midas.splitobs.Abs --aide
```

Note that this program processes SQLite and BURP files transparently.

### Splitting a file to prepare a unit test input

In the case, a user wants to split an observation file to prepare a
unit test input, this program can be called like this:
```bash
midas.splitobs.Abs -obsin ${OBS_FILE_IN} -obsout ${OBS_FILE_OUTPUT_PREFIX} -round-robin -npex ${NPEX} -npey ${NPEY}
```
where `${OBS_FILE_IN}` is the input observation file,
`${OBS_FILE_OUTPUT_PREFIX}` is the prefix of the output files,
`${NPEX}` is the number of parts in the X direction and `${NPEY}` is
the number of parts in the Y direction.  Here the `${NPEX}` and
`${NPEY}` must be the same as the MPI topology of the `run` task of
the test.

For example, if you have `cpu="2x40x4"` as the CPUs setting in the
`UnitTest/run.xml` resource file of the test, then `NPEX=2` and
`NPEY=40` and the program will generate $`2\times 40=80`$ files.

For example, if you call the command
```bash
midas.splitobs.Abs -obsin obs_input -obsout obs_split -round-robin -npex 2 -npey 3
```
it will split the input observation file `obs_input` and generates 6 files with the name:
```
obs_split_0001_0001
obs_split_0001_0002
obs_split_0001_0003
obs_split_0002_0001
obs_split_0002_0002
obs_split_0002_0003
```

# Tools

Several tools related to MIDAS are included in the codebase.  Those
tools have a code separated from the main code in MIDAS.

## `midas_scripts`

Those are the helper scripts which launch MIDAS programs.

Refer to the [`midas_scripts/README.md`](tools/midas_scripts/README.md) for more details.

## `midas.splitobs`

This program is used to split the observations into several files
according to one of the following strategy:
 * round-robin
 * lat-lon tiles of a grid

It can also select the observations that lies in a domain defined by a
RPN grid.

Refer to the [`splitobs/README.md`](tools/splitobs/README.md) for more details.

## `midas.monitor`

This program monitors a file to react to its content.

See [`monitor/README.md`](tools/monitor/README.md) for more details.

## `midas.findTrials`

This scripts finds the trial name extensions in an assimilation window.

See [`findTrials/README.md`](tools/findTrials/README.md) for more details.

