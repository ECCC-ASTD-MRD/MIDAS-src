# Advanced Unit Testing Topics

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

### Debugging tools

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

### Preparing interactive mode without changing any configuration (advanced)

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


## Activating the code coverage report

You can activate the code coverage reporting by setting
```bash
export MIDAS_COMPILE_ADD_CODECOVERAGE_OPTIONS=yes
``̀
before compiling the program.


## Updating Test Results

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


## Splitting an observation file to prepare a unit test input

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

### Gaining access to `midas.splitobs.Abs`

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
