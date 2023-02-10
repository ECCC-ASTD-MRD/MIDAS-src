# Advanced Unit Testing Topics

## Test Driven Development in MIDAS

Test-driven development (TDD) is a development process relying on the use of test
cases being developed at the start of a development cycle (or at least before
its end).

This has many advantages.  First it allows to clearly define a priori the
objectives of a given development, but maybe more importantly it implies that we
develop and maintain a suite of tests that reduce the risk of breaking something
while adding a new feature or modifying an existing feature.

It also allows for verification after each incremental modification or refactoring
step allowing for easier correction of errors before they become entangled in
important modifications in many parts of the code.

The [MIDAS test suite](maestro/suites/midas_system_tests/config/Tests) is based on
this development paradigm.

### MIDAS Test Suite

The MIDAS test suite is based on maestro module environment.  It instantiates the
[`UnitTest` module](maestro/suites/midas_system_tests/modules/UnitTest) for each
test case through a sequence of 4 tasks:
  1. `get`: fetch the inputs from the reference path provided in the test
     configuration
  2. `run`: submit the program to be tested with the provided configuration
  3. `check`: fetch the expected results from the reference path and compare
     them with the program outputs (aborting if there is a difference)
  4. `clean`: remove the work directories 
If a task aborts, the following one won't launch allowing the user to inspect the
test state.

## Building a Test

### Appending to the Flow

Edit [`maestro/suites/midas_system_tests/modules/Tests/flow.xml`](maestro/suites/midas_system_tests/modules/Tests/flow.xml)
and add a `SUBMIT` for the new test
```xml
  <SUBMITS sub_name="${yourNewTest}"/>
```
Then add below the `FAMILY` description.
If there is only a single test, the description will contain only a single `UnitTest`
instance `SUBMITS` element
```xml
  <FAMILY  name="${yourNewTest}">
    <SUBMITS sub_name="UnitTest"/>
    <MODULE  name="UnitTest"/>
  </FAMILY>
```
If you are creating more than one test (or adding a new one), then there will be
many sub-`FAMILY` element:
```xml
  <FAMILY  name="${yourNewTest}">
    <FAMILY name="testA">
      <SUBMITS sub_name="UnitTest"/>
      <MODULE  name="UnitTest"/>
    </FAMILY>
    <FAMILY name="testB">
      <SUBMITS sub_name="UnitTest"/>
      <MODULE  name="UnitTest"/>
    </FAMILY>
  </FAMILY>
```

### Test Case Configuration

In `maestro/suites/midas_system_tests/config/Tests/`, create a directory for your
test configuration.
If it is a single test, create the test configuration at that same level, for instance
```
yourNewTest/
yourNewTest.cfg
```
Then edit this configuration file and provide the basic `UnitTest` variable definitions
(all the variables and default values can be consulted in the [`UnitTest` module definition](maestro/suites/midas_system_tests/modules/UnitTest/container.cfg))
```sh
UnitTest_run_namelist=${SEQ_EXP_HOME}/config/Tests/${yourNewTest}/nml
UnitTest_run_exe=${ABS_DIR}/midas-${yourTestProgram}_${ORDENV_PLAT}-${MIDAS_version}.Abs

UnitTest_mpiscript=var.sh 
# most probably unless a test specific launch script is provided

UnitTest_reference=${pathToReference}
#    * by convention the leaf directory should be version number represented as
#      four digits padded with 0 (such as `0001`)
#    * when your TTD contribution will be merged back, this directory will be moved
#      to a protected account (such as `~sanl000`).
UnitTest_reference_update=${pathForUpdate}
# the path where the _updated_ inputs and expected 
# results would be placed if the task `update` was launched.

UnitTest_maximum_execution_time=${expectedTime}
```

### Preparing Inputs

The `get` task will fetch the contents of all `inputs*.ca` (`cmcarc` archives) 
from the directory pointed by `${UnitTest_reference}`.
When creating a test (or modifying its inputs), assemble all the required 
inputs and group them in different archives, for instance separating observations
from trials or constants and so forth.
To archive a group of files, the `cmcarc` program:
```sh
cmcarc -f ${UnitTest_reference}/inputs_group.ca -a file1 file2 ...
```
Note that some files (such as ensemble members and observation files) are expected
to be within subdirectories in the working directory and therefore must also be
stored in the same subdirectories within the archive. For example:
```sh
cmcarc -f ${UnitTest_reference}/inputs_ensemble.ca -a ensemble/${Date}_*
```

### Allocating Resources

Maestro resources must also be provided for each task in the test.
Theses files are in `maestro/suites/midas_system_tests/resources/Tests/${yourNewTest}`
and you should create (or copy from another test and adapt them):
```
${yourNewTest}/
├── container.xml
└── UnitTest
    ├── check.xml
    ├── clean.xml
    ├── container.xml
    ├── get.xml
    ├── run.xml
    └── update.xml
```
If there are multiple tests in the group, then each one of them will need that
file tree of resource description.


### Providing Expected Results

Similar to inputs, results provided at the same path (`${UnitTest_reference}`) as
archives named `results*.ca`.  They will be used by the task `check` to compare
the outputs produced by the program tested and the ones expected by the test.


---

## Interactive debugging

You can debug interactively a MIDAS program by launching a job and
login into it.  For that, you can set
```bash
UnitTest_stop_for_interactive_work=yes
```
in `maestro/suites/midas_system_tests/experiment.cfg`
to prepare the working directory for interactive debugging.

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

The script `launch_interactive.sh` launches an interactive job and
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
script `launch_program.sh` supports two tools for debugging:
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
export MIDAS_COMPILE_CODECOVERAGE_DATAPATH='an absolute path where the code coverage data can be saved'
```
before compiling the program.  Keep in mind that this directory will
contain more than 10GB of data so make sure you specify a directory
with enough space.

After running all the tests with those binaries, you can obtain a code
coverage report using [`codecov.sh`](src/codecov.sh):
```
>>> cd src
>>> ./codecov.sh -h

 *** SEQUENCE D'APPEL ***

 [positionnels]
 IN       -cov_dir [:] code coverage directory for files generated by the execution (default: value of env variable 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH'
 IN       -src_dir [${PWD}/../src:${PWD}/../src]  (default: result of 'git rev-parse --show-toplevel')
 IN       -web_dir [${HOME}/public_html/midas/codecoverage-${revnum}:${HOME}/public_html/midas/codecoverage-${revnum}] directory for storing html pages (default: ${HOME}/public_html/midas/codecoverage-${revnum})
 IN       -setmx [no:yes] Insert 'set -x' in the script (default: no)
          [-- positionnels]

```

A [code coverage report has been generated for version
`v_3.8.1-516-g746c074` if you want to have a
look](http://goc-dx.science.gc.ca/~erv000/midas/codecoverage-v_3.8.1-516-g746c074/CODE_COVERAGE.HTML)

## Updating Test Results

The results can be updated by running the task `UnitTest/update` for
the wanted test.  The new results will be stored in the directory specified
with the variable `UnitTest_reference_update` in the test
configuration file (the task will abort if the directory already exists).

The listings will be collected at the same time.

Once the new results have been collected, you can update the variable
`UnitTest_reference` with the path to the new results. However, it is 
recommended to replace `${USER}` with your actual username to allow other
users to run the test using your updated results.

You can then commit the changes to the configuration.

When you will ask a merge-request, the new results will be copied in a
safe directory with all other reference results by one of the
git repository maintainers.

## Splitting an observation file to prepare a unit test input

In the case a user wants to split an observation file to prepare a
unit test input, this program can be called like this (after loading an 
appropriate MIDAS ssm package):
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
. ssmuse-sh -d eccc/mrd/rpn/anl/midas/3.8.1
```
or point directly to the program in the SSM domain:
```
/fs/ssm/eccc/mrd/rpn/anl/midas/3.8.1/${ORDENV_PLAT}/bin/midas.splitobs.Abs
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
