# Breaking change!

If you compile MIDAS code after tag `v_3.6.0-a1`, then you must update your
login profile to version `1.11.0`:
```bash
ln -svi /fs/ssm/eccc/mrd/ordenv/profile/1.11.0 ~/.profile_1.11.0
rm -v ~/.profile && ln -svi .profile_1.11.0 ~/.profile
```

This change is backward compatible for your suites but you absolutely
need to update your profile to compile any MIDAS code after version
`v_3.6.0-a1`.

To know if your code is after `v_3.6.0-a1`, you can execute:
```bash
git describe
```

If the output is containing the string `v_3.6.0`, then you need to
update your profile.  If not, then you have to use the previous
version of the profile which is
`/fs/ssm/eccc/mrd/ordenv/profile/1.10`.

# MIDAS Fortran coding standards:

* [List of standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards)
* [Automatic documentation standards](https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Documentation_Standards)

# MIDAS releases general and code documentations:

The documentation for officially supported branches is available:
* `master` branch
  * [General documentation (`README.md`) - this page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/master/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-master)
* `v_3.6` branch 
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.6/README.md) 
  * [Fortran code documentation](http://hpfx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.6) (IC-3)

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
or if one is interested in the latest version of the master branch
```bash
clone_projet --no-central -c master git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-master
```

## Getting code related to IC-3 system

```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
clone_projet --no-central -c v_3.6 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.6
```

If you created a new branch with the GitLab web UI, then the branch
has been created using the default branch which is `master`.  One must
reset it to the release branch.  One can simply do:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
clone_projet --no-central -c v_3.6 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
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

# SSM

The [CI](CI.md) has been configured to produce a SSM domain under
```bash
/fs/ssm/eccc/mrd/rpn/anl/midas
```
automatically when a tag is pushed.

## Procedure to follow when creating a new version

Once a release is decided to be published, identify the version name
by following the [semantic
versioning](http://semver.org/spec/v2.0.0.html).  The tag name
will be prepend by `v_`.  So if the version is `3.6.6`, the tag name
will be `v_3.6.6`.

### Warning!

First, avoid to create a SSM domain on the last business day of a week
(for example, a Friday).  Although we took many precautions, we are
creating files under `/fs/ssm/eccc/mrd/rpn/anl/midas` a directory
directly used by the CMC Operations.  We do not want to interrupt the
operational system by doing a mistake in R&D!

Now, let's detail the procedure to publish a new MIDAS version.

### Update CHANGELOG

When the version name is set, modify the [`CHANGELOG`](CHANGELOG.md)
by replacing `[Unreleased]` by the version name and reintroduce the
`[Unreleased]` section with empty subsections.  You can take example
on the commit 6136c4241b5016f5241bf868f73a10d2b84d3504 which did this
change for version `v_3.6.6`.

Once this changelog is done, commit with the command
```bash
git add -v CHANGELOG.md
git commit -m "Prepare CHANGELOG for version 'v_${VERSION}'"
```
and push this change and wait for the [CI automatic
tests](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/pipelines)
to be finished.  This is **very** **very** important since if you push
the tag immediately, the current CI pipeline actived by this commit
will use the new tag in the program names and the next CI pipeline,
when you will push the tag, won't work because it cannot overwrite any
programs already generated.

Do not ask to avoid running the CI by including some string like
`[skip CI]` because then, when the tag will be pushed, the CI pipeline
will not be triggered.

### Create the tag

When the CI pipeline for the CHANGELOG commit is done, you can create
an annotated tag by prepending `v_` in front of the version name.  You
can use this command to create the tag:
```bash
git tag -a v_${VERSION} -F - <<EOF
This version is available in the SSM domain:
    /fs/ssm/eccc/mrd/rpn/anl/midas/${VERSION}

See CHANGELOG for more details.
EOF
```

### Create the SSM domain

Then you push the tag:
```bash
git push origin v_${VERSION}
```
and you can monitor the [CI
pipeline](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/pipelines)
to check if
 * all the programs are compiled (`build` stage),
 * all tests are running correctly (`test` stage),
 * the documentation is generated (`doc` stage) and
 * the SSM domain is published under `/fs/ssm/eccc/mrd/rpn/anl/midas/${VERSION}` (`deploy` stage).

### Merge changes from release branch to `master`

Once a version is published, there are probably some changes (like
bugfixes) that needs to be also made in the `master` branch which is
our main development branch.

We suggest to open a merge request using the title "Merge tag
'v_${VERSION}'" which describes the changes that will be introduced.
See for example, what has been done in the merge request !476.

## Updating the scripts under `sanl000`

For security considerations, the scripts that the user `sanl000` is
using are not coming directly from the MIDAS depot itself but reviewed
copies under his control.  So when scripts under directory `ssm` are
modified, we must update them manually.

If the version of MIDAS is `${VERSION}` in the directory
`${MIDAS_SOURCE_CODE}`, then here are the commands the user `sanl000`
has to do for this update:
```bash
cd ${HOME}/ssm

mkdir -v midas/${VERSION}
cd midas/${VERSION}
cp -vi ${MIDAS_SOURCE_CODE}/ssm/publish .
cp -vi ${MIDAS_SOURCE_CODE}/ssm/ssm_publish .
cp -vi ${MIDAS_SOURCE_CODE}/ssm/post-install .

cd ..  ## current directory is now '${HOME}/ssm/midas'

## update the script '${HOME}/ssm/midas/post-install'
echo "Removing ${PWD}/post-install which is now pointing to $(true_path post-install)"
rm -v post-install
ln -svi ${VERSION}/post-install .

## update the script '${HOME}/ssm/midas/publish'
echo "Removing ${PWD}/publish which is now pointing to $(true_path publish)"
rm -v ssm_publish
ln -svi ${VERSION}/publish .

cd ..  ## current directory is now '${HOME}/ssm'

## update the script '${HOME}/ssm/ssm_publish'
echo "Removing ${PWD}/ssm_publish which is now pointing to $(true_path ssm_publish)"
rm -v ssm_publish
ln -svi midas/${VERSION}/ssm_publish .
```

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

# Automatic Testing using GitLab-CI

An automatic system of tests has been developed.  For each push in the
`master` branch the system tests are launched to guarantee that the
all the tests pass for the `master` branch.  The [instructions for
automatic testing using GitLab-CI are available in a separate
file](CI.md).

## Cleaning of the programs directory

The programs compiled in the CI pipeline are copied into the directory:
```
/home/sanl888/data_maestro/ords/midas/gitlab-ci/abs
```

If everything goes as planned they should be moved automatically to
```
/home/sanl888/data_maestro/ppp5/midas/gitlab-ci/abs
```

But, if the files do not get moved automatically, they may fill the
`ords` directory.

Once in a while, we must clean this directory and move the programs
elsewhere to avoid filing the `ords` directory of user `sanl888`.

You can use the script
```
/home/sanl888/data_maestro/ords/midas/gitlab-ci/abs/move_abs.sh
```
which moves the MIDAS programs to
```
/home/sanl888/data_maestro/ppp5/midas/gitlab-ci/abs
```
