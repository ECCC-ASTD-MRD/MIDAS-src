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
* `v_3.4` branch
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.4/README.md) 
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.4) (IC-2)

# Contributing

We strongly suggest anyone considering to contribute to the MIDAS
 project, to follow the workflow documented in the [contributing
 guide](CONTRIBUTING.md).

# Getting a local copy of the code

To simply get a local copy of the code from an existing branch
associated with an issue, we suggest the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.3
clone_projet --no-central -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
or if one is interested in the latest version of the master branch
```bash
clone_projet --no-central -c master git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-master
```

## Getting code related to IC-3 system
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.3
clone_projet --no-central -c v_3.6 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.6
```

## Getting code related to IC-2 (operational) system
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.3
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
Running the script without argument defaults to running it sequentially for all 
programs, but that can be too long for the PPPs and it often gets killed.
As a workaround, one can either launch it only for the modified (or new) 
programs or launch it in parallel for all programs, from `src/programs/`:
```sh
for pgm in $(ls *.f90); do ./make_src_files.sh $pgm & done
```

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
`ubuntu-18.04-amd64-64` and `sles-11-broadwell-64-xc40`.  The
compiling is done in parallel.  You can call it with:
```bash
cd src/programs
./compile_all_plat.sh
```

# Testing the new compilation using `make`
All the information and instructions can be found in the [src](./src) directory, by clicking on that folder in gitlab and browsing down the files, you will see [`src/README.md`](./src/README.md) content appended there.

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

The programs compiled in the CI pipeline are moved into the directory:
```
/home/sanl888/data_maestro/ords/midas/gitlab-ci/abs
```

Once in a while, we must clean this directory to save them elsewhere
to avoid filing the `ords` directory of user `sanl888`.

The script
```
/home/sanl888/data_maestro/ords/midas/gitlab-ci/abs/move_abs.sh
```
moves the MIDAS programs to
```
/home/sanl888/data_maestro/eccc-ppp4/midas/gitlab-ci/abs
```
and links them to keep a reference.
