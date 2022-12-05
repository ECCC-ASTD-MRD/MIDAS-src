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

* [List of standards](docs/codingStandards.md)
* [Automatic documentation standards](docs/documentationStandards.md)

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

## `dumpBmatrix`

This standalone program dump 1D covariance values from the binary file `Bmatrix.bin`.

See [`scripts/convenient_tools/dumpBmatrix/README.md`](scripts/convenient_tools/dumpBmatrix/README.md) for more details.
