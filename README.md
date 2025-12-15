
# Contributing

We strongly suggest anyone considering to contribute to the MIDAS
 project, to follow the workflow documented in the [contributing
 guide](CONTRIBUTING.md).

# MIDAS Fortran coding standards:

* [List of standards](docs/codingStandards.md)
* [Automatic documentation standards](docs/documentationStandards.md)

# MIDAS releases general and code documentations:

The documentation for officially supported branches is available:
* `main` branch
  * [General documentation (`README.md`) - this page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/main/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-main)
* `v_4.1` branch: includes non-backward compatible changes with respect to `v_3.9` and `v_3.10`
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_4.1/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_4.1)
* `v_3.10` branch
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.10/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.10)
* `v_3.10-RandD` branch: The code in this branch validates with `v_3.10` but contains new features for testing using IC-4 final cycles as the reference.
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.10-RandD/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.10-RandD)
* `v_3.9` branch: related to IC-4 release
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.9/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.9)
* `v_3.9-RandD` branch: The code in this branch validates with `v_3.9` but contains new features for testing using IC-4 final cycles as the reference.
  * [General documentation (`README.md`)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/v_3.9-RandD/README.md)
  * [Fortran code documentation](http://goc-dx.science.gc.ca/~sanl888/midas-sphinx-doc/latest-v_3.9-RandD)

# Getting a local copy of the code

To simply get a local copy of the code from an existing branch
associated with an issue, we suggest the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
or if one is interested in the latest version of the `main` branch
```bash
clone_projet --no-central -c main git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-main
```

To use this version, it is necessary to load LibRMN 20 utilities
starting from version `rpn/utils/20231219`.

## Getting the code for tag `v_4.0.0`

The code related to the `v_4.0.0` release is available with the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_4.0.0 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-4.0.0
```

To use this version, it is necessary to load LibRMN 20 utilities
starting from version `rpn/utils/20231219`.

## Getting code related to IC-4 implementation on HPCR-U2

The official code related to IC-4 implementation is available with the command:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_3.9 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.9
```

One can also use the branch `v_3.9-RandD`.  The code in this branch
validates with `v_3.9` but contains new features for testing using
IC-4 final cycles as the reference.
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_3.9-RandD git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.9-RandD
```


If you created a new branch with the GitLab web UI without choosing the origin branch, then it
has been created using the default branch which is `main`.  One must
reset it to the release branch.  One can simply do:
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_3.9 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
cd midas-${ISSUE_NUMBER}
git checkout -b ${ISSUE_NUMBER}-complete-the-name-of-the-branch-as-on-GitLab
git push origin ${ISSUE_NUMBER}-complete-the-name-of-the-branch-as-on-GitLab --force-with-lease
```

## Getting code related to IC-4 implementation on HPCR-U3

```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_3.10 git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.10
```

One can also use the branch `v_3.10-RandD`.  The code in this branch
validates with `v_3.10` but contains new features for testing using
IC-4 HPCR-U3 final cycles as the reference.
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.6
clone_projet --no-central -c v_3.10-RandD git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-3.10-RandD
```

# Compiling MIDAS

[`src/midas_build`](./src/README.md)
is the official compilation tool to build MIDAS.

To compile all programs (`src/programs/*.f90` as well as
[`splitobs`](./src/README.md#splitobs-an-external-program)),
simply do:
```bash
cd ${where_your_code_is}
cd src
./midas_build
```

## Complete documentation on using `midas_build`, `cmake` and `make`

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
