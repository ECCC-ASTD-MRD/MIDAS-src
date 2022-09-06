# SSM

The [CI](CI.md) has been configured to produce a SSM domain under
```bash
/fs/ssm/eccc/mrd/rpn/anl/midas
```
automatically when a tag is pushed.

## Procedure to follow when creating a new version

### Our tag convention
Once a release is decided to be published, identify the version name
by following the [semantic versioning](http://semver.org/spec/v2.0.0.html).

Tags will use the prefix `v_` and contain 3 fields (for example `v_A.B.C`):

* the first field (`A`) will change after a version is delivered to Operations
* the second field (`B`) will be updated when one of those cases happen:
    * the results are changing
    * the API is modified in a non-backward-compatible way
        * a namelist variable is changed
        * an input file is changed
        * a new mandatory input file is added
        * an output file is added
* the third digit (`C`) is for backward-compatible changes:
    * bugfix
    * new functionality that does not affect other programs (for example, adding SQLite observations files).

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

### Merge changes from release branch to `main`

Once a version is published, there are probably some changes (like
bugfixes) that needs to be also made in the `main` branch which is our
main development branch.

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

