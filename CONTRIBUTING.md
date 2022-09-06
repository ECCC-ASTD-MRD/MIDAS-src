# Contributing Guide

The general contribution workflow described in this guide can be summarized in these steps

1. [Create an issue](#create-an-issue)
2. [Make modifications in a local branch](#make-modifications-in-a-local-branch)
3. [Merge back your branch to the main one](#merge-back-your-branch-to-the-main-one)

## Create an issue

Anyone can open an issue, even if they don't intend to actually address it 
themself; the author of an issue does not have to be the assignee.
Opening an issue to signify and document a bug or a feature request is an
important part of the workflow.

To create an issue:

1. Go to the **GitLab** [*issue* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues)
  and press [`New issue`](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/new) button.
2. Choose a *template* from the drop-down menu.
3. Add a description and fill out the template.
4. Then press `Submit issue`.



## Make modifications in a local branch

Once assigned to someone, the assignee should create an issue branch and pull
it locally on their workstation to start working on it.

### Get a local copy of the code

**Create a branch**  
in the issue page, just under the issue description there is green drop down green menu
choose `Create a branch`. 
That button will create automatically a new branch at the same point as the `main` branch

**Get a local copy of that branch**  
```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
. setcentraldepot.dot unset  ## needed only if 'setcentraldepot.dot' has already been called
clone_projet -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```
(For this to work, you need to create a file `~/.gitlab-private-token`
like explained in the [Git documentation (see section "Interaction avec GitLab en utilisant des scripts")](https://wiki.cmc.ec.gc.ca/wiki/Git/Doc#Outil_de_partage_de_code:_gitlab.science.gc.ca))

When these steps have been completed make sure to move the issue to the ~Doing
column on the [issue board](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/boards)

### Make modifications to your local branch


**In your filesystem, go to the local directory associated with the issue**  
```bash
cd midas-${ISSUE_NUMBER}
```
**edit some files**  
following the MIDAS [coding standards](docs/codingStd_top10.md)
```
${EDITOR} $file1 $file2
```
**stage your modifications**  
```bash
git add $file1 $file2
```
**when you are satisfied with a group of modifications, you can commit them**  
```bash
git commit
```
Write a commit message.  The editor will open and a template for commit message will appear.

This is a very simplified contribution workflow.
You should also use `git status` and `git diff [--cached]` when staging and before committing your changes.

## Merge back your branch to the main one

### Push your code to the [GitLab project](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas):

Once the code is ready to share (not necessarily final) or ready to merge :
```bash
git push origin HEAD
```
(Note: if you are using `set_centraldepot`, the remote depot is named `central` instead of `origin`)

  - The `git push` command can be done many times during the course of development.
  - to get a status of the work done, first run
    - `git fetch origin`, and then
    - `git log --oneline origin/main..HEAD` or
    - `git diff origin/main`.

### Synchronise your branch with the `main` branch:

Once the code related to the proposed change is ready, make sure your
local branch is up-to-date with the `main` branch on **GitLab** :

```bash
git fetch origin
git log --oneline ..origin/main
```

If the command returns something, that means the branch `main` has
evolved since the time you created your branch.  To synchronize your
branch, you can *rebase* your changes in your branch on top of the
latest commit in the `main` branch in the GitLab project.

It is often a good idea to create a new branch to work on (which
starts out as a copy of your existing branch). It is mandatory to do
this if your branch has already been pushed to the central repository.

```bash
## ...not mandatory, if the branch has not yet been pushed
##                or if you dont mind to change the working branch in-place
git checkout -b <name of branch>_v1
git pull --rebase origin main
```

Resolve any possible conflicts and do a final :
```bash
git push origin HEAD
```

### Testing your modifications
An automatic system of tests has been developed that garantee that modifications
do not brake existing features or change expected results.

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

### Create and assign the `Merge Request`

  - Go to **gitlab** [*branch* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/branches) and press on `Merge Request` button appearing next to the appropriate *branch* (`${ISSUE_NUMBER}-...`).
  - Verify the information, in particular the source branch should be the branch you worked on and the target branch should be `main`.
  - Select the template, then follow the guidelines when filling in the Description box.
  - Assign the `Merge Request` to a colleague.
  - And press `Submit merge request`.
  - Move the issue to the ~"Under Review" column on the [issue board](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/boards)
  - The assignee is expected to review the code, discuss/comment and finally accept the `Merge Request`, which will:
    - merge the code for the branch into `main`,
    - close the associated *issue* and
    - delete the *branch*.
  - If problems or better solutions come up from the review, the `Merge Request` can be closed without merging the code.

A merge request does not need to be a final review step.  You can use
it as a development process to share code with colleagues.  To prevent
the branch to be merged accidentally, you can prefix the title of the
merge request with [`WIP: ` ("Work In
Progress")](https://docs.gitlab.com/ce/user/project/merge_requests/work_in_progress_merge_requests.html).
You will see instructions about this feature in the GitLab merge
request page.

Once the `Merge Request` is accepted with with the `main` branch as target (as 
is generally the case), the modifications are pushed and the system tests are
automatically launched to guarantee that all the tests pass for the `main` branch.  The [instructions for automatic
testing using GitLab-CI are available in a separate file](docs/CI.md).


# Advanced Topics

* [Creating SSM packages](docs/ssm.md)
* [Managing Continuous Integration](docs/CI.md) 
* [Advanced Unit Testing Topics](docs/unitTests.md) 
    * [Updating Test Results](docs/unitTests.md#updating-test-results)
    * [Interactive debugging](docs/unitTests.md#interactive-debugging)
* [Adding New Observations](docs/newObs.md)
