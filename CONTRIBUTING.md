## Contributing Guide
----

### Create an issue

 - Go to the **GitLab** [*issue* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues) and press `New issue` button.
 - Choose a *template* from the drop-down menu.
 - Add a description and fill out the template.
 - Then press `Submit issue`.

### Obtain a copy of the code and modify in a local branch :

#### Get a local copy of the code

 - in the issue page, there is a button `New Branch`
   - that button will create automatically a new branch at the same point as the `master` branch

```bash
. ssmuse-sh -d eccc/cmd/cmdi/utils/2.5
. setcentraldepot.dot unset  ## needed only if 'setcentraldepot.dot' has already been called
clone_projet -c ${ISSUE_NUMBER} git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git midas-${ISSUE_NUMBER}
```

#### Make modifications to your local branch :

Edit and commit your code as follows :
```bash
cd midas-${ISSUE_NUMBER}
```
(edit some files)
```bash
git add $file1 $file2
git commit
```
(Write a commit message.  The editor will open and a template for commit message will appear)

### Push your code to the [GitLab project](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas):

Once the code is ready to share (not necessarily final) or ready to merge :
```bash
git push origin HEAD
```

  - The `git push` command can be done many times during the course of development.
  - to get a status of the work done, first run
    - `git fetch origin`, and then
    - `git log --oneline origin/master..HEAD` or
    - `git diff origin/master`.

### Synchronise your branch with the master branch : 

Once the code related to the proposed change is ready, make sure your
local branch is up-to-date with the `master` branch on **GitLab** :

```bash
git fetch origin
git log --oneline ..origin/master
```

If the command returns something, that means the branch `master` has
evolved since the time you created your branch.  To synchronize your
branch, you can *rebase* your changes in your branch on top of the
latest commit in the `master` branch in the GitLab project.

It is often a good idea to create a new branch to work on (which
starts out as a copy of your existing branch). It is mandatory to do
this if your branch has already been pushed to the central repository.

```bash
## ...not mandatory, if the branch has not yet been pushed
##                or if you dont mind to change the working branch in-place
git checkout -b <name of branch>_v1
git pull --rebase origin master
```

Resolve any possible conflicts and do a final :
```bash
git push origin HEAD
```


### Create and assign the `Merge Request`

  - Go to **gitlab** [*branch* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/branches) and press on `Merge Request` button appearing next to the appropriate *branch* (`${ISSUE_NUMBER}-...`).
  - Verify the information, in particular the source branch should be the branch you worked on and the target branch should be `master`.
  - Select the template, then follow the guidelines when filling in the Description box.
  - Assign the `Merge Request` to a colleague.
  - And press `Submit merge request`.
  - The assignee is expected to review the code, discuss/comment and finally accept the `Merge Request`, which will:
    - merge the code for the branch into `master`,
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
