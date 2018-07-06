## Contributing Guide
----

### Create an issue

 - Go to the **gitlab** [*issue* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/derivate/issues) and press `New issue` button.
 - If available, choose a *template* from the drop-down menu.
 - Add a description and fill out the template.
 - Then press `Submit issue`.

### Obtain a copy of the code and modify in a local branch :

#### Get a local copy of the master branch :
```bash
. setcentraldepot.dot unset;. setcentraldepot.dot git@gitlab.science.gc.ca:atmospheric-data-assimilation/derivate.git
clone_suite -c master derivate_master_issue${ISSUE_NUMBER}
```

 - `derivate_master_issue${ISSUE_NUMBER}` is the name of the suite directory in `$HOME/.suites`
 - a local branch is created named : `$USER.derivate_master_issue${ISSUE_NUMBER}`

#### Make modifications to your local branch :
Edit and commit your code as follows :
```bash
cd $HOME/.suites/derivate_master_issue${ISSUE_NUMBER}
```
(edit some files)
```bash
git add $file1 $file2
git commit
```
(Write a commit message.  Please refer to [the wiki for guidelines](https://wiki.cmc.ec.gc.ca/wiki/Assimilation/env#Faire_des_commits_avec_Git).)

### Push your code to the [gitlab project](https://gitlab.science.gc.ca/atmospheric-data-assimilation/derivate) : 

Once the code is ready to share (not necesserily final) or ready to merge :
```bash
git push central HEAD
```

  - The `git push` command can be done many times during the course of development.
  - to get a status of the work done, first run `git fetch central`, and then `git log --oneline central/master..HEAD` 
 or 
`git diff central/master`.

### Synchronise your branch with the master branch : 

Once the code related to the proposed change is ready, make sure your local branch is up-to-date with the `master` branch on **gitlab** :

```bash
git fetch central
git log --oneline ..central/master
```
If the command returns something

```bash
git branch <name of branch_v1>
git rebase central/master $USER.derivate_master_issue${ISSUE_NUMBER}
```

Resolve any possible conflicts and do a final :
```bash
git push central HEAD
```

### Create and assign the `Merge Request`

  - Go to **gitlab** [*branch* page](https://gitlab.science.gc.ca/atmospheric-data-assimilation/derivate/branches) and press on `Merge Request` button appearing next to the appropriate *branch* (`$USER.derivate_master_issue${ISSUE_NUMBER}`).
  - Verify the information, in particular the source should be `$USER.derivate_master_issue${ISSUE_NUMBER}` and the target branch should be `master`.
  - Assign the `Merge Request` to colleague.
  - And press `Submit merge request`.
  - The assignee is expected to review the code, discuss/comment and finally accept the `Merge Request`, which will close the *issue* and delete the *branch*.
  - If problems or better solutions come up from the review, the `Merge Request` can be closed without merging the code.
