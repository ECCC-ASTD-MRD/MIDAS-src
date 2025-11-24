# Automatic Testing using GitLab-CI

Here are the instructions to install and maintain the automatic launch of
the MIDAS test suite and automatic generation of the documentation for
this code.

It is based on the [GitLab CI](https://docs.gitlab.com/ce/ci)
functionalities which are available with our code repository.

It is configured such that, for each push in the `main` branch, the
test suite is launched.  If the test suite compiles, runs the programs
with error, the documentation is generated.  Then, since the link in
the main [README](README.md) to the documentation is generic, the end
user will access the latest documentation for the `main` branch.

## Cleaning of the programs directory

The programs compiled in the CI pipeline are copied into the directory:
```
/home/sanl888/data_maestro/ords/midas/gitlab-ci/abs
```

If everything goes as planned they should be moved automatically to
```
/home/sanl888/data_maestro/midas/midas/gitlab-ci/abs
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
/home/sanl888/data_maestro/midas/midas/gitlab-ci/abs
```

## Install and Register the GitLab Runner

Here are the steps to install the a `gitlab-runner`.  For now, the
`gitlab-runner` will be used for automatic generation of the
documentation as it is configured in the file
[`.gitlab-ci.yml`](.gitlab-ci.yml).  But, it can be used for any other
automatic task.

First, one must download a program, called `gitlab-runner`, which
listens to the GitLab server for triggers to CI.  This program must be
comptatible with the GitLab API version where the code is hosted.  For
the original CI install, we got the program with the command:
```bash
wget https://gitlab-ci-multi-runner-downloads.s3.amazonaws.com/v1.11.5/binaries/gitlab-ci-multi-runner-linux-amd64
chmod +x gitlab-ci-multi-runner-linux-amd64
```

But since, RPN-SI has been developping a specialized version for our
own HPC environment.

### The CMC-owned running

@phc001 worked on a Gitlab runner which can submit jobs.  The path to
that runner is:
```
/home/sici000/bin/gitlab-runner-science-15.13.0
```
It allows to specify resources for a ̀ord_soumet` job submission
through variables `ORD_SOUMET_*`.

Then a runner has to be registered to the GitLab server.  You must
execute that command:
```bash
/home/sici000/bin/gitlab-runner-science-15.13.0 register      \
         --non-interactive                                    \
         --url https://gitlab.science.gc.ca                   \
         --registration-token ${GITLAB_CI_TOKEN}              \
         --description "GitLab runner running under user '${USER}' on '${TRUE_HOST}' using 'ordsoumet' executor."    \
         --tag-list  ${runner_tag}                            \
         --executor  ordsoumet                                \
         --builds-dir   ${HOME}/data_maestro/ords/midas/gitlab-ci/builds \
         --cache-dir    ${HOME}/data_maestro/ords/midas/gitlab-ci/cache
```

where the `${GITLAB_CI_TOKEN}` is the token found in the [CI
settings](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/pipelines/settings)
of the project (put the `Runner token`) and `${runner_tag}` is a tag
to identify which runner to use.  Here `${runner_tag}` can be
`hpcr-u2` or `hpcr-u3`.

At the end, the `${HOME}/.gitlab-runner/config.toml` should contain
the information above.

Congratulations, you have registered your GitLab runner.  You can
confirm that it has been configured correctly by looking at the
[runners page of the
project](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/runners).

### Concurrency

For several pipeline jobs can be run concurrently by the runner, one
must edit the file `${HOME}/.gitlab-runner/config.toml' to set:
```
concurrent = 10
```
We only test `10` as a value.

### Registering, step by step

You can do a step by step configuration by following the instructions
below.

Then, one has to [register the
runner](https://docs.gitlab.com/runner/register) with the command:
```bash
/home/sici000/bin/gitlab-runner-science-15.13.0 register
```
The program will ask for the `gitlab-ci coordinator URL' which is in
our case:
```
https://gitlab.science.gc.ca/ci
```

Then it will ask for the `gitlab-ci` token which can be found in the [CI
settings](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/pipelines/settings)
of the project (put the `Runner token`).  You also have to put a
description.  I put this:
```
Runner for '${USER}' connected to 'gitlab.science.gc.ca:atmospheric-data-assimilation/midas'
```

It will ask for tags which you can ignore.  The next question is the
`executor` for which we want a `ssh` and you put
`ppp7.science.gc.ca` as the `SSH server address` when asked, then
the script will be executed by doing a SSH connection to
`ppp7.science.gc.ca`.

The last questions are the port of the SSH server which is `22` and DO
NOT ENTER YOUR PASSWORD, just do return and it will ask you the path
to the `SSH identity file` which is `${HOME}/.ssh/id_rsa`.

## Launch the GitLab Runner

We suggest to launch the GitLab Runner in the daemon queue of one of the PPPs.  To do that, we must create a job:
```bash
[ ! -d ~/bin ] && mkdir -v ~/bin
cat > ~/bin/gitlab_runner.sh <<EOF
#!/bin/bash

set -ex

runhost=\${1:-ppp7}
qname=dev_daemon

gitlabrunner_exists=true
jobst -c \${runhost} -u \${USER} -q \${qname} | grep gitlab || gitlabrunner_exists=false

if [ "\${gitlabrunner_exists}" != true ]; then
    cat > \${TMPDIR}/gitlab_runner <<ENDOFGITLABRUNNER
#!/bin/bash

set -ex

/home/sici000/bin/gitlab-runner-science-15.13.0 --log-level debug run
ENDOFGITLABRUNNER

    ord_soumet \${TMPDIR}/gitlab_runner -mach \${runhost} -queue \${qname} -cpus 1 -w \$((90*24*60))

    rm \${TMPDIR}/gitlab_runner
fi
EOF
chmod +x ~/bin/gitlab_runner.sh
~/bin/gitlab_runner.sh
```

This script will launch a job on the queue `dev_daemon` (which has no time limit) on `ppp7`.

### Maintain the runner with `hcron`

To install a `hcron` rule to check if the gitlab runner is running, do this
```bash
mkdir -pv ~/.hcron/hcron-dev7.science.gc.ca/events/ppp7
cat > ~/.hcron/hcron-dev7.science.gc.ca/events/ppp7/gitlab-runner <<EOF
as_user=
host=ppp7.science.gc.ca
command=echo ~/bin/gitlab_runner.sh | bash --login
notify_email=
notify_message=
when_month=*
when_day=*
when_hour=*
when_minute=$((RANDOM % 60 ))
when_dow=*
EOF
ssh hcron-dev7 hcron reload
```
