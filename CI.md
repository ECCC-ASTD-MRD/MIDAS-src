# Automatic Testing using GitLab-CI

Here are the instructions to install and maintain the automatic launch of
the MIDAS test suite and automatic generation of the documentation for
this code.

It is based on the [GitLab CI](https://docs.gitlab.com/ce/ci)
functionalities which are available with our code repository.

It is configured such that, for each push in the `master` branch, the
test suite is launched.  If the test suite compiles, runs the programs
with error, the documentation is generated.  Then, since the link in
the main [README](README.md) to the documentation is generic, the end
user will access the latest documentation for the `master` branch.

Here are the steps to install the a `gitlab-runner`.  For now, the
`gitlab-runner` will be used for automatic generation of the
documentation as it is configured in the file
[`.gitlab-ci.yml`](.gitlab-ci.yml).  But, it can be used for any other
automatic task.

## Register the GitLab Runner

First, one must download a program, called `gitlab-runner`, which
listens to the GitLab server for triggers to CI.  This program must be
comptatible with the GitLab API version where the code is hosted.  In our
case, we must download the program with the command:
```bash
wget https://gitlab-ci-multi-runner-downloads.s3.amazonaws.com/v1.11.5/binaries/gitlab-ci-multi-runner-linux-amd64
chmod +x gitlab-ci-multi-runner-linux-amd64
```
But, a copy has already been installed here:
```
/home/sidr000/bin/gitlab-ci-multi-runner-linux-amd64
```

Then a runner has to be registered to the GitLab server.  You must
execute that command:
```bash
/home/sidr000/bin/gitlab-ci-multi-runner-linux-amd64 register \
         --non-interactive                                    \
         --url https://gitlab.science.gc.ca/ci                \
         --registration-token ${GITLAB_CI_TOKEN}              \
         --description "Some phrase describing the runner"    \
         --run-untagged                                       \
         --executor shell                                     \
         --builds-dir   ${HOME}/data_maestro/ords/midas/gitlab-ci/builds \
         --cache-dir    ${HOME}/data_maestro/ords/midas/gitlab-ci/cache
```

where the `${GITLAB_CI_TOKEN}` is the token found in the [CI
settings](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/pipelines/settings)
of the project (put the `Runner token`) and the `description` is a
phrase describing the runner.

At the end, the `${HOME}/.gitlab-runner/config.toml` should contain
the information above.

Congratulations, you have registered your GitLab runner.  You can
confirm that it has been configured correctly by looking at the
[runners page of the
project](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/runners).

### Registering, step by step

You can do a step by step configuration by following the instructions
below.

Then, one has to [register the
runner](https://docs.gitlab.com/runner/register) with the command:
```bash
/home/sidr000/bin/gitlab-ci-multi-runner-linux-amd64 register
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
Runner for 'erv000' connected to 'gitlab.science.gc.ca:atmospheric-data-assimilation/midas'
```

It will ask for tags which you can ignore.  The next question is the
`executor` for which we want a `ssh` and you put
`eccc-ppp1.science.gc.ca` as the `SSH server address` when asked, then
the script will be executed by doing a SSH connection to
`eccc-ppp1.science.gc.ca`.

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

runhost=\${1:-ppp1}
qname=dev_daemon

gitlabrunner_exists=true
jobst -c \${runhost} -q \${qname} | grep "gitlab_runner *\${USER}" || gitlabrunner_exists=false

if [ "\${gitlabrunner_exists}" != true ]; then
    cat > \${TMPDIR}/gitlab_runner <<ENDOFGITLABRUNNER
#!/bin/bash
set -ex

/home/sidr000/bin/gitlab-ci-multi-runner-linux-amd64 run
ENDOFGITLABRUNNER

    ord_soumet \${TMPDIR}/gitlab_runner -mach eccc-\${runhost} -queue \${qname} -cpus 1 -w \$((90*24*60))
    rm \${TMPDIR}/gitlab_runner
fi
EOF
chmod +x ~/bin/gitlab_runner.sh
~/bin/gitlab_runner.sh
```

This script will launch a job on the queue `dev_daemon` (which has no time limit) on `eccc-ppp1`.

### Maintain the runner with `hcron`

To install a `hcron` rule to check if the gitlab runner is running, do this
```bash
mkdir -pv ~/.hcron/hcron1.science.gc.ca/events/eccc-ppp1
cat > ~/.hcron/hcron1.science.gc.ca/events/eccc-ppp1/gitlab-runner <<EOF
as_user=
host=\$HCRON_EVENT_NAME[1]
command=echo ~/bin/gitlab_runner.sh | bash --login
notify_email=
notify_message=
when_month=*
when_day=*
when_hour=*
when_minute=$((RANDOM % 60 ))
when_dow=*
EOF
ssh hcron1 hcron-reload
```