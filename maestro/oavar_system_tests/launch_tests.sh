#!/bin/bash

. ssmuse-sh -d arma/maestro/dev/2.8
. ssmuse-sh -d isst/maestro/1.4.3-rc4

if [ -a ~/.suites/oavar_system_tests ]; then
  echo ' '
  echo 'The maestro suite oavar_system_tests already exists...'
  echo 'Therefore simply launching xflow without doing any setup.'
  echo ' '
else
  ln -s $PWD ~/.suites/
  make_links oavar_system_tests
fi

export SEQ_EXP_HOME=~/.suites/oavar_system_tests
xflow &
