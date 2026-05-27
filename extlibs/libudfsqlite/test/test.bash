#!/bin/bash

set -e

. d.compile

set -x

echo +++
export LD_LIBRARY_PATH=../src/:$LD_LIBRARY_PATH
echo ===

../bin/d.sqlite junk.sqlite "select concat('cmda', 'est', 'cool')" > out2.test
grep . out2.test | grep -v 'DSQLITE.*=' > out.test
rm out2.test
rm junk.sqlite

diff -B out.expected out.test
