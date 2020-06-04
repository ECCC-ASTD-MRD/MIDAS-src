## `midas.monitor`

This program contains a loop to check for the content of a file (given as the first argument):
 * if it contains `VAR3D_STATUS=REBM_DONE\n` then execute a script (given as the second argument)
 * if it contains `VAR3D_STATUS=VAR3D_END\n` then exits.

This is used in script `midas.launch` which monitors the content of
the file `VAR3D_STATUS.dot` to launch other task as soon as the
analysis is done prior to the update of BURP files which is done
after.

## UnitTests

You can run the unit tests very easily with the command `./test.sh` or `make test` which should give you an output like this.
```
Testing program './midas.monitor.Abs'
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 0 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 11 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
Could only read 22 bytes out of 23 bytes in file 'VAR3D_STATUS.dot'
Continue...
The 'midas.monitor.Abs' tests are successull!!!
```

When a developer want to make changes, he should always follow those steps:
 1. create a test which represents the functionality or the bug to be fixed
 2. run the tests as is and all the tests should pass except the new one
 3. implement the new feature or fix the bug
 4. run the tests again
 5. then choose between these two cases:
   * if the tests pass, go to step 6.
   * if the tests do not pass, repeat the process by going to step 3.
 6. All the tests pass, so you can commit your changes
 7. Push to changes to the central depot
