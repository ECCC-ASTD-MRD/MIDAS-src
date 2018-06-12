## `envar.findTrials`

This scripts finds the trials file in an assimilation window.

Here is the `help` page:
```
./envar.findTrials -h
usage: envar.findTrials [-h] [--trialfrequency TRIALFREQUENCY]
                        [--trialoutputfrequency TRIALOUTPUTFREQUENCY]
                        [--assimilationwindowwidth WIDTH] [--date DATE]
                        [--unittest] [--verbose] [--version]

plot timings extracted by 'getTimings.py'

optional arguments:
  -h, --help            show this help message and exit
  --trialfrequency TRIALFREQUENCY, -t TRIALFREQUENCY
                        The frequency of trials in hours (may be fractional)
  --trialoutputfrequency TRIALOUTPUTFREQUENCY, -o TRIALOUTPUTFREQUENCY
                        The frequency of trials in minutes (must be an
                        integer, default is 15 minutes)
  --assimilationwindowwidth WIDTH, -w WIDTH
                        Width of the assimilation window in hours (may be
                        fractional)
  --date DATE, -d DATE  Date of the analysis in format 'YYYYMMDDHH'
  --unittest, -u        Ignore any other arguments and run the UnitTests
  --verbose, -v         Explain what is being done
  --version, -V         Output version information and exit
```

### Examples

For example, when calling the command `./envar.findTrials -t 6 -o 15 -w 6 -d 2018050212`, we obtain the output:
```
2018050206_180m
2018050206_195m
2018050206_210m
2018050206_225m
2018050206_240m
2018050206_255m
2018050206_270m
2018050206_285m
2018050206_300m
2018050206_315m
2018050206_330m
2018050206_345m
2018050206_360m
2018050206_375m
2018050206_390m
2018050206_405m
2018050206_420m
2018050206_435m
2018050206_450m
2018050206_465m
2018050206_480m
2018050206_495m
2018050206_510m
2018050206_525m
2018050206_540m
```

For one hour assimilation window with trials at each 6 hours for which
output each 15 minutes, we call `./envar.findTrials -t 6 -o 15 -w 1 -d 2018050210`
which gives:
```
2018050206_210m
2018050206_225m
2018050206_240m
2018050206_255m
2018050206_270m
```

## UnitTests

You can run the unit tests very easily with the command `./envar.findTrials -u` which should give you an output like this.
```
test_1hr (__main__.Test_findTrials) ... ok
test_1hr_1bin (__main__.Test_findTrials) ... ok
test_6hrs (__main__.Test_findTrials) ... ok
test_6hrs_15m (__main__.Test_findTrials) ... ok
test_error (__main__.Test_findTrials) ... ok
test1 (__main__.Test_trialDate) ... ok
test2 (__main__.Test_trialDate) ... ok
test3 (__main__.Test_trialDate) ... ok
test4 (__main__.Test_trialDate) ... ok

----------------------------------------------------------------------
Ran 9 tests in 0.004s

OK
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
