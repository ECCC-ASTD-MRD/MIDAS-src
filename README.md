## `midas.findTrials`

This scripts finds the trials file in an assimilation window.

Here is the `help` page:
```
./midas.findTrials -h
usage: midas.findTrials [-h] [--trialfrequency TRIALFREQUENCY]
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

### UnitTests

You can run the unit tests very easily with the command `./midas.findTrials -u` which should give you an output like this.
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
