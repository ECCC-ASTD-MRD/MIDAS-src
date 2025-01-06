## Timings

The script [`tools/timingTool/midas.timingTool`](midas.timingTool)
allows to extract timing information from MIDAS task listings.  It parses a
listing file and looks for the "`TMG:`" where timing information is printed out.

The output will look like this:
```
LABEL                                     minTime  maxTime meanTime   minCount   maxCount  meanCount
=====                                     =======  ======= ========   ========   ========  =========
Main....................................   409.73   412.30   409.97          1          1          1
--ReadEnsemble..........................    69.30    69.33    69.31          1          1          1
--WriteEnsMeanRms.......................    40.75    59.77    43.23          4          4          4
--StateToColumn.........................   121.00   124.12   121.91         66         66         66
----s2c_NL..............................   116.51   121.97   119.47         68         68         68
--LETKFanalysis.........................    80.22    80.32    80.26          1          1          1
----CommWeights.........................    20.26    63.18    26.24        340       6370       5411
----CommWeights-signals.................     0.13    53.36     0.98         85       3100       2621
low-level--gio_writeToFile..............     5.10    68.83    17.05          7         10          8
low-level--gsv_tilesToVarsLevs_r4.......    73.22    82.41    79.64         66         66         66
low-level--gio_writeToFile-gather.......     3.61    55.97    12.02       2053       2053       2053
```
where the call hierarchy is made explicit.

The script can be called like this to simply extract the timings:
```sh
./midas.timingTool ${LISTING}
```

### Comparing timings

It can also be useful to use it to compare two runs.  To do so, it is important
to make sure timing lists extracted from both listings contain the same elements
in the same order for the comparison to be useful. To achieve this, one will
call `midas.timingTool` on each listing, but using the other listing as _reference_
and specifying the order for the timing outputs be ordered properly for the
comparison:
```sh
./midas.timingTool ${LISTING_1} -r ${LISTING_2} -o 1 > timings_1.dat
./midas.timingTool ${LISTING_2} -r ${LISTING_1} -o 2 > timings_2.dat
```
then one can use their prefered `diff` tool to compare both run:
```
diff -y -W 210 timings_1.dat timings_2.dat
```
`xxdiff` is GUI alternative to this command line.

### Filtering timings

`midas.timingTool` can be called with the argument `-h` or `--help` to get some
help on it's calling interface.
There is the option to specify the minimal timing threshold for an element be
retained in the list:
```sh
./midas.timingTool ${LISTING} -l 50.0
```
for instance, will only list elements that took 50 seconds or more to run.
