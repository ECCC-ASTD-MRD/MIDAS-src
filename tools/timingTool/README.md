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


## Calling synopsis
`./midas.timingTool -h` will list all calling options.
The script can be called like this to simply extract the timings:
```sh
./midas.timingTool ${LISTING}
```

One can feed a listing directly from standard input:
```sh
cat ${LISTING} | ./midas.timingTool
```
Note that this allows to use it in conjunction with `nodelister` to get the timings directly from a suite:
```sh
SEQ_EXP_HOME=${path_to_suite} nodelister -n ${maestro_node} -d ${datetime} | ./midas.timingTool
```
It also allows to process a gzipped listing directly,
```sh
./midas.timingTool ${COMPRESSED_LISTING}
```
however, it cannot process a gzipped stream (such as `cat ${COMPRESSED_LISTING} | ./midas.timingTool`).  You could however do
```sh
zcat ${COMPRESSED_LISTING} | ./midas.timingTool
```


### Comparing timings

It can also be useful to use it to compare two runs.  To do so, it is important
to make sure timing lists extracted from both listings contain the same elements
in the same order for the comparison to be useful. To achieve this, one will
call `midas.timingTool` on each listing, but using the other listing as _reference_
for the timing outputs be ordered properly for the comparison:
```sh
./midas.timingTool ${LISTING_1} -d -r ${LISTING_2}
```
Note, that it is also possible to use `nodelister` for the référence using process
substitution (`<( ${command} )`), for example to compare two listings directly from two suites, you
could use either of these commands:
```sh
SEQ_EXP_HOME=${path_to_suite1} nodelister -n ${maestro_node1} -d ${datetime1} | ./midas.timingTool -d -r <(SEQ_EXP_HOME=${path_to_suite2} nodelister -n ${maestro_node2} -d ${datetime2})
```
```sh
./midas.timingTool <(SEQ_EXP_HOME=${path_to_suite1} nodelister -n ${maestro_node1} -d ${datetime1}) -d -r <(SEQ_EXP_HOME=${path_to_suite2} nodelister -n ${maestro_node2} -d ${datetime2})
```
The reference listing can also be gzipped.

Another example is to use `cmcarc` to extract a listing from an archive like this:
```bash
archive=~sanl888/data/ppp5/UnitTests/midas/listings/letkf/glb_15km/listings_v_4.1.0-b1-7-g88dba0fad.ca
archive_reference=~sanl888/data/ppp5/UnitTests/midas/listings/letkf/glb_15km/listings_v_4.1.0-b1-6-ge9675fe6c.ca

./midas.timingTool -l 0 <(cmcarc -f $archive -x Tests.letkf.glb_15km.UnitTest.run -O)           \
                     -r <(cmcarc -f $archive_reference -x Tests.letkf.glb_15km.UnitTest.run -O) \
                     -d --title with-blas -2 no-blas
```


To maintain retro-compatibility, it is also possible to do it manually:
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

### Unit Tests
The script contains a self testing feature; simply call it with `-u`:
```sh
./midas.timingTool -u
```
