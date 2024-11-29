# Building MIDAS using make

The present `README` **assumes you are in the `src` directory** (the directory in which is this `README`), meaning onward `./` points to the `src` directory.

## Using midas_build - for most use cases

Although the present build strategy is based on [cmake](https://cmake.org)
which in turn produces makefiles for [GNU make](https://www.gnu.org/software/make/),
we provide a fully automated build wrapper script that should be used in most use cases:
**`midas_build`**.  It builds MIDAS executables using multicore compilation and does some
error checking.

### Configuring the compilation and linking process

The compilation and linking is configured through some environment variables
declared in [`./config.dot.sh`](./src/config.dot.sh), but it is suggested **not to modify directly that file**
since it is part of the versioned repository.
You can either modify and export them in your shell or add them to your profile.

Their default values (in parentheses), **should be good for most users**.

* `MIDAS_COMPILE_DIR_MAIN (unset)`: directory where build directories and the
  executables directory will be.  Each git version will have its build
  directory `${MIDAS_COMPILE_DIR_MAIN}/midas_bld-${VERSION}`, but
  executables **of all versions** will be in
  `${MIDAS_COMPILE_DIR_MAIN}/midas_abs` with the version included in
  the absolute name.  If `${MIDAS_COMPILE_DIR_MAIN}` is equal to
  `build_directory_local_to_the_repository` then this variable is set
  to the `compiledir` directory at the root of the Git repository.  If
  it is not defined, then it will be set to
  `${HOME}/data_maestro/ords/midas-bld/${leaf}` where `${leaf}` is
  `basename` of the Git repository where the code is put.
* `MIDAS_COMPILE_ADD_DEBUG_OPTIONS (no)` : activate the debug flag for
  the compilation if set to `yes`.  You can activate this option by
  using [`--debug` when calling
  `midas_build`](#activate-the-debug-options-for-compilation).  Note
  that enabling debug options **may subtly alter the results and
  therefore cause some unit tests to fail**.  So, to avoid that, we
  set `CHECK_RESULTS_CATCHUP` and `CLEAN_UNITTEST_CATCHUP` to `9` in
  the resources file which prevent the `check` and `clean` tasks to
  fail.  We also increase the memory request for the task
  `/Tests/letkf/glb_15km/UnitTest/run` in debug mode.  Also, make sure
  to **`make clean` before recompiling when you change that variable
  value** if `MIDAS_COMPILE_CLEAN=false`.  Otherwise, some already
  compiled object will keep the impact of the debug options and may
  result in inconsistencies.
* `MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR (true)`: append the version
  identifier to the build directory.  It can be `true` (default) or
  `false`.
* `MIDAS_COMPILE_CODECOVERAGE_DATAPATH` : path to store the code
   coverage diagnostics files.  Same as for the debug options, this
   may subtly alter the results.  So, we avoid running the tasks
   `check` and `clean` for each test.  You can activate this option by
   using `--codecov ${datapath}` when calling `midas_build`.
* `MIDAS_COMPILE_OPTIMIZE_REPORT` : if `yes`, the compiler will
  produce optimization reports.  This will produce files with the
  `.optrpt` extension in the compilation directory with lots of
  information about compiler code optimization.  You can activate this
  option by using `--opt-report` when calling `midas_build`.
* `MIDAS_COMPILE_FRONTEND (ppp5)` : cluster on which to proceed with the compilation
* `MIDAS_COMPILE_CLEAN (true)` : if `true`, remove the build directory after a
  successful installation of the absolutes (if applicable)
* `MIDAS_COMPILE_KEEP_LISTING (false)` : if `false`, remove listings on
  successful compilation (and linking if applicable)
* `MIDAS_COMPILE_HEADNODE_FRONTEND (false)` : if `true`, frontend multicore
  compilation is done directly on headnode, this should only be used on a
  dedicated node obtained through `jobsubi`
  (see [this section](#calling-make-in-parallel) for instructions.)
* `MIDAS_COMPILE_COMPF_GLOBAL ()` : additional user-specified compilation options
* `MIDAS_COMPILE_JOBNAME (midasCompilation)` : name for the job submission
  and the prefix for listings
* `MIDAS_COMPILE_NCORES (8)` : numbers of cores to be used on each machine to
  compile (more than 8 provide no significant improvement).
* `MIDAS_COMPILE_VERBOSE (2)` : verbosity level

Notice that we are **uniformizing the MIDAS compilation environment variable
naming convention**, please consult
[this section](#new-environment-variable-convention) if you used to define
compilation variable in your profile.

Note also that the values of `CHECK_RESULTS_CATCHUP` and
`CLEAN_UNITTEST_CATCHUP` in
`maestro/suites/midas_system_tests/resources/resources.def` are set
according to the values of `MIDAS_COMPILE_ADD_DEBUG_OPTIONS` and
`MIDAS_COMPILE_CODECOVERAGE_DATAPATH` because we want to avoid running
the task `check` and `clean` from the `UnitTest` module since it is
expected that activating any of these features will change the results
of the programs.  We also increase the memory request for the task
`/Tests/letkf/glb_15km/UnitTest/run`.

### Building all

Simply execute **`./midas_build`** from a frontend machine.

Successful compilation, linking and installing will be confirmed with the
display
```
+------------------------------+
|                              |
| MIDAS INSTALLATION COMPLETED |
|   ALL PROGRAMS INSTALLED     |
|                              |
+------------------------------+
```

It will
1. build MIDAS executables using the number of cores specified
2. install (copy) the absolutes in the directory
   `${MIDAS_COMPILE_DIR_MAIN}/midas_abs` using the format
   `midas-${program}_${ORDENV_PLAT}-${VERSION}.Abs`
3. if `MIDAS_COMPILE_CLEAN=true`, will delete the build directory if the
   compilation was successful
   (directory `midas_bld-${VERSION}`)
4. if `MIDAS_COMPILE_KEEP_LISTING=false`, will delete the listings if the
   compilation was successful


Please remember to remove listings (`./midasCompilitation.*`) that are in the
`src` directory.


### Activate the debug options for compilation

If the environment variable `MIDAS_COMPILE_ADD_DEBUG_OPTIONS` is set
to `yes`, the debug options will be enabled at compilation.

The `midas_build` options `--debug`, `-debug` or `-d` will also enable
the debug options avoiding to set the environment variable prior to
compilation.

### New environment variable convention
In the spirit of uniformizing environment variable convention across our
different tools, we decided to change some variable names used in the previous
compilation strategy.
All environment variables now **start** with the prefix `MIDAS_`.
These former variables have been renamed:

* `COMPILEDIR_MIDAS_MAIN` is now `MIDAS_COMPILE_DIR_MAIN`
* `COMPILE_MIDAS_ADD_DEBUG_OPTIONS` is now `MIDAS_COMPILE_ADD_DEBUG_OPTIONS`

If any of those are defined in your profile, they won't be taken into account
and you should remove them.

## `splitobs` an *external* program

The program `splitobs` is built by default with the other programs as described
in the [section Building all](#building-all).  It can also be built as a specific
program in the same manner as described in the
[advanced use cases section](#using-make-advanced-use-cases).

The sources can be found in [`../tools/splitobs`](./tools/splitobs).

## Adding a new program or changing external dependencies

Internal dependencies (the `use` statements in the programs and modules) are
dealt with automatically by `cmake`.

External dependencies however, that are needed at link time by the programs,
need to be specified by the contributor in the file
[`./programs/CMakeLists.txt`](programs/CMakeLists.txt).

So **when a new program is added** or when **external libraries change for an
existing program**, edit the [`./programs/CMakeLists.txt`](programs/CMakeLists.txt)
file and list all external libraries (previously in `programs/programs.mk`)
as prerequisite of the absolute target, such as:
```cmake
add_executable(var var.f90)
target_link_libraries(var PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)
```

### New external dependencies in a module
When new external dependencies are added in a module, it will potentially impact
multiple programs (for which the dependencies will have to be added as described
above).

<!-- TODO: is it still an issue-->

### Addressing circular dependency issues
A circular dependency happens when a module uses another module
which also use the first one (or use a module that use the first one and so
forth).
One can detect these by analysing the `use` statements and the
[dependency tree](#automatic-dependencies).  but it can get tricky.

`make` will detect a circular dependency issues and print a message such as
```
make[1]: Circular earthconstants_mod.o <- mathphysconstants_mod.o dependency dropped.
```
`midas_build` will in turn detect this as well.

Still to help debugging or finding the proper hierarchy of `use` statements, it
can be useful to know the order of compilation that will result from the
dependency tree.
[`analyzeDep.py`](analyzeDep.py) can help you with that; calling it with `-c`
will return, in order, all modules that will need to be compiled *prior* to
compile a given target object (module or program).

For instance, if we want to compile `innovation_mod.o`,
```bash
$ source config.dot.sh
$ ./analyzeDep.py -c innovation
...
Building innovation will result in compiling these objects:
1 : randomnumber
2 : clib_interfaces
3 : utilities
4 : obsFamilyList
5 : mathphysconstants
6 : mpi
...
65 : bmatrixchem
66 : obsOperatorsChem
67 : obsoperators
```
`innovation_mod` depends on all these objects; they will be compiled first and
no `use innovation_mod` can be found in any of those modules - it would
result in a circular dependency error.

Consult `./analyzeDep.py -h` for more information.

--------



## Using make - advanced use cases


If you want to have a more fine grained control, you can call `make` directly,
but for **most users,** `midas_build` **should do just fine**.
I invite you to read its man pages (short and straight to the point).

You'll first need to source the compilation environment (if you want to modify
environment variable values, do it in your profile or through explicit `export`
in the shell - don't modify `config.dot.sh`):
```
$ source ./config.dot.sh
```
Otherwise, only a few targets will be available: `clean`, `cleanabs`,
`cleanall`, `cleaninst`, `cleandep`, `cleanobj`  and `help`.

You can call `make` to build any target; target may be an object file, a
specific program or a label (or *phony* target such as `all`, `objects` or other
label that are not a file *per se*), you can always use autocompletion by
pressing `<TAB>`:
```
$ make <TAB><TAB>
Display all 156 possibilities? (y or n) <y>
absolutes                       install
adjointTest.Abs                 kdtree2_mod.o
adjointTest.o                   lambmatrixhi_mod.o
advection_mod.o                 lamspectraltransform_mod.o
advector.Abs                    letkf.Abs
advector.o                      letkf.o
all                             localizationfunction_mod.o
...
$ make var<TAB>
var.Abs            var.o              varnamelist_mod.o  varqc_mod.o
```
If `<TAB>` does not work (or just show you `help` and `clean`), it is probably
because did not `source ./config.dot.sh`.

When you ask `make` to build a target, it will determine everything that needs
to be done to achieve that goal, for instance if want to build `var.Abs`:
```
$ make var.Abs --dry-run
Preprocessing codeprecision_mod.f90 inplace
Preprocessing clib_interfaces_mod.f90 inplace
Preprocessing rttov_interfaces_mod.f90 inplace
Generating object dependencies > dep.obj.inc
Generating executables dependencies > dep.abs.inc
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c .../midas/src/modules/clib_interfaces_mod.f90 -o clib_interfaces_mod.o > /dev/null
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c .../midas/src/modules/utilities_mod.f90 -o utilities_mod.o > /dev/null
...
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c .../midas/src/programs/var.f90 -o var.o > /dev/null
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4\
    -lf90sqlite ... obsfilter_mod.o -o var.Abs
```
The option `--dry-run` (or `-n`) only prints what it would be doing without actually doing it.


Some frequently used phony targets are:
* `help` : print a short synopsis and important targets
* `all` : compile all programs on current architecture
* `info` : print information about compilation setup
* `install` : install all compiled programs
  Copy them in `${MIDAS_COMPILE_DIR_MAIN}/midas_abs/` and rename them with
  version number:
  `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the
  `../midas.version` script.
* `clean` : remove the build directory for the current version
* `cleanabs` : remove programs in the current build directory
* `cleanobj` : remove objects in the current build directory
* `cleandep` : remove dependencies files in the current build directory
  ([see Automatic dependencies below](#automatic-dependencies))
* `cleanall` : remove **all** build directories
* `cleaninst` : remove **all** installed programs

Omitting the target defaults to `all`.




### The install target
Calling `make install` **after** `make [all| target_1 [target_2 ... ]]` will copy the target absolutes **on the present architecture** to the binaries directory at `${MIDAS_COMPILE_DIR_MAIN}/midas_abs`.  All binaries are copied at the same place with the naming convention `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the `../midas.version` script.


A complete install is then
```
(source ./config.dot.sh && make && make install)
```

Note that when `make` is called to compile and link an absolute (`*.Abs`),
it will first **remove the corresponding installed program** in the `midas_abs`
directory.  This is to address issue #513; always do `make install` after
calling `make` on a program **before** launching the test maestro suite.


### Calling make in parallel

To compile a target using multiple cores use
```
make -j ${MIDAS_COMPILE_NCORES} [-O] [<target>]
```
The `-O` ensure outputs are collected together rather that interspersed with output from other jobs (more in `man make` on that).

Careful to not overload the head node.
One can specify a load average maximum when calling `make` in parallel using `-l`, for instance
```
$ make -j 10 -l 8.5
```
will not start new jobs unless the load average on the node is a below 8.5.
You can also use `ord_soumet` or if you are having a compile-open-house, ask for a interactive party node:
```
$ r.load /fs/ssm/main/opt/jobsubi/jobsubi-0.3
$ jobsubi --show-request  -r ncores=${nCores} ppp4
...
$ cd ${YOUR_MIDAS_PROJECT}/src
$ make -j ${nCores} -O
```
(`nCores=8` is enough.)






### Incremental builds


If you modify a single source file, `make` will determine which targets (objects and absolutes) depend on it and will reprocess only those.
For example, let say you have checked out and compiled all programs, then you work on `modules/varqc_mod.f90`.  When you are done, only `minimization_mod.o`, `var.o` will have to be recompiled and `var.Abs` relinked.  This is done using file time stamps and [automatically generated dependency files](#automatic-dependencies).

What if you just modify the implementation of a function or subroutine without touching its interface?
You should not have to recompile other modules.
But `make` won't know what you did, it will just look at the time stamp and decide to recompile everything that depends on the modified file.  Actually just touching the file, will make `make` thinks it needs to reprocess all dependent files:
```
$ make -n
...
make[2]: Nothing to be done for 'absolutes'
...
$ touch modules/varqc_mod.f90
$ make -n
s.f90 ... -c .../midas/src/modules/varqc_mod.f90 -o varqc_mod.o
s.f90 ... -c .../midas/src/modules/minimization_mod.f90 -o minimization_mod.o
s.f90 ... -c .../midas/src/programs/var.f90 -o var.o
s.f90 ... -o var.Abs
```

One solution to that is to first recompile the modified object (here
`varqc_mod.o`)
```
$ make varqc_mod.o
...
$ make -n | grep '^\-o'
-o minimization_mod.o
-o var.o
-o var.Abs
```
and `touch`  **the other intermediate targets** (using `make --touch`) , making
them newer that the dependencies.
Now because we are compiling out of the source tree, using `make --touch`
directly will work, but would also create empty spurious files (#444)
in the `src` directory.
So you should instead use `make touch-objects`, that will mark **every objects**
as up-to-date.
```
$ make touch-objects
touch minimization_mod.o
touch var.o
$ make -n | grep '^\-o'
-o var.Abs
```
Now we see that only the linking will be done.
```
$ make
...
```
You should be **very careful** with the `touch-%` targets as they override
`make` in determining what should be compiled.



### Automatic dependencies

`make` determine dependencies at build time.
It is thought a better practice to find out about these dependencies automatically.


When `make` is called on a `clean` (or `cleandep`) state, it will determine all objects dependencies using `makedepf90` and include it dynamically in the `Makefile` and then launch the build.

Making only the dependencies:
```
$ make depend
$ tree ../compiledir
../compiledir
└── midas_bld-v_3.5.2-133-g111551f_M
    └── ubuntu-18.04-skylake-64
        └── intel-19.0.3.199
            ├── clib_interfaces_mod.f90
            ├── codeprecision_mod.f90
            ├── dep.abs.inc
            ├── dep.obj.inc
            └── rttov_interfaces_mod.f90
```
And this is what it looks like inside of `dep.obj.inc`:
```make
obstimeinterp_mod.o : obstimeinterp_mod.f90 obsspacedata_mod.o timecoord_mod.o utilities_mod.o mpivar_mod.o mpi_mod.o
thinning_mod.o : thinning_mod.f90 utilities_mod.o obsspacedata_mod.o bufr_mod.o codeprecision_mod.o
increment_mod.o : increment_mod.f90 varnamelist_mod.o bmatrix_mod.o gridVariableTransforms_mod.o utilities_mod.o humiditylimits_mod.o verticalcoord_mod.o horizontalcoord_mod.o gridstatevector_mod.o timecoord_mod.o mpivar_mod.o mpi_mod.o
...
bgckAtms.o : bgckAtms.f90 bgckmicrowave_mod.o
```
`dep.obj.inc` are generated by `makedepf90` and produce *superficial* dependencies, those needed to compile **objects**, but **not enough to link absolutes**.

To proceed with linking we need the fully recursive dependencies.
We parse (with a simple python script, [`recursiveDep.py`](./recursiveDep.py)) the superficial dependency file and deduce the fully recursive ones needed at link time; they are in `dep.abs.inc` files.


## SSM packaging

To publish the absolutes in a SSM domain, one have to
1. make sure to keep the build directory by exporting
   `MIDAS_COMPILE_CLEAN=false` in your shell and build:
   ```
   (export MIDAS_COMPILE_CLEAN=false ; ./midas_build)
   ```
2. update `MIDAS_SSM_*` variables in `./config.dot.sh` or export them in the
   shell (making sure you have write privilege to `${MIDAS_SSM_TARGET}`)
3. for **each architecture**
   ```
   (source ./config.dot.sh && make ssm)
   ```
4. protect the domain
   ```
   (source ./config.dot.sh && make ssm_protect)
   ```


## What is left to do


* automated `doc` building and `diagrams`, etc.

But most of all... taking into account your input.
Don't hesitate to contact-me for your input or for some guidance: @mad001


Many thanks to @phb001 and [mad-scientist](http://make.mad-scientist.net/)
