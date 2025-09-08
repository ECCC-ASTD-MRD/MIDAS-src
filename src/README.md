# Building MIDAS using make

The present `README` **assumes you are in the `src` directory** (the directory
in which is this `README`), meaning onward `./` points to the `src` directory.
When it's not the case, it will be pointed out explicitely.

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
* `MIDAS_COMPILE_CLEAN (true)` : if `true`, clean the build directory after a
  successful installation of the absolutes (if applicable)
* `MIDAS_COMPILE_KEEP_LISTING (false)` : if `false`, remove listings on
  successful compilation (and linking if applicable)
* `MIDAS_COMPILE_HEADNODE_FRONTEND (false)` : if `true`, frontend
  multicore compilation is done directly on headnode (activated when
  `--local` is specified to `midas_build`)
  (see [this section](#calling-make-in-parallel) for instructions.)
* `MIDAS_COMPILE_COMPF_GLOBAL ()` : additional user-specified compilation options
* `MIDAS_COMPILE_JOBNAME (midasCompilation)` : name for the job submission
  and the prefix for listings
* `MIDAS_COMPILE_NCORES (8)` : numbers of cores to be used on each machine to
  compile (more than 8 provide no significant improvement).
* `MIDAS_COMPILE_VERBOSE` : Defaults to `TRUE` when using `midas_build` otherwise it is `FALSE`.


Note that the values of `CHECK_RESULTS_CATCHUP` and
`CLEAN_UNITTEST_CATCHUP` in
`maestro/suites/midas_system_tests/resources/resources.def` are set
according to the values of `MIDAS_COMPILE_ADD_DEBUG_OPTIONS` and
`MIDAS_COMPILE_CODECOVERAGE_DATAPATH` because we want to avoid running
the task `check` and `clean` from the `UnitTest` module since it is
expected that activating any of these features will change the results
of the programs.  We also increase the memory request for the task
`/Tests/letkf/glb_15km/UnitTest/run`.

#### Environment variable convention

In the spirit of uniformizing environment variable convention across our
different tools, all environment variables **start** with the prefix `MIDAS_`.

These former variables have been renamed:
* `COMPILEDIR_MIDAS_MAIN` is now `MIDAS_COMPILE_DIR_MAIN`
* `COMPILE_MIDAS_ADD_DEBUG_OPTIONS` is now `MIDAS_COMPILE_ADD_DEBUG_OPTIONS`

If any of those are still defined in your profile, they won't be taken into account
and you should remove them.


### Building all

Simply execute **`./midas_build`** from a frontend machine.

Successful compilation, linking and installing will be confirmed with the display
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
3. if `MIDAS_COMPILE_CLEAN=true`, will clean the build directory if the
   compilation was successful (but keep installed binaries)
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

## `splitobs` an *external* program

The program `splitobs` is built by default with the other programs as described
in the [section Building all](#building-all).  It can also be built as a specific
program in the same manner as described in the
[advanced use cases section](#using-make-advanced-use-cases).

The sources can be found in [`../tools/splitobs`](./tools/splitobs).

## Adding a new program or changing external dependencies

Internal dependencies (defined through the `use` statements in the programs and modules) are
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
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module random hpcoperf)
```
and add the program to the list `PROJECT_PROGRAMS` in that same file
[`./programs/CMakeLists.txt`](programs/CMakeLists.txt).

### New external dependencies in a module

When new external dependencies are added in a module, it will potentially impact
multiple programs (for which the dependencies will have to be added as described
above).  One can analyze the `use` statements in the source code or `grep` the
files `depend.make` in the `programs/CMakeFiles/` sub-directories of the build
directory after a build (or an attempt at building) the code.  These files are
part of [the `cmake` build directory strcuture](#what-does-cmake-do).


### Addressing circular dependency issues

A circular dependency happens when a module uses another module
which also uses the first one (or use a module that use the first one and so
forth).

`make` will detect a circular dependency issues and print a message such as
```
make[1]: Circular earthconstants_mod.o <- mathphysconstants_mod.o dependency dropped.
```
`midas_build` will in turn detect this as well.

To fix this situation, the code will potentially need to be modified and `use` statements
modified such that no two modules uses one another.
It also possible to analyze the dependencies of any module or program by looking
at the files `depend.make` in the corresponding build directory after a build
(or an attempt at building) the code.  These files are part of
[the `cmake` build directory strcuture](#what-does-cmake-do).


--------

## Why do we now use `cmake`?
We had a build strategy, efficient and specifically designed to build MIDAS, so
why do we need to change it to use `cmake`?  The motivation is twofold:
standardization and portability.

RPN-SI is trying to unify build methods such that all major components of our
operational products be built using the same standard.  This obviously makes it
easier to support, but it also allows to integrate all products in a global build
chained process through nested `cmake` configurations.  `cmake` is also an
universal industrial standard, such that using it allow for more portable code.
There are efforts toward a better portability and distribution of MIDAS (such as
on [github.com](https://github.com/ECCC-ASTD-MRD/MIDAS-src)) in order that MIDAS
be properly referenced for scientific publications and for academic and other
partners access.  There is still work to be done to acheive real portability, but
using a standardized build system is an important stepping stone; we previously used many
in-house scripts - we still do (for instance SSM packages), but diminishing
their number is necessary to evolve toward this goal.


## Using make - advanced use cases


If you want to have a more fine grained control, you can call `cmake` and `make`
directly, but for **most users,** `midas_build` **should do just fine**. I invite
you to read `make` man pages (short and straight to the point).

Onward, we'll note the root of your git project `${topLevel}` and your build
directory `${buildDir}`.  Note that `${topLevel}/src/config.dot.sh`, configuring
the compilation environment, will create a build directory at `${MIDAS_COMPILE_DIR_MAIN}`
and link it to `${topLevel}/compildir`. You could in fact launch `cmake` from any
directory of your choice to define it as a build directory.

The development workflow will be:
1. export some specific [`MIDAS_` environment variables](#configuring-the-compilation-and-linking-process)
   or define them in your profile
2. `. ${topLevel}/src/config.dot.sh`: source the compilation environment and
   create the default build directory from `${MIDAS_COMPILE_DIR_MAIN}`
3. `cd ${buildDir}` : the default build directory is linked at ` ${topLevel}/compiledir`
4. `cmake ${MIDAS_SOURCE_DIR}` : [initialize and configure the build directory](#what-does-cmake-do)
5. `make -j [${MIDAS_COMPILE_NCORES}]` : build everything using [multiple cores](#calling-make-in-parallel)
6. iterative development loop which we cover in detail in the [section on incremental builds](#incremental-builds)

One can add `VERBOSE=1` to the `make` command to have the exact
compilation command done.

### What does `cmake` do?
To better understand the process, it is relevant to understand some elements of
the `cmake` compilation strategy.  **`cmake` does not do the build**, it
prepares and organises the build _recipes_,  a bunch of `Makefile` for each
parts to be build.

It is meant to be launched from the desired destination for the binaries
(objects, modules and program binairies), the build directory (which could be
anywhere and does not need to be linked to the source directory).
`cmake` will create a directory for every real (non phony) target with `make`
recipes in `.make` files.  Let's look at `${buildDir}/src/` directory structure
(up to 3 depth levels):
```
├── modules
│   ├── CMakeFiles
│   │   ├── CMakeDirectoryInformation.cmake
│   │   ├── midas.dir
│   │   └── progress.marks
│   ├── cmake_install.cmake
│   ├── CTestTestfile.cmake
│   ├── Makefile
│   └── modules
└── programs
    ├── CMakeFiles
    │   ├── adjointTest.dir
    │   ├── analysisErrorOI.dir
   ... ...
    │   └── var.dir
    ├── cmake_install.cmake
    ├── CTestTestfile.cmake
    └── Makefile
```
As can be seen, every programs got his own directory in `programs/CMakeFiles/` each containing
various files and in particular a `build.make` file containing the main `cmake`-generated
recipes (not meant to be modified).  One can see all the prerequisites needed to
build the absolute binary.  One can appreciate that no MIDAS module object is
present as a prerequisite.  That is because they are archived in a static
library, `src/modules/libmidas.a`, listed in the prerequisites.

The modules structure is different, there is a single `midas.dir`.  This
directory is the one containing the recipe to build `libmidas.a` from all the
module objects.

This is different from our previous `make` build strategy where each object was
built and linked individually to programs and, has explained in the [incremental
builds section](#incremental-builds), it has somewhat important consequences.

So in summary, calling `cmake` configure the different `Makefile` and calling
`make` from `${buildDir}` will:
1. build all module objects
2. link them together in `libmidas.a`
3. build program objects
4. link each program with `libmidas.a`

Worth mentionning, for each program and for the static library, there is a
`depend.make` file generated at compile time and located in
`${buildDir}/src/programs/CMakeFiles/${programBaseName}.dir/depend.make` for
each program and at `${buildDir}/src/modules/CMakeFiles/midas.dir/depend.make` for
the static library and containing the dependency tree for all modules in a
single file.  Analyzing these files is another way to determine the dependency
structure.


### Incremental builds
Once you set up the build environment and got a fresh build (presumably from a
clean git working tree) you are ready to make and test some modifications.  This
section assumes you are in `${buildDir}`.

The point of using the method described here is to save a **couple of tens of
seconds** per build, so if you plan on making one or two modifications and you
don't intend many iterations, you should definitely stick to the
[integrated `midas_build` method](#building-all).  If your objective is to make
a modification, build, test, make adjustments and iterate many times like this
and you are comfortable with `make`, this section is for you.

The compromise you do to gain these seconds, is to do `cmake` dependancy
analysis on your own and evaluate what need to be recompiled; error in this
analysis **could result in compilation error, unexpected behavior, faulty results
or silent errors**.  Also, when the iterative incremental development is completed,
make sure to **do a final build using `midas_build`**.  This is a requirement before
any merge request, but also a good practice to ensure you did not oversee some dependancy
during the incremental process.

Now that warnings have been stated - here we go!
In any case, you first have to go through the first build steps
[described previously](#using-make-advanced-use-cases).

#### Incremental build for programs
If your modifications are exclusively done in programs, then the process is really
simple since you are modifying only the leaves of the dependancy tree.  You can
then simply follow these steps:
1. modify the program(s)
2. adjust, if necessary, the [external dependancies](#adding-a-new-program-or-changing-external-dependencies)
3. `make [-j] [ ${programBaseName} [ ${anotherProgram} ...]]` : rebuild and link all modified programs.  Alternatively,
   you can specify which program(s) to build by passing them as arguments to `make`.
   In such case, `${programBaseName}` is a program source file name without the
   extension. Note that autocompletion will suggest you available targets by pressing `<TAB>`.
4. test, modify and reiterate again as long as no modules are modified, if a module is
   modified, you will have to follow the [incremental build method for modules](#incremental-build-for-modules)


#### Incremental build for modules
Remember from [previous section on `cmake`](#what-does-cmake-do), that modules are
compiled and archived together in a static library, `libmidas.a`, which is in turn
prerequisite to all programs.  This is important, because it means that **any modification**
to a module will trigger the reconstruction of the static library, which in turn
will inform `make` that all programs will have to be recompiled and linked since
it is a prerequisite of each of them and it has been modified.

When a module is modified, every module that `use` it will have to be recompiled
as well as all modules down the dependency tree.  One can consult the
`${buildDir}/src/modules/CMakeFiles/midas.dir/depend.make`
(see [the `cmake` build directory strcuture](#what-does-cmake-do)) or analyze the
 `use` statements in the module source files to determine which are those
 modules.  Thus, if you simply call `make -j midas`, all those modules will be
 recompiled.

If your modification __does not impact any subroutine interface__, then you know
depending modules need not be recompiled.  Yet `make` does not know this,
because it uses the last modification timestamp and some very shallow code analysis
of the change made; if you actually only `touch` a source file, it will recompile
it, but not its dependent modules.  If, however, you make an actual modification
**in a subroutine** (not only in the module header for instance), then it will recompile
all dependent modules even if no interface have been modified.  In such circumstances
and to save a few seconds, you could trick `make` in not recompiling the depending
modules.  To do so you'll need to call a somewhat lengthy `make` incantation from
`${buildDir}`:
```sh
make -f src/modules/CMakeFiles/midas.dir/build.make src/modules/CMakeFiles/midas.dir/${modifiedModule}_mod.f90.o
make -f src/modules/CMakeFiles/midas.dir/build.make --touch src/modules/CMakeFiles/midas.dir/*.o
make -j midas
```
What is going on here?
First let's observe the similitude between the first two lines; in each case
`make` is told to use a specific makefile through the argument `-f`,
`src/modules/CMakeFiles/midas.dir/build.make`, which contains the recipes for
`libmidas.a`, including all module objects.  The first line tells `make` to build
`${modifiedModule}_mod.f90.o` given the recipes for the library.
The second line is the trick _per se_, it tells `make` to `touch` all the
elements in the library, in effect brigning their timestamp newer than the
recompiled module such that when, in the third command, the library is rebuilt,
only linking will be done.


The static library is also rerequisite to all programs, a simple call to `make` will
recompile all the programs.  It may or may not be what you need depending on the tests you
will be conducting.  If you intend to test only a single program, then you should
relink the static library and compile only the specific program.  Such workflow
will then be:
1. modify module(s) and possibly program(s)
2. adjust, if necessary, the [external dependancies](#adding-a-new-program-or-changing-external-dependencies)
3. follow the preceding steps to rebuild the MIDAS static library
4. `make ${programBaseName}` : rebuild a specific program.


#### Is this worth your time?

* Building all carelessly from scratch but in parallel on the headnode using `make -j`
  is around 120 seconds.
* Building and linking the library is less than a second if a single module has
  been modified.
* Building a single program is around 2 seconds.
* If only a couple of modules need to be recompiled, building only the shared library
  is less than 20 seconds in most cases (depending on how deep is the modification).

So if you make a modification, build and test once or twice, a simple
`make -j [ ${programBaseName} ]` is what you need.  If you intend on doing this
development iteratively many times with short tests in between, then maybe the
seconds saving is worth the more complex incantation.

In any case, doing the build in parallel on the frontend using `make` instead of
`midas_build` is a certainly a gain in time.

A complete compilation from `${topLevel}` is then
```sh
## This command will run 'cmake' and change directory to build directory
. ./src/config.dot.sh
make -j all
```

### The `install` target

The `install` target is quite a standard among compilation process.
It does exist in this project but the programs generated by this
command are not used anywhere in the development flow.

The MIDAS project is not using it because it has several drawbacks:
 * The name of the files generated must be decided when `cmake` is
   called.  In MIDAS context, it is not possible because the version
   can evolve quite a bit during the development process.  So to have
   the expected names, `cmake` would need to be called each time and
   it would take to much time.
 * The installing process manipulates the binaries by adding `rpath`
   so the installed programs are not the ones that are tested.
   Preferring to keep the binaries that have been tested, the
   installed programs will never been used.

### The `prepare_test` target

A special target is available to update the testing environment which
creates the files
`maestro/suites/midas_system_tests/resources/resources.def` and
`maestro/suites/midas_system_tests/abs.dot` that are necessary for the
testing `maestro` suite.

### Calling make in parallel

To compile a target using multiple cores use
```
make -j [${MIDAS_COMPILE_NCORES}] [-O] [<target>]
```
The `-O` ensure outputs are collected together rather that interspersed with output from other jobs (more in `man make` on that).

Careful to not overload the head node.
One can specify a load average maximum when calling `make` in parallel using `-l`, for instance
```
$ make -j 8 -l 10
```
will not start new jobs unless the load average on the node is a below 10.
Also, 8 cores is enough to build and exploit parallelization to its maximum.
You can also leave `-j` without a value to use as many cores as needed.

--------

Many thanks to @cpi001 and @phb001 for their contributions!
