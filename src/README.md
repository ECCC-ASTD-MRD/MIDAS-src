# Building MIDAS using make

The present `README` **assumes you are in the `src` directory** (the directory in which is this `README`), meaning onward `./` points to the `src` directory.

The former compilation solution is still functionning, although there is a 
change in the [environment variable naming convention](#new-environment-variable-convention).


## Using midas_build - for most use cases

Although the present build strategy is based on [GNU make](https://www.gnu.org/software/make/), we provide a fully automated build wrapper script that should be used in most use cases: **`midas_build`**.  It builds MIDAS executables on both front and backend using multicore compilation and does some error checking.

### Configuring the compilation and linking process

The compilation and linking is configured through some environment variables.
You can either modify and export them in your shell or add them to your profile.
Here are these variables.

Their default values (in parentheses), **should be good for most users**.

* `MIDAS_COMPILE_BACKEND (daley)` and `MIDAS_COMPILE_FRONTEND (eccc-ppp4)`: 
  machines used for each architecture
* `MIDAS_COMPILE_JOBNAME (midasCompilation)`: name for the job submission
  and the prefix for listings
* `MIDAS_COMPILE_DIR_MAIN (${HOME}/data_maestro/ords/midas-bld)`: 
  directory where build directories and the executables directory will be.
  Each git version will have its build directory 
  `${MIDAS_COMPILE_DIR_MAIN}/midas_bld-${VERSION}`, but executables 
  **of all versions** will be in `${MIDAS_COMPILE_DIR_MAIN}/midas_abs` with the
  version included in the absolute name.
* `MIDAS_COMPILE_NCORES (8)`: numbers of cores to be used on each machine to 
  compile (more than 8 provide no significant improvement).
* `MIDAS_COMPILE_VERBOSE (2)`: verbosity level
* `MIDAS_COMPILE_CLEAN (true)`: if `true`, remove the build directory after the
  installation of the absolutes
* `MIDAS_COMPILE_HEADNODE_FRONTEND (false)`: if `true`, frontend multicore 
  compilation is done directly on headnode, this should only be used on a 
  dedicated node obtained through `jobsubi` 
  (see [this section](#calling-make-in-parallel) for instructions.)
* `MIDAS_COMPILE_KEEP_LISTING (false)`: if `false`, remove listings on 
  successful compilation (and linking if applicable)
* `MIDAS_COMPILE_ADD_DEBUG_OPTIONS (no)` : add debug options to compiler if
  set to `yes`.

These variables are declared in `./config.dot.sh`, but it is suggested **not to 
modify directly that file** since it is part of the versioned repository.

Notice that we are **uniformizing the MIDAS compilation environment variable 
naming convention**, please consult 
[this section](#new-environment-variable-convention) if you used to define 
compilation variable in your profile.  


### Building all

Simply execute **`./midas_build`** from a frontend machine.

Successful compilation, linking and installing will be confirmed with the 
display
```
╔═════════════════════════════╗
║                             ║
║ MIDAS COMPILATION COMPLETED ║
║                             ║
╚═════════════════════════════╝
```

It will 
1. build MIDAS executables on both architectures using the number of cores 
   specified 
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

### Using midas_build for specific targets
`midas_build` is a wrapper around `make`; it defaults to compiling, 
linking and installing all the absolutes on both architectures, but it can also
be used to build specific targets by passing it as arguments:
```
$ ./midas_build obsSelection.Abs var.Abs
```

A *target* is something (often a file) to build; you can get information on
available targets by calling `make help`.
See [this section](#using-make-advanced-use-cases) for more on targets.


### Auto-completion

`midas_build` comes with a bash auto-completion feature, such that argument 
passing can be auto-completed by pressing `<TAB>`:
```
$ midas_build <TAB><TAB>
Display all 143 possibilities? (y or n)
$ midas_build obsImpact.<TAB><TAB>
obsImpact.Abs  obsImpact.o
```

However, this feature needs to be installed:
```
./install_build_completion.sh
...
Auto completion for midas_build installed

To use it directly (in the present shell):
   `source /home/${USER}/.profile.d/interactive/post`
   `source /home/${USER}/.bash_completion`
(in any case, it will be automatically loaded on next shells)

```
This will create a file `~/.bash_completion` and a directory 
`~/.bash_completion.d` in your home, and append the current directory to your
`${PATH}` by adding a line to your `~/.profile.d/interactive/post`.
For it to be functionnal in the present shell, you'll have to source the two
files (it will be automatic in future shells).

### New environment variable convention
In the spirit of uniformizing environment variable convention across our 
different tools, we decided to change some variable names used in the previous 
compilation strategy.
All environment variables now **start** with the prefix `MIDAS_`.
These former variables have been renamed:

* `COMPILEDIR_MIDAS_MAIN` is now `MIDAS_COMPILE_DIR_MAIN`
* `COMPILE_MIDAS_ADD_DEBUG_OPTIONS` is now `MIDAS_COMPILE_ADD_DEBUG_OPTIONS`

If any of those are defined in your profile, you will see a corresponding 
warning and should change them to respect this new convention in order to
obtain the expected result.


## Adding a new program or changing external dependencies

In the previous solution, when a new program is added, two files needed to be changed in `./programs/src_files/`:
* `src_files_${PGM}.sh`
* `compile_setup_${PGM}.sh`

The first contains dependencies information and this is dealt with automatically ([see this section](#automatic-dependencies)).
The second one contain external libraries that are needed at link time by the programs.
The information contained in **all** `compile_setup_*.sh` files is now found in [`./programs/programs.mk`](programs/programs.mk).

So **when a new program is added** or when **external libraries change for an existing program**, edit the [`./programs/programs.mk`](programs/programs.mk) file and list all external libraries (previously in `compile_setup_${PGM}.sh`) as prerequisite of the absolute target, such as:
```
var.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
         rttov_parallel rttov_main rttov_emis_atlas rttov_other\
         $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random
 ```



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
`cleanall`, `cleandep`, `cleanobj`  and `help`.

There is however [a pending bug (#453)](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/453) that is **triggered by sourcing 
`src/programs/commons/compile_setup.sh`** (sourced in `config.dot.sh`).
After the config sourcing, the shell becomes unstable with respect to some 
command and/or auto-completion features.

You can call make to build any target.

A target may be an object file, a specific program or a label (or *phony* 
target such as `all`, `clean` or other label that are not a file *per se*), you
can always use autocompletion by pressing `<TAB>`:
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
Auto-completion does not work on the backends.

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
  `../midas.version.sh` script.
* `clean` : remove the build directory for the current version 
* `cleanabs` : remove programs in the current build directory
* `cleanobj` : remove objects in the current build directory
* `cleandep` : remove dependencies files in the current build directory
  ([see Automatic dependencies below](#automatic-dependencies))
* `cleanall` : remove **all** build directories

Omitting the target defaults to `all`.




### The install target
Calling `make install` **after** `make [all]` will copy the absolute **on the present architecture** to the binaries directory at `${MIDAS_COMPILE_DIR_MAIN}/midas_abs`.  All binaries are copied at the same place with the naming convention `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the `../midas.version.sh` script.


A complete install is then 
```
(source ./config.dot.sh && make && make install)
```
launched from all platforms (what is done by `./midas_build` without argument).


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
`varqc_mod.o`) and `touch`  **the other intermediate targets** (using
`make --touch` or `-t`) , making them newer that the dependencies:
```
$ make varqc_mod.o
...
$ make -n | grep '^\-o'
-o minimization_mod.o
-o var.o
-o var.Abs
$ make --touch minimization_mod.o var.o
touch minimization_mod.o
touch var.o
touch minimization_mod.o
touch var.o
$ make -n | grep '^\-o'
-o var.Abs
```
(one can also use the phony target `objects` refering to all `.o` files.)
Now we see that only the linking will be done.
```
$ make
...
```
This has an anoying drawback, it creates empty files (here `minimization_mod.o` and  `var.o`, the `make --touch` targets) in the `src` directory; they can be deleted or ignored (it is related to the [out-of-tree compilation](#out-of-tree-compilation)).
It is an issue (#444) we are aware of.




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
   (export MIDAS_COMPILE_CLEAN=false ; midas_build)
   ```
2. update `MIDAS_SSM_*` variables in `./config.dot.sh` or export them in the 
   shell (making sure you have write privilege to `${MIDAS_SSM_TARGET}`)
3. for **each architecture**
   ```
   (source ./config.dot.sh && make ssm)
   ```
4. **once all architectures have been published**, protect the domain
   (only need to be done once, from either front or backend):
   ```
   (source ./config.dot.sh && make ssm_protect)
   ```


## What is left to do


* understand the shell instability triggered by the sourcing of 
  `programs/commons/compile_setup.sh` (#453) that hinders the direct use of 
  `make`
* address the `make --touch` spurious empty file bug (#444)
* automated `doc` building and `diagrams`, etc.

But most of all... taking into account your input.  
Don't hesitate to contact-me for your input or for some guidance: @mad001


Many thanks to @phb001 and [mad-scientist](http://make.mad-scientist.net/)
