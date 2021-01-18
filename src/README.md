# Building MIDAS using make

The present `README` **assumes you are in the `src` directory** (the directory in which is this `README`), meaning onward `./` points to the `src` directory.


## I just want to build

From `PPP[34]`, you can edit `./config.dot.sh` (for instance to modify `MIDAS_COMPILE_BACKEND`, `MIDAS_COMPILE_FRONTEND` or `MIDAS_COMPILE_DIR_MAIN`).
Notice that we are **uniformizing the MIDAS compilation environment variable naming convention**, please consult [this section](#new-environment-variable-convention) if you used to define compilation variable in your profile.
then
```
$ ./midas_build
...
#####################################
... Launching direct compilation on daley
    > listing_daley
######################################
... Launching compilation on eccc-ppp3
...
=================== ord_soumet version 1.27 =================
...
######################################
#
#  MIDAS COMPILATION COMPLETED
#
#  > /home/mad001/data_maestro/ords/midas-bld
#
######################################
```
It will build MIDAS executables on both `PPP4` (or `PPP3`) and `Daley` (or 
`Banting`) submitting jobs with the number of cores specified in 
`config.dot.sh` (8 seems to be optimal).
It will then install (copy) the absolutes in the directory 
`${MIDAS_COMPILE_DIR_MAIN}/midas-bld` using the format `midas-$(program)_$(ORDENV_PLAT)-$(VERSION).Abs`. 

## Using `build-midas` for specific targets
`build-midas` is a wrapper around 
[GNU `make`](https://www.gnu.org/software/make/); it defaults to compiling, 
linking and installing all the absolutes on both architectures, but it can also
be used to build specific targets by passing it as arguments:
```
$ ./midas_build obsSelection.Abs var.Abs
```

A *target* is something (often a file) to build; you can get information on
available targets by calling
```
$ make help
USAGE:
    source ./config.dot.sh
    make [-j ${MIDAS_COMPILE_NCORES} -O] [OPTIONS] [TARGETS] [VERBOSE=(1|2)]
OPTIONS:
    consult make manual: man make
TARGETS:
    %.Abs                          link an absolute 
    %.f90                          preprocess an ftn90 file
    %.o                            compile an object 
    all                            compile all programs 
    clean                          delete all build directories 
    cleanabs                       delete all absolutes
    cleandep                       delete all dependency file
    cleanobj                       delete all objects
    depend                         generate all dependency files
    diagrams                       build diagrams (not available yet)
    doc                            build documentation (not available yet)
    help                           print this help
    install                        install all programs
    ssm                            build SSM package 
```
See [next section](#using-make) for more on targets.


### Auto-completion

`midas_build` comes with a bash auto-completion feature, such that argument 
passing can be auto-completed by pressing `<TAB>`:
```
$ midas_build <TAB><TAB>
Display all 141 possibilities? (y or n)
$ midas_build obsImpact.<TAB><TAB>
obsImpact.Abs  obsImpact.o
```

However, this feature needs to be installed:
```
./install_build_completion.sh
...
Auto completion for midas_build installed

To use it directly (in the present shell):
   `source /home/mad001/.profile.d/interactive/post`
   `source /home/mad001/.bash_completion`
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
* `COMPILE_MIDAS_COMPF_GLOBAL` is now `MIDAS_COMPILE_COMPF_GLOBAL`
* `COMPILE_MIDAS_ADD_DEBUG_OPTIONS` is now `MIDAS_COMPILE_ADD_DEBUG_OPTIONS`

If any of those are defined in your profile, you should change them to respect 
this new convention in order to obtain the expected result.

## Using `make`

If you want to have a more fine grained control, you can call `make` directly,
I invite you to read its man pages (short and straight to the point).  

You'll first need to source (and edit if you want) the compilation environment:
```
$ source ./config.dot.sh
```
Otherwise, only two targets will be available: `clean` and `help`.

Then you're good to go!
You can call make to build any target.

A target may be an object file, a specific program or a label (or *phony* 
target such as `all`, `clean` or other label that are not a file *per se*), you
can always use autocompletion by pressing `<TAB>`:
```
$ make <TAB><TAB>
Display all 147 possibilities? (y or n) <y>
Makefile                       install
absolutes                      kdtree2_mod.o
addIncrement.Abs               lambmatrixhi_mod.o
addIncrement.o                 lamspectraltransform_mod.o
adjointTest.Abs                letkf.Abs
adjointTest.o                  letkf.o
advection_mod.o                localization_mod.o
advector.Abs                   localizationfunction_mod.o
advector.o                     localizationspectral_mod.o
all                            mathphysconstants_mod.o
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
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/clib_interfaces_mod.f90 -o clib_interfaces_mod.o > /dev/null 
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/utilities_mod.f90 -o utilities_mod.o > /dev/null 
...
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/programs/var.f90 -o var.o > /dev/null 
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4\
    -lf90sqlite ... obsfilter_mod.o -o var.Abs
```
The option `--dry-run` (or `-n`) only prints what it would be doing without actually doing it.


Some frequently used phony targets are:
* `help` : print a short synopsis and important targets
* `all` : compile all programs on current architecture 
* `info` : print information to stdout
* `install` : install all compiled programs   
  Copy them in `${MIDAS_COMPILE_DIR_MAIN}/midas_abs/` and rename them with version number:  
  `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the
  `../midas.version.sh` script.
* `info` : print information to stdout
* `clean` : remove all objects, programs, intermediate files, everything that was produced by `make` **from all versions of MIDAS** in the build directory.
* `cleanabs` : remove all but installed programs in `${MIDAS_COMPILE_DIR_MAIN}/midas_abs`
* `cleanobj` : remove objects and dependencies
* `cleandep` : remove dependencies files ([see Automatic dependencies below](#automatic-dependencies))

Omitting the target defaults to `all`.

### The `install` target
Calling `make install` **after** `make [all]` will copy the absolute **on the present architecture** to the binaries directory at `${MIDAS_COMPILE_DIR_MAIN}/midas_abs`.  All binaries are copied at the same place with the naming convention `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the `../midas.version.sh` script.


A complete install is then 
```
(source ./config.dot.sh && make && make install)
```
launched from all platforms (what is done by `./midas_build` without argument).


## Calling make in parallel

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




## Out-of-tree compilation

By default `make` will build in `${HOME}/data_maestro/ords/midas-bld` and will symlink it to `../compiledir`.  
You may modify that default in `./config.dot.sh`.  



## Incremental builds


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
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/varqc_mod.f90 -o varqc_mod.o
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/minimization_mod.f90 -o minimization_mod.o
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/programs/var.f90 -o var.o
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




## Automatic dependencies

`make` determine dependencies at build time.
It is thought a better practice to find out about these dependencies automatically.


When `make` is called on a `clean` (or `cleandep`) state, it will determine all objects dependencies using `makedepf90` and include it dynamically in the `Makefile` and then launch the build.

Making only the dependencies: 
```
$ make depend 
$ tree ../compiledir
../compiledir
└── v_3.5.2-133-g111551f_M
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


This is all well, but there is (at least) one case where a user could make that strategy fail.  **Once dependencies are evaluated, they are considered static**.  That means that if after having compiled, someone modify a module's dependencies (adding a `use` statement somewhere) and want to use the incremental awesome feature of `make` it will most probably fail, because the dependencies won't be automatically updated.
There is no simple way to do that in Fortran 2003 (it is possible using Fortran 2008 submodules, but that would require a lot of refactoring.)

So if you find yourself in such a situation, rebuild explicitly the dependencies **before** rebuilding, using `cleandep`
```
$ make cleandep [depend|all|...]
```


## Adding a new program or changing external dependencies

In the previous solution, when a new program is added, two files needed to be changed in `./programs/src_files/`:
* `src_files_${PGM}.sh`
* `compile_setup_${PGM}.sh`

The first contains dependencies information and this is dealt with automatically ([see previous section](#automatic-dependencies)).
The second one contain external libraries that are needed at link time by the programs.
The information contained in **all** `compile_setup_*.sh` files is now found in [`./programs/programs.mk`](programs/programs.mk).

So **when a new program is added** or when **external libraries change for an existing program**, edit the [`./programs/programs.mk`](programs/programs.mk) file and list all external libraries (previously in `compile_setup_${PGM}.sh`) as prerequisite of the absolute target, such as:
```
var.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
         rttov_parallel rttov_main rttov_emis_atlas rttov_other\
         $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random
 ```


## SSM packaging

To publish the absolutes in a SSM domain, one have to
1. update `./config.dot.sh` `MIDAS_SSM_*` variables  
   (making sure they have write privilege to `${MIDAS_SSM_TARGET}`)
2. for **each architecture** `(source ./config.dot.sh && make ssm)`
3. once all architectures have been published, protect the domain:  
   `(source ./config.dot.sh && make ssm_protect)`


## What is left to do

Candies.  

* automated `doc` building and `diagrams`, etc.

But most of all... taking into account your input.  
Don't hesitate to contact-me for your input or for some guidance: @mad001


Many thanks to @phb001 and [mad-scientist](http://make.mad-scientist.net/)
