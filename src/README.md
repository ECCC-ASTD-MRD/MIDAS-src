Building MIDAS using make
=========================

The present `README` **assumes you are in the `./src` directory**.


I just want to build!
---------------------

From `PPP[34]`, edit `./config.dot.sh`, then
```
> ./build.sh
```
It will build MIDAS executables on both `PPP4` and `Daley` submitting jobs with the requested number of cores.

Be aware that it is dependent on `makedepf90`, which [is installed but not *persistent* yet](https://gitlab.science.gc.ca/hpc_migrations/hpcr_upgrade_1/issues/980).
It should be solved soon, but if you get an error, it could well be that.




Using `make`
------------
If you want to have a more fine grained control, you can call `make` directly, I invite you to read its man pages (short and straight to the point).  

You'll first need to source the compilation environment:
```
> source ./programs/commons/compile_setup.sh
```
Then you're good to go!

You can call make to build any *target*.
A target may be an object file, a specific program or a label (or *phony* target such as `all`, `clean` or other label that are not a file *per se*), you can always use autocompletion by pressing `<TAB>`:
```
> make <TAB>
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
> make var <TAB>
var.Abs            var.o              varnamelist_mod.o  varqc_mod.o
```
If `<TAB>` does not work (or just show you `all` and `clean*`), it is probably because either you don't have `makedepf90` installed or did not `source ./programs/commons/compile_setup.sh`.
Autocompletion does not work on the backends.

When you ask `make` to build a target, it will determine everything that needs to be done to achieve that goal, for instance if want to build `var.Abs`:
```
> make var.Abs --dry-run
Preprocessing codeprecision_mod.f90 inplace
Preprocessing clib_interfaces_mod.f90 inplace
Preprocessing rttov_interfaces_mod.f90 inplace
Generating object dependencies > dep.real8.inc
Generating executables dependencies > dep.real8.abs.inc
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/clib_interfaces_mod.f90 -o clib_interfaces_mod.o > /dev/null 
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/utilities_mod.f90 -o utilities_mod.o > /dev/null 
...
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4 -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/programs/var.f90 -o var.o > /dev/null 
s.f90 -openmp -mpi -mkl -check noarg_temp_created -no-wrap-margin -O 4\
    -lf90sqlite -ludfsqlite -lrttov_coef_io -lrttov_hdf -lrttov_parallel -lrttov_main -lrttov_emis_atlas -lrttov_other -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lburp_module -lvgrid -lirc -lrpn_comm -lrandom -lhpcoperf -lsqlite3 -lrmnMP var.o randomnumber_mod.o tt2phi_mod.o burpfiles_mod.o ... obsfilter_mod.o -o var.Abs
```
The option `--dry-run` (or `-n`) only prints what it would be doing without actually doing it.


Some frequently used phony targets are:
* `all` : compile and install all programs on current architecture 
  (in both `REAL` precisions, [see Multi-precision problematic below](#multi-precision-problematic))  
* `info` : print information to stdout
* `install` : install all compiled programs for both precisions
  Copy them in `${DIR_BLD_ROOT}/midas_abs/` and rename them with version number:  
  `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the
  `../midas.version.sh` script.
* `real4` : compile `REAL 4` programs (`obsIO.Abs` and `prepcma.Abs` for the EnKF)
* `real8` : compile `REAL 8` programs (most MIDAS programs) 
* `info` : print information to stdout
* `install` : synonym for `all`
* `install4` and `install8` : install only specific `REAL` precision programs
* `clean` : remove all objects, programs, intermediate files, everything that was produced
  by `make`
* `cleanabs` : remove all but installed programs in `${DIR_BLD_ROOT}/midas_abs`
* `cleanobj` : remove objects and dependencies
* `cleandep` : remove dependencies files ([see Automatic dependencies below](#automatic-dependencies))

Omitting the target defaults to `all`.

### The `install` target
Calling `make install` **after** `make [all]` will copy the absolute **on the present architecture** to the binaries directory at `${DIR_BLD_ROOT}/midas_abs`.  All binaries are copied at the same place with the naming convention `midas-_${ORDENV_PLAT}-${VERSION}.Abs` where `${VERSION}` is obtained by the `../midas.version.sh` script.


A complete install is then 
```
make all && make install
```
launched from all platforms (what is done by `./build.sh`).

To compile and install only single precision:
```
make real4 && make install4
```
(`make install` would fail)



Calling make in parallel
------------------------

To compile a target using multiple cores use
```
make -j${NCORES} [<target>]
```

Careful to not overlaod the head node.  You can either use `ord_soumet` or if you are having a compile-open-house, ask for a interactive party node:
```
> r.load /fs/ssm/main/opt/jobsubi/jobsubi-0.3
> jobsubi --show-request  -r ncores=20 ppp4
```
However, `makedepf90` is [not installed yet on compute nodes](https://gitlab.science.gc.ca/hpc_migrations/hpcr_upgrade_1/issues/980), but could be available through conda([see #255](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/255#note_143881), or ask-me: @mad001 )




Out-of-tree compilation
-----------------------

By default `make` will build in `./build/`.  You can explicitly pass `${DIR_BLD_ROOT}` to `make` to do otherwise, it will `mkdir -p` the directory and build from there:
```
> make DIR_BLD_ROOT=../compile
```
Note that when you call `make clean` or `make wipe` **you have also to specify the build directory**, otherwise, it will use the default.
You can always `export DIR_BLD_ROOT`.



Incremental builds
------------------

If you modify a single source file, `make` will determine which targets (objects and absolutes) depend on it and will reprocess only those.
For example, let say you have checked out and compiled all programs, then you work on `modules/varqc_mod.f90`.  When you are done, only `minimization_mod.o`, `var.o` will have to be recompiled and `var.Abs` relinked.  This is done using file time stamps and dependency tree.

What if you just modify the implementation of a function without touching it's interface?
You should not have to recompile other modules.
But `make` won't know what you did, it will just look at the time stamp and decide to recompile everything that depends on the modified file.  Actually just touching the file, will make `make` thinks it needs to reprocess all dependent files:
```
> make -n
...
make[2]: Nothing to be done for 'absolutes'
...
> touch modules/varqc_mod.f90
> make -n
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/varqc_mod.f90 -o varqc_mod.o
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/modules/minimization_mod.f90 -o minimization_mod.o
s.f90 ... -c /fs/homeu1/eccc/aq/arqi/mad001/code/midas/src/programs/var.f90 -o var.o
s.f90 ... -o var.Abs
```

One solution to that is to `touch`  **targets** (using `make --touch` or `-t`) , making them newer that the dependencies:
```
> make --touch var.Abs
touch varqc_mod.o
touch minimization_mod.o
touch var.o
touch var.Abs
touch var.Abs
> make -n
...
make[2]: Nothing to be done for 'absolutes'.
...
```



Multi-precision problematic
---------------------------

Two programs (`obsIO.Abs` and `prepcma.Abs` used in the EnKF) need to be compiled in single precision, `REAL 4`, while the others are compiled in `REAL 8`. **But, they all share the rest of the code base**, such that some objects need to be compiled in both precisions.

In the former strategy, is that there is a separate compilation of every objects for each program and palatalization is done at that level.  The advantage is that there is no special care to be done with objects `REAL` precision.  On the down side, multiple copies of the exact same objects are recompiled and time is wasted.

The strategy here is different.  In the build directory, the structure is:
```
./
└── ${ORDENV_PLAT}
    └── ${EC_ARCH}
        ├── real4
        └── real8
```
And all programs for a given `(${ORDENV_PLAT},${EC_ARCH})` and a given `REAL` precision will share their objects which are compiled in parallel.  So: much more efficient parallelization.

If you just want to build `REAL 4`:
```
> make real4
```
and `make real8` to build only `REAL 8`; `make all` does both.

If you want to build a single non-phony target, say `bmatrixhi_mod.o`, **you'll need to explicitely specify you are building with `REAL 4`** by passing `MIDAS_COMPILE_REAL_SINGLE=true` to `make`:
```
make bmatrixhi_mod.o MIDAS_COMPILE_REAL_SINGLE=true
```
otherwise, you'll get an error.
This is necessary for the compilation to happen in the same directory (`.../real4/`) as the objects compiled with the same `REAL` precision.
(remember that the object `bmatrixhi_mod.o` could also be compiled in `REAL 8`, so make cannot know which precision is wanted)

However, I added an *ad hoc* fix for both simple precision programs, `obsIO.*` and `prepcma.*`, since these are always compiled in single precision, `make` will implicitely append `MIDAS_COMPILE_REAL_SINGLE=true` when called, such that
```
> make obsIO.Abs
```
is equivalent to 
```
> make obsIO.Abs MIDAS_COMPILE_REAL_SINGLE=true
```
Be careful however with that, it is preferable to append it explicitely, I could not garantee what would be the result of `make bmatrixhi_mod.o obsIO.Abs` for instance...



Automatic dependencies
----------------------

`make` determine dependencies at build time.
It is thought a better practice to find out about these dependencies automatically.


When `make` is called on a `clean` (or `cleandep`) state, it will determine it's dependency trees using `makedepf90` and include it dynamically in the `Makefile` and then launch the build.

`makedepf90` is not yet [permanently installed on the network yet](https://gitlab.science.gc.ca/hpc_migrations/hpcr_upgrade_1/issues/980), but is installed on `PPP4` head node and is available through conda([see #255](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/255#note_143881)) for those who would like to try out on interractive nodes (with `jobsubi`).

Making only the dependencies: 
```
> make depend DIR_BLD_ROOT=./compile                          
> tree compile
compile/
└── ubuntu-18.04-skylake-64
    └── intel-19.0.3.199
        ├── real4
        │   ├── dep.real4.abs.inc
        │   └── dep.real4.inc
        └── real8
            ├── dep.real8.abs.inc
            └── dep.real8.inc
```
And this is what it looks like inside of `dep.real4.inc`:
```make
obstimeinterp_mod.o : obstimeinterp_mod.f90 obsspacedata_mod.o timecoord_mod.o utilities_mod.o mpivar_mod.o mpi_mod.o 
thinning_mod.o : thinning_mod.f90 utilities_mod.o obsspacedata_mod.o bufr_mod.o codeprecision_mod.o 
increment_mod.o : increment_mod.f90 varnamelist_mod.o bmatrix_mod.o gridVariableTransforms_mod.o utilities_mod.o humiditylimits_mod.o verticalcoord_mod.o horizontalcoord_mod.o gridstatevector_mod.o timecoord_mod.o mpivar_mod.o mpi_mod.o
...
bgckAtms.o : bgckAtms.f90 bgckmicrowave_mod.o
```
`dep.real[48].inc` are generated by `makedepf90` and produce *superficial* dependencies, those needed to compile **objects**, but **not enough to link absolutes**.

To proceed with linking we need the fully recursive dependencies.
We parse (with a simple python script, [`recursiveDep.py`](./recursiveDep.py)) the superficial dependency files and deduce the fully recursive ones needed at link time; they are in `dep.real[48].abs.inc` files.

Notice also that dependencies are **only generated for the frontend** ([`makedepf90` is not available on the backend for now](https://gitlab.science.gc.ca/hpc_migrations/hpcr_upgrade_1/issues/980)).
This is not a real problem since the dependencies are the same; we will create the relevant directory structure and copy the dependencies there.
This is done in [copy_depend_backend.sh](./copy_depend_backend.sh) (and also [build.sh](./build.sh) that does the whole multiplateform compilation).

This is all well, but there is (at least) one case where a user could make that strategy fail.  **Once dependencies are evaluated, they are considered static**.  That means that if after having compiled, someone modify a module's dependencies (adding a `use` statement somewhere) and want to use the incremental awesome feature of `make` it will most probably fail, because the dependencies won't be automatically updated.
There is no simple way to do that in Fortran 2003 (it is possible using Fortran 2008 submodules, but that would require a lot of refactoring.)

So if you find yourself in such a situation, rebuild explicitly the dependencies **before** rebuilding.
```
> make cleandep depend DIR_BLD_ROOT=./compile
```

### Building dependencies on backend
Backends don't have `makedepf90` yet and probably won't have it.
But dependencies are the same, so to build dependency files you need to do it on the frontend and copy them at the right place.

This is done automatically using the script `./copy_depend_backend.sh` (which call a the same bash function that `./build.sh` does).


Adding a new program or changing external dependencies
------------------------------------------------------
In the previous solution, when a new program is added, two files needed to be changed in `./programs/src_files/`:
* `src_files_${PGM}.sh`
* `compile_setup_${PGM}.sh`

The first contains dependencies information and this is dealt with automatically ([see previous section](#automatic-dependencies)).
The second one contain external libraries that are needed at link time by the program.
The information contained in **all** `compile_setup_*.sh` files is now found in `./programs/programs.mk`.
This file is separated in two sections, one for each `REAL` precision group of programs (all but two programs for the EnKF are double precisions, **no new program should be single precision**).
Each programs need to be listed in `PGM_DBL`

> If LETKF replace EnKF, this whole two precisions thing is going to be removed and it is going to simplify the structure **a lot**.  Amongst other things, the `PGM_DBL` list will then be automatic, using `$(wildcard *.f90)`.

So **when a new program is added** or when **external libraries change for an existing program** two things need to be done in the `./programs/programs.mk` file:
1. if needed, add the program name in the `PGM_DBL` list variable
2. list all external libraries (previously in `compile_setup_${PGM}.sh`)
   as prerequisite of the absolute target, such as:
   ```
   var.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
           rttov_parallel rttov_main rttov_emis_atlas rttov_other\
           $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random
   ```

What is left to do
------------------

Candies.  

* automated `ssm` packaging
* automated `doc` building and `diagrams`, etc.
* some more control in the `build.sh` script

But most of all... taking into account your input.  
Don't hesitate to contact-me for your input or for some guidance: @mad001


Many thanks to @phb001 and [mad-scientist](http://make.mad-scientist.net/)
