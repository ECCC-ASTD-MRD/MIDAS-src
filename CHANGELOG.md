# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
 * Minor fix in tt2phi_mod.f90 (slightly affects results) : now setting near-sfc temperature and momentum altitude levels to their known height offset (#180 and !212).
 * Change '!' to '!!' in those places that bothered Sphinx.
	
### Added

 * The environment variable `MIDAS_MAKE_LINKS_MACHINE_LIST` can be
   used to control the hosts on which links will be created by
   `install_suite.sh` in the maestro test suite.  By default, only the
   links on which the suite will run are created. (#231 and !216)
 * The logical namelist variable ltopofilt has been removed. You probably MUST update your namelist. The new namelist variable is called list_topoFilt. This string array variable allow to activate the topographic rejection criteria for selected observation families. See the namelist in the unit tests from examples (#225 and !211).
 * Add ability to define a local domain and control inclusion of each variable for energy norm (#207 and !204)
 * Add the program `prepcma` to reproduce the similar program in the EnKF codebase (#189 and !198)
 * 3DVar analysis of SST data (family 'TM') can now be computed without MPI (#203 and !195)
 * The scripts to build the MIDAS SSM domain are now in the MIDAS
   depot (#187 and !186).  See the [README](README.md) for more information.
 * A column `OBS_CRPS` has been added to `obsSpaceData` (#185 and !188).
 * `Bmatrixensemble_mod` can now read an ensemble of perturbations like lagged forecast differences (#193 and !190)
 * Variance smoothing is now possible in `bMatrixEnsemble_mod` (#193 and !190)
 * ScaleFactor added to `ensManip` for vertically scaling ensemble perturbations when recentering (#186 and !179)
 * Minor changes to extend capability of generating OmP (and OmA) diagnostics via the CH obs family (#184 and !177)
 * Visibility and wind gust near the surface can now be assimilated (#173 and !176)
 * A sea ice concentration analysis can now be done with CIS daily ice charts (other obs types to come) (#163 and !175)
 * A horizontal land/sea mask can now be included in the analysisgrid file (#163 and !175)
 * The program, midas_obsSelection, was created (so far only for aladin HLOS obs), comprising O-P computations, background check, and thinning (#113 and !174)
 * Enable `vcode=5005` for ensemble B matrix (#188 and !173)
 * The aladin HLOS wind observations can now be adjusted to compensate for the observations having been calculated at another meteorolgical facility (#139 and !172)
   * A new namelist, NAMALADIN_OBS, is required when there are height-level observations (i.e. aladin) data to be treated.
 * Able to write contents of obsSpaceData to simplified sqlite files (useful when using another input file format and sqlite wanted for diagnostics) (#167 and !166)
 * Add observation operator and background check for aladin (HLOS winds) (#114 and !163)
 * Level-dependant steering flow scalefactor capability for advection (#168 and !146)
 * Add checks on humidity limits in ensManip (#164 and !143)
 * New script to automatically generate module dependencies: `make_src_files.sh` (#149 and !136)
 * Removal of constraints on spectral truncation and number of levels relative to the MPI topology (#135 and !135)
 * New functionality: now able to read various types of sea ice data (family =*GL*) (#127 and !131)

### Changed

 * The height/pressure are computed on the grid, before horizontal interpolation to observation locations, to prepare for using slanted columns and footprint operators. (#124 and !220)
    * The height/pressure are part of statevector/columndata main data storage arrays (gd_r4/gd_r8/all) and are calculated for the trial fields and the increments.
    * Allocation of height/pressure is set to true, by default, and it is done if the necessary variables for their calculation are available in the statevector/columndata.
    * Z_M/Z_T and P_M/P_T are the height and pressure on grid on the TH and MM levels in varNameList.
    * dPdPsfc is no longer used in any observation operators since the increment of pressure is calculated on the grid and is interpolated to the observation location.
    * Change namelist variable addGZsfcOffset to addHeightSfcOffset.
    * Variable/function/subroutine names that include `gz` are changed to `height` to reflect the fact that geometric altitude/height is now the primary variable instead of geopotential.
    * Memory requirements are higher for some programs and configurations (but not gdps and rdps configurations). 
    * The execution time is also increased for some (e.g. gdps takes ~100 seconds longer, but this can be reduced by increasing number of nodes to 30 or 36).
 * The namelist variable `scaleFactor` in NAMBHI must now be specified in all 3DVar configurations because default value was changed from 1.0 to 0.0. (#224 and !209)
 * CalcStats in LAM mode was made MPI compatible (#158 and !202)
 * Replacing the old numerical recipe for generating gaussian random
   values by a much more efficient method (#82 and !192)
   [`random_tools`](https://gitlab.com/mfvalin/random_tools).
    * This is changing the results only for the program `randomPert.Abs`.
 * Change to the height (GZ) calculation within MIDAS so that GPS-RO and GPS-gb now use same heights as other obs types (#141 and !191)
    * The variable was also changed from geopotential to altitude
    * New namelist variable to allow using static GPS-RO observation error variance
    * Small change to H(x) results due to these modifications for GPS-RO, GPS-gb, GZ and possibly sfc obs
 * When compiling with `COMPILE_MIDAS_ADD_DEBUG_OPTIONS=yes`, the options `-debug -check all -O 0` are added to the compile command (#182 and !187)
 * The body of the program `oMinusF` was extracted into a new module, `ominusf_mod` (#113 and !174)
 * **Major change:** New approach for horizontal/temporal interpolation of background state and increment to observation locations/times (#80 and !147)
    * Changes to how the background state is read; now relies on gridstatevector_mod
    * File copy to ramdisk for all files can be done by fortran code (simplifies scripts)
    * Linear time interpolation can be applied to the background state when computing the innovations (controlled by namelist variable, default is nearest neighbour)
    * Small change to results for most applications and some increases in time and memory requirements
 * Use constituent BUFR elements from official tableburp file released 30 Sept 2018 (#150 and !137)
 * The input file trlp is no longer necessary in varbc mode (Slight unsignificant change in the results of varbc) (#145 and !129)

### Fixed

 * Correction of the conversion factor used to compute air mass predictors in biascorrection_mod.f90. Problem introduced in !191. (#219 and !219)
 * Fixed some potential bugs detected while compiling with `-check all` (#182 and !187)
 * Improved efficiency of ensemble amplitude memory access and writing of `rehm` and `anlm` files (#170 and !153)
 * Fix the selection of GPSRO-bending angle observations (#151 and !145)

### Removed

 * (Nothing yet)

## [3.3.2]:

### Changed

 * Initialize 4 variables in routine write_info in the module burpread_mod (#172 and !182)
 * The positional records will now have the same etiket at the fields they represent (#190 and !178)

## [3.3.1]:

### Changed

* `ensManip` now support humidity adjustments in recentering mode (#174 and !157)
* The program `obsIO` is using the value of `OBS_REAL` to work at
  single precision when working with `obsSpaceData_mod`
  (#175 and !161).
   * The module `burpread_mod` now includes the function
     `brpr_getTypeResume` which returns the module private variable
     `TYPE_RESUME`
* Improved efficiency of ensemble amplitude memory access and writing
  of `rehm` and `anlm` files (#170 and !151)

## [3.3.0]:

### Added

* The program `ensManip` can now compute the standard deviation of an
  ensemble of forecasts and recenter the ensemble forecasts around a
  specified mean (#55, #65, #104, #131, !48, !52, !94 and !116).
* The program `ensManip` can now read a file `targetgrid` in the
  working directory on which grid all fields will be interpolated
  (#138 and !128).  By default, all files will be interpolated on the
  grid of the ensemble.
* A program `addIncrement` has been added (#38, #53, #54, #123, #126, !102, !41, !47, !100, !106 and !115)
  * The namelist `NAMADDINC` has been renamed `NAMINC`.  The namelist
    variable `CETIKINC` in namelist `NAMMIN` has been moved to
    `NAMINC` with name `ETIKET_REBM`. The namelist variable
    `WRITEANALYSIS` was moved from `NAMMIN` to `NAMCT0`.
* Variational bias correction functionality have been added (#41 and !43)
* Add a 2D mode (#32, #51, !45)
* Add assimilation of SST observations (code is `22042`) (#111 and !123)
* Add a `thinning_mod` fortran module (#110 and !120)
* First step towards controlling precision (#47, !46)
* Benjamin Menetrier's localization lengthscale diagnostics is now
  available in global and lam mode in calcstats (#31 and !34)
* Add a new type of global static B matrix with some latitudinal heterogeneity of the
  correlations (#39, #48, !37 and !40)
* Add a new program, `ensembleH`, to apply H to ensemble (#60 and !51)
* Add automatic increment normalization in diagBmatrix
* Allow computation of local horizontal correlation in calcstats for LAM
* Add namelist variables `nelems_altDiffMax`, `list_altDiffMax` and
  `value_altDiffMax` in `NAMFILT` for maximum difference of altitude
  for surface data (#71 and !57)
* Add a program `obsImpact` for FSO (#56, #58, #105, !65, !86 and !103)
* Add a program `adjointTest`
* Add a program `diagHBHt` to compute HBHt using a randomization approach (#72 and !69)
* The documentation is automatically generated for each commit in the `master` branch (#78, #125, !71 and !110)
* Use GitLab-CI to run automatic tests (#88 and !74)
* Add `scaleFactorCC` to `NAMBHI` namelist
* Exclude ensembles for specified variables (#98, !85 and !93)
* Add an implicit version of the diffusion operator (#89 and !88)
* Add new functionality to read and write observations in sqlite format (#64, #117, !98 and !102)
* Add a new program `obsIO` for testing observation I/O routines (#118 and !108)
* Introduce module for computing slant path positions (#116 and !109)
* Include Yin-Yang support in `horizontalCoord_mod` (#134 and !119)
* Include a custom `r.run_in_parallel` based on the one available in
  `rpn/utils/16.2.2` but with changes to use `/bin/bash` instead of
  `/bin/ksh` in the script launched in parallel (#136 and !122)
   * This is to be removed when `r.run_in_parallel` will be officially
     released.

### Changed

* Using RTTOV-12 v1.1 from which `lapack.o` has been removed to use
  the system library which is faster
* The directories have been reorganized (#50, !42)
* Relax constraints on MPI topology (#49, !44)
* Make the Bnmc-LAM 2D mode fullly functionnal (#62 and !50)
* Make `controlVector_mod` more general (#73 and !58)
* Refactor all system tests to test directly all the programs (#59 and !60)
* Improve `tim_getDateStampFromFile`: using trial/analysis input file
  to determine the date (#74, #76 and !63)
* Add chemical consituents capacity (#98, !95 and !89)
* Make advection modular and flexible (!87 and !99)
* Improve advection code (#87, #119, !99 and !105)
* `HU` rather than `LQ` in `gridStateVector` outside $`B`$ matrix modules. This
  has significant impact on the results and also requires changes to
  the use of the randomPert program within the EnKF (a background
  state must now be supplied) (#67 and !111).
* The copy of files to RAMDisk is now done directly in the fortran
  code and no longer in the scripts (#133 and !117)
* Removed namelist variable write_mpi from NAMENKF (#134 and !119)
* Centralize unit conversion and convert to Kelvin (#134 and !119)
* Simplify calculations in `windRotation_mod` (#134 and !119)
* Decrease a threshold for vertical interpolation for TOVS (#134 and !119)
* Write the hessian after 'rebm', 'rehm' and 'anlm' files (#142 and !126)

### Fixed

* Reactivation of the Scale-Dependent Localization (#40 and !35)
* Fix and a modification to the LQ to HU tangent-linear transform (#61 and !49)
* Fix for the vertical interpolation of LAM ensembles (#46 and !53)
* Fix bugs in gsv_writeToFile related to tic toc records (#70 and !62)
* Fix the I/O for hessian in LAM mode (#128 and !112)
* Fix to statetocolumn_mod in the extremely rare case where a
  processor had some observations before the load balancing but none
  after (#129 and !113)
* Fix 'get_avhrr_emiss' when some channels are missing (#140 and !125)
* A fix was done to control the minimum value for 'HU' after
  interpolating profiles from background state levels to analysis
  levels (#144 and !127).  This is affecting very weakly the results of
  most test for program `midas-var`.
* Reject observations with unrealistic lat-lon values (#137 and !128)

## [3.2.2] - 2018-05-09

No change to the fortran code so it is equivalent to `v_3.2.1`.

## [3.2.1] - 2018-05-03

### Fixed

- Fix a bug when IASI are missing at the background check step.  It
was causing an aborting in the following analysis.  We also correct a
non-initialized variable.  (#94 and #96)

## [3.2.0] - 2018-03-02

This release comes from `v_3.0.4` to ignore changes with `lapack.o` in
RTTOV-12.

### Fixed

- Fix a bug which affects `BURP_update` (#66 and #79)

## [3.1.0] - 2018-01-28

This version will be ignored.

Since the release `v_3.0.5` has changed the results, we should have
tagged it `v_3.1.0` to follow semantic versioning.

## [3.0.5] - 2018-01-12

This version was not implemented in operations.  This release has been
ignored.

### Changed

Using `rttov/12v1.1` which does not contain `lapack.o`.  We want to
use the system librairies.  This does not impact the timings but the
results are very sligthly changed.

## [3.0.4] - 2017-11-28

### Fixed

- Fix a bug where an assumption was incorrectly made that the number
of vertical levels of the trial field was the same as in the analysis
vertical grid.

## [3.0.3] - 2017-11-28

### Fixed

- Fix the situation where the grid for TG is not the same as for other
     fields like TT.  This is the case where the trial field comes
     from a LAM (#45).

## [3.0.2] - 2017-11-16

### Fixed

- Fix MPI related problems in `tvs_rttov_read_coefs` (!39 and #44)

## [3.0.1] - 2017-11-08

### Fixed

- Fix the rare case where all obs on one mpi task are sent to another
mpi task during the redistribution step (!38 and #42)

## [3.0.0] - 2017-10-30

This is the initial version delivered in final cycles for the GDPS 6.1 project in 2018.

### Added
- Includes backward comptatible changes to conventional observations by Stéphane Laroche

### Changed
- Introducing the use of RTTOV-12 library (release `1.0`) (non backward compatible) (Sylvain Heilliette)
- Using `cmda/libs/16.2-6`
- Reduce memory usage for IR bgck by a factor of 5

## [2.2.0] - 2016-09

This is the first version published and use on the HPCR platforms
`ubuntu-14.04-amd64-64` and `sles-11-amd64-64` on the `science.gc.ca`
network.

Some other `v_2.2.*` subsequent versions have been published but we
are not documenting them here.

[Unreleased]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.2...HEAD
[3.3.2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.1...v_3.3.2
[3.3.1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.0...v_3.3.1
[3.3.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.2.2...v_3.3.0
[3.2.2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.2.1...v_3.2.2
[3.2.1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.2.0...v_3.2.1
[3.2.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.4...v_3.2.0
[3.1.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.5...v_3.1.0
[3.0.5]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.4...v_3.0.5
[3.0.4]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.3...v_3.0.4
[3.0.3]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.2...v_3.0.3
[3.0.2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.1...v_3.0.2
[3.0.1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.0.0...v_3.0.1
[3.0.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_2.2.0...v_3.0.0
