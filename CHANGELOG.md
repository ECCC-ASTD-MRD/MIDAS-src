# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

 * (Nothing yet)

### Changed

 * (Nothing yet)

### Fixed

 * (Nothing yet)

### Removed

 * (Nothing yet)

## [3.3.0-rc3]: Unreleased

### Added

* The program `ensManip` can now read a file `targetgrid` in the
  working directory on which grid all fields will be interpolated
  (#138 and !128).  By default, all files will be interpolated on the
  grid of the ensemble.

### Changed

* Write the hessian after 'rebm', 'rehm' and 'anlm' files (#142 and !126)

### Fixed

* Fix 'get_avhrr_emiss' when some channels are missing (#140 and !125)
* A fix was done to control the minimum value for 'HU' after
  interpolating profiles from background state levels to analysis
  levels (#144 and !127).  This is affecting very weakly the results of
  most test for program `midas-var`.
* Reject observations with unrealistic lat-lon values (#137 and !128)

## [3.3.0-rc2]

This version is identical to [3.3.0-rc1].

## [3.3.0-rc1]

### Added

* The program `ensManip` can now compute the standard deviation of an
  ensemble of forecasts and recenter the ensemble forecasts around a
  specified mean (#55, #65, #104, #131, !48, !52, !94 and !116).
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

### Fixed

* Reactivation of the Scale-Dependent Localization (#40 and !35)
* Fix and a modification to the LQ to HU tangent-linear transform (#61 and !49)
* Fix for the vertical interpolation of LAM ensembles (#46 and !53)
* Fix bugs in gsv_writeToFile related to tic toc records (#70 and !62)
* Fix the I/O for hessian in LAM mode (#128 and !112)
* Fix to statetocolumn_mod in the extremely rare case where a
  processor had some observations before the load balancing but none
  after (#129 and !113)

### Removed

* (Nothing yet)

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

[Unreleased]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.0-rc3...HEAD
[3.3.0-rc3]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.0-rc2...v_3.3.0-rc3
[3.3.0-rc2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.0-rc1...v_3.3.0-rc2
[3.3.0-rc1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.2.2...v_3.3.0-rc1
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
