# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

 * When `midas.splitobs` is splitting an observation file into equal
   parts we process not only the tables `header` and `data` but all
   the tables containing `id_obs` as a key and we split them in the
   same way as `header` and `data` (#376 and !371).
 * The scripts `midas.*` in the directory `tools/midas_scripts` have
   been taught to process the files `obsfiles_${fam}.beforeThinning`
   produced by the program `midas-obsSelection` which are equivalent
   to the `evalalt` files (#375 and !370).
 * Added NAMTOV namelist variables (#386 and !365)
   * doAzimuthCorrection(:) and userDefinedDoAzimuthCorrection to control correction of Satellite Azimuth Angle
   * isAzimuthValid(:) and userDefinedIsAzimuthValid to specify if the azimuth in observation files is valid
 * Added temporal thinning algorithm for surface observations (#381 and !369)
 * Added namelist variable to include year in random seed in `ensPostProcess_mod` (#389 and !368)
 * Added ability to perform background check for SAR winds (wind speed retrieval) (#299 and !364)
 * Added thinning algorithm for radiosondes for use in `obsSelection` program (#373 and !361)
 * Added `crisCloudFractionThreshold` to namelist section `NAMBGCKIR` (#327 and !360)
 * Added thinning for most obs types in `obsSelection` program (#367 and !357)
 * Added thinning algorithms for most obs types (#357 and !355)
 * Incorporated in `tools/splitobs` the program `midas.splitobs` (#362 and !352)
   * extensive tests can be run with `tools/splitobs/unittest` (see [`README`](tools/splitobs/README.md) for more details)
   * those tests will be run automatically by the CI system only when a tag is pushed
 * Added the arguments `-analysis`, `-forecast_a` and `-forecast_b` to `tools/midas_scripts/midas.launch` (#369 and !351)
   * The program `midas-obsImpact` needs them.
 * Incorporated in `tools/midas_scripts` all the helper scripts that was before in `oavar_scripts` (#364 and !349)
   * renamed them from `oavar.*` to `midas.*`
     * renamed `midas.var_mpi` to `midas.mpi` to remove any reference to the VAR program
   * renamed the environment variables from `OAVAR_` to `MIDAS_*`
   * renamed the arguments in the script starting with `oavar_` to `midas_*`
   * modified the SSM publication process to take the scripts directly in the MIDAS code base
 * Added two tools that were maintained separately before: (#363 and !346)
   * `midas.monitor` which monitors a file to react to its content and
   * `midas.findTrials` which finds the trial name extensions in an assimilation window.

### Changed

 * Make diffusion operator B matrix mpi compatible (#355 and !354)
 * Reduce memory requirements for bgck and other programs (in s2c_nl) (#371 and !353)
   * also speed up burp update by using kdtree
 * Move background procedure from `var` to `obsSelection` program (#359 and !347)
 * Improve documentation about SSM publishing for a single user (#363 and !346)
 * RTTOV radiative transfer performed directly on model levels - NOTE: results affected and not backward compatible (#253 and !323)

### Fixed

 * Fixed bugs in thinning for radiosondes, hyperIR, CSR and aladin observations - changes results (#390 and !366)
 * Fixed bug in `burpClean` that was causing errors for integer valued elements in the data blocks (#373 and !361)
 * Fixed channel indexing bug in hyperspectral IR background check (#379 and !359)
 * Fixed midas.reunir_obs_mpi so as to preserve the content of RESUME and RDB4_SCHEMA tables (#374 and !356)
 * Fixed `ramdisk_mod.copyFile` when file size is bigger than maximum integer (#366 and !350)
 * Fixed bugs introduced during !323 preventing compilation and execution of system tests in debug mode (#365 and !348)

### Removed

 * Remove the programs `write_subdomains` and `reunir_obs` from SSM domain publication (#360, #361 and !345)
 * Compile using only one precision and remove the publication of the MIDAS library (#358 and !344)
   * remove the program `midas-obsIO`
   * convert the program `midas-prepcma` to `CODEPRECISION_OBS_REAL_DOUBLE`
      * this affects the results of the `midas-prepcma`system test
   * remove the automatic generation of the MIDAS library when building the MIDAS SSM domain
   * remove the script `midas.compile.sh` from the SSM domain
      * this script was designed to help users of the MIDAS librairies to compile and link but it has never been used
   * remove the description of the library in the automatic documentation

## [3.5.2]

### Added

 * Add ability to bias correct AI and GP observations in Midas (#305 and !335)
 * Introduction of vertical correlation ansatz functions, with coefficients (one set per satid) in GPSRO namelist (#296 and !332)
 * New namelist variable HTPMAXER in NAMGPSRO namelist section (#343 and !328)
 * Add ability to assimilate ASCAT data (backscatter anisotropy and open water retrievals) for sea ice concentration analysis (#332 and !326)
 * Add ability to generate random perturbations with `LQ` humidity variable in `randomPert` program (#313 and !325 )
 * Add handling of variables needed for land surface analysis, so far only for `ensManip` (#206 and !324)

### Changed

 * Using `oavar_scripts` to version `2.2.7` which includes the following changes since version `2.2.6`:
   * `oavar.launch`:
      - On rend facultatif l'etape `reunir_obs` si on utilise
       `oavar.launch ... -oavar_reunir_obs no`
   * La variable `OAVAR_OBS_MPI_ORDERING` est mise a `regular` par
     défaut
      - Cela est cohérent avec le code de MIDAS depuis la version
       `v_3.5.0`.
   * Adaptation des scripts pour pouvoir tourner `midas-genCoeff`
   * `oavar.mpi_barrier`:
      - Le script est beaucoup moins verbose qu'auparavant.  On peut
       utiliser la variable d'environnement
       `OAVAR_MPI_BARRIER_VERBOSE=yes` pour réactiver le `set -x` dans ce
       script.
   * `oavar.launch` et `oavar.var_mpi`:
      - Ajout du mode `distribute` pour `-splitobs_mode` pour distribuer
      les fichiers sur chacune des tuiles MPI plutôt que d'utiliser le
      programme `splitobs.Abs` pour ce faire
   * `oavar.check_ensemble`:
      - Adaptation de la manipulation du namelist pour transformer un
        EnVar en 3D-Var pour les versions après `v_3.4.2`
 * Changed IR quality control and background check to add protection against missing values for angles (#349 and !341)
  * Move RTPP ensemble inflation and it's namelist variable from `letkf` to `ensPostProcess` (#352 and !339)
  * Efficiency improvements (mostly for global EnVar) (#235 and !337):
    * Allow single precision in parts of the code, controlled by environment variables
    * Compute height and pressure increments on the column instead of the grid, controlled by namelist variable
  * Improvement of the filtering functionality of module biascorrection_mod.f90 (#341 and !330)
  * Program `ensPostProcess` can now be used to just do recentering or computing trial mean and spread (#334 and !327)
    * **Note:** changes must be made to namelist block names for the `ensPostProcess` and `letkf` programs!

### Fixed

 * Fix bug in the rejection filter of GPSRO data such that bending angle with missing azimuth is now correctly rejected (#353 and !343)
 * Fix bug in the way radiosonde weights are read and interpolated in genCoeff (#354 and !342)
 * Fix minor bug in `midas-letkf` when no obs near analysis grid point (#352 and !339)

### Removed

 * Remove the programs `addIncrement` and `seaIce` because they are no longer needed (#345 and !334)

## [3.5.1]

### Added

 * Add ability to output unperturbed subspace ensemble (#338 and !322)
 * Add program ensPostProcess for recentering and inflating letkf analysis ensemble (#329 and !314)
 * Add ability to recenter LETKF analysis ensemble on a supplied EnVar analysis (#328 and !313)
 * Addition of the OmP error std dev OBS_OMPE element and its possible use during background check (#246 and !311)
 * Add vertically varying horizontal localization length scale for LETKF (#322 and !305)
 * Add option `fullyUseExtremeTimeBins` (default is `.false.`) to namelist section NAMTIME (#323 and !304)
 * Add ability to "clean" a burp file by removing observations flagged for thinning (#319 and !302)
 * The variable names O3L, CH4L and N2OL included in GEM as of v5.1.a9 are now recognized (#321 and !300)
 * Add check on mpi imbalance of `gridStateVector` object, abort if too imbalanced or allow imbalance with
   `abortOnMpiImbalance` namelist variable in `namstate` (#312 and !298)

### Changed

 * Using `oavar_scripts` to version `2.2.6` which includes the following changes since version `2.2.4`:
       * On a generalise les scripts pour tourner les programmes du LETKF.
         * Ces changements sont compatibles arriere.
       * `oavar.mpirun`: set `TBB_MALLOC_USE_HUGE_PAGES=1` on `sles-15-*`
       * `oavar.var_mpi`: On corrige le mode `fasttmp=no` pour éviter que les
         fichiers complets se retrouvent dans le meme répertoire que les
         fichiers splittés.
       * `oavar.launch`: Ajout d'une cle `-analinc_mask`
 * Update to `rpn/utils/19.5.1` and `cmda/utils/19.5-3` (#339 and !321)
 * Unify the `bgckMW` program for AMSUA and ATMS QC (#308 and !315)
 * Improve memory usage in `stateToColumn_mod` (#306 and !310)
 * The `midas-prepcma` program is now mpi (#325 and !309)
 * Diag sqlite files now produce unique id_obs, id_data values across mpi tasks (#318 and !303)
 * For consistency with EnKF: (#300 and !293)
   * allow vertical interpolation of HU instead of log(HU) for radiance computation
   * optionally change raobs topo filter
   * read AMV obs error from obsfile
   * apply humidity adjustment before random pert

### Fixed

 * Bug fix to avoid overflow in the computation of secant of sat zenith angle when variable is missing in the obs file (#337 and !318)
 * Fix bugs for the `genCoeff` program and application of bias correction in `var` and `oMinusF` (#330 and !316)
 * Fix global mode for `calcStats` program, which was not working, and make it mpi (#307 and !308)
 * Bug fix for the `genCoeff` program: radiosonde weighting was not working properly (#323 and !304)
 * Fix recently introduced bug in adjoint of HU vertical interp for radiance computation (#314 and !299)

## [3.5.0]

### Highlights

 * Satellite radiance bias correction: apply bias correction and estimate bias correction coefficients
 * Slant-path interpolation operator for use with radiance observations
 * Implementation of LETKF, both standard approach and original approach with cross-validation: `midas-letkf`
 * New observations can be assimilated: MWHS2 radiances, radar-derived precipitation, SAR wind speed
 * Ability to use model ozone, instead of climatology, for RTTOV
 * 2D analyses now possible: near surface atmosphere, sea-ice and SST

### Added

 * Add ability to apply radiance bias correction and estimate bias correction coefficients in the MIDAS framework (#210 and !279)
 * Add missing functionality to `midas-prepcma` program (#260 and !289)
 * Add stochastic LETKF with cross validation algorithm (#276 and !282)
 * An SSM domain is created if a tag is pushed (#292 and !281)
 * Add ability to assimilate of SAR wind speed (#218 and !272)
 * Add ability to assimilate MWHS2 data (#287 and !274)
 * Add ability to read and use 2D-fields of correlation lenth scale and background STD for diffusion B matrix (#274 and !270)
 * Add lake operator (so far only for CIS lake ice obs) for horizontal interpolation from the grid to the observation location (#271 and !256)
 * Add ability to assimilate (log-transformed) precipitation in EnVar and LETKF (#267 and !252)
 * Add ability to use ozone profiles from trial field instead of climatology. Controlled by namelist variable `useO3Climatology`. (#195 and !246)
 * Add ability to use multiple instances of `bMatrixEnsemble_mod` are now possible which enables e.g. scale-dependent localization with spectral localization (SDLwSL) (#198 and !242)
 * Add many new features to the `midas-letkf` program, including a new cross-validation algorithm (#249 and #262, !241)
   * Also includes Yin-Yang grid compatibility and additional procedures for quality control, data selection and modification of obs error to facilitate comparison with the current EnKF
 * Add new variables to varnamelist_mod useful for `midas-ensManip` to compute mean and stddev (HR,TD,PN,PR,I2,I3,I4,I5,I6,I8,DN,FB,FI) (#236 and !240)
 * Implementation of the slant-path radiative transfer for the radiance observations, on background and analysis increment states (#243 and !238).
   * Two new namelists are added: `nams2c` activates the slant-path calculation for background and/or analysis increment states, `namSlantPath` defines the parameters for iterations to resolve the slant line-of-sight.
 * Add first implementation of the Local Ensemble Transform Kalman Filter in MIDAS (#245 and !233)
 * Add footprint operator (so far only for sea ice obs) for horizontal interpolation from the grid to the observation location (#197 and !222)
 * The environment variable `MIDAS_MAKE_LINKS_MACHINE_LIST` can be
   used to control the hosts on which links will be created by
   `install_suite.sh` in the maestro test suite.  By default, only the
   links on which the suite will run are created. (#231 and !216)
 * The logical namelist variable `ltopofilt` has been removed. Note: you probably **must update your namelist** (#225 and !211)
   * The new namelist variable is called `list_topoFilt`. This string array variable allow to activate the topographic rejection criteria for selected observation families. See the namelist in the unit tests from examples.
 * Add ability to define a local domain and control inclusion of each variable for energy norm (#207 and !204)
 * Add the program `midas-prepcma` to reproduce the similar program in the EnKF codebase (#189 and !198)
 * 3DVar analysis of SST data (family 'TM') can now be computed without MPI (#203 and !195)
 * The scripts to build the MIDAS SSM domain are now in the MIDAS
   depot (#187 and !186).  See the [README](README.md) for more information.
 * A column `OBS_CRPS` has been added to `obsSpaceData` (#185 and !188).
 * `bMatrixEnsemble_mod` can now read an ensemble of perturbations like lagged forecast differences (#193 and !190)
 * Variance smoothing is now possible in `bMatrixEnsemble_mod` (#193 and !190)
 * ScaleFactor added to `midas-ensManip` for vertically scaling ensemble perturbations when recentering (#186 and !179)
 * Minor changes to extend capability of generating OmP (and OmA) diagnostics via the CH obs family (#184 and !177)
 * Visibility and wind gust near the surface can now be assimilated (#173 and !176)
 * A sea ice concentration analysis can now be done with CIS daily ice charts (other obs types to come) (#163 and !175)
 * A horizontal land/sea mask can now be included in the analysisgrid file (#163 and !175)
 * The program, `midas_obsSelection`, was created (so far only for aladin HLOS obs), comprising O-P computations, background check, and thinning (#113 and !174)
 * Enable `vcode=5005` for ensemble B matrix (#188 and !173)
 * The aladin HLOS wind observations can now be adjusted to compensate for the observations having been calculated at another meteorolgical facility (#139 and !172)
   * A new namelist, `NAMALADIN_OBS`, is required when there are height-level observations (i.e. aladin) data to be treated.
 * Able to write contents of `obsSpaceData` to simplified sqlite files (useful when using another input file format and sqlite wanted for diagnostics) (#167 and !166)
 * Add observation operator and background check for aladin (HLOS winds) (#114 and !163)
 * Level-dependant steering flow scalefactor capability for advection (#168 and !146)
 * Add checks on humidity limits in `midas-ensManip` (#164 and !143)
 * New script to automatically generate module dependencies: `make_src_files.sh` (#149 and !136)
 * Removal of constraints on spectral truncation and number of levels relative to the MPI topology (#135 and !135)
 * New functionality: now able to read various types of sea ice data (family =*GL*) (#127 and !131)
 
### Changed

 * Update the expected execution timings for the ones on HPCR-U1 (#298 and !291)
   * On XC50, we use `TBB_MALLOC_USE_HUGE_PAGES=1` to save 7% in execution freely.w
   * Some users are notified when a timing outlier is found.
 * The GitLab runner which runs the continuous integration process is now executed under user `sanl888` (#288 and !277)
 * Also using `eccc/mrd/rpn/anl/rttov/12v1.4` which have been compiled
    with `code-tools/01.3` (#275 and !257)
 * The automated (CI) system tests now runs on both available platforms: HPCR-U0 and HPCR-U1 (#270 and !250)
 * The observation variable transforms was generalized and modernized (#247 and !239).
    * Previously only wind (speed,direction) -> (u,v) was supported. Now users can add (relatively) easily any type of variable transform. For now, only visibility -> log(visibility) was added.
    * An observation transform is activated when an assimilated observation is not found in the observations read from burp or sqlite files.
    * On outputs, the o-p and/or o-a of the transformed variables are converted and added to the original/source variable. In the output burp files, only the source variable appears because it was found too difficult to modify `brpr_updateBurp`. However, both source and transformed variable info are written to the sqlite files.
 * The height/pressure are computed on the grid, before horizontal interpolation to observation locations, to prepare for using slanted columns and footprint operators. (#124 and !220)
    * The height/pressure are part of statevector/columndata main data storage arrays (gd_r4/gd_r8/all) and are calculated for the trial fields and the increments.
    * Allocation of height/pressure is set to true, by default, and it is done if the necessary variables for their calculation are available in the statevector/columndata.
    * Z_M/Z_T and P_M/P_T are the height and pressure on grid on the TH and MM levels in `varNameList_mod`.
    * dPdPsfc is no longer used in any observation operators since the increment of pressure is calculated on the grid and is interpolated to the observation location.
    * Change namelist variable `addGZsfcOffset` to `addHeightSfcOffset`.
    * Variable/function/subroutine names that include `gz` are changed to `height` to reflect the fact that geometric altitude/height is now the primary variable instead of geopotential.
    * Memory requirements are higher for some programs and configurations (but not gdps and rdps configurations). 
    * The execution time is also increased for some (e.g. gdps takes ~100 seconds longer, but this can be reduced by increasing number of nodes to 30 or 36).
 * Minor change in `tt2phi_mod` (slightly affects results): now setting near-sfc temperature and momentum altitude levels to their known height offset (#180 and !212).
 * The namelist variable `scaleFactor` in `NAMBHI` must now be specified in all 3DVar configurations because default value was changed from 1.0 to 0.0. (#224 and !209)
 * CalcStats in LAM mode was made MPI compatible (#158 and !202)
 * Replacing the old numerical recipe for generating gaussian random values by a much more efficient method (#82 and !192)
   [`random_tools`](https://gitlab.com/mfvalin/random_tools).
    * This is changing the results only for the program `midas-randomPert`.
 * Changes to the height (GZ) calculation within MIDAS so that GPS-RO and GPS-gb now use same heights as other obs types (#141 and !191)
    * The variable was also changed from geopotential to altitude
    * New namelist variable to allow using static GPS-RO observation error variance
    * Small change to H(x) results due to these modifications for GPS-RO, GPS-gb, GZ and possibly sfc obs
 * When compiling with `COMPILE_MIDAS_ADD_DEBUG_OPTIONS=yes`, the options `-debug -check all -O 0` are added to the compile command (#182 and !187)
 * The body of the program `midas-oMinusF` was extracted into a new module, `ominusf_mod` (#113 and !174)
 * **Major change:** New approach for horizontal/temporal interpolation of background state and increment to observation locations/times (#80 and !147)
    * Changes to how the background state is read; now relies on `gridstatevector_mod`
    * File copy to ramdisk for all files can be done by fortran code (simplifies scripts)
    * Linear time interpolation can be applied to the background state when computing the innovations (controlled by namelist variable, default is nearest neighbour)
    * Small change to results for most applications and some increases in time and memory requirements
 * Use constituent BUFR elements from official tableburp file released 30 Sept 2018 (#150 and !137)
 * The input file trlp is no longer necessary in varbc mode (Slight unsignificant change in the results of varbc) (#145 and !129)

### Fixed

 * Add missing O-P background check for surface (2m) dewpoint depression ES (#294 and !288)
 * Fix bug triggered by GB-GPS reports with missing ZTD in background check of conventional observations (#295 and !287)
 * Correction to the surface humidity (BUFR element 13214) writen in BURP files info block during hyperspectral IR background check. Due to the error introduced in !111 exp(HU) was written instead of HU. (#257 and !247)
 * Fixed bugs, compilation procedure and system tests to allow compatibility on new machines (ppp3/4, banting/daley) - note, only the latest version of the reference results are available on the new machines (#258 and !243)
 * Correction of the conversion factor used to compute air mass predictors in `biascorrection_mod`. Problem introduced in !191. (#219 and !219)
 * Fixed some potential bugs detected while compiling with `-check all` (#182 and !187)
 * Improved efficiency of ensemble amplitude memory access and writing of `rehm` and `anlm` files (#170 and !153)
 * Fix the selection of GPSRO-bending angle observations (#151 and !145)

## [3.4.2]

 * Using `rpn/libs/19.5` instead of `rpn/libs/19.4` (no impact of the results) (#272 and !253)

## [3.4.1]

 * Using `rpn/libs/19.4` instead of `rpn/libs/19.2` (no impact of the results) (#269 and !251)

## [3.4.0]

### Fixed

 * Fixed bugs, compilation procedure and system tests to allow
   compatibility on new machines (ppp3/4, banting/daley) - note, only
   the latest version of the reference results are available on the
   new machines (#234 and !244)
	* This work is aimed to be introduced in a new release branch
	`v_3.4` which origins from `v_3.3.5`.

## [3.3.5]

### Added

 * Added a test 'var/EnVar/geps' to check the configuration used in the Operational ENKF (#241 and !229)

### Fixed

 * Skip update instead of abort in the case where one of the BURP input files contains no valid data (#244 and !231)

## [3.3.4]

### Fixed

 * Reject unknown satellites instead of aborting (#221 et !221)
 * Correction of a bug in the update of cloud parameters and
   emissivity in IR bgcheck mode that was affecting CrIS FSR. (#240 and !223)
 * Adding two files for instrument CrIS FSR for observations inputs in
   case we receive more than 2 satellites.  We already introduced:
    * `obscrisfsr1` and `obscrisfsr2` recently and we add
    * `obscrisfsr3` and `obscrisfsr4`.
 * Bug fix for the BURP cloud parameters and emissivity update bugfix
   above to handle properly missing data cases (#240 and !225)

## [3.3.3]

### Added

 * The analysis increment can now be masked in the blending zone like in the former program `addAnalInc` (#213 and !200)
 * Adding support for CrIS FSR radiances (#205 and !193)
   * It needs a modified version of RTTOV-12.
   * New input files were added to work around the file size limit of
     BURP:
     * `brpcrisfsr`, `brpcrisfsr1`, `brpcrisfsr2`
     * `obscrisfsr`, `obscrisfsr1`, `obscrisfsr2`
 * With GOES-R, it will be ABI instrument instead of `goesimager`.
   The code was modified to account for that. (#211 and !196)

### Removed

 * A quality control test specific to CrIS was forgottten for CrIS FSR (#232 and !219)

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
* First step towards controlling precision (#47 and !46)
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
* The directories have been reorganized (#50 and !42)
* Relax constraints on MPI topology (#49 and !44)
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

[Unreleased]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.5.2...HEAD
[3.5.2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.5.1...v_3.5.2
[3.5.1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.5.0...v_3.5.1
[3.5.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.4.2...v_3.5.0
[3.4.2]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.4.1...v_3.4.2
[3.4.1]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.4.0...v_3.4.1
[3.4.0]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.5...v_3.4.0
[3.3.5]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.4...v_3.3.5
[3.3.4]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.3...v_3.3.4
[3.3.3]: https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/compare/v_3.3.2...v_3.3.3
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
