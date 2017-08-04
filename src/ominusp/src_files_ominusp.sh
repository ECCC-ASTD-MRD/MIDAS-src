#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="utilities_mod.ftn90 mpi_mod.ftn90 codtyp_mod.ftn90 toplevelcontrol_mod.ftn90 ramdisk_mod.ftn90 bufr_mod.ftn90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpivar_mod.ftn90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90 obsspacedata_mod.ftn90 analysisgrid_mod.ftn90"
SRC_FILES="$SRC_FILES localizationfunction_mod.ftn90 controlvector_mod.ftn90 varnamelist_mod.ftn90 timecoord_mod.ftn90"
SRC_FILES="$SRC_FILES obsspacedata_mod.ftn90 verticalcoord_mod.ftn90 rmatrix_mod.ftn90 hirchannels_mod.ftn90 tovs_nl_mod.ftn90"
SRC_FILES="$SRC_FILES gridstatevector_mod.ftn90 lqtoes_mod.ftn90 presprofileoperators_mod.ftn90"
SRC_FILES="$SRC_FILES windrotation_mod.ftn90 columndata_mod.ftn90 tt2phi_mod.ftn90 tovs_lin_mod.ftn90 multi_ir_bgck_mod.ftn90 randomnumber_mod.ftn90"
SRC_FILES="$SRC_FILES statetocolumn_mod.ftn90 ozoneclim_mod.ftn90 tovs_extrap_mod.ftn90 burpread_mod.ftn90 gps_mod.ftn90 globalspectraltransform_mod.ftn90"
SRC_FILES="$SRC_FILES obssubspacedata_mod.ftn90 burpfiles_mod.ftn90 chem_obserrors_mod.ftn90 obserrors_mod.ftn90 chem_setup_mod.ftn90 bmatrixchem_mod.ftn90 chem_obsoperators_mod.ftn90 obsfilter_mod.ftn90 obsoperators_mod.ftn90" # obsspacediag_mod.ftn90"
SRC_FILES="$SRC_FILES innovation_mod.ftn90 residuals.ftn90 emissivities_mod.ftn90"
SRC_FILES="$SRC_FILES tovs_allocate_transmission.ftn90 tovs_fill_profiles_tl.ftn90 tovs_rttov.ftn90 tovs_calc_jo.ftn90 tovs_fill_profiles.ftn90 tovs_rttov_ad.ftn90 tovs_deallocate_transmission.ftn90 tovs_fill_profiles_ad.ftn90  tovs_rttov_tl.ftn90 compute_HBHT.ftn90 compute_HBHT_ensemble.ftn90 compute_HBHT_static.ftn90 compute_HBHT_static_chem.ftn90"
