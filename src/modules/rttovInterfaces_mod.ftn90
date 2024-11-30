
module rttovInterfaces_mod
  ! MODULE rttovInterfaces_mod (prefix='' category='9. Global interfaces')
  !
  ! :Purpose: To *include* all of the needed files containing the interfaces
  !           for the rttov library subroutines. This allows all other 
  !           files in MIDAS to be of the *f90* type (instead of *ftn90*),
  !           thus avoiding a pass through the pre-processor.
  !
implicit none
public

#include "rttov_coeffname.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_alloc_tl.interface"
#include "rttov_alloc_ad.interface"
#include "rttov_alloc_k.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_opts.interface"
#include "rttov_get_emis.interface"
#include "rttov_calcemis_ir.interface"
#include "rttov_setup_emis_atlas.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_parallel_ad.interface"
#include "rttov_parallel_k.interface"
#include "rttov_direct.interface"

#include "rttov_scatt.interface"
#include "rttov_scatt_tl.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"

#include "rttov_init_prof.interface"
#include "rttov_init_transmission.interface"
#include "rttov_init_rad.interface"
#include "rttov_iniscatt.interface"
#include "rttov_eddington.interface"
#include "rttov_errorreport.interface"
#include "rttov_calcbt.interface"
#include "rttov_scatt_emis_terms.interface"

end module rttovInterfaces_mod
