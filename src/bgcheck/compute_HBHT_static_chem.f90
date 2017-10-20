SUBROUTINE COMPUTE_HBHT_STATIC_CHEM(lcolumng,lobsSpaceData,active)
!
! Author  : M. Sitwell, May 2015
!
! Purpose: To compute the background error standard deviations in
!          observation space, sqrt(diag(H*B_static*H^T)).
!
! Arguments:
!
!  Input
!
!           lcolumng             column at observation location
!
!  Inout:
!
!           lobsSpaceData        observation space data, output saved in OBS_HPHT column
!
!  Output:
!           active               flag to indicate if chemical consituents are to be used
!
! Revision:
!
!-----------------------------------------------------------------------------------------

  use mpivar_mod
  use codtyp_mod
  use obsSpaceData_mod
  use columnData_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use bmatrixchem_mod
  use chem_obsoperators_mod

  implicit none
  
  type(struct_hco), pointer :: hco_anl
  type(struct_vco), pointer :: vco_anl
  type(struct_obs)        :: lobsSpaceData
  type(struct_columnData) :: lcolumng
  logical                 :: active
      
  integer :: cvdim

  
  !- Get the appropriate Horizontal and Vertical Coordinate
  hco_anl => agd_getHco('ComputationalGrid')
  vco_anl => col_getVco(lcolumng)
  
  call bchm_setup( hco_anl,vco_anl, &  ! IN
                   cvdim, &            ! OUT
                  'BackgroundCheck' )  ! IN

  active = bchm_is_initialized()
  
  if (active) then
     write(*,*)
     write(*,*) 'Computing H*B*H^T using B_static_chm - Start'
  else
     if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT_static_chem: option NOT ACTIVATED'
     return
  end if
          
  call chm_observation_operators(lcolumng,lobsSpaceData,kmode=1) ! kmode=1 for background check to compute HBH^T
  
  write(*,*)
  write(*,*) 'Computing H*B*H^T using B_static_chm - End'
  
  RETURN
END SUBROUTINE COMPUTE_HBHT_STATIC_CHEM