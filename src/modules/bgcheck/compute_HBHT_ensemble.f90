!--------------------------------------- LICENCE BEGIN -----------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

!--------------------------------------------------------------------------
!! *Purpose*: Compute background error stddev in observation space using
!!            ensemble-based statistics.
!!
!--------------------------------------------------------------------------
subroutine compute_HBHT_ensemble(columng,columnhr,obsSpaceData,active)
  use mpivar_mod
  use obsSpaceData_mod
  use columnData_mod
  use timeCoord_mod
  use bMatrixEnsemble_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use gridStateVector_mod
  use stateToColumn_mod
  use obsOperators_mod
  implicit none

  type(struct_obs)        :: obsSpaceData ! Observation-related data
  type(struct_columnData) :: columng      ! Columns of the background interpolated 
                                          ! to analysis levels and to obs horizontal locations
  type(struct_columnData) :: columnhr     ! Columns of the background interpolated 
                                          ! to obs horizontal locations
  logical                 :: active

  type(struct_columnData) :: column
  type(struct_gsv)        :: statevector

  type(struct_hco), pointer :: hco_anl
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable :: HBHT_ens(:)

  integer :: memberIndex, index_body, cvdim

  !
  !- 1.  Initialization
  !

  !- 1.1 Get vertical and horizontal analysis grid attributes
  vco_anl => col_getVco(columng)
  hco_anl => agd_getHco('ComputationalGrid')

  !- 1.2 Initialize/Read the flow-dependent ensemble perturbations
  call ben_Setup( hco_anl,             & ! IN
                  vco_anl,             & ! IN
                  cvdim,               & ! OUT
                  'BackgroundCheck' )    ! IN

  if ( cvdim > 0 ) then
     write(*,*)
     write(*,*) 'Computing HBHT from ensemble perturbations - START'
     active = .true.
  else
     if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT_ensemble: option NOT ACTIVATED'
     active = .false.
     return
  end if

  !- 1.3 Create a gridstatevector to store the ensemble perturbations
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
       mpi_local_opt=.true.)

  !- 1.4 Create column vectors to store the ens perturbation interpolated to obs horizontal locations
  call col_setVco(column,vco_anl)
  call col_allocate(column,col_getNumCol(columng),mpiLocal_opt=.true.)
  call col_copyLatLon(columng,column)

  !- 1.5 Create a working a array to sum H ensPert HT
  allocate(HBHT_ens(obs_numBody(obsSpaceData)))
  HBHT_ens(:) = 0.d0

  !- 1.6
  call tim_timeBinning(obsSpaceData,tim_nstepobsinc)
  call tim_sutimeinterp(obsSpaceData,tim_nstepobsinc)

  !
  !- 2.  Compute HBHT from the ensemble perturbations
  !
  do memberIndex = 1, ben_getnEns()

     !- 2.1 Extract perturbations from the current memberIndex
     write(*,*)
     write(*,*) 'Reading ensemble perturbation from member = ', memberIndex
     call ben_getPerturbation( statevector,    & ! OUT
                               memberIndex,    & ! IN
                               'ConstantValue' ) ! IN

     !- 2.2 Interpolation to the observation horizontal locations
     call s2c_tl( statevector,           & ! IN
                  column,                & ! OUT (H_horiz EnsPert)
                  columng, obsSpaceData )  ! IN

     !- 2.3 Interpolation to observation space
     call oop_Htl( column, columng, & ! IN
                   obsSpaceData,    & ! OUT (Save as OBS_WORK: H_vert H_horiz EnsPert = H EnsPert)
                   1 )                ! IN

     !- 2.4 alpha * HBH^T = sum(OBS_WORK^2)
     do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_ens(index_body) = HBHT_ens(index_body) + &
                               (obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body))**2
     end do

  end do

  !- 2.5 Insert the standard deviations in OBS_WORK
  do index_body = 1, obs_numBody(obsSpaceData)
     call obs_bodySet_r(obsSpaceData,OBS_WORK,index_body,sqrt(HBHT_ens(index_body)))
  end do

  !
  !- 3.  Ending/Deallocation
  !
  deallocate(HBHT_ens)
  call col_deallocate(column)
  call gsv_deallocate(statevector)

  write(*,*)
  write(*,*) 'Computing HBHT from ensemble perturbations - END'

end subroutine compute_HBHT_ensemble
