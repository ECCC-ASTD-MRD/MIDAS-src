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
!! *Purpose*: Compute background error stddev in observation space.
!!
!! Revision:   
!!v            M. Sitwell (ARQI/AQRD)  May 2015
!!v               - Added call to compute_HBHT_static_chem for computation
!!v                 with B_chm for chemical constituents
!!v            Y.J. Rochon (ARQI/AQDR) March 2016
!!v               - Allowed static=ensemble=.false. when static_chm=.true.
!!v               - Allowed possibility of static_chm=.true. and ensemble=.true.
!!v                 with and without static=.true.
!--------------------------------------------------------------------------
subroutine compute_HBHT(columng,columnhr,obsSpaceData)

  use mpivar_mod
  use obsSpaceData_mod
  use columnData_mod
  use utilities_mod
  implicit none

  type(struct_obs)        :: obsSpaceData ! Observation-related data
  type(struct_columnData) :: columng      ! Columns of the background interpolated 
                                          ! to analysis levels and to obs horizontal locations
  type(struct_columnData) :: columnhr     ! Columns of the background interpolated 
                                          ! to obs horizontal locations

  real(8) :: HBHT_static, HBHT_ensemble, HBHT_hybrid

  logical :: static   = .false.
  logical :: ensemble = .false.
  logical :: static_chm = .false.

  integer :: index_body
  integer :: fnom, fclos, ierr, nulnam

  character(len=12) :: hybrid_mode
  character(len=2)  :: fam

  !namelist
  NAMELIST /NAMHBHT/hybrid_mode

  !
  !- 1.  Compute HBHT (sigma_b in observation space)
  !

  !- 1.1 HBHT from the Bnmc matrix
  call compute_HBHT_static( columng,      & ! IN
                            columnhr,     & ! IN
                            obsSpaceData, & ! INOUT (HBnmcHT std. dev. outputted in OBS_HPHT)
                            static        ) ! OUT   (Active if TRUE)

  !- 1.2 HBHT from the Bens
  call compute_HBHT_ensemble( columng,      & ! IN
                              columnhr,     & ! IN
                              obsSpaceData, & ! INOUT (HBensHT std. dev. outputted in OBS_WORK)
                              ensemble )      ! OUT   (Active if TRUE)

  !- 1.3 HBHT from the B_static matrix for chemistry
  call compute_HBHT_static_chem( columng,      & ! IN
                                 obsSpaceData, & ! INOUT ( sqrt(diag(H*B*H^T)) with B_static_chm outputted in OBS_HPHT )
                                 static_chm )    ! OUT   (Active if TRUE)

  !
  !- 2. Select/Blend HBHT
  !
  if ( (static.or.static_chm) .and. .not. ensemble ) then
     ! Bnmc only
     write(*,*)
     write(*,*) 'compute_HBHT: Using B_static ONLY'
     ! HBnmcHT std. dev. already in OBS_HPHT
  else if ( .not. static .and. ensemble .and. .not. static_chm ) then
     write(*,*)
     write(*,*) 'compute_HBHT: Using B_ensemble ONLY'
     ! Transfer HBensHT std. dev. values in OBS_WORK to OBS_HPHT
     do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_ensemble)
     end do
  else if ( (static.or.static_chm) .and. ensemble ) then
     ! Read Namelist first
     hybrid_mode = 'WEIGHTED_SUM' ! default value
     nulnam = 0
     ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
     read(nulnam,nml=namhbht,iostat=ierr)
     if ( ierr /= 0) call utl_abort('compute_HBHT: Error reading namelist')
     if ( mpi_myid == 0 ) write(*,nml=namhbht)
     ierr = fclos(nulnam)

     write(*,*)
     write(*,*) 'compute_HBHT: Using hybrid approach (blend of B_static and B_ensemble) in mode = ', trim(hybrid_mode)

     do index_body = 1, obs_numBody(obsSpaceData)
        fam = obs_getFamily(obsSpaceData,obs_bodyElem_i(obsSpaceData,OBS_HIND,index_body))
        if ( (trim(fam).eq.'CH'.and..not.static_chm) .or. (trim(fam).ne.'CH'.and..not.static) ) then
           HBHT_static = 0.0D0
        else
           HBHT_static = obs_bodyElem_r(obsSpaceData,OBS_HPHT,index_body)
        end if
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        select case ( trim(hybrid_mode) )
        case ('WEIGHTED_SUM')
           HBHT_hybrid = sqrt(HBHT_static**2 + HBHT_ensemble**2)
        case ('MAX_VALUE')
           HBHT_hybrid = max(HBHT_static,HBHT_ensemble)
        case default
           write(*,*)
           write(*,*) 'compute_HBHT: Unknown hybrid_mode ', trim(hybrid_mode)
           call utl_abort('compute_HBHT')
        end select
        call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_hybrid)
     end do

  else
     call utl_abort('compute_HBHT: no B matrix was initialized')
  end if

end subroutine compute_HBHT