!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!! MODULE chem_obsoperators_mod (prefix='chm' category='4. Observation operators')
!!
!! *Purpose*: Provides observation operator routines.
!!
!! @author Mike Sitwell and Yves Rochon (ARQI/AQRD)
!!
!! Public routines:
!!v       - "chm_obsoperators": Applies observation operators.
!--------------------------------------------------------------------------
module chem_obsoperators_mod

  use utilities_mod
  use chem_setup_mod
  use obsSubSpaceData_mod
  use bmatrixchem_mod
  use physicsfunctions_mod
  use MathPhysConstants_mod
  use earthconstants_mod
  use bufr_mod

  ! Following modules needed only by chm_observation_operators 
  ! and or chm_obsoper_init
  
  use mpi_mod, only: mpi_myid
  use codtyp_mod
  use obsSpaceData_mod
  use columnData_mod
  use gridStateVector_mod
  use stateToColumn_mod
  use verticalCoord_mod
  use HorizontalCoord_mod
  use timeCoord_mod
  use varNameList_mod
   
  implicit none
  private

! public procedures
! -----------------

  public :: chm_observation_operators

! Arrays for integration upper boundary of retrieved total column measurements 
  type(struct_oss_obsdata) :: chm_column_boundary

! Arrays for background error std. dev.
  type(struct_oss_obsdata) :: chm_sigma_trial


contains

!------------------ CONTROL ROUTINES FOR OBSERVATION OPERATORS ---------------------
!-----------------------------------------------------------------------------------


!--------------------------------------------------------------------------
!! *Purpose*: Apply the observation operators for chemical constituents.
!!            Mode of operator set by kmode.
!!
!! @author M. Sitwell, Aug 2015
!!
!! Revisions: 
!!v          M. Sitwell, ARQI/AQRD Nov 2015
!!v          - Modified the flagging of observations after chm_obsoperators is called
!!v          M. Sitwell, ARQI/AQRD April 2016
!!v          - Moved most of the input to chm_obsoperators into obsoper
!!v          Y. Rochon, ARQI/AQRD May 2016
!!v          - Modified skipping over processing when not needed
!!
!! Input
!!
!!v   kmode           Mode of observation operator, with values
!!v                     0 for general (non-linear and linear) simulation operator
!!v                     1 for determination of sqrt(diag(H*B*H^T))
!!v                     2 for tangent linear operator
!!v                     3 for adjoint of tangent linear operator
!!
!!v   column_bkgrnd   Column of x_background interpolated to observation location. Can
!!v                   have the same vertical levels as the trial field (columnhr) or
!!v                   as the increment field (columng)
!!
!!v   obsSpaceData    Observation space data structure
!!
!! Output
!!
!!v   For kmode values of 0) OmP and total Jo(x_background) for CH. 
!!v                          OmP saved in OBS_OMP of obsSpaceData
!!v                       1) background error standard deviations in
!!v                          observation space, sqrt(diag(H*B_static*H^T)).
!!v                          Saved in OBS_HPHT of obsSpaceData
!!v                       2) Hdx saved in OBS_WORK of obsSpaceData
!!v                       3) H^T * R^-1 (OmP-Hdx) in columnInc_opt
!!
!!v   jobs            Optional output of total Jo(x_background) for chemical constituents. 
!!v                   Required for kmode=0 and not provided otherwise.
!!
!! Inout
!!
!!v   columnInc_opt   Optional argument for input/output of column of increment (column).
!!v                   For kmode=2, used as input for increment H_horiz dx interpolated
!!v                   to observation location. For kmode=3, used as output for
!!v                   H^T * R^-1 (OmP-Hdx). Required for kmode=2,3.
!!
!!  Comments:
!!v     - See type struct_chm_obsoperators in chem_setup_mod.ftn90 for 
!!v       description of obsoper elements.
!!v     - Currently can only handle the case when nlev_bkgrnd == nlev_inc
!!v     - Two equivalent methods for looping over a report body.
!!
!!v       Method 1:
!!
!!v            call obs_set_current_body_list(obsSpaceData,headerIndex)
!!v            BODY: do
!!
!!v               bodyIndex = obs_getBodyIndex(obsSpaceData)
!!v               if (bodyIndex < 0) exit BODY1
!!
!!v               ... obs_bodyElem_r(obsSpaceData, ... ,bodyIndex)
!!  
!!v            enddo BODY
!!
!!v       Method 2:
!!
!!v            bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
!!v            bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1
!!v            do  bodyIndex=bodyIndex_start,bodyIndex_end
!!v               ... obs_bodyElem_r(obsSpaceData, ... ,bodyIndex)   
!!v            end do
!!
!--------------------------------------------------------------------------
  subroutine chm_observation_operators(column_bkgrnd,obsSpaceData,kmode,columnInc_opt,jobs_opt)

    implicit none
    
    ! Subroutine arguments

    type(struct_columnData), intent(inout) :: column_bkgrnd
    type(struct_columnData), intent(inout), optional :: columnInc_opt
    type(struct_obs), intent(inout) :: obsSpaceData
    integer, intent(in) :: kmode
    real(8), intent(out), optional :: jobs_opt

    ! Local variables
    
    real(8) :: zomp,zinc
    integer :: unit,ier
    integer, external :: fclos

    ! Obs space local variables

    integer :: headerIndex,bodyIndex,bodyIndex_start,bodyIndex_end
    integer :: icodtyp,iobslev,nobslev,varno
    character(len=12) :: stnid

    integer, allocatable :: ixtr(:),iass(:),flag(:)
    logical, allocatable :: success(:),process_obs(:)

    ! Model space profile local variables

    real(8), allocatable :: obs_col(:)
    real(8), pointer :: col(:),model_col(:)
    integer :: nlev_bkgrnd,nlev_inc,imodlev
    character(len=2), parameter :: varLevel = 'TH'
    
    type(struct_chm_obsoperators) :: obsoper

    if ((kmode.eq.2.or.kmode.eq.3) .and. (.not.present(columnInc_opt))) then
       write(*,*) "chm_observation_operators: columnInc_opt must be specified for kmode = ",kmode
       call utl_abort("chm_observation_operators")
    end if
    
    ! Initializations
        
    if (present(jobs_opt)) jobs_opt = 0.d0

    nlev_bkgrnd = col_getNumLev(column_bkgrnd,varLevel)
    
    ! Allocate memory for model_col. Not necessary for kmode=0 since model_col points to obsoper%trial.
    select case(kmode)
    case(2)
       nlev_inc = col_getNumLev(columnInc_opt,varLevel)
       allocate(model_col(nlev_inc))
    case(1,3)
       allocate(model_col(nlev_bkgrnd))
    end select

    ! Allocations outside chm_obsoper_init since this can be done outside the HEADER loop.
    ! See chm_obsoper_init for assignment of array content.
    
    ! Model obs background, GZ, TT, and HU profiles.
    allocate(obsoper%trial(nlev_bkgrnd),obsoper%gz(nlev_bkgrnd),obsoper%tt(nlev_bkgrnd),obsoper%hu(nlev_bkgrnd))
    ! Model PP and pressure model layer boundaries taken as the middle between model levels.
    allocate(obsoper%pp(nlev_bkgrnd),obsoper%vmodpress(nlev_bkgrnd+1))
    ! Work array: derivative for any variable transform to be applied to a profile
    allocate(obsoper%dtransform(nlev_bkgrnd))
    ! Work array: integration weights array associated to the second order Lagrangian interpolation - see routines
    ! chm_vertintg* in chem_obsoperators_mod.ftn90
    allocate(obsoper%vweights(nlev_bkgrnd,nlev_bkgrnd))
    
    ! Loop over all header indices of the 'CH' family:
    
    call obs_set_current_header_list(obsSpaceData,'CH')
    HEADER: do

      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
  
      icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (icodtyp.ne.codtyp_get_codtyp('CHEMREMOTE').and.icodtyp.ne.codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
      stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
      bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1

      ! Set number of obs profile elements by removing count of BUFR_SCALE_EXPONENT elements
      nobslev = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)  
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).eq.BUFR_SCALE_EXPONENT) nobslev = nobslev-1
      end do

      ! varno is expected to be the same for all profile points where OBS_VNM value .ne. BUFR_SCALE_EXPONENT
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
            varno = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
            exit
         end if
      end do

      ! Allocate memory for remaining profile data not in obsoper
      allocate(obs_col(nobslev),success(nobslev),ixtr(nobslev),iass(nobslev),process_obs(nobslev),flag(nobslev))

      ! Check to see if background error variances available
      if (kmode.eq.1) process_obs(:) = bchm_StatsExistForVarName(vnl_varnameFromVarnum(varno,obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex)))
 
      ! Prepare for checking if any processing is needed according to initial flag values     
      iobslev=0
 
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then

            iobslev=iobslev+1

            ixtr(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) ! indicates if obs extends outside model profile vertical range
            iass(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) ! indicates if obs is to be assimilated
            flag(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex) ! observation integer flag
               
            ! Indicates if this obs should be processed by chm_obsoperators
            if (kmode.eq.1) then
               process_obs(iobslev) = ixtr(iobslev).eq.0.and.iass(iobslev).eq.obs_assimilated.and.process_obs(iobslev)
            else
               process_obs(iobslev) = ixtr(iobslev).eq.0.and.iass(iobslev).eq.obs_assimilated
            end if

         end if
      end do

      ! Initialize processing success flag
      success(1:nobslev) = process_obs(1:nobslev)
      
      if (all(.not.process_obs)) then

         ! All observations in the profile flagged so can skip obs operator for current measurement

         if (kmode.eq.3) then
            model_col(:) = 0.0D0
            obsoper%varName = vnl_varnameFromVarnum(varno,obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex))
         end if

      else  

         ! Initialize obsoper variables and allocate arrays
         call chm_obsoper_init(obsoper,obsSpaceData,headerIndex,column_bkgrnd,nlev_bkgrnd,nobslev,kmode,varno,stnid)
 
         ! Initialize model_col, dependent on kmode. Used for input for kmode=0,2, output for kmode=3.
         ! model_col represents for kmode 0) the horizontally interpolated background H_horiz(x_b)
         !                                1) not used
         !                                2) the analysis increment H_horiz dx
         !                                3) the result of applying the adjoint of H_vert 

         select case(kmode)
         case(0)
            model_col => obsoper%trial
         case(2)
            do imodlev=1,nlev_inc
               model_col(imodlev) = col_getElem(columnInc_opt,imodlev,headerIndex,obsoper%varName)
            end do
         case(1,3)
            model_col(:) = 0.0D0
         end select
               
         ! Loop over all body indices (profile elements) to aquire remaining data
      
         iobslev=0
         call obs_set_current_body_list(obsSpaceData,headerIndex)
         BODY1: do

            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY1
            
            ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements

            if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
               iobslev=iobslev+1
            else
               cycle BODY1
            end if

            ! Get vertical coordinate data. Valid for point data values in profiles.
            ! For layer data values, vertical coordinate data will instead be assigned within chm_obsoperators.

            obsoper%obslev(iobslev) = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

            ! Get normalized increment
            if (kmode.eq.3) then
               if (iass(iobslev).eq.1) then
                  obs_col(iobslev) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
               else
                  obs_col(iobslev) = 0.0D0
               end if
            end if

         enddo BODY1
      
         ! Apply observation operator. model_col,obs_col are the inputs/outputs in model,observation
         ! space, respectively. Other required inputs are in obsoper. Input/output is as follows:
         !
         !    kmode      model_col      obs_col
         !    -----      ---------      -------
         !      0           in            out
         !      1         not used        out
         !      2           in            out
         !      3           out           in

         call chm_obsoperators(obsoper,model_col,obs_col,kmode,ixtr,success)
 
      end if

      ! Output results
      
      if (kmode.eq.3) then
         ! Store H^T * R^-1 (OmP-Hdx) in columnInc
                     
         col => col_getColumn(columnInc_opt,headerIndex,obsoper%varName)
         col(1:nlev_bkgrnd) = model_col(1:nlev_bkgrnd)

      else
         ! Store results in obsSpaceData

         iobslev=0
         call obs_set_current_body_list(obsSpaceData,headerIndex)
         BODY2: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY2

            ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements
            if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
               iobslev=iobslev+1
            else
               cycle BODY2
            end if

            ! Check for success in calculations
            if (process_obs(iobslev).and..not.success(iobslev)) then
               ! Observation was flagged within this call of chm_observation_operators
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,0)
               call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
               call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),9) )
               cycle BODY2
            else if (iass(iobslev).eq.0) then
               ! Observation was flagged previous to this call of chm_observation_operators
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,0)
               call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
               cycle BODY2
            else if (.not.process_obs(iobslev) .and. kmode == 1) then
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               cycle BODY2
            end if

            ! Store result in appropriate location in obsSpaceData
            select case(kmode)
            case(0)
            
               ! Store OmP in OBS_OMP and add to Jo(x_background) of CH.   
                           
               zomp = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex) - obs_col(iobslev)
               call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,zomp)

               if (chm_diagn_only('CH',stnid,varno,nobslev,flag(iobslev))) then
                  ! Observation is for diagnostics and is not to be assimilated
                  call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,3)
               else if (present(jobs_opt).and.iass(iobslev).eq.1) then
                  ! Add to Jo contribution (factor of 0.5 to be applied outside report loop)
                  zinc = zomp/obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
                  jobs_opt = jobs_opt + zinc**2
               end if

            case(1)
            
               ! Background error standard deviations in
               ! observation space, sqrt(diag(H*B_static*H^T)), 
               ! saved in OBS_HPHT of obsSpaceData

               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,obs_col(iobslev))

            case(2)
            
               !   Store Hdx in OBS_WORK of obsSpaceData               
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,obs_col(iobslev))
            end select

         enddo BODY2

      end if

      ! Deallocate profile data
      deallocate(process_obs,success,ixtr,iass,obs_col,flag)
      call chm_obsoper_dealloc(obsoper)
  
    enddo HEADER

    deallocate(obsoper%trial,obsoper%pp,obsoper%tt,obsoper%dtransform)
    deallocate(obsoper%gz,obsoper%vmodpress,obsoper%vweights,obsoper%hu)
    if (kmode.ne.0) deallocate(model_col)
  
    if (present(jobs_opt)) jobs_opt = 0.5d0*jobs_opt

    ! Close message file
    ier=fclos(chm_setup_get_int('message_unit'))
    
  end subroutine chm_observation_operators

!--------------------------------------------------------------------------
!! *Purpose*: Initializes struct_chm_obsoperators variables and allocates arrays.
!!
!! @author M. Sitwell and Y. Rochon, April 2016
!!
!! Revisions:
!!v           J-F Caron, ARMA/MRD, Jan. 2018
!!v           - Account for change from LQ to Q for 'HU'
!!v 
!! Input
!!
!!v     column_bkgrnd   Column of x_background interpolated to observation location. Can
!!v                     have the same vertical levels as the trial field (columnhr) or
!!v                     as the increment field (columng)
!!v     obsSpaceData    Obs database struture
!!v     headerIndex     Measurement index in obsSpaceData (see main)
!!v     bodyIndex       Measurement element index in obsSpaceDate (see chm_obsoper_proceed)
!!v     kmode           0 for non-linear/linear model in assimilation (all models included are currently linear)
!!v                     1 for determination of sqrt(diag(H*B*H^T))
!!v                     2 for tangent linear model
!!v                     3 for adjoint model 
!!v     nobslev         Number of obs elements (see chm_obsoper_proceed)
!!v     nmodlev         Number of background field (model) levels
!!v     varno           BUFR number
!!v     stnid           Station ID
!!
!! InOut
!!
!!v     obsoper        Structure for constituents associated to obs (see struct_chm_obsoperators 
!!v                    in chem_obsoperators_mod.ftn90)
!!
!! Comments: 
!!v          - Allocation of arrays that are only dependent on nlev_bkgrd (nmodlev) have been moved
!!v            outside this subroutine so that they are only allocated once.
!--------------------------------------------------------------------------
  subroutine chm_obsoper_init(obsoper,obsSpaceData,headerIndex,column_bkgrnd,nmodlev,nobslev,kmode,varno,stnid)

    implicit none

    integer, intent(in) :: headerIndex,nmodlev,nobslev,kmode,varno
    character(len=12), intent(in) :: stnid
    type(struct_chm_obsoperators), intent(inout) :: obsoper
    type(struct_obs), intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: column_bkgrnd

    integer :: bodyIndex,jl,nmodlev_uv
    real(8), pointer :: col_ptr_gzb(:)
    real(8), allocatable :: uu(:),vv(:)
    character(len=2), parameter :: varLevel = 'TH'

    obsoper%nmodlev = nmodlev
    obsoper%nobslev = nobslev
    obsoper%obs_index = headerIndex
    obsoper%varno = varno
    obsoper%stnid = stnid

    ! Get obs space info that are part of the profile header
    obsoper%date  = obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex)
    obsoper%hhmm  = obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)
    ! Constituent identifyer following local version of WMO GRIB Table 08046 (similar to BUFR Table 08043)
    obsoper%constituent_id = obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex)
    if (obsoper%constituent_id.gt.chm_var_maxnumber().and.(obsoper%constituent_id.lt.7000.and.obsoper%constituent_id.ge.0)) then
        ! chm_constituents_size=chm_var_maxnumber() <7000 as values >=7000 (and <0) restricted to NWP fields.
        write(*,*) 'chm_obsoper_init: chm_constituents_size less than ',obsoper%constituent_id,' for STNID ',obsoper%stnid
        write(*,*) 'chm_obsoper_init: may need to increae chm_constituents_size'
        write(*,*) 'chm_obsoper_init: or edit BURP file values for Table 08046 element.'
        call utl_abort('chm_obsoper_init')                  
    end if

    obsoper%lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) 
    obsoper%lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) 

    ! Body info that we only need for first point in the profile
    bodyIndex       = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex) 
    ! Index of vertical coordinate type    
    obsoper%vco     = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)  
    ! Model field name (NOMVAR value)
    obsoper%varName = vnl_varnameFromVarnum(obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex),obsoper%constituent_id)

    ! Allocate arrays

    allocate(obsoper%obslev(nobslev))       ! Reference vertical levels
    allocate(obsoper%vlayertop(nobslev))    ! Layer tops for layer measurements
    allocate(obsoper%vlayerbottom(nobslev)) ! Layer bottoms for layer measurements
    allocate(obsoper%imodlev_top(nobslev))  ! Index of highest model level (lowest index) involved with obs element
    allocate(obsoper%imodlev_bot(nobslev))  ! Index of lowest model level (highest index) involved with obs element

    allocate(obsoper%zh(nobslev,nmodlev))   ! Local model operator H (excluding conversion constants and horizontal interpolation)
    allocate(obsoper%zhp(nobslev,nmodlev))  ! Part of zh that excludes aspects related to vertical resolition

    obsoper%vweights(:,:)=0.0D0            
    obsoper%zh(:,:)=0.0D0
    obsoper%zhp(:,:)=0.0D0
    obsoper%imodlev_top(:)=1
    obsoper%imodlev_bot(:)=nmodlev
    obsoper%dtransform(:)=1.0D0

    if (.not.col_varExist('TT')) then
       if (chm_required_field('TT',obsoper%varno)) then
          write(*,*) "chm_observation_operators: TT required for BUFR code ",obsoper%varno 
          call utl_abort("chm_observation_operators")
       end if
    end if

    ! Get background profiles at observation location
    do jl=1,nmodlev
       obsoper%pp(jl) = col_getPressure(column_bkgrnd,jl,headerIndex,varLevel)
       obsoper%trial(jl) = col_getElem(column_bkgrnd,jl,headerIndex,obsoper%varName)
       obsoper%tt(jl) = col_getElem(column_bkgrnd,jl,headerIndex,'TT')
    enddo

    if (col_varExist('TT').and.col_varExist('HU').and.col_varExist('P0')) then     
       ! GZ would have been generated in the call to sugomobs. 
       ! Convert from geopotential to geopotential height.
       col_ptr_gzb => col_getColumn(column_bkgrnd,headerIndex,'GZ','TH')
       obsoper%gz(1:nmodlev) = col_ptr_gzb(1:nmodlev)/RG
    else
       obsoper%gz(:) = -1.
    end if

    ! Get specific humidity if available
    if (col_varExist('HU')) then
       do jl=1,nmodlev
         obsoper%hu(jl) = col_getElem(column_bkgrnd,jl,headerIndex,'HU')       ! lnq was replaced by q
       enddo
    else
       obsoper%hu(:)=-1
    end if
    
    ! If applicable, get column upper boundaries for use with total column measurements when the related increment profile
    ! is to be restricted to the lower atmosphere (e.g. troposphere or PBL; 
    ! when m_setup_get_int('tropo_bound',:)>0 )
    if (obsoper%vco.eq.4.and.nobslev.eq.1.and.kmode.ne.1) then
       if (kmode.eq.0) then
          
          if (col_varExist('HU')) then
             if (col_varExist('UU').and.col_varExist('VV')) then
                nmodlev_uv=col_getNumLev(column_bkgrnd,'MM')
                allocate(uu(nmodlev_uv),vv(nmodlev_uv))
                do jl=1,nmodlev_uv
                   uu(jl) = col_getElem(column_bkgrnd,jl,headerIndex,'UU')
                   vv(jl) = col_getElem(column_bkgrnd,jl,headerIndex,'VV')
                enddo
                obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%gz,hu_opt=obsoper%hu,uu_opt=uu,vv_opt=vv)
                deallocate(uu,vv)
             else 
                obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%gz,hu_opt=obsoper%hu)   
             end if

          else
              obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%gz)
          end if

          call chm_add_col_boundary(headerIndex,obsoper%column_bound)  ! save boundary for kmode>0 calls using headerIndex

       else
          obsoper%column_bound = chm_retrieve_col_boundary(headerIndex)
       end if
    else
       obsoper%column_bound = -1.
    end if

  end subroutine chm_obsoper_init

!--------------------------------------------------------------------------
!! *Purpose*: Deallocate arrays for struct_chm_obsoperators.
!!
!! @author M. Sitwell, April 2016
!!
!! Comments:
!!v          - Deallocation of arrays that are only dependent on nmodlev have been moved
!!v            outside this subroutine so that they are only deallocated once.
!--------------------------------------------------------------------------
  subroutine chm_obsoper_dealloc(obsoper)

    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper

    if (allocated(obsoper%obslev))       deallocate(obsoper%obslev)
    if (allocated(obsoper%vlayertop))    deallocate(obsoper%vlayertop)
    if (allocated(obsoper%vlayerbottom)) deallocate(obsoper%vlayerbottom)
    if (allocated(obsoper%zh))           deallocate(obsoper%zh)
    if (allocated(obsoper%zhp))          deallocate(obsoper%zhp)
    if (allocated(obsoper%imodlev_top))  deallocate(obsoper%imodlev_top)
    if (allocated(obsoper%imodlev_bot))  deallocate(obsoper%imodlev_bot)

  end subroutine chm_obsoper_dealloc

!--------------------------- Routines for observation operators ----------------------------
!-------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------
!! *Purpose*: Apply observation operator for indicated observation data and condition.
!!
!! @author Y. Rochon, Feb 2015
!!
!! Input
!!
!!v     kmode  mode of observation observational operator
!!v              0) general (potentially non-linear) simulation operator
!!v              1) determination of sqrt(diag(H*B*H^T))
!!v              2) tangent linear operator
!!v              3) linear adjoint operator
!!
!! Inout
!!
!!v     obsoper    Structure that holds information needed by the observation operator
!!
!!v     obs_col    Observation space input/output profile
!!v                  kmode     input/output       profile
!!v                  -----     ------------       -------
!!v                    0           out            H(xb)
!!v                    1           out            sqrt(diag(H*B*H^T))
!!v                    2           out            H*dx
!!v                    3           in             R**-1 (Hdx-d)
!!
!!v     model_col  Model space input/output profile
!!v                  kmode     input/output       profile
!!v                  -----     ------------       -------
!!v                    0           in             xb
!!v                    1           not used       not used
!!v                    2           in             dx at obs location
!!v                    3           out            adjoint product H^T(...)
!!
!!v     ixtr       Flag indicating if obs within the model vertical coord range (.ne.0 for no) 
!!v                Can be modified internally - hence intent(inout) - even though
!!v                these changes will not be needed outside this routine.
!!
!!v     success   Indicates if the observation was successfully assimilated
!!
!! Revisions:
!!v           Y. Rochon, ARQI/AQRD, Feb. 2015
!!v           - Modifications of adaptations and update of BUFR elements
!!v           Ping Du, Mar 2015
!!v           - Finalization of Ht*grad contribution (case(3))
!!v           M. Sitwell, ARQI/AQRD, Mar 2015
!!v           - Modified to calculate whole profile within a single call
!!v           M. Sitwell, ARQI/AQRD, May 2015
!!v           - Added vertical interpolation operator
!!v           M. Sitwell, ARQI/AQRD, June 2015
!!v           - Changed calculation of HBH^T to forgo calculating off-diagonal elements
!!v             of the final product
!!v           M. Sitwell, ARQI/AQRD, April 2016
!!v           - Modified input arguments so that most inputs are passed through obsoper
!!
!! Further changes required for generalization:
!!
!!v 1) Add layer average operators.
!!v 2) Add AOD operators (summation over model layers).
!!v 3) Add option to include use of obs error correlation matrix for kmode=2,3 
!!v    (This may/will need to be done in oop_Hchm and oop_HTchm where the division 
!!v     by sigma_obs is applied. A new routine will be needed for this
!!v     operation - and others for reading the correlation matrices similarly to the
!!v     averaging kernels.)
!!
!! Comments:
!!v     - When kmode=0, call from chm_observation_operators passes model_col as a pointer to
!!v       obsoper%trial.
!!v     - Does not yet account for potential future applications of obs 
!!v       vertical correlation matrices.
!!v     - Potential specification of background error std. dev. (sigma_trial(:,2)) and correlation matrices 
!!v       for the ensemble-based and lam cases to be done when stats for these become in use with constituents.
!!v     - 'struct_chm_obsoperators' defined in chem_setup_mod.ftn.
!--------------------------------------------------------------------------
  subroutine chm_obsoperators(obsoper,model_col,obs_col,kmode,ixtr,success)
  
  implicit none

! Declarations

! Structure to hold observation operator information
    
  type(struct_chm_obsoperators), intent(inout) :: obsoper

! I/O arguments: obs space variables
    
  integer, intent(in) :: kmode 
  integer, intent(inout) :: ixtr(obsoper%nobslev)
  logical, intent(inout) :: success(obsoper%nobslev)

! I/O arguments: model space profile data and others

  real(8), intent(inout) :: model_col(obsoper%nmodlev), obs_col(obsoper%nobslev)

! Local variables
  
  real(8) :: press_obs(obsoper%nobslev),zwork(obsoper%nmodlev),unit_conversion(obsoper%nmodlev)
  real(8) :: sigma_trial(obsoper%nmodlev,2),temp_eff(1),zhmin
  integer :: iobslev,imodlev,jvar
  integer, parameter :: code_len=90
  character(len=code_len) :: code    ! Must be at least as large as chm_code_len
  real(8), allocatable :: avg_kern(:,:)

  if (code_len.lt.oss_obsdata_code_len()) call utl_abort('chm_obsoperators: Length of code string' &
                                          // ' needs to be increased to ' // trim(utl_str(oss_obsdata_code_len())))
       
! Determine if layer boundaries are assigned to this data source.
! If so, obtain them for use in this routine. 
! Routine provides obsoper%layer_identified,
!                  obsoper%vlayertop(nobslev), 
!                  obsoper%vlayerbottom(nobslev) 
 
  call chm_get_layer_boundaries(obsoper%stnid,obsoper%varno,obsoper%vco, &
       obsoper%nobslev,obsoper%pp(1),obsoper%pp(obsoper%nmodlev), &
       obsoper%layer_identified,obsoper%vlayertop,obsoper%vlayerbottom)
            
! Identify observation operator based on observation units and presence or
! not of layer boundaries

  if (bufr_IsIntegral(obsoper%varno)) then
     if (.not.obsoper%layer_identified) then
        write(*,*)   '----------------------------------------------------------'
        write(*,*)   'STNID, BUFR index, nobslev: ',obsoper%stnid,' ',obsoper%varno,obsoper%nobslev
        call utl_abort('chm_obsoperators: Required layer boundaries not available!')
     else
        ! Vertical integration operator
        obsoper%modelIndex=3
     end if
  else if (obsoper%layer_identified) then
      ! Layer averaging operator
     obsoper%modelIndex=4
  else if (obsoper%vco.eq.5.and.obsoper%nobslev.eq.1) then  
      ! Surface point (in-situ) measurement
      obsoper%modelIndex=2
  else    
      ! Vertical interpolation operator
     obsoper%modelIndex=1
  end if  

! Indicate if the generalized innovation operator is to be applied.

  obsoper%apply_genoper=.false.
  if (obsoper%constituent_id.ge.0) then
    if (chm_setup_get_int('genoper',obsoper%constituent_id).gt.0 .and. kmode.ne.1 .and. &
        (obsoper%modelIndex.eq.3 .or. obsoper%modelIndex.eq.4)) then
       
      if (kmode.eq.0) then
          ! Set reference profile for use with generalized innovation operator when kmode>=2
          call chm_set_reference_obsdata(obsoper)
  
          ! Get background error std dev profile at obs locations
          sigma_trial(:,:)=0.D0
          call bchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!          call blamchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!          call benschm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,2)) 

          call chm_add_sigma_trial(obsoper%obs_index,sigma_trial)

      else if (kmode.ge.2) then
          obsoper%apply_genoper = .true.
      end if
    end if
  end if

! Apply unit conversion and any necessary (nonlinear) transformations (apply unit conversion later for kmode=3)

  call chm_transform_profile(obsoper, computeDtransform_opt=kmode.ne.0) ! transformations on obsoper%trial

  select case(kmode)
  case(0)
     ! Perform transformation and unit conversion on obsoper%trial.
     ! Note that model_col => obsoper%trial for kmode=0.
     call chm_convert_units(obsoper,obsoper%trial)
  case(1)
     ! Save the conversion factor in <unit_conversion>
     unit_conversion(:) = 1.0
     call chm_transform_profile(obsoper,unit_conversion)
     call chm_convert_units(obsoper,unit_conversion)
  case(2)
     call chm_transform_profile(obsoper,model_col)
     call chm_convert_units(obsoper,model_col)
  end select

  if (obsoper%apply_genoper) then
     ! Perform unit conversion on obsoper%trial when applying the generalized
     ! obs operator for kmode=2,3. Keep obsoper%trial in ug/kg in this case.
     call chm_convert_units(obsoper,obsoper%trial,ppb_opt=.true.)
  end if

! Convert observation vertical coordinate value(s) to pressure if needed

  select case(obsoper%vco)
  case(1)
     ! Convert altitude to pressure
     select case(obsoper%modelIndex)
     case(1)
        press_obs = phf_convert_z_to_pressure(obsoper%obslev,obsoper%gz,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
     case(2,3)
        obsoper%vlayertop = phf_convert_z_to_pressure(obsoper%vlayertop,obsoper%gz,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
        obsoper%vlayerbottom = phf_convert_z_to_pressure(obsoper%vlayerbottom,obsoper%gz,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
     end select
  case(2)
     ! Pressure, no conversion needed
     if (obsoper%modelIndex.eq.1) press_obs = obsoper%obslev
  case(4,5)
     ! No actions taken
  case default
     call utl_abort("chm_obsoperators: vertical coordinate type vco = " // trim(utl_str(obsoper%vco)) // " not available for this operator.")
  end select

! Determine if averaging kernel is to be applied

  obsoper%iavgkern = chm_find_avgkern(obsoper%stnid,obsoper%varno,obsoper%nobslev)

! Apply appropriate core observation operator
   
  select case(obsoper%modelIndex)
  case(1)

!    Vertical interpolation operator

     call chm_vert_interp_operator(obsoper,press_obs,ixtr,success)

  case(2)

!    Surface point measurement

!    Set weight to unity for lowest model level.
     obsoper%zh(1,1:obsoper%nmodlev-1)=0.0
     obsoper%zh(1,obsoper%nmodlev) = 1.0

!    Set range of elements for model vertical levels
     obsoper%imodlev_top(1) = obsoper%nmodlev
     obsoper%imodlev_bot(1) = obsoper%nmodlev 
    
  case(3)
  
!    Layer integration operator

     call chm_layer_integ_operator(obsoper,ixtr,success,kmode)

!  case(4)

!    Layer averaging operator

!    call chm_layer_avg_operator  ! see 3dvar_chem routine ch_vavg

  end select
  
! Apply averaging kernels if requested

  if (obsoper%iavgkern.gt.0) then

     allocate(avg_kern(obsoper%nobslev,obsoper%nobslev))
     
     call chm_get_avgkern(obsoper%iavgkern,obsoper%stnid,obsoper%nobslev,obsoper%lat,obsoper%lon,obsoper%date,obsoper%hhmm,avg_kern)
     do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then

!          Apply averaging kernels to observation operator(s)
           obsoper%zh(iobslev,:)= matmul(avg_kern(iobslev,:),obsoper%zh(:,:))
           if (obsoper%apply_genoper) obsoper%zhp(iobslev,:) = matmul(avg_kern(iobslev,:),obsoper%zhp(:,:))

!          Extend vertical range of obs operator according to the influence of
!          the averaging kernel. Either extend to the entire model vertical range
!          (commented out below) or to the vertical range with non-negligable values.

!          obsoper%imodlev_top(iobslev) = 1
!          obsoper%imodlev_bot(iobslev) = obsoper%nmodlev

           zhmin=1.0D-10*maxval(abs(obsoper%zh(iobslev,:)))
           do imodlev=1,obsoper%imodlev_top(iobslev)
              if (abs(obsoper%zh(iobslev,imodlev)).gt.zhmin) exit
           end do
           if (imodlev.gt.obsoper%imodlev_top(iobslev)) imodlev=obsoper%imodlev_top(iobslev)
           obsoper%imodlev_top(iobslev) = imodlev
           do imodlev=obsoper%nmodlev,obsoper%imodlev_bot(iobslev),-1
              if (abs(obsoper%zh(iobslev,imodlev)).gt.zhmin) exit
           end do
           if (imodlev.lt.obsoper%imodlev_bot(iobslev)) imodlev=obsoper%imodlev_bot(iobslev)
           obsoper%imodlev_bot(iobslev) = imodlev

        end if
     end do

     deallocate(avg_kern)

  end if

! Apply generalized innovation operator if requested

  if (obsoper%apply_genoper) call chm_genoper(obsoper,kmode,success)

! Finalize required quantities depending on kmode
   
  select case(kmode)

  case(0,2)
!
!     Finalize non-linear/linear operator step
!
      do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
           obs_col(iobslev)=dot_product(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)), &
                                 model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))
        end if
      end do

!     Calculate concentration-weighted effective temperature (for output purpose)

      if (kmode.eq.0) then 
         if ((obsoper%constituent_id.ge.0.and.obsoper%constituent_id.le.chm_var_maxnumber()).and. &
            obsoper%modelIndex.eq.3.and.obsoper%nobslev.eq.1.and.obsoper%vco.eq.4.and.success(1)) then
            if (all(obsoper%tt.gt.0.0).and.obs_col(1).gt.0.0) then
                temp_eff(1)=dot_product(obsoper%zh(1,obsoper%imodlev_top(1):obsoper%imodlev_bot(1)), &
                        obsoper%tt(obsoper%imodlev_top(1):obsoper%imodlev_bot(1))*model_col(obsoper%imodlev_top(1):obsoper%imodlev_bot(1))) &
                        /obs_col(1)
                code=oss_obsdata_get_header_code(obsoper%lon,obsoper%lat,obsoper%date,obsoper%hhmm,obsoper%stnid)
                call chm_add_efftemp_obsdata(code,temp_eff)
             end if
          end if
       end if
                       
  case(1)
!
!    Compute sqrt(diag(H*B*H^T))
!     
!    Apply unit conversion to observation operator
     do iobslev=1,obsoper%nobslev
       if (success(iobslev)) then
          obsoper%zh(iobslev,:) = unit_conversion * obsoper%zh(iobslev,:)
       end if
     end do

!    Get background error std dev profile at obs location
     sigma_trial(:,:)=0
     call bchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!      call blamchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!      call benschm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,2)) 

!
!    Identify variable position index in background error correlation matrices
!
     jvar=1
     do while (trim(bchm_varnamelist(jvar)).ne.'') 
        if (trim(bchm_varnamelist(jvar)).eq.trim(obsoper%varName)) exit
        jvar=jvar+1
     end do
     if (trim(bchm_varnamelist(jvar)).eq.'') call utl_abort('chm_genoper: Correlation matrix not found for ' // trim(obsoper%varName) )

     do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
!!           call chm_corvert_mult(obsoper%varName,obsoper%zh(iobslev,:),obs_col(iobslev), &
!!                obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
!!                1,.true.,3,sigma_trial) ! get h*B*h^T

           do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
               zwork(imodlev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                         *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
                         *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
           end do
           obs_col(iobslev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                *zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))

           obs_col(iobslev) = sqrt(obs_col(iobslev))  ! save as sqrt(h*B*h^T)
        else
           obs_col(iobslev) = 0.0
        end if
     end do

  case(3)
!
!     H^T*grad contribution from adjoint of tangent linear model.
!
      model_col(:) = 0.0

      do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
           zwork(:)=0.0D0
            zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) = &
                 obs_col(iobslev)*obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))
                
            call chm_transform_profile(obsoper,zwork)
            call chm_convert_units(obsoper,zwork)
           
            model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) = &
                 model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) + &
                 zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))
        end if
      end do

  end select

  end subroutine chm_obsoperators

!--------------------------------------------------------------------------
!! *Purpose*: Perform any required nonlinear transformations on a model-space profile.
!!
!! @author Y. Rochon, April 2016
!!
!! Revisions:
!!v           J-Caron (ARMA/MRD) and Y Rochon (ARQI/ARQD), Jan. 2018
!!v           - Account for change from LQ to Q for 'HU' (see  lines with !! for original code)
!!
!! Inout
!!   
!!v    obsoper               Contains basic information related to the observation operator
!!v    incrementCol_opt      Model-space increment profile to transform (optional)
!!v    computeDtransform_opt Computes obsoper%dtransform (optional, only used when incrementCol
!!v                          is not provided and default is true in this case)
!!
!!  Comments:
!!v    - If incrementCol is not provided, transformations will be done on obsoper%trial. If
!!v      incrementCol is provided, then transformations will done on incrementCol only and will
!!v      not be done on obsoper%trial.
!!v    - When called with incrementCol provided, it is assumed that obsoper%dtransform has
!!v      already been computed in a previous call.
!!v    - At some point may be merged with chm_apply_transform.
!!
!--------------------------------------------------------------------------
  subroutine chm_transform_profile(obsoper,incrementCol_opt,computeDtransform_opt)
    
    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
    real(8), intent(inout), optional  :: incrementCol_opt(obsoper%nmodlev)
    logical, intent(in), optional :: computeDtransform_opt
    
    logical :: comp_dtransform
    integer :: transform_id

    if (obsoper%constituent_id.lt.0) return
    transform_id=chm_setup_get_int('transform',obsoper%constituent_id)
    if (transform_id.le.0.and.trim(obsoper%varname).ne.'HU') return
    
    if (present(incrementCol_opt)) then

       incrementCol_opt = obsoper%dtransform * incrementCol_opt

    else

       if (present(computeDtransform_opt)) then
          comp_dtransform = computeDtransform_opt
       else
          comp_dtransform = .true.
       end if
       
       if (obsoper%constituent_id.eq.1.and.trim(obsoper%varname).eq.'HU') then

          ! Transformations for water vapor, converting from
          ! specific humidity HU to volume mixing ratio
          !! (replaced) specific humidity LQ to volume mixing ratio
          
          !! Converts LQ=ln(q) to HU=q
          !!obsoper%trial = exp(obsoper%trial)
          !!where(obsoper%trial.gt.0.8) obsoper%trial = 0.8
          !!where(obsoper%trial.lt.-1.D-5) obsoper%trial = -1.D-5
          
          ! Converts specific humidity (q) to mass mixing ratio (rm) rm = q/(1-q)
          ! then from mass mixing ratio (rm) to volume mixing ratio (r) r = m_a/m_H2O rm
          obsoper%trial = obsoper%trial / (1.0-obsoper%trial) &
               * MPC_MOLAR_MASS_DRY_AIR_R8 &
               /chm_setup_get_float('amu',obsoper%constituent_id)
     
          !! For conversion of LQ increment (dlnq) to volume mixing ratio increment (dr) via mass mixing ratio (drm)
          !! drm = dq/(1-q)^2 = q*dlnq/(1-q)^2 = rm(rm+1) dlnq, dr = m_a/m_H2O drm = r(r*m_H2O/m_a + 1) dlnq          
          !!if (comp_dtransform) obsoper%dtransform = obsoper%trial * &
          !!     (obsoper%trial*chm_setup_get_float('amu',obsoper%constituent_id)/MPC_MOLAR_MASS_DRY_AIR_R8 + 1.0)

          ! For conversion of q increment (dq) to volume mixing ratio increment (dr) via mass mixing ratio (drm)
          ! drm = dq/(1-q)^2 = (rm+1)^2 dq, dr = m_a/m_H2O drm = m_a/m_H2O (r*m_H2O/m_a + 1)^2 dq          
          if (comp_dtransform) obsoper%dtransform = MPC_MOLAR_MASS_DRY_AIR_R8/chm_setup_get_float('amu',obsoper%constituent_id) * &
               (obsoper%trial*chm_setup_get_float('amu',obsoper%constituent_id)/MPC_MOLAR_MASS_DRY_AIR_R8 + 1.0)**2
     
       else if (transform_id.eq.1) then

          ! Transformations for ln(x) to x
          
          ! Converts ln(x) to x
          obsoper%trial = exp(obsoper%trial)
          
          ! For conversion of dlnx to dx, dx = x dlnx
          if (comp_dtransform) obsoper%dtransform = obsoper%trial
          
       else
          
          call utl_abort('chm_transform_profile: Transformation #' &
               // trim(utl_str(chm_setup_get_int('transform',obsoper%constituent_id))) // &
               ' for constituent ' // trim(utl_str(obsoper%constituent_id)) // &
               " and variable name " // trim(obsoper%varname) // &
               ' is not defined.')
          
       end if

    end if

  end subroutine chm_transform_profile

!--------------------------------------------------------------------------
!! *Purpose*: Set unit conversion factor for consistency of Hx units with obs units. 
!!
!! @author Y. Rochon and Y. Yang, June 2005 to June 2011
!!
!! Input
!!
!!v     obsoper         Contains basic information related to the observation operator
!!v     ppb             Logical indicating if model_col should be kept in ug/kg instead of
!!v                     the units dictated by the BUFR number (optional, .false. by default)
!!
!! InOut
!!
!!v     model_col       Model space profile to have its units changed
!!
!! Other
!!
!!v     chm_setup_get_float('amu',*)  Molecular mass of constituent (g/mol).
!!
!! Revisions:
!!v           Ping Du, CDMA,Jan 2015
!!v           - Adapted for the EnVar
!!v           Y. Rochon, ARQI/AQRD, Feb. 2015
!!v           - Modifications of adaptations and update of BUFR elements
!!v           M. Sitwell, ARQI/AQRD, April 2016
!!v           - Moved most of the inputs into obsoper
!!v           Y. Rochon, ARQI/AQRD, April 2016
!!v           - Added use of varName and consideration of ug/kg for model fields
!!v           M. Sitwell, ARQI/AQRD, April 2016
!!v           - Moved nonlinear transformations into chm_transform_profile
!!v           Y. Rochon, ARQI/AQRD, Oct 2018
!!v           - Updated to latest set of BUFR_UNIT_* in bufr_mod.f90
!!
!! Further changes required for generalization 
!!
!!v 1) May ultimately have two versions, with one called outside the minimization section. 
!!v 2) Conversion for surface emissions not included as yet (if any is needed)
!!
!! Comments:
!!
!!v     A. Standard model/analysis species field provided as mass mixing 
!!v        ratio in ug/kg (ppb). Conversion to ppb is applied when this is 
!!v        not the case except for AOD and surface emissions.
!!v        As this is hardcoded, any changes in analysis variable must
!!v        be reflected by correspondingly modifying this module.
!!
!!v     B. Unit conversion factor is calculated in chm_convert_units
!!v        from the following factors:
!!v          (1) physical constants
!!v          (2) parameters related to a particular species such as molecular
!!v              mass
!!v          (3) variables such as T and P from background field at each
!!v              iteration
!!
!!v     C. The baseline integral observation operator can be interpreted as being 
!!v        integrals of the gas partial pressure, giving products in kg/m^2, 
!!v        e.g. with sample discretized layer integrals
!!v
!!v               (mass density) * dz = - d(gas partial pressure)/g 
!!v                                   = - [rho(gas)/rho(air)]*dP/g
!!v                                   = - 1E-9 * [mass mixing ratio in parts per billion (ppb)]*dP/g 
!!v
!!v        The actual integration in pressure (in Pascal) is performed outside this routine.
!!v        For integral products in kg/m^2, the output of this routine is to be in mmr/(m/s^2) 
!!v        (mmr=mass mixing ratio), which is equivalent to (1E-9 ug/kg)/(m/s^2) and kg/(m^2*Pa). 
!!v        Therefore, the input value in ug/kg has to be multiplied by 1E-9/g (g=RG below).
!!v
!!v        For integral products in other units, additional conversion factors are to also be applied.
!!
!!v     D. Coefficients related to unit conversion
!!
!!v        rho_stp=1.293                     Air density at STP (1.293 kg/m^3)
!!v        RG=9.807 (=g)                     Acceleration due to gravity (m/s^2)
!!v        MPC_AVOGADRO_R8 = Na              Avogadro's number. 6.023E23 molecules/mole
!!v        MPC_MOLAR_MASS_DRY_AIR_R8 (m_air) Dry air molecular mass. 28.9644 g/mole
!!v        MPC_RGAS_IDEAL_R8 = R             Ideal gas constant. 8.341 J/mole/K  (J=kg m^2/s^2)
!!
!!v                                            PV = nRT (n=number of moles)
!!
!!v        MPC_RGAS_DRY_AIR_R8  = Rd         Dry air constant. 287.1 J/kg/K  (J=kg m^2/s^2)
!!v                                          = MPC_RGAS_IDEAL_R8 * 1000 g/km / MPC_MOLAR_MASS_DRY_AIR_R8
!!
!!v                                            P=rho*Rd*T = [n*m_air*0.001 kg/g]*Rd*T
!!v                                                       = n*[m_air*0.001*Rd]*T
!!v                                                       = nRT
!!
!!v      E. List should be revised following changes to the 'tableburp' file.
!!
!--------------------------------------------------------------------------
  subroutine chm_convert_units(obsoper,model_col,ppb_opt)
  
    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
    real(8), intent(inout) :: model_col(obsoper%nmodlev)
    logical, intent(in), optional :: ppb_opt
    
    ! Declaration of local variables
    real(8) :: zcoef
    integer :: exp_P,exp_T
    logical :: ppb_out
    real(8), parameter :: rho_stp=1.293  ! kg/m^3
    
    if (obsoper%constituent_id.lt.0) return
    
    ! No conversion necessary for these BUFR numbers
    if (any( obsoper%varno .eq. (/ BUFR_UNIT_OptDepth,BUFR_UNIT_OptDepth2, &
            BUFR_UNIT_OptDepth3, BUFR_UNIT_MR_NVaerosol, BUFR_NETT /)  )) return
    
    if (present(ppb_opt)) then
       ppb_out = ppb_opt
    else
       ppb_out = .false.
    end if
      
    zcoef = 1.
    exp_T = 0   ! exponent of multiplicative factor obsoper%tt
    exp_P = 0   ! exponent of multiplicative factor obsoper%pp
    
    ! Convert to ug/kg if not already in those units

    if (obsoper%varName(1:2).eq.'AF' .or. obsoper%varName(1:2).eq.'AC') then

       ! PM2.5 or PM10
       
       if (any(obsoper%varno .eq. (/ BUFR_UNIT_VMR, BUFR_UNIT_VMR2,  &
                    BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2,   &
                    BUFR_UNIT_NumberDensity, BUFR_UNIT_MolarDensity, & 
                    BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2, &
                    BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                    BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2 /) )) &
            call utl_abort("chm_convert_units: BUFR # " // trim(utl_str(obsoper%varno)) // " is not valid for PM" )
       
       ! Conversion from ug/m^3 to ug/kg  (scaling by Rd*T/P)
       zcoef = zcoef * MPC_RGAS_DRY_AIR_R8
       exp_T = exp_T+1 ! multiply by T
       exp_P = exp_P-1 ! divide by P
       
    end if
  
    ! Convert from ug/kg to desired unit if ppb_out = .false.

    if (.not.ppb_out) then
       select case (obsoper%varno)
       
       ! The first four cases below are for integral observations which
       ! require a conversion, in this routine, from ug/kg to the units of 
       ! the integrand values for integrals in pressure (Pascal). Comment C above.
       ! Note: 1 ug/kg = 1 ppb = 1E9 mmr (mass mixing ratio)
       
       case(BUFR_UNIT_IntegDens, BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3) 
        
          ! For conversion from ug/kg to integrand values in kg/(m^2*Pa) 
          ! Note: 1 kg/(m^2**Pa) = = 1 mmr / RG = 1E-9 ug/kg / RG 
          !
          zcoef = zcoef * 1.E-9 / RG
          
       case(BUFR_UNIT_IntegMolarDens) 
        
          ! For conversion from ug/kg to integrand values in moles/(m^2*Pa)
          ! Note: 1 moles/(m^2*Pa) = 1E-9 ug/kg / [RG * (1E-3 kg/g * m_gas)]
          !
          ! To convert from kg/m^2 for the gas to moles/m^2, one must 
          ! divide by the molar mass of the gas in kg/mole.
          !
          ! Note: One u or Da (unified atomic mass unit or dalton) is numerically equivalent to 1 g/mole.
          ! So 1 kg is equivalent to (1E3/(atomic mass)) moles
          
          zcoef = zcoef * 1.E-6 / (RG*chm_setup_get_float('amu',obsoper%constituent_id))
          
       case(BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2)
          
          ! For conversion from ug/kg to integrand values in molecules/(m^2*Pa)
          ! Note: 1 molecule/(m^2*Pa) = 1E-9 ug/kg * Na / [RG * (1E-3 kg/g * m_gas)]
          ! 
          ! To convert from kg/m^2 for the gas to molecules/m^2, one must 
          ! divide by the gas molar mass (kg/mole) and multiply by the Avogrado number          

          zcoef = zcoef * 1.E-6 * MPC_AVOGADRO_R8 &
               / (RG*chm_setup_get_float('amu',obsoper%constituent_id))

       case(BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4) 
      
          ! For conversion from ug/kg to integrand values in DU/Pa
          ! Note: 1 DU/Pa =  1E-9 ug/kg / [RG * m_gas * rho_stp/m_air * 1E-5 m] 
          !
          ! 1 DU = 0.01 mm of gas at STP = 1E-5 m of gas at STP
          !      = integral of gas number density at STP over 1E-5 m.
          !      = Na*P/(RT) * 1E-5   at STP  (Na=Avogadro's number)
          !      = Na*(molar density) * 1E-5   at STP
          !      = Na*rho(STP)/(molar mass) * 1E-5 m
          !      = integral of air number density at STP over 1E-5 m.
          !      = Na*rho(air,STP)/m_air * 1.E-5 m
          !               
          ! Hence 1 DU equivalent to Na*rho_stp/m_air * 1E-5 m (= 2.69E20 molecules/m^2)
          !
          ! To convert from kg/m^2 for the gas to molecules/m^2, one must 
          ! divide by the gas molar mass (kg/mole) and multiply by the Avogrado number 
          !
          ! To convert from molecules/m^2 to DU, one must divide by 2.69E20 or (Na*rho_stp/m_air * 1E-5).
          ! So for conversion from kg/m^2 to DU, must divide by (m_gas*rho_stp/m_air * 1E-5)       
          
          zcoef = zcoef * 1.E-4 * MPC_MOLAR_MASS_DRY_AIR_R8 &
                /(chm_setup_get_float('amu',obsoper%constituent_id)*RG*rho_stp)
        
       case(BUFR_UNIT_Density, BUFR_UNIT_Density2, &
            BUFR_UNIT_AirDensity, BUFR_UNIT_PMDensity)

          ! For conversion from ug/kg to kg/m^3
          !
          ! rho(gas) = mass mixing ratio * rho(air) = mass mixing ratio * P/Rd/T
          
          zcoef = zcoef * 1.E-9 / MPC_RGAS_DRY_AIR_R8
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P
          
       case(BUFR_UNIT_MMR, BUFR_UNIT_MMR2) 

          ! For conversion from ug/kg to kg/kg
     
          zcoef = zcoef * 1.E-9 
     
       case(BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2) 
     
          ! For conversion from ug/kg to partial pressure (PA)
          !
          ! parial pressure = P * vmr
          !                 = P * m_air/m_gas * mass mixing ratio
          
          zcoef = zcoef * 1.E-9 * MPC_MOLAR_MASS_DRY_AIR_R8 &
               /chm_setup_get_float('amu',obsoper%constituent_id)
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_NumberDensity)
          
          ! For conversion from ug/kg to molecules/m^3
          !
          ! Number density of gas = Na*rho(gas)/m_gas = Na*rho(air) * mass mixing ratio /m_gas
          !                       = Na * P/Rd/T * mass mixing ratio /m_gas
        
          zcoef = zcoef * 1.E-6 * MPC_AVOGADRO_R8/MPC_RGAS_DRY_AIR_R8 &
               /chm_setup_get_float('amu',obsoper%constituent_id)
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_MolarDensity)
          
          ! For conversion from ug/kg to moles/m^3
          !
          ! Mole density of gas = rho(gas)/m_gas = rho(air) * mass mixing ratio /m_gas
          !                       = P/Rd/T * mass mixing ratio /m_gas

          zcoef = zcoef * 1.E-6 /MPC_RGAS_DRY_AIR_R8 &
               /chm_setup_get_float('amu',obsoper%constituent_id)
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_VMR, BUFR_UNIT_VMR2, BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2)
          
          ! For conversion from ug/kg to vmr (or moles/mole)
          
          zcoef = zcoef * 1.E-9 * MPC_MOLAR_MASS_DRY_AIR_R8 &
               /chm_setup_get_float('amu',obsoper%constituent_id)

       !case(15192,15011) 
           
          ! Code to be revised when actually applied for the first time
          ! according to model field units.
           
       case default 
        
          call utl_abort('CHM_CONVERT_UNITS: Unknown obs units for varno = ' // trim(utl_str(obsoper%varno)) )
         
       end select
    end if
  
    ! Apply constant scaling
    model_col = model_col * zcoef
    
    if (exp_T.ne.0) then
       if (any(obsoper%tt.le.0.)) call utl_abort("CHM_CONVERT_UNITS: Missing valid temperature for conversion.")
       model_col = model_col * obsoper%tt**exp_T
    end if
    
    if (exp_P.ne.0) model_col = model_col * obsoper%pp**exp_P
    
  end subroutine chm_convert_units
  
!--------------------------------------------------------------------------
!! *Purpose*: Determines if the specifed field name is required somewhere
!!            in the observation operators for a particular observation type.
!!
!! @author M. Sitwell Dec 2015
!!
!! Input
!!
!!v     varName   Name of field 
!!v     varno     BUFR descriptor element
!--------------------------------------------------------------------------
  function chm_required_field(varName,varno) result(needed)

    implicit none

    character(len=*), intent(in) ::varName
    integer, intent(in) :: varno
    logical :: needed
    
    select case(trim(varName))
    case('TT')
 
       select case (varno)
       case(BUFR_UNIT_Density,BUFR_UNIT_Density2,BUFR_UNIT_AirDensity, &
            BUFR_UNIT_PMDensity,BUFR_UNIT_NumberDensity,BUFR_UNIT_MolarDensity)
          needed = .true.
       case default
          needed = .false.
       end select
         
    case default
       needed = .false.
    end select

  end function chm_required_field

!--------------------------------------------------------------------------
!! *Purpose*: Perform layer integration and calculations.
!!
!! @author Y. Rochon, Feb 2015
!!
!! Input 
!!
!!v     ixtr         Flag indicating if obs outside model vertical range
!!v                  0 for no.
!!v     kmode        Observation model stage used to allow option of
!!v                  tropo increment determination from total column data when
!!v                  obsoper%column_bound > obsoper%vlayertop for
!!v                  nobslev=1. For kmode=3, calc only for region between
!!v                  obsoper%column_bound and surface. For kmode=2, split
!!v                  calc for model top to obsoper%column_bound and 
!!v                  obsoper%column_bound and surface.
!!v     iconstituent_id  BUFR constituent ID.
!!
!! InOut
!!
!!v     obsoper  observation operator object
!!v     success      success of integration
!!v    ixtr         Modified ixtr as needed.
!!
!! Revisions:
!!v           M. Sitwell, ARQI/AQRD, Mar 2015
!!v           - Modified to calculate whole profile within a single call
!!v           Y. Rochon ARQOI/AQRD, May 2015
!!v           - Added input of ixtr and iavgkern
!!v           Y. Rochon ARQOI/AQRD, Sept 2015
!!v           - Added iconstituent_id, kmode, ... to varName and related option of producing
!!v             tropo increments from total column data.
!!v           M. Sitwell, ARQI/AQRD, April 2016
!!v           - Many of the input arguments moved into obsoper
!--------------------------------------------------------------------------
  subroutine chm_layer_integ_operator(obsoper,ixtr,success,kmode)

    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(in) :: kmode
    integer, intent(inout) :: ixtr(obsoper%nobslev)
    logical, intent(inout) :: success(obsoper%nobslev)
    
    integer :: ij,iobslev,tropo_mode
    real(8) :: vlayertop_ref,vlayerbottom_ref,imodlev_bot_ref

    ! Conduct initial setup for vertical integration components

    call chm_vertintg_setup(obsoper)

    ! Ensure that each layer is within model vertical range.
    
    do iobslev=1,obsoper%nobslev
     
       if (obsoper%vlayerbottom(iobslev).lt.obsoper%vlayertop(iobslev)) then
          success(1:obsoper%nobslev)=.false.
          write(chm_setup_get_int('message_unit'),*) 'chm_layer_integ_operator: WARNING. Layer top/bot value problem.', &
               obsoper%vlayertop(iobslev), obsoper%vlayerbottom(iobslev), &
               '. Entire profile skipped over.'
          return
       else if (obsoper%vlayerbottom(iobslev).lt.obsoper%pp(1)*1.01 .or. &  
            obsoper%vlayertop(iobslev).gt.obsoper%pp(obsoper%nmodlev)*0.99) then
          success(iobslev)=.false.
          if (obsoper%vlayerbottom(iobslev).lt.obsoper%pp(1)*1.01) then
             ixtr(iobslev)=1
          else
             ixtr(iobslev)=2
          end if
          write(chm_setup_get_int('message_unit'),*) 'chm_layer_integ_operator: WARNING. Layer top/bot value problem.', &
               obsoper%vlayertop(iobslev), obsoper%vlayerbottom(iobslev)
          cycle
       end if
       if (obsoper%vlayerbottom(iobslev).gt.obsoper%pp(obsoper%nmodlev)*0.999) obsoper%vlayerbottom(iobslev)=obsoper%pp(obsoper%nmodlev)*0.999
       if (obsoper%vlayertop(iobslev).lt.obsoper%pp(1)*1.001) obsoper%vlayertop(iobslev)=obsoper%pp(1)*1.001
      
    end do

    ! Check for special treatment if tropo_mode>=1, kmode=2,3, and nobslev=1 for
    ! column observations (obsoper%vco=4)  that extend to the surface.

    tropo_mode = chm_setup_get_int('tropo_mode',obsoper%constituent_id)
    
    if (obsoper%nobslev.eq.1.and.kmode.ge.2.and.obsoper%vlayerbottom(1).gt.obsoper%pp(obsoper%nmodlev)*0.99.and. &
         obsoper%constituent_id.ge.0.and.obsoper%vco.eq.4) then
       
       if (obsoper%constituent_id.gt.chm_var_maxnumber()) &
            call utl_abort("chm_layer_integ_operator: Invalid constituent ID with value " // trim(utl_str(obsoper%constituent_id)))
       
       vlayerbottom_ref=obsoper%vlayerbottom(1)
       
       if ( tropo_mode.ge.1 .and. obsoper%column_bound.gt.obsoper%vlayertop(1) ) then
          
          if (obsoper%iavgkern.ne.0) &
               call utl_abort("chm_layer_integ_operator: Use of averaging kernels not possible with reduced range of increment profile.")
          
          if (kmode.eq.2.and.tropo_mode.eq.1) then
             
             ! When kmode=2, split calc in two. This is done due to difference in 
             ! calc at the interface region when producing zh and zhp. The tangent linear
             ! model in the lower region for kmode=2 must be consistent with 
             ! that associated to kmode=3.
             
             ! Start with bottom region in order to use correct zhp with chm_genoper
             ! when use of this operator is requested.
             
             vlayertop_ref=obsoper%vlayertop(1)
             
             obsoper%vlayertop(1) = obsoper%column_bound
             
             call chm_vertintg(obsoper,ixtr,success)
             
             ! Apply generalized innovation operator if requested
             
             if (obsoper%apply_genoper) call chm_genoper(obsoper,kmode,success)
             
             obsoper%apply_genoper=.false.
             
             imodlev_bot_ref=obsoper%imodlev_bot(1)
             
             ! Reset top and bottom values for integration of the remaining region.
             ! The second integration provides the change in upper level contributions to the
             ! total column from the assimilation of other observations.
             
             obsoper%vlayertop(1)=vlayertop_ref
             obsoper%vlayerbottom(1)=obsoper%column_bound
          else
             ! Reset top new value. Restricts adjoint/tangent linear calcs to this reduced region.            
             obsoper%vlayertop(1)=obsoper%column_bound
          end if
          
       end if
    else
       vlayerbottom_ref=0.0
    end if
    
    ! Calculate vertical integration components for specified obs layer.
    
    call chm_vertintg(obsoper,ixtr,success)
    
    ! If tropo_mode=1, reset original vertical range 
    ! for the tangent linear operator 
    if (obsoper%nobslev.eq.1.and.kmode.eq.2.and.obsoper%constituent_id.ge.0.and.obsoper%vco.eq.4) then
       if (tropo_mode.eq.1.and.  &
            obsoper%column_bound.gt.obsoper%vlayertop(1).and. &
            vlayerbottom_ref.gt.obsoper%pp(obsoper%nmodlev)*0.99) then    
          
          obsoper%vlayerbottom(1)=vlayerbottom_ref
          obsoper%imodlev_bot(1)=imodlev_bot_ref
       end if
    end if
   
  end subroutine chm_layer_integ_operator

!--------------------------------------------------------------------------
!! *Purpose*: Preliminary calculations for producing components required for vertical 
!!            integration w.r.t. pressure to calculate partial (or total)
!!            column value of model state profile or used for adjoint calc.
!!
!! @author Y. Rochon, P. Du, and M. Sitwell, Feb 2015 (based on pre-Envar version by
!!         Y. Rochon, Y. Yang and S. Ren, Nov 2004 to Dec 2012)
!!
!! Output  
!!v        obsoper%vmodpress(nmodlev+1)       -- Model layer boundaries given that pressmod are
!!v                                              taken at mid-layer values.
!!v        obsoper%vweights(nmodlev,nmodlev)  -- Second order Lagrangian interp integration weights
!!
!! Comments:
!!v        - This subroutine does the following:
!!v           - Setting of model layer boundaries
!!v           - Determining integration weights associated to second order Lagrangian interpolation.
!!v        - Layer boundaries are taken as mid-point between eta levels in lnP coordinate. Layer
!!v          values are set to be the values interpolated to the mid-point in P within the various
!!v          layers. Interpolation in P is done quadratically. 
!--------------------------------------------------------------------------
  subroutine chm_vertintg_setup(obsoper)

    implicit none
      
    type(struct_chm_obsoperators), intent(inout) :: obsoper

    ! Declaration of local variables

    integer   :: jk
    real(8)   :: zp, zp1, zp2, zp3, zr1, zr2, zr3

    ! Determine P boundaries of analysis layers and save weights for
    ! use in setting innovation operator array.
    ! N.B.: Boundaries of layers set to mid-point of model levels
      
    ! Calculate layer boundaries

    obsoper%vmodpress(1)=obsoper%pp(1)
    obsoper%vmodpress(obsoper%nmodlev+1)= obsoper%pp(obsoper%nmodlev)

    DO JK = 2, obsoper%nmodlev
       obsoper%vmodpress(jk)=sqrt(obsoper%pp(jk-1)*obsoper%pp(jk))
    END DO

    ! Interpolation to mid-layer level in P using
    ! second degree Lagrangian interpolator.
    ! N.B.: Integration is w.r.t. P
    
    ! Calculating for jk=1
    
    zp1= obsoper%pp(1)
    zp2= obsoper%pp(2)
    zp3= obsoper%pp(3)
    zp = (obsoper%vmodpress(2)+obsoper%vmodpress(1))/2.0
    zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
    zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
    zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
    obsoper%vweights(1,1)=zr1
    obsoper%vweights(2,1)=zr2
    obsoper%vweights(3,1)=zr3

    DO JK=2,obsoper%nmodlev-1
       zp1=obsoper%pp(jk-1)
       zp2=obsoper%pp(jk)
       zp3=obsoper%pp(jk+1)
       zp=(obsoper%vmodpress(jk+1)+obsoper%vmodpress(jk))/2.0
       zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
       zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
       zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
       obsoper%vweights(jk-1,jk)=zr1
       obsoper%vweights(jk,jk)=zr2
       obsoper%vweights(jk+1,jk)=zr3
    ENDDO
    
    ! Calculating  for jk=obsoper%nmodlev
    
    zp1= obsoper%pp(obsoper%nmodlev-2)
    zp2= obsoper%pp(obsoper%nmodlev-1)
    zp3= obsoper%pp(obsoper%nmodlev)
    zp = (obsoper%vmodpress(obsoper%nmodlev+1)+obsoper%vmodpress(obsoper%nmodlev))/2.0
    zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
    zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
    zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
    obsoper%vweights(obsoper%nmodlev-2,obsoper%nmodlev)=zr1
    obsoper%vweights(obsoper%nmodlev-1,obsoper%nmodlev)=zr2
    obsoper%vweights(obsoper%nmodlev,obsoper%nmodlev)=zr3

  end subroutine chm_vertintg_setup

!--------------------------------------------------------------------------
!! *Purpose*: Calculate components required for vertical integration w.r.t.
!!            pressure to calculate partial (or total) column value of model
!!            state profile or used for adjoint calc.
!!
!! @author Y. Rochon, Y. Yang and S. Ren, Nov 2004 to Dec 2012
!!
!! Input   
!!v        obsoper%nmodlev       -- # of model vertical levels
!!v        obsoper%nobslev       -- # of obs vertical levels
!!v        obsoper%vweights      -- See routine chm_vertintg_setup
!!v        obsoper%vmodpress
!!v        success               -- Logical indicating if calc are to be performed.
!!v        ixtr                  -- Flag indicating if obs outside model vertical range (0 for no)
!!
!! Output  
!!v        obsoper%zh(obsoper%nobslev,obsoper%nmodlev)  -- Initial innovation model array 
!!v                                                        (other than conversion constants)
!!v        obsoper%zhp(obsoper%nobslev,obsoper%nmodlev) -- Part of innovation operator not 
!!v                                                        related to resolution
!!
!! Revisions:
!!v        Ping Du and Y. Rochon, Jan-Feb 2015
!!v         - Adapted for the EnVar
!!v        M. Sitwell, ARQI/AQRD, Mar 2015
!!v         - Modified to calculate whole profile within a single call
!!v        Y. Rochon ARQOI/AQRD, May 2015
!!v         - Added input and use of ixtr and iavgkern
!!v        M. Sitwell, ARQI/AQRD, April 2016
!!v         - Some input arguments moved into obsoper
!--------------------------------------------------------------------------
  subroutine chm_vertintg(obsoper,ixtr,success)

    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(in) ::ixtr(obsoper%nobslev)
    logical, intent(in) :: success(obsoper%nobslev)
    
    integer, parameter :: ivweights=2  ! Order of Lagrangian interpolation.

    ! Declaration of local variables
    
    integer   :: J,JK,ILMAX2,ILMIN2
    integer   :: ILMIN, ILMAX, iobslev
    real(8)   :: zp, zp1, zp2, zp3, zr1, zr2, zr3, ptop, pbtm

    do iobslev=1,obsoper%nobslev

       if (success(iobslev).or.(ixtr(iobslev).eq.0.and.obsoper%iavgkern.ne.0)) then

          ptop = obsoper%vlayertop(iobslev)
          pbtm = obsoper%vlayerbottom(iobslev)
         
          ! Find the range of vertical levels over which to perform the integration
          ! and set innovation operator ZH over this range.
          
          ilmin=1
          ilmax=obsoper%nmodlev
          if (ptop.le.obsoper%vmodpress(1)*1.01.and.pbtm.ge.obsoper%vmodpress(obsoper%nmodlev+1)*0.99) then

             ! Total column integration part
             
             do jk = 1,obsoper%nmodlev
                do j=max(1,jk-ivweights),min(obsoper%nmodlev,jk+ivweights)
                   obsoper%zh(iobslev,jk)=obsoper%zh(iobslev,jk)+(obsoper%vmodpress(j+1) &
                        -obsoper%vmodpress(j))*obsoper%vweights(jk,j)
                   obsoper%zhp(iobslev,jk)=obsoper%zhp(iobslev,jk)+obsoper%vweights(jk,j)
                end do
             end do
             
          else

             ! Partial column integration part (special treatment at boundaries)
     
             ! Identify analysis layer boundaries just within obs layer.
             
             ilmin = chm_igetmodlev(ptop, obsoper%vmodpress, 'top', obsoper%nmodlev+1)
             ilmax = chm_igetmodlev(pbtm, obsoper%vmodpress, 'btm', obsoper%nmodlev+1)
               
             if (ilmin.eq.ilmax+1) then

                ! Entire obs layer within one analysis layer
                
                j=ilmin
                if (j.lt.3) j=3
                if (j.gt.obsoper%nmodlev) j=obsoper%nmodlev
                zp1=obsoper%vmodpress(j-2)
                zp2=obsoper%vmodpress(j-1)
                zp3=obsoper%vmodpress(j)
                zp=(ptop+pbtm)/2.0
                zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                
                obsoper%zh(iobslev,j-2)=(pbtm-ptop)*zr1
                obsoper%zh(iobslev,j-1)=(pbtm-ptop)*zr2
                obsoper%zh(iobslev,j)=(pbtm-ptop)*zr3
                obsoper%zhp(iobslev,j-2)=zr1
                obsoper%zhp(iobslev,j-1)=zr2
                obsoper%zhp(iobslev,j)=zr3
                ilmin=j-2
                ilmax=j
                  
             else
                
                ! Determine terms from the inner layers (excluding the lower and upper
                ! boundary layers when these layers not covering entire analyses layers)
                
                if (pbtm.ge.obsoper%vmodpress(obsoper%nmodlev)*0.99) then
                   ilmax2=obsoper%nmodlev
                else
                   ilmax2=ilmax-1
                end if
                if (ptop.le.obsoper%vmodpress(1)*1.01) then
                   ilmin=1
                   ilmin2=ilmin
                else
                   ilmin2=ilmin
                end if
                if (ilmin2.le.ilmax2) then
                   do jk = ilmin2,ilmax2
                      do j=max(1,jk-ivweights),min(obsoper%nmodlev,jk+ivweights)
                         obsoper%zh(iobslev,jk)=obsoper%zh(iobslev,jk)+(obsoper%vmodpress(j+1) &
                               -obsoper%vmodpress(j))*obsoper%vweights(jk,j)
                         obsoper%zhp(iobslev,jk)=obsoper%zhp(iobslev,jk)+obsoper%vweights(jk,j)
                      end do
                   end do
                end if
                
                ! Determine terms from the lower and upper boundary layers
                ! when these layers do not cover entire analyses layers.
                
                if (pbtm.lt.obsoper%vmodpress(obsoper%nmodlev)*0.99) then
                     
                   j=ilmax+1
                   if (j.gt.obsoper%nmodlev) j=obsoper%nmodlev
                   if (j.lt.3) j=3
                   zp1=obsoper%vmodpress(j-2)
                   zp2=obsoper%vmodpress(j-1)
                   zp3=obsoper%vmodpress(j)
                   zp=(obsoper%vmodpress(ilmax)+pbtm)/2.0
                   zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                   zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                   zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
                   obsoper%zh(iobslev,j-2)=obsoper%zh(iobslev,j-2)+(pbtm - obsoper%vmodpress(ilmax))*zr1
                   obsoper%zh(iobslev,j-1)=obsoper%zh(iobslev,j-1)+(pbtm - obsoper%vmodpress(ilmax))*zr2
                   obsoper%zh(iobslev,j)=obsoper%zh(iobslev,j)+(pbtm - obsoper%vmodpress(ilmax))*zr3
                   obsoper%zhp(iobslev,j-2)=obsoper%zhp(iobslev,j-2)+zr1
                   obsoper%zhp(iobslev,j-1)=obsoper%zhp(iobslev,j-1)+zr2
                   obsoper%zhp(iobslev,j)=obsoper%zhp(iobslev,j)+zr3
                   ilmax=j
                  
                end if
                  
                if (ptop.gt.obsoper%vmodpress(1)*1.01) then
                     
                   j=ilmin-1
                   if (j.lt.1) j=1
                   if (j.gt.obsoper%nmodlev-2) j=obsoper%nmodlev-2
                   zp1= obsoper%vmodpress(j)
                   zp2= obsoper%vmodpress(j+1)
                   zp3= obsoper%vmodpress(j+2)
                   zp = (obsoper%vmodpress(ilmin)+ptop)/2.0
                   zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                   zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                   zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
                   obsoper%zh(iobslev,j)=obsoper%zh(iobslev,j)+(obsoper%vmodpress(ilmin)-ptop)*zr1
                   obsoper%zh(iobslev,j+1)=obsoper%zh(iobslev,j+1)+(obsoper%vmodpress(ilmin)-ptop)*zr2
                   obsoper%zh(iobslev,j+2)=obsoper%zh(iobslev,j+2)+(obsoper%vmodpress(ilmin)-ptop)*zr3
                   obsoper%zhp(iobslev,j)=obsoper%zhp(iobslev,j)+zr1
                   obsoper%zhp(iobslev,j+1)=obsoper%zhp(iobslev,j+1)+zr2
                   obsoper%zhp(iobslev,j+2)=obsoper%zhp(iobslev,j+2)+zr3
                   ilmin=j
                   if (ilmax.lt.j+2) ilmax=j+2
                   
                end if
                if (ilmin.gt.ilmax-2) ilmin=ilmax-2
             end if
          end if

          obsoper%imodlev_top(iobslev)=ilmin
          obsoper%imodlev_bot(iobslev)=ilmax
            
       else
          obsoper%zh(iobslev,:) = 0.0D0
          obsoper%zhp(iobslev,:) = 0.0D0
          
          obsoper%imodlev_top(iobslev)=1
          obsoper%imodlev_bot(iobslev)=1
       end if
       
    end do

  end subroutine chm_vertintg

!--------------------------------------------------------------------------
!! *Purpose*: Interpolation to point in profile. Uses piecewise linear vertical
!!            interpolation in log(Pressure).
!!
!! @author M. Sitwell, March 2015 (based on ch_vprof from 3DVAR-CHEM by Y. Rochon July 2005)
!!
!! Revisions:
!!v        Y. Rochon, May 2015
!!v         - Added input and use of ixtr and iavgkern
!!v        M. Sitwell, April 2015
!!v         - Some input arguments moved into obsoper
!! Input
!!
!!v       obsoper%pp      pressure on model levels, assumed to be in ascending order
!!v       pres_obs        pressure on observation levels
!!v       ixtr            Flag indicating if obs outside model vertical range (0 for no)
!!
!! Output
!!
!!v       obsoper%zh      interpolation coefficients
!!v       success         success of interpolation
!!v       ixtr            Modified ixtr as needed.
!!
!! Comments:
!!v      - Current implementation searches for index of nearest model level. This step
!!v        is redundant since this information is already saved in obsSpaceData in OBS_LYR.
!!v        This step is repeated so that the routines in chem_obsoperators_mod are more independent
!!v        of the rest of the EnVar code. If it is desired to skip this redundant step,
!!v        the content of the OBS_LYR column could be passed to chm_obsoperators and
!!v        subsequentially to this subroutine.
!--------------------------------------------------------------------------
  subroutine chm_vert_interp_operator(obsoper,pres_obs,ixtr,success)

    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(inout) :: ixtr(obsoper%nobslev)
    real(8), intent(in) :: pres_obs(obsoper%nobslev)
    logical, intent(inout) :: success(obsoper%nobslev)

    integer :: iobslev,jmodlev

    do iobslev=1,obsoper%nobslev

       ! check if obs is above or below model boundaries
       if ( pres_obs(iobslev).lt.obsoper%pp(1) .or. &
            pres_obs(iobslev).gt.obsoper%pp(obsoper%nmodlev) ) then
          success(iobslev)=.false.
          if (pres_obs(iobslev).lt.obsoper%pp(1)) then
              ixtr(iobslev)=1
          else
              ixtr(iobslev)=2
          end if 
       end if

       if (success(iobslev).or.(ixtr(iobslev).eq.0.and.obsoper%iavgkern.ne.0)) then

          ! Find model layers directly above and below obs.
          ! After exit of loop, the obs will be between model
          ! levels jmodlev and jmodlev+1.
          do jmodlev=1,obsoper%nmodlev-1
             if ( pres_obs(iobslev).ge.obsoper%pp(jmodlev) .and. &
                  pres_obs(iobslev).lt.obsoper%pp(jmodlev+1) ) then
                exit
             end if
          end do

          ! Set interpolation weights
          obsoper%zh(iobslev,jmodlev+1) = LOG(pres_obs(iobslev)/obsoper%pp(jmodlev)) &
                                        / LOG(obsoper%pp(jmodlev+1)/obsoper%pp(jmodlev))
          obsoper%zh(iobslev,jmodlev) = 1.0D0 - obsoper%zh(iobslev,jmodlev+1)

          ! set range of nonzero elements for model vertical levels
          obsoper%imodlev_top(iobslev) = jmodlev
          obsoper%imodlev_bot(iobslev) = jmodlev+1

       else
          obsoper%imodlev_top(iobslev) = 1
          obsoper%imodlev_bot(iobslev) = 1
       end if

    end do

  end subroutine chm_vert_interp_operator

!--------------------------------------------------------------------------
!! *Purpose*: Get the vertical level index for the pressure in rppobs
!!            within obs layer and nearest specified obs layer boundary.
!!
!! @author Y. Yang, May 2004
!!
!! Arguments
!!
!!v       rpress        pressure value in Pascal
!!v       rppobs        profile of pressure at obs. location
!!v       topbtm        indicating whether we are looking for top or bottom pressure
!!v       ntotlev       total number of levels of rppobs
!--------------------------------------------------------------------------
  integer function chm_igetmodlev(rpress, rppobs, topbtm, ntotlev)

    implicit none
 
    integer, intent(in) :: ntotlev
    real(8), intent(in) :: rpress, rppobs(ntotlev)
    character(len=*), intent(in) :: topbtm
    
    integer     :: ilev1, ilev2
    integer     :: jk

    ! Find the model levels adjacent to pressure level rpress

    ! Default values
    
    if (rpress .lt. 0.) then
       if ((topbtm .eq. 'btm') .or. (topbtm .eq. 'BTM')) then
          chm_igetmodlev = ntotlev
       endif
       if ((topbtm .eq. 'top') .or. (topbtm .eq. 'TOP')) then
          chm_igetmodlev = 1
       endif                                                  
    endif
      
    ilev1=0
    ilev2=1
    do jk=1,ntotlev
       if (rpress.gt.rppobs(jk)) then
          ilev1=jk
          ilev2=jk+1
       else
          exit
       endif
    enddo

    ! Find the model level index

    ! If we are looking for top level, the index is the level immediately 
    ! below. if looking for bottom level, the index is the one immediately 
    ! above.
    
    if ((topbtm .eq. 'btm') .or. (topbtm .eq. 'BTM')) then
       chm_igetmodlev=ilev1
    else if ((topbtm .eq. 'top') .or. (topbtm .eq. 'TOP')) then
       chm_igetmodlev=ilev2
    endif

    if (chm_igetmodlev .lt. 1) chm_igetmodlev=1
    if (chm_igetmodlev .gt. ntotlev) chm_igetmodlev=ntotlev
  
  end function chm_igetmodlev

!--------------------------------------------------------------------------
!! *Purpose*: Set generalized innovation operator for integral or layer avg obs.
!!            Relevant only for 3D incremental fields. This version is intended
!!            to vertically distribute the obs increments proportionally to the
!!            background state.
!!
!! @author Y. Rochon, April 2015 (from pre-EnVar version of Nov 2004)
!!
!! Input   
!!
!!v      obsoper%zh,zhp       see routine chm_vertintg_setup
!!v      kmode                index specifying if content to be applied (i.e. if kmode>1)
!!v      success              logical indicating if calc are to be performed.
!!
!! Output  
!!v      obsoper%zh(obsoper%nmodlev)     a*w: Final innovation model array (other than conversion constants)
!!v      obsoper%zhp(obsoper%nmodlev)    w (see comments section)
!!
!! Revisions:
!!v         M. Sitwell, April 2016
!!v          - Some input arguments moved into obsoper
!!v         Y. Rochon, May 2016
!!v          - Added weighting cases 
!!
!! Comments:
!!
!!v     (1) This routine prepares an alternative innovation operator g, called
!!v     the generalized innovation operator, to take the place of the
!!v     innovation (TLM) operator h (row of zh). The operator g is
!!v     specified as:
!!v
!!v            g = a*w
!!v
!!v     where the modified innovation operator 'w' can be set as:
!!v
!!v            w = P[ (h'x)^T ] *  B^{-1}     
!!v
!!v     with  h' is the part of h which excludes resolution dependence
!!v              (only/mostly contains the physics part of h; zhp),
!!v           x is the state profile rval
!!v           P is a window cutoff operator (sets small values to zero), and
!!v           B is the original/initial total "vertical" covariance matrix (in 2D)
!!v         
!!v
!!v     and 'a' is a proportionality constant ensuring that the innovation
!!v             increment remains unchanged for the 1D case in the absence
!!v             of other obs., i.e.,
!!v
!!v                 a^2 = (h*B*h^T)(w*B*w^T)^{-1},
!!v
!!v     Application of the state profile x (rval) is to make the
!!v     increment profile be more proportional to the state profile.
!!v
!!v     The presence of B^{-1} is to negate the weight re-distribution from the later 
!!v     application of B in grad(Jo).
!!v
!!v     While dx is provided to the obs operator, the minimization is done for
!!v     dx/sigma where sigma is the background error std. dev. in B and so C (correlation matrix)
!!v     is used instead of B in the minimization. Moreover, the transformation from dx/sigma 
!!v     to dx is done outside the forward model operators (at the spectral to physical space 
!!v     transformation step). For this reason, the expression for w should technically be replaced by
!!v 
!!v            w = P[ (h'x/sigma^2)^T ] *  C^{-1}             (Option 1 below)
!!v
!!v     still with
!!v
!!v            a^2 = (h*B*h^T)(w*B*w^T)^{-1}
!!v 
!!v     The presence of C^{-1} does/can give difficulty to the iterative variational 
!!v     minimization. It can result in oscillations in the increment profile depending on
!!v     where the iterations are stopped. Moreover, if the spectral space C matrix is for 
!!v     non-seperable vertical and horizontal correlations, there will be oscillations 
!!v     due to inconstencies between the total inverse vertical correlation C^{-1} in physical space 
!!v     and the inverse vertical correlation matrix for each spectral wavenumber.
!!v
!!v     As alternatives, one can completely omit the role of  C^{-1} from w, i.e.
!!v
!!v            w = P[ (h'x/sigma^2)^T ]                        (Option 3)
!!v
!!v      or use the following substitute to approximate the role of C^{-1} in approximatily
!!v      negating the weight re-distribution from later application of C in grad(Jo), i.e.
!!v
!!v            w(i) = P[ (h'x/sigma^2)^T ]_i / sum(C(:,i))     (Option 2 - preferred)
!!v
!!v     (2) The matrices B and B^{-1} are the total error covariance matrix (in physical space)
!!v     and its inverse with the related error correlation matrices 'corvert' and 'corverti' 
!!v     provided from 'bmatrixchem_mod.ftn90'.
!!v
!!v     (3) In the presence of both (a) neighbouring measurements and (b) horizontal background error 
!!v     correlation lengths that vary in the vertical, the increments will also be subject to the
!!v     the latter, displaying larger increments in vertical regions with larger horizontal
!!v     error correlation lengths - this distorting the vertical increment distribution stemming
!!v     from chm_genoper alone. This stems from w accounting for vertical correlations via C^{^-1}
!!v     and not the horizontal correlations in B.
!!v
!!v     (4) NOTE: Cases with ensemble-based and or lam-based background covariances
!!v     are not taken into account in this version.
!--------------------------------------------------------------------------
  subroutine chm_genoper(obsoper,kmode,success)
      
    implicit none
      
    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(in) :: kmode
    logical, intent(in) :: success(obsoper%nobslev)
    
    ! Declaration of local variables
    
    real(8), parameter :: pwin=0.01
    integer  :: iobslev,imodlev,irmse,jvar
    logical  :: lrgsig
    real(8)  :: zwbw(1),zhbh(1),za,work(obsoper%nmodlev),sigma_trial(obsoper%nmodlev,4)
    real(8), parameter :: threshold=1.D-20
    real(8)  :: zmin,rvalw(obsoper%nmodlev),rvalr(obsoper%nmodlev)
    real(8)  :: rvalc(obsoper%nmodlev),rmse,w1 
    character(len=22) :: code
   
    if (kmode.le.1) return

    ! Retrieve from stored background error std dev [elemements (:,1-2)] at obs location [and inverses at elements (:,3-4)]    
    sigma_trial=chm_retrieve_sigma_trial(obsoper%obs_index)               

    ! Identify variable position index in background error correlation matrices
    
    jvar=1
    do while (trim(bchm_varnamelist(jvar)).ne.'') 
       if (trim(bchm_varnamelist(jvar)).eq.trim(obsoper%varName)) exit
       jvar=jvar+1
    end do
    if (trim(bchm_varnamelist(jvar)).eq.'') call utl_abort('chm_genoper: Correlation matrix not found for ' // trim(obsoper%varName) )
    
    ! Initialize reference mass (mixing ratio) weighting profile as profile from trial field
     
    rvalr(1:obsoper%nmodlev)=obsoper%trial(1:obsoper%nmodlev)
    
    ! Check on rvalr
            
    if (any(abs(rvalr(:)).le.threshold)) then
       ! write(*,*) 'chm_genoper: Trial field profile segment is ~zero. Could cause an abort in this routine.'
       ! write(*,*) 'chm_genoper: To prevent abort and still be effective, corresponding elements set to a non-zero constant.'
       ! write(*,*) 'chm_genoper: ',varName,rlat,rlong,nobslev
       if (all(abs(rvalr(:)).le.threshold)) then
          rvalr(:)=1.0 
       else  
          zmin=minval(abs(rvalr(:)), &
               mask=abs(rvalr(:)).gt.threshold)
          where (abs(rvalr(:)).le.threshold) &
               rvalr(:)=zmin  
       end if
    end if
    
    ! Loop over obs elements
    
    lrgsig=.true. 
    write(code,'(I22)') obsoper%obs_index
    
    do iobslev=1,obsoper%nobslev
       
       if (.not.success(iobslev)) cycle   
       
       ! Set reference mass (mixing ratio) weighting profile
       
       rvalw(1:obsoper%nmodlev)=rvalr(1:obsoper%nmodlev)
       if (chm_setup_get_int('genoper',obsoper%constituent_id).gt.1) then

          ! Set reference mass weighting profile according to the difference between 
          ! an external reference (such as a climatology) and trial field profiles.
          !
          ! This is a mechanism to force the solution profile shape somewhat towards that 
          ! of the external reference when the reliability of the vertical structure of the 
          ! trial field is not high.
          !
          ! This can be used at the beginning of long assimilation periods if the initial
          ! trial field is not as realistic as may be desired or somewhat mitigate
          ! gradual biasing of vertical structures that might otherwise occur from assimilation 
          ! of integrated quantities (when there is insufficient data from other observation types).
          !
          ! The larger the rms difference of xc (external reference) with xb (trial field profile), 
          ! the greater is the influence of this difference in the weighting. If there is 
          ! little difference, then either xb or xc can be directly used as the weighting 
          ! profile; xb is used below.
          
          rmse=0.0
          irmse=0
          rvalc(1:obsoper%nmodlev)=chm_get_ref_column(code)
            
          zmin=pwin*maxval(abs(obsoper%zhp(iobslev,1:obsoper%nmodlev)))
          do imodlev=1,obsoper%nmodlev
             if (obsoper%zhp(iobslev,imodlev).gt.zmin) then
                irmse=irmse+1
                rmse=rmse+((rvalc(imodlev)-obsoper%trial(imodlev))*sigma_trial(imodlev,3))**2
             end if
          end do
          rmse=rmse*2./irmse
          if (rmse.lt.1.0) then
             rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)
          else 
             rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)-obsoper%trial(1:obsoper%nmodlev)
             if (rmse.lt.5.0) then
                w1=(5.0-rmse)/(5.0-1.0)
                rvalw(1:obsoper%nmodlev)=(1.0-w1)*rvalw(1:obsoper%nmodlev)+w1*obsoper%trial(1:obsoper%nmodlev)
             end if
             where (abs(rvalw(1:obsoper%nmodlev)).lt.0.01*rvalc(1:obsoper%nmodlev)) &
                  rvalw(1:obsoper%nmodlev)=sign(0.01*rvalc(1:obsoper%nmodlev),rvalw(1:obsoper%nmodlev))
          end if
       end if
       
       ! Begin preparation of the new innovation operator w (=new zhp)
       
       if (obsoper%nobslev.eq.1.and.obsoper%imodlev_top(iobslev).eq.1.and. &
            obsoper%imodlev_bot(iobslev).eq.obsoper%nmodlev) then
          
          ! Treat as total column obs. Here, zhp would be approx. equal
          ! to 1 except for the near-end points of the model vertical domain,
          ! the latter due to the discretized domain. Not using zhp avoids this
          ! discretization issue from weakly affecting results at the boundaries.
          
          work(1:obsoper%nmodlev)=obsoper%trial(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3) 
       else 

          ! Account for localized obs function (e.g. partial columns, Jacobians. For Jacobians,
          ! zhp must also be independent of the model layer thicknesses.)
          work(1:obsoper%nmodlev)=obsoper%zhp(iobslev,1:obsoper%nmodlev)* &
               obsoper%trial(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3) 
       end if

       ! Apply cutoff (apply to zhp*obsoper%trial/sigma_trial(:,1) instead of the resultant zh)
       ! Resultant outside cutoff region should be zhp*obsoper%trial/sigma_trial(:,1)/sigma_trial(:,1).
       ! Note: sigma_trial(:,3)=1.0/sigma_trial(:,1)
 
       zmin=pwin*maxval(abs(work(1:obsoper%nmodlev)))
       where (abs(work(1:obsoper%nmodlev)).lt.zmin) 
          work(1:obsoper%nmodlev)=0.0D0 
       elsewhere        
          work(1:obsoper%nmodlev)=work(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3)
       endwhere

       ! Application of C^{-1} or substitute (for negating the weight impact of the
       ! later application of C in finalizing grad(Jo). Option 2 is favoured.
       !  
       ! Option 1: Application of C^{-1}
       ! call chm_corvert_mult(obsoper%varName,work(1:obsoper%nmodlev),obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                       obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev,obsoper%nmodlev, &
       !                       .false.,-1)
       !
       ! Option 2: Application of 1/sum(C(:,j)) to approximately negate the weight re-distribution from C
       !           in the calc of grad(Jo).
       !
       ! call chm_corvert_mult(obsoper%varName,work(1:obsoper%nmodlev),obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                      obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev,obsoper%nmodlev, &
       !                      .false.,0)
       obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))= &
            obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
            *bchm_invsum(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar)

       ! Option 3: Just skip over consideration attempt at negating the weight re-distribution from C.
       ! chm_obsoper%zhp(iobslev,1:obsoper%nmodlev)=work(1:obsoper%nmodlev)
       !
       ! Determine proportionality factor 'a' = (h*B*h^T)(w*B*w^T)^{-1}
       !
       ! Determine/estimate w*B*w^T (zwbw(1))
       !
       ! call chm_corvert_mult(obsoper%varName,obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                       zwbw,obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
       !                       1,lrgsig,3,sigma_trial))
       do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
          work(imodlev)=sum(obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
               *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
               *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
       end do
       zwbw(1)=sum(obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
            *work(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))
       
       ! Determine/estimate h*B*h^T (zhbh(1))
       !
       !  call chm_corvert_mult(obsoper%varName,obsoper%zh(iobslev,1:obsoper%nmodlev), &
       !                        zhbh,obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
       !                        1,lrgsig,3,sigma_trial)
         do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
            work(imodlev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                         *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
                         *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
         end do
         zhbh(1)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
              *work(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))

         ! Set proportionality factor 'a'
         
         za=sqrt(zhbh(1)/zwbw(1))
         ! if (abs(obsoper%lat*180./3.1415-78.).lt.2.0.and.abs(rlon*180./3.1415-185.).lt.2.0) then
         !     write(6,*) 'ZA  ',obsoper%lat*180.0/3.1415,rlon*180.0/3.1415,za,zhbh(1),zwbw(1)
         !     write(6,*) 'obsoper%trial',obsoper%trial(1:obsoper%nmodlev)
         !     write(6,*) 'sigma_trial',sigma_trial(1:obsoper%nmodlev,1)
         !     write(6,*) 'ZH  ',chm_obsoper%zh(iobslev,1:obsoper%nmodlev)
         !     write(6,*) 'ZHP ',chm_obsoper%zhp(iobslev,1:obsoper%nmodlev)*za
         ! end if
         
         ! Set final innovation operator
         
         obsoper%zh(iobslev,1:obsoper%nmodlev)=obsoper%zhp(iobslev,1:obsoper%nmodlev)*za
         
      end do
      
  end subroutine chm_genoper

!--------------------------------------------------------------------------
!! *Purpose*: Multiplys a matrix by the covariance matrix or its inverse
!!            (or combinations of these).
!!
!! @author Y. Rochon, April 2015
!!
!! Input
!!
!!v    itype                    type of multiplication, for an input matrix A the output matrix is
!!v                               0) D(i,j)=A(i,j)/sum(C(1:n,i))
!!v                               1) A*C
!!v                               2) C*A
!!v                               3) A*C*A^T
!!v                              -1) A*CI
!!v                              -2) CI*A
!!v                              -3) A*CI*A^T
!!v                             where C is the covariance matrix and CI is its inverse
!!v    varName                  variable name
!!v    rmat_in(ndim1,ndim2)     input matrix/vector A (see comments sections)
!!v    imodlev_top(ndim1)       top level of non-zero values in rmat_in
!!v    imodlev_bot(ndim1)       bottom level of non-zero values in rmat_in
!!v    ndim1,ndim2              matrix dimensions (ndim1 to be 1 one 1D input vectors)
!!v    ndim3                    expected output dimension
!!v                                   =ndim1 for itype= +/-3
!!v                                   =ndim2 otherwise
!!v    lrgsig                   index to indicate if rgsig to be included as part of C or CI below.
!!v    rsig(ndim2,2)            background error std. dev. at obs locations (must be provided when lrgsig=.true.)
!!
!! Output
!!
!!v    rmat_out(ndim1,ndim2)    output matrix/vector
!!
!! Comments:
!!
!!v   - If rmat_in is a 1-D vector, then
!!v        for cases +/- 2, one should have set ndim2=1 and ndim1=vector-length.
!!v        for cases +/- 1,3, one should have set ndim1=1 and ndim2=vector-length.
!!v
!!v   - Revisions required whem LAM and ensembles cases become available.
!--------------------------------------------------------------------------
  subroutine chm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot, &
                              ndim1,ndim2,ndim3,lrgsig,itype,rsig_opt)
 
      implicit none

      character(len=*), intent(in) :: varName
      logical, intent(in)    :: lrgsig
      integer, intent(in)    :: ndim1,ndim2,ndim3,itype
      integer, intent(in)    :: imodlev_top(ndim1),imodlev_bot(ndim1)
      real(8), intent(in)    :: rmat_in(ndim1,ndim2)
      real(8), intent(out)   :: rmat_out(ndim1,ndim3)
      real(8), intent(in), optional :: rsig_opt(ndim2,2)
   
      integer :: nsize
      real(8) :: rsig(ndim2,2)
      
      rmat_out(:,:)=0.0D0
      rsig=0.0
      if (present(rsig_opt)) rsig=rsig_opt
        
      ! Apply operation related to static background error covariance/correlation matrix.
      ! Applicability tests within b*chm_corvert_mult.

      call bchm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
                             lrgsig,itype,rsig(:,1))
!      call blamchm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
!                             lrgsig,itype,rsig(:,1))

      ! Apply operation related to ensemble-based background error covariance/correlation matrix.
      ! Applicability test within benschm_corvert_mult.

!      call benschm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
!                                lrgsig,itype,rsig(:,2))

  end subroutine chm_corvert_mult
  
!--------------------------------------------------------------------------
!! *Purpose*: Determine and store the boundary (e.g. tropopause or PBL) pressure levels if needed by
!!            the observation operators.
!!
!! @author Y. Rochon, Oct 2015
!!
!! Input
!!
!!v      iconstituent_id   BUFR code element of Table 08046 identifying the constituent.
!!v      nmodlev           number of model levels for variables other than uu and vv
!!v      pressmod          model pressure array
!!v      tt                model temperature (Kelvin)
!!v      gz                model geopotential height (meters)
!!v      hu                specific humidity 
!!v      uu                model zonal wind component (m/s)
!!v      vv                model meridional wind component (m/s)
!!
!! Output
!!
!!v      bound_press       pressure level of boundary to be imposed
!--------------------------------------------------------------------------
  function chm_get_col_boundary(iconstituent_id,nmodlev,pressmod,tt,gz,hu_opt,uu_opt,vv_opt) result(bound_press)

    implicit none

    integer, intent(in) :: nmodlev,iconstituent_id
    real(8), intent(in) :: pressmod(nmodlev),tt(nmodlev),gz(nmodlev)
    real(8), optional, intent(in) :: uu_opt(:),vv_opt(:),hu_opt(nmodlev)
   
    real(8) :: bound_press
    integer :: tropo_bound
    
    bound_press = -1.
    
    if (chm_setup_get_int('tropo_mode',iconstituent_id).eq.0) return
   
    tropo_bound=chm_setup_get_int('tropo_bound',iconstituent_id)
    
    if (tropo_bound.gt.0) then
       if (.not.all(tt.lt.0.) .and. .not.all(gz.lt.0.) ) then
    
          select case(tropo_bound)
          case(1)
    
             ! Get tropopause pressure level
      
             bound_press = phf_get_tropopause(nmodlev,pressmod,tt,gz,hu_opt=hu_opt)
    
          case(2)
 
             ! Get PBL pressure level
      
             bound_press = phf_get_pbl(nmodlev,pressmod,tt,gz,hu_opt=hu_opt,uu_opt=uu_opt,vv_opt=vv_opt) 
      
          case default
             call utl_abort("chm_get_col_boundary: Unrecognized value for tropo_bound of " &
                  // trim(utl_str(tropo_bound)) )
          end select
                                     
      end if     
    end if
      
    ! Use tropo_column_top value if tropo_bound=0 or model derived boundary was unsuccessful
    if (bound_press.lt.0.0) &
         bound_press = chm_setup_get_float('tropo_column_top',iconstituent_id)
      
  end function chm_get_col_boundary
       
!--------------------------------------------------------------------------
!! *Purpose*: Adds column boundary data to chm_column_boundary which can be retrieved later
!!            using a header index.
!!
!! @author M. Sitwell, April 2016
!--------------------------------------------------------------------------
  subroutine chm_add_col_boundary(headerIndex,bound_press)

    implicit none 
    
    integer, intent(in) :: headerIndex
    real(8), intent(in) :: bound_press
    
    integer obsdata_maxsize
    
    obsdata_maxsize = chm_setup_get_int('obsdata_maxsize')   
     
    if (.not.associated(chm_column_boundary%data1d)) then
       call oss_obsdata_alloc(chm_column_boundary,obsdata_maxsize, dim1=1)
       chm_column_boundary%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    chm_column_boundary%nrep = chm_column_boundary%nrep+1 

    if (chm_column_boundary%nrep.gt.obsdata_maxsize) &
         call utl_abort('chm_add_col_boundary: Reach max size of array ' // trim(utl_str(obsdata_maxsize)) )
  
    ! Use the header number as the unique code for this obs data
    write(chm_column_boundary%code(chm_column_boundary%nrep),'(I22)') headerIndex

    chm_column_boundary%data1d(1,chm_column_boundary%nrep) = bound_press

  end subroutine chm_add_col_boundary

!--------------------------------------------------------------------------
!! *Purpose*: Retrieves previously saved column boundary data in chm_column_boundary from
!!            the header index.
!!
!! @author M. Sitwell, April 2016
!--------------------------------------------------------------------------
  function chm_retrieve_col_boundary(headerIndex) result(bound_press)

    implicit none

    integer, intent(in) :: headerIndex
    real(8) :: bound_press
    character(len=22) :: code

    write(code,'(I22)') headerIndex
    
    bound_press = oss_obsdata_get_element(chm_column_boundary,code,1)

  end function chm_retrieve_col_boundary

!--------------------------------------------------------------------------
!! *Purpose*: Adds background sigma profiles (and inverse) to chm_sigma_trial which can be retrieved later
!!            using a header index.
!!
!! @author Y. Rochon, May 2016 (based on chm_add_col_boundary by M. Sitwell)
!--------------------------------------------------------------------------
  subroutine chm_add_sigma_trial(headerIndex,sigma)

    implicit none 
    
    integer, intent(in) :: headerIndex
    real(8), intent(in) :: sigma(:,:)
   
    integer :: obsdata_maxsize
    
    obsdata_maxsize = chm_setup_get_int('obsdata_maxsize')
     
    if (.not.associated(chm_sigma_trial%data2d)) then
       call oss_obsdata_alloc(chm_sigma_trial, obsdata_maxsize, dim1=size(sigma,dim=1), dim2_opt=max(4,size(sigma,dim=2)))
       chm_sigma_trial%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    chm_sigma_trial%nrep = chm_sigma_trial%nrep+1 

    if (chm_sigma_trial%nrep.gt.obsdata_maxsize) &
         call utl_abort('chm_sigma_trial: Reached max size of array ' // trim(utl_str(obsdata_maxsize)) )
  
    ! Use the header number as the unique code for this obs data
    write(chm_sigma_trial%code(chm_sigma_trial%nrep),'(I22)') headerIndex

    chm_sigma_trial%data2d(:,1:2,chm_sigma_trial%nrep) = sigma(:,1:2)

    where (sigma(:,1).gt.0.0D0)
       chm_sigma_trial%data2d(:,3,chm_sigma_trial%nrep) = 1.0D0/sigma(:,1)
    elsewhere
       chm_sigma_trial%data2d(:,3,chm_sigma_trial%nrep) = 0.0D0
    end where

    where (sigma(:,2).gt.0.0D0)
       chm_sigma_trial%data2d(:,4,chm_sigma_trial%nrep) = 1.0D0/sigma(:,2)
    elsewhere
       chm_sigma_trial%data2d(:,4,chm_sigma_trial%nrep) = 0.0D0
    end where

  end subroutine chm_add_sigma_trial

!--------------------------------------------------------------------------
!! *Purpose*: Retrieves previously saved background sigma profiles chm_sigma_trial from
!!            the header index.
!!
!! @author Y. Rochon, May 2016 (based on chm_retrieve_col_boundary by M. Sitwell)
!--------------------------------------------------------------------------
  function chm_retrieve_sigma_trial(headerIndex) result(sigma)

    implicit none

    integer, intent(in) :: headerIndex
    real(8) :: sigma(chm_sigma_trial%dim1,chm_sigma_trial%dim2)
    character(len=22) :: code

    write(code,'(I22)') headerIndex
    
    sigma = oss_obsdata_get_array2d(chm_sigma_trial,code)

  end function chm_retrieve_sigma_trial

end module chem_obsoperators_mod
