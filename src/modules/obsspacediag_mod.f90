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

module obsSpaceDiag_mod
  ! MODULE obsSpaceDiag_mod (prefix="osd" category='1. High-level functionality')
  !
  ! :Purpose: Some experimental procedures for computing various diagnostics in
  !           observation space.
  !
  use codePrecision_mod
  use midasMpi_mod
  use bufr_mod
  use codtyp_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use controlVector_mod
  use obsSpaceData_mod
  use columnData_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use bMatrix_mod
  use bMatrixHi_mod
  use bMatrixEnsemble_mod
  use bCovarSetupChem_mod
  use varNameList_mod
  use stateToColumn_mod
  use randomNumber_mod
  use obsOperators_mod
  use utilities_mod
  use physicsfunctions_mod
  use obsSubSpaceData_mod
  use obsfiles_mod
  use obsOperatorsChem_mod
  use obsFamilyList_mod
  
  implicit none
  save
  private
  
  ! public procedures
  !------------------
  public :: osd_ObsSpaceDiag

  type :: struct_osd_diagn

     ! Structure for storing observation-space diagnostic arrays
     !
     ! Real arrays indexed by (lat,lon,lev,stat), where stat=1 for RMS and stat=2 for mean for all except Jo_stats.
     ! For Jo_stats, stat=1 is Jo for x=x_analysis and stat=2 is Jo for x=x_background.
     ! The array counts is indexed by (lat,lon,lev).
     ! The array status_count is indexed (lat,lon,lev,status), where the status number ranges from 0 to 2.
     !     
     !  Variable               Description
     !  --------               -----------
     !  OmP_stats              obs - background statistics
     !  OmA_stats              obs - analysis statistics
     !  obs_stats              observation statistics
     !  Jo_stats               cost function statistics 
     !                         - First four elements (of last dimension) apply only prescribed obs error variances 
     !                           as normalizing-scaling denomicators.
     !                         - Fifth element (of last dimension) applies the sum of the prescribed obs error variances 
     !                           and the diag(HPHT) as normalizing-scaling denominator. This does not include consideration  
     !                           of spatial correlations of (O-P) between obs points 
     !                           associated to HPHT in the normalizing denominator.
     !  Jpa_stats              statistics (P-A) in observation space.
     !                         - Not exactly equivalent to Jb of the cost function 
     !                         - Applies the prescribed diag(HPHT) as normalizing-scaling denomiator.     
     !                         - Does not include consideration of spatial correlations of (P-A) between obs points 
     !                           associated to HPHT in normalizing denominator
     !                         - Vert. coordinate binning included but not currently output.
     !  diagR_stats            Elements for the calc of mean{[(O-P)-mean(O-P)][(O-A)-mean(O-A)]} 
     !                         (with each O-P and O-A difference divided by sigma_obs)
     !                         for scaling factor adjustments of obs std. dev via the Desroziers approach.
     !  diagHPHT_stats         Elements for the calc of mean{[(O-P)-mean(O-P)][[(O-P)-mean(O-P)]-[(O-A)-mean(O-A)]]}
     !                         = mean{[(O-P)-mean(O-P)][(A-P)-mean(A-P)]}(with each O-P and O-A difference divided by sqrtHPHT)
     !                         for scaling factor adjustments of background std. dev. in obs space via the Desroziers approach.
     !  counts                 number of observations in a (lat,lon,lev) bin
     !  nlat                   number of latitude levels
     !  nlon                   number of longitude levels
     !  nlev                   number of vertical levels
     !  nbin                   total number of (lat,lon,lev) bins
     !  nstat                  number of different statistics types
     !  nstatus                indicates the number of observations with a certain status,
     !                         with status values denoting:
     !                           0 - observation has been rejected and not included for diagnostics
     !                           1 - observation has been assimilated
     !                           2 - observation has been used for diagnostics only (not assimilated)
     !  deltaLat               latitude bin size (deg)
     !  deltaLon               longitude bin size (deg)
     !  deltaLogPressure       vertical bin size in log(pressure)
     !  allow_print_summary    indicates fs printing of summary diagnostics is allowed
     !  assim_mode             indicates if assimilation was performed for this dataset

     real(8), allocatable :: OmP_stats(:,:,:,:),OmA_stats(:,:,:,:),obs_stats(:,:,:,:),Jo_stats(:,:,:,:)
     real(8), allocatable :: diagR_stats(:,:,:,:),diagHPHT_stats(:,:,:,:), Jpa_stats(:,:,:)
     integer, allocatable :: counts(:,:,:),nstatus(:,:,:,:)
     
     integer :: nlev,nlat,nlon,nbin,nstat
     real(8) :: deltaLat,deltaLon,deltaLogPressure
     logical :: allow_print_summary,assim_mode

  end type struct_osd_diagn

  ! namelist variables
  integer, parameter :: max_cfg_size=100
  real(8) :: deltaLat,deltaLon,deltaPressure,deltaHeight
  integer :: numFamily
  character(len=2) :: familyList(ofl_numFamily)
  integer :: numElement
  integer :: elementList(ofl_numFamily) 
  integer :: nrandseed
  logical :: lrandom  

  !  Parameters identifying obs sets and related actions for diagnostic calcs of each family
  ! 
  !  diagn_num           Prescribed (starting) number of (stnid, bufr, nlev) for the diagnostics calc
  !
  !  diagn_stnid         Prescribed (starting) list of stnid (with *s as needed) for the diagnostics calc
  !                      with '*' denoting wild cards  
  !
  !  diagn_varno         Prescribed (starting) list of bufr elements for the diagnostics calc  
  !
  !  diagn_unilev        Prescribed (starting) list of logicals indicating uni-level obs for the diagnostics calc
  !
  !  diagn_pressmin      Bottom of top layer for diagnostics (in Pa). 
  !
  !  diagn_save          Logical indicating gridded diagnostics are to be saved
  !                      in an ascii file in addition to overall diagnostics. 
  !
  !  diagn_nset          Integer indicating grouping of diagnostics with
  !                      1: group by stnid
  !                      2: group by (stnid,bufr)
  !                      3: group by (stnid,bufr,nlev)
  ! 
  !  diagn_all           Logical indicating if all combinations specified by diagn_nset are to be
  !                      used in diagnostics or only those specified by the diagn_* arrays
  ! 
  !  obsspace_diagn_filename 
  !                     File name for file containing obs space diagnostics related to chemical constituents.

  integer :: diagn_num(ofl_numFamily),diagn_nset(ofl_numFamily)
  integer :: diagn_varno(ofl_numFamily,max_cfg_size)
  character(len=9) :: diagn_stnid(ofl_numFamily,max_cfg_size)
  logical :: diagn_save(ofl_numFamily),diagn_all(ofl_numFamily),diagn_unilev(ofl_numFamily,max_cfg_size)
  real(8) :: diagn_pressmin(ofl_numFamily)
  character(len=100) :: obsspace_diagn_filename(ofl_numFamily)

contains

  !--------------------------------------------------------------------------
  ! osd_ObsSpaceDiag
  !--------------------------------------------------------------------------
  subroutine osd_ObsSpaceDiag( obsSpaceData, columnTrlOnAnlIncLev, hco_anl, analysisMode_opt )
    !           
    ! :Purpose: Calls routines to perform observation-space diagnostic tasks
    !
    ! :Arguments:
    !           :obsSpaceData: Obs space data structure
    !           :columnTrlOnAnlIncLev: Structure of vertical columns at obs locations.
    !                                  Expected to be for analysis vertical levels if to be used.
    !           :analysisMode: logical indicating if following analysis mode or not (optional)
    !                          Assumed .true. if not present.
    !

    implicit none
    
    !Arguments:
    type(struct_obs)              :: obsSpaceData
    type(struct_columnData)       :: columnTrlOnAnlIncLev
    type(struct_hco), pointer     :: hco_anl
    logical, intent(in), optional :: analysisMode_opt
    
    !Locals:
    logical :: nmlExists,anlm_mod    
    integer :: ierr
    integer :: dateprnt,timeprnt,newdate
    
    ! write(*,*) 'osd_ObsSpaceDiag: Starting'

    if (present(analysisMode_opt)) then
       anlm_mod = analysisMode_opt
    else
       anlm_mod = .true.
    end if
    
    call osd_setup(nmlExists)
    ierr = newdate(tim_getDatestamp(),dateprnt,timeprnt,-3)
    dateprnt=dateprnt*100+timeprnt/1000000
    
    ! Perform diagnostics based on OmP (and OmA if available)
 
    call osd_obsPostProc(obsSpaceData,columnTrlOnAnlIncLev,deltaLat,deltaLon,deltaPressure,anlm_mod)
    
    if ((.not. anlm_mod) .or. (.not.lrandom) .or. (.not.nmlExists)) return

    ! Perform diagnostics from random perturbations
    
    call osd_calcInflation(obsSpaceData,columnTrlOnAnlIncLev,hco_anl,dateprnt)

   ! write(*,*) 'osd_obsspace_diag: Completed'

  end subroutine osd_ObsSpaceDiag
  
  !--------------------------------------------------------------------------
  ! osd_calcInflation
  !--------------------------------------------------------------------------
  subroutine osd_calcInflation( obsSpaceData, columnTrlOnAnlIncLev, hco_anl, dateprnt )
    !      
    ! :Purpose: Calculates observation-space diagnostics from random perturbations
    !

    implicit none
 
    ! Arguments:
    type(struct_obs)            :: obsSpaceData
    type(struct_columnData)     :: columnTrlOnAnlIncLev
    type(struct_hco), pointer   :: hco_anl
    integer                     :: dateprnt

    ! Locals:
    type(struct_gsv)            :: statevector
    type(struct_columnData)     :: column
    type(struct_vco), pointer   :: vco_anl
    real(8), allocatable,target :: controlVector(:)

    integer :: familyIndex,elementIndex,bodyIndex,headerIndex,latIndex,lonIndex,verticalIndex
    integer :: maxLat,maxLon,maxVertical
    real(8), allocatable :: innovStd(:,:,:),bmatHiStd(:,:,:),bmatEnStd(:,:,:)
    integer, allocatable :: counts(:,:,:)

    real(8), allocatable :: my_innovStd(:,:,:),my_bmatHiStd(:,:,:),my_bmatEnStd(:,:,:)
    integer, allocatable :: my_counts(:,:,:)
    
    integer :: ierr,nulinnov,nulBmatHi,nulBmatEn,nulcount,fnom,fclos,ivco,ivco_recv,iseed,jj,jlev,jvar
    integer :: ivar_count,nlev_max
    logical :: lpert_static, lpert_ens
    real(8), pointer         :: cvBhi(:), cvBen(:), cvBchm(:)
    real(pre_incrReal), pointer :: field(:,:,:,:)
    real(8), allocatable     :: HxBhi(:), HxBen(:)
    real(8), allocatable     :: scaleFactor(:),scaleFactorChm(:,:)
    character(len=128) :: innovFileName,bmatHiFileName,bmatEnFileName,countFileName
    character(len=6)   :: elementStr
    character(len=10)  :: dateStr
    
    write(*,*) 'osd_calcInflation: Starting'

    if( nrandseed == 999 ) nrandseed=dateprnt ! if seed not set by namelist, use valid date/time
    write(*,*) 'osd_calcInflation: random seed set to ',nrandseed

    maxLat = nint(180.0d0/deltaLat)
    maxLon = nint(360.0d0/deltaLon)
    maxVertical = max(1+nint(110000.0d0/deltaPressure),1+nint(80000.0d0/deltaHeight),200)

    write(*,*) 'osd_calcInflation: Compute random realization of background error'

    ! allocate vectors to store Hx for the static and ensemble-based covariance matrices
    allocate(HxBhi(obs_numbody(obsSpaceData)))
    allocate(HxBen(obs_numbody(obsSpaceData)))

    ! initialize columnData object for increment
    call col_setVco(column,col_getVco(columnTrlOnAnlIncLev))
    call col_allocate(column,col_getNumCol(columnTrlOnAnlIncLev),mpiLocal_opt=.true.)

    ! initialize gridStateVector object for increment
    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, mpi_local_opt=.true.)

    nlev_max=max(col_getNumLev(columnTrlOnAnlIncLev,'MM'),col_getNumLev(columnTrlOnAnlIncLev,'TH'))

    allocate(controlVector(cvm_nvadim))
    allocate(scaleFactor(nlev_max))
    allocate(scaleFactorChm(100,nlev_max))

    ! COMPUTE BMATRIX PERTURBATION FOR THE STATIC COVARIANCES CASE; from Bhi and or BChm 

    if (all(familyList(1:numFamily).eq.'CH').or.(.not.cvm_subVectorExists('B_HI'))) then
       cvBhi => null()
    else
       cvBhi => cvm_getSubVector(controlVector,'B_HI')
    end if
    if (any(familyList(1:numFamily).eq.'CH').and.obs_famExist(obsSpaceData,'CH').and.cvm_subVectorExists('B_CHM')) then
       cvBChm => cvm_getSubVector(controlVector,'B_CHM')
    else
       cvBChm => null()
    end if

    HxBhi(:) = 0.0d0
    
    iseed = abs(nrandseed)
    call rng_setup(iseed)

    if (cvm_subVectorExists('B_HI').or.cvm_subVectorExists('B_CHM')) then
       
       ! compute random control vector

       controlVector(:) = 0.0d0

       if (cvm_subVectorExists('B_HI')) then
          do jj = 1,size(cvBhi)
             cvBhi(jj)=rng_gaussian()
          enddo
          ! initialize vector of scaleFactors
          call bhi_getScaleFactor(scaleFactor)
       else
          scaleFactor(:)=1.0      
       end if

       if (cvm_subVectorExists('B_CHM')) then
          do jj = 1,size(cvBChm)
             cvBChm(jj)=rng_gaussian()
          enddo
          ! initialize vector of scaleFactors
          call bcsc_getScaleFactor(scaleFactorChm)
       else
          scaleFactorChm(:,:)=1.0
       end if
 
       ! multiply vector by B^1/2
       call bmat_sqrtB(controlVector,cvm_nvadim,statevector)

       ! undo the scaleFactor (THIS IS NOT CORRECT FOR 2D VARIABLES!!! (P0 and TG) ) 

       ivar_count=0
       do jvar=1,vnl_numvarmax 
          if(gsv_varExist(statevector,vnl_varNameList(jvar))) then

             call gsv_getField(statevector,field,vnl_varNameList(jvar))

             if (vnl_varKindFromVarname(vnl_varNameList(jvar)) == 'MT') then
                do jlev = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))   
                   if(scaleFactor(jlev) > 0.0d0) field(:,:,jlev,:)=field(:,:,jlev,:)/scaleFactor(jlev)
                enddo
             else if (vnl_varKindFromVarname(vnl_varNameList(jvar)) == 'CH') then
                ivar_count=ivar_count+1
                do jlev = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))   
                   if(scaleFactorChm(ivar_count,jlev).gt.0.0d0) field(:,:,jlev,:)=field(:,:,jlev,:) &
                         /scaleFactorChm(ivar_count,jlev)
                end do
             end if
          endif
       enddo

       ! multiply by H
       call s2c_tl(statevector,column,columnTrlOnAnlIncLev,obsSpaceData)  ! put in column H_horiz dx
       call oop_Htl(column,columnTrlOnAnlIncLev,obsSpaceData,1)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
       do bodyIndex=1,obs_numBody(obsSpaceData)
          HxBhi(bodyIndex) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
       enddo

    end if

    ! COMPUTE BMATRIX PERTURBATION FOR THE ENSEMBLE COVARIANCES CASE; from Ben

    if (cvm_subVectorExists('B_ENS')) then

       cvBen => cvm_getSubVector(controlVector,'B_ENS')
       HxBen(:) = 0.0d0

       ! compute random control vector

       controlVector(:) = 0.0d0

       do jj = 1,size(cvBen)
          cvBen(jj)=rng_gaussian()
       enddo

       ! initialize vector of scaleFactors
       call ben_getScaleFactor(scaleFactor)
       
       scaleFactorChm(:,:)=1.0

       ! multiply vector by B^1/2
       call bmat_sqrtB(controlVector,cvm_nvadim,statevector)

       ! undo the scaleFactor
       
       ivar_count=0
       do jvar=1,vnl_numvarmax 
          if(gsv_varExist(statevector,vnl_varNameList(jvar))) then

             call gsv_getField(statevector,field,vnl_varNameList(jvar))

             if (vnl_varKindFromVarname(vnl_varNameList(jvar)) == 'MT') then
                do jlev = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))   
                   if(scaleFactor(jlev) > 0.0d0) field(:,:,jlev,:)=field(:,:,jlev,:)/scaleFactor(jlev)
                enddo
             else if (vnl_varKindFromVarname(vnl_varNameList(jvar)) == 'CH') then
                ivar_count=ivar_count+1
                do jlev = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))   
                   if(scaleFactorChm(ivar_count,jlev) > 0.0d0) field(:,:,jlev,:)=field(:,:,jlev,:) &
                         /scaleFactorChm(ivar_count,jlev)
                end do
             end if
          endif
       enddo
              
       ! multiply vector by H
       call s2c_tl(statevector,column,columnTrlOnAnlIncLev,obsSpaceData)  ! put in column H_horiz dx
       call oop_Htl(column,columnTrlOnAnlIncLev,obsSpaceData,1)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
       do bodyIndex=1,obs_numBody(obsSpaceData)
          HxBen(bodyIndex) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
       enddo
    else
       cvBen => null()
       HxBen(:) = 0.0d0
    end if

    deallocate(scaleFactor)
    deallocate(scaleFactorChm)
    call col_deallocate(column)
    call gsv_deallocate(statevector)
       
    allocate(my_innovStd(maxLat,maxLon,maxVertical))
    allocate(my_bmatHiStd(maxLat,maxLon,maxVertical))
    allocate(my_bmatEnStd(maxLat,maxLon,maxVertical))
    allocate(my_counts(maxLat,maxLon,maxVertical))

    allocate(innovStd(maxLat,maxLon,maxVertical))
    allocate(bmatHiStd(maxLat,maxLon,maxVertical))
    allocate(bmatEnStd(maxLat,maxLon,maxVertical))
    allocate(counts(maxLat,maxLon,maxVertical))

    FAMILY: do familyIndex = 1, numFamily

      ELEMENT: do elementIndex = 1, numElement  

        ! Initialize logicals for calc of perturbation diagnostics.

        lpert_static=.false.
        lpert_ens=.false.
        
        if (familyList(familyIndex) /= 'CH') then
           if (cvm_subVectorExists('B_HI')) lpert_static=.true.
        else        
           if (cvm_subVectorExists('B_CHM')) lpert_static=.true.
        end if
        if (cvm_subVectorExists('B_ENS')) lpert_ens=.true.
      
        ivco = -999
        my_innovStd(:,:,:) = 0.0d0
        my_bmatHiStd(:,:,:) = 0.0d0
        my_bmatEnStd(:,:,:) = 0.0d0
        my_counts(:,:,:) = 0

        call obs_set_current_body_list(obsSpaceData,familyList(familyIndex))
        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          if(obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == elementList(elementIndex) .and. &
             obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then

            call osd_getIndices(obsSpaceData,bodyIndex,latIndex,lonIndex,verticalIndex)
            
            if(verticalIndex == -1) then
               ! skip this obs for whatever reason
               cycle BODY
            else if(latIndex > maxLat .or. lonIndex > maxLon .or. verticalIndex > maxVertical) then
               write(*,*) 'osd_calcInflation: index too big: lat,lon,vertical=',latIndex,lonIndex,verticalIndex, &
                          ' lat_max,lon_max,vertical_max=',maxlat,maxlon,maxvertical
               call utl_abort('osd_calcInflation')
            endif

            ivco = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)
            my_counts(latIndex,lonIndex,verticalIndex) = my_counts(latIndex,lonIndex,verticalIndex) + 1
            my_innovStd(latIndex,lonIndex,verticalIndex) = my_innovStd(latIndex,lonIndex,verticalIndex) +     &
                                                        obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)* &
                                                        obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
            if (lpert_static) my_bmatHiStd(latIndex,lonIndex,verticalIndex)  = my_bmatHiStd(latIndex,lonIndex,verticalIndex) +     &
                                                          HxBhi(bodyIndex)*HxBhi(bodyIndex)
            if (lpert_ens) my_bmatEnStd(latIndex,lonIndex,verticalIndex)  = my_bmatEnStd(latIndex,lonIndex,verticalIndex) +     &
                                                          HxBen(bodyIndex)*HxBen(bodyIndex)

            headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)

          endif
        enddo BODY

        call rpn_comm_allreduce(ivco,ivco_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ierr)
        ivco = ivco_recv

        call rpn_comm_allreduce(my_counts,counts,maxLat*maxLon*maxVertical,"MPI_INTEGER","MPI_SUM","GRID",ierr)
        call rpn_comm_allreduce(my_innovStd,innovStd,maxLat*maxLon*maxVertical,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
        if (lpert_static) call rpn_comm_allreduce(my_bmatHiStd,bmatHiStd,maxLat*maxLon*maxVertical,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
        if (lpert_ens) call rpn_comm_allreduce(my_bmatEnStd,bmatEnStd,maxLat*maxLon*maxVertical,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)

        where (counts > 0) innovStd = sqrt(innovStd/counts)
        if (lpert_static) then
           where (counts > 0) bmatHiStd = sqrt(bmatHiStd/counts)
        end if
        if (lpert_ens) then
           where (counts > 0) bmatEnStd = sqrt(bmatEnStd/counts)
        end if 
        
        if(mmpi_myid == 0 .and. sum(counts(:,:,:)) > 0) then

         ! determine file names
          write(dateStr,'(i10.10)') dateprnt
          write(elementStr,'(i6.6)') elementList(elementIndex)
          innovFileName = 'innov' // dateStr // '_'  // trim(familyList(familyIndex)) // '_' // trim(elementStr) // '.dat'
          if (lpert_static) bmatHiFileName =  'bmathi'  // dateStr // '_'  // trim(familyList(familyIndex)) // '_' // trim(elementStr) // '.dat'
          if (lpert_ens) bmatEnFileName =  'bmaten'  // dateStr // '_'  // trim(familyList(familyIndex)) // '_' // trim(elementStr) // '.dat'
          countFileName = 'count' // dateStr // '_'  // trim(familyList(familyIndex)) // '_' // trim(elementStr) // '.dat'

          ! open files
          nulinnov=0
          nulBmatHi =0
          nulBmatEn =0
          nulcount=0
          ierr = fnom(nulinnov,innovFileName,'FMT+R/W',0)
          if (lpert_static) ierr = fnom(nulBmatHi ,bmatHiFileName ,'FMT+R/W',0)
          if (lpert_ens) ierr = fnom(nulBmatEn ,bmatEnFileName ,'FMT+R/W',0)
          ierr = fnom(nulcount,countFileName,'FMT+R/W',0)

          ! write data for this family/element
          write(nulinnov,*) '***maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight='
          write(nulinnov,*) maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight
          
          if (lpert_static) then 
             write(nulBmatHi,*)  '***maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight='
             write(nulBmatHi,*)  maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight
          end if
          if (lpert_ens) then 
             write(nulBmatEn,*)  '***maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight='
             write(nulBmatEn,*)  maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight
          end if
          
          write(nulcount,*) '***maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight='
          write(nulcount,*) maxLon,maxLat,deltaLon,deltaLat,deltaPressure,deltaHeight
          do verticalIndex = 1,maxVertical
            if(sum(counts(:,:,verticalIndex)).gt.0) then
              write(nulinnov,*) '***verticalIndex,vco='
              write(nulinnov,*) verticalIndex,ivco
              
              if (lpert_static) then 
                 write(nulBmatHi,*)  '***verticalIndex,vco='
                 write(nulBmatHi,*)  verticalIndex,ivco
              end if
              if (lpert_ens) then 
                 write(nulBmatEn,*)  '***verticalIndex,vco='
                 write(nulBmatEn,*)  verticalIndex,ivco
              end if
              
              write(nulcount,*) '***verticalIndex,vco='
              write(nulcount,*) verticalIndex,ivco
              do latIndex = 1,maxLat
                write(nulinnov,*) innovStd(latIndex,:,verticalIndex)
                write(nulcount,*) counts(latIndex,:,verticalIndex)
              enddo
 
              if (lpert_static) then 
                 do latIndex = 1,maxLat
                   write(nulBmatHi ,*) bmatHiStd(latIndex,:,verticalIndex)
                 enddo
              end if
              if (lpert_ens) then 
                 do latIndex = 1,maxLat
                    write(nulBmatEn ,*) bmatEnStd(latIndex,:,verticalIndex)
                 enddo
              end if

            endif
          enddo

          ! close files
          ierr = fclos(nulinnov)
          if (lpert_static) ierr = fclos(nulBmatHi)
          if (lpert_ens) ierr = fclos(nulBmatEn)
          ierr = fclos(nulcount)
        endif
      enddo ELEMENT

    enddo FAMILY

    deallocate(my_counts) 
    deallocate(my_innovStd)  
    deallocate(my_bmatHiStd)  
    deallocate(my_bmatEnStd)  

    deallocate(innovStd)
    deallocate(bmatHiStd)
    deallocate(bmatEnStd)
    deallocate(counts)

    deallocate(HxBhi)
    deallocate(HxBen)

    write(*,*) 'osd_calcInflation: Finished'

  end subroutine osd_calcInflation

  !--------------------------------------------------------------------------
  ! osd_getIndices
  !--------------------------------------------------------------------------
  subroutine osd_getIndices( obsSpaceData, bodyIndex, latIndex, lonIndex, verticalIndex )
    !
    implicit none

    !Arguments:
    type(struct_obs) :: obsSpaceData
    integer          :: bodyIndex
    integer          :: latIndex
    integer          :: lonIndex
    integer          :: verticalIndex

    !Locals:
    real(8), parameter :: epsilon=0.001
    integer            :: headerIndex, bodyElem_i

    ! codtypes for tovs: 164(AMSUA) 168 180 181 182 183 185 186 192 193

    ! epsilon is added below to handle case where lon/dlon~1 or lat/dlat~1
    headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
    latIndex = 1 + floor((90.0d0 + obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)*MPC_DEGREES_PER_RADIAN_R8)/deltaLat - epsilon)
    if (latIndex == 0) latIndex=1
    lonIndex = 1 + floor(obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)*MPC_DEGREES_PER_RADIAN_R8/deltaLon - epsilon)
    if (lonIndex == 0) lonIndex=1

    select case(obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex))
      case(1)
        ! height coordinate
        verticalIndex = 1 + nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)/deltaHeight)
      case(2)
        ! pressure coordinate
        verticalIndex = 1 + nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)/deltaPressure)
      case(3)
        ! channel number
        verticalIndex = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
        if(obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == codtyp_get_codtyp("AMSUA")) then
          ! amsu-a
          verticalIndex = verticalIndex - 27
        else
          ! ignore other types of TOVS for now
          verticalIndex = -1
        endif
      case(4,5)
         ! Integrated column or surface value - assign to first level
         verticalIndex = 1
      case default
        ! unknown vertical coordinate
        bodyElem_i = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)
        write(*,*) 'osd_getIndices: Unknown VCO! ', bodyElem_i
        verticalIndex = -1
    end select
 
  end subroutine osd_getIndices

  !--------------------------------------------------------------------------
  ! osd_setup
  !--------------------------------------------------------------------------
  subroutine osd_setup(nmlExists) 
    !
    implicit none

    !Arguments:
    logical :: nmlExists

    !Locals:
    integer :: nulnam,ierr,fnom,fclos,j
    
    namelist /namosd/nrandseed,deltaLat,deltaLon,deltaPressure,deltaHeight, &
        numFamily,familyList,numElement,elementList,lrandom, &
        diagn_save,diagn_all,diagn_num,diagn_stnid,diagn_varno,diagn_unilev,     &
        diagn_nset,diagn_pressmin        

    ! set default values for namelist variables
    nrandseed = 999
    lrandom=.true.
    deltaLat = 10.0d0
    deltaLon = 10.0d0
    deltaPressure = 10000.0d0
    deltaHeight = 5000.0d0

    !numFamily = ofl_numFamily
    !familyList(:) = ofl_familyList(:)
    numFamily = 7
    familyList(:) = '  '
    familyList(1) = 'UA'
    familyList(2) = 'AI'
    familyList(3) = 'SC'
    familyList(4) = 'RO'
    familyList(5) = 'TO'
    familyList(6) = 'SW'
    familyList(7) = 'SF'
    
    numElement = 11
    elementList(:) = -1
    elementList(1) = BUFR_NETT
    elementList(2) = BUFR_NEUU
    elementList(3) = BUFR_NEVV
    elementList(4) = BUFR_NEES
    elementList(5) = BUFR_NEUS
    elementList(6) = BUFR_NEVS
    elementList(7) = BUFR_NBT1
    elementList(8) = BUFR_NBT2
    elementList(9) = BUFR_NBT3
    elementList(10)= BUFR_NERF
    elementList(11)= BUFR_NEPS

    diagn_save(:)=.false.
    diagn_all(:)=.true. 
    diagn_pressmin(:)=10.  !  0.1 hPa
    diagn_nset(:)=2
    diagn_num(:)=0
    diagn_stnid(:,:)='*********'
    diagn_varno(:,:)=0
    diagn_unilev(:,:)=.false.
    obsspace_diagn_filename(:)='obsspace_diag_'
    do j=1,numFamily
       obsspace_diagn_filename(j) ='obsspace_diag_'//familyList(j)//'_'
    end do

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namosd,iostat=ierr)
    if(ierr.ne.0) then
      write(*,*) 'osd_setup: No valid namelist NAMOSD found, skipping some diagnostics'
      nmlExists = .false.
      ierr = fclos(nulnam)
      return
    else
      nmlExists = .true.
    endif
    if(mmpi_myid.eq.0) write(*,nml=namosd)
    ierr = fclos(nulnam)

    do j=1,numFamily
       if ( .not. ofl_isFamilyTypeInList(familyList(j)) ) &
         call utl_abort('osd_setup: Specified family '//familyList(j)//' was not recognized')
       obsspace_diagn_filename(j) ='obsspace_diag_'//familyList(j)//'_'
       if (diagn_num(j).gt.max_cfg_size) call utl_abort('osd_setup: Number exceeds allowed size of max_cfg_size')
    end do
    
  end subroutine osd_setup

  !--------------------------------------------------------------------------
  ! osd_update_obsfile
  !--------------------------------------------------------------------------
  subroutine osd_update_obsfile( obsSpaceData )
    !
    ! :Purpose: Update of obs file(s) for content other
    !           than OBS,OMA,OMP,OER,FGE,MRK in obsSpaceData
    !           Content can be augmented as needed.
    !

    implicit none
    
    !Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData

    ! If needed, add effective temperature values in CH family obs file 
    ! for total column measurements

    if (obs_famExist(obsSpaceData,'CH')) call oopc_addEfftempObsfile()

  end subroutine osd_update_obsfile

  !-----------------------------------------------------------------------------------------
  !------------------- Observation-space diagnostic functions and routines -----------------

  !--------------------------------------------------------------------------
  ! osd_obsPostProc
  !--------------------------------------------------------------------------
  subroutine osd_obsPostProc( obsSpaceData, columnTrlOnAnlIncLev, deltaLat, deltaLon, deltaPressure, anlm_mode )
    !
    ! :Purpose: Interface for observation-space post-processing procedures.
    !
    ! :Arguments:
    !       :obsSpaceData:    Obs space data structure
    !       :obsfam:          Target obs family (e.g. CH)
    !       :codtyplist:      Code type list asscoiated to obsfam.
    !       :columnTrlOnAnlIncLev: Columns from analysis vertical coordinate in obs space (at obs location)
    !       :date:            YYYYMMDDHH
    !       :deltaLat:        Size of latitude bins for diagnostics (degrees)
    !       :deltaLon:        Size of longitude bins for diagnostics (degrees)
    !       :deltaPressure:   Size of vertical bins for diagnostics (Pascal)
    !       :anlm_mode:       Logical indicating if OmA (and Jo) diagnostics to be generated.
    ! 

    implicit none

    !Arguments:
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columnTrlOnAnlIncLev
    real(8), intent(in)     :: deltaLat
    real(8), intent(in)     :: deltaLon
    real(8), intent(in)     :: deltaPressure
    logical, intent(in)     :: anlm_mode
    
    !Locals:
    integer, allocatable :: codtyplist(:)
    integer :: jelm,ifam
    
    if (mmpi_myid == 0) then
       write(*,*)
       write(*,*) "Enter osd_obsPostProc: Observation-space post-processing tasks for chemical constituents"
       write(*,*)
    end if

    ! Generate and output cost function, OmP, and OmA diagnostics 

    if (obs_famExist(obsSpaceData,'CH')) then
  
       ifam=0
       do jelm=1,numFamily
          if (familyList(jelm) == 'CH') then
             ifam=jelm
             exit
          end if
       end do
       if (ifam == 0) then
       
          write(*,*) 'osd_obsPostProc: Warning - No post-processing requested for CH family.'
                     
       else
         
          ! Initialize oss_comboIDlist and add (stnid,varno) pairs from the namelist
          ! Sets list of identifiers for observations to be processed in osd_obsDiagnostics within the CH family
               
          call oss_comboIdlist(initialize_opt=.true., nset_opt=diagn_nset(ifam), all_combos_opt=diagn_all(ifam))
          do jelm=1,diagn_num(ifam)
             call oss_comboIdlist(stnid_add_opt=diagn_stnid(ifam,jelm), varno_add_opt=diagn_varno(ifam,jelm), unilev_add_opt=diagn_unilev(ifam,jelm))
          end do

          ! Diagnostics for retrievd chemical constituents (CH family)
    
          allocate(codtyplist(2))

          codtyplist(1)=codtyp_get_codtyp('CHEMREMOTE')
          codtyplist(2)=codtyp_get_codtyp('CHEMINSITU')
                       
          call osd_obsDiagnostics(obsSpaceData,columnTrlOnAnlIncLev,'CH',codtyplist,trim(obsspace_diagn_filename(ifam)), &
                      diagn_save(ifam),deltaLat,deltaLon,deltaPressure,diagn_pressmin(ifam),anlm_mode)
       
          deallocate(codtyplist)
          
       end if
        
    end if
      
    ! Generate and output cost function, OmP, and OmA diagnostics for 
    ! channels/instruments of the TO family (when processed with accompanying
    ! CH obs). 
    
    ! call osd_TO_obsDiagnostics(obsSpaceData,columnTrlOnAnlIncLev,date,deltaLat,deltaLon,deltaPressure,anlm_mode)
    
    ! Apply any required obs file update
    
    call osd_update_obsfile(obsSpaceData)

    if (mmpi_myid == 0) then
       write(*,*)
       write(*,*) "Exit osd_obsPostProc"
       write(*,*)
    end if

  end subroutine osd_obsPostProc

  !--------------------------------------------------------------------------
  ! osd_obsDiagnostics
  !--------------------------------------------------------------------------
  subroutine osd_obsDiagnostics( obsSpaceData, columnTrlOnAnlIncLev, obsfam, codtyplist, filename, save_diagn, &
                                 deltaLat, deltaLon, deltaPressure, pressmin, anlm_mode )
    !       
    ! :Purpose: Calculates and prints observation-space diagnostics for chemical constituents
    !
    ! :Arguments:
    !       :obsSpaceData:    Obs space data structure
    !       :columnTrlOnAnlIncLev: Columns from analysis vertical coordinate in obs space (at obs location)
    !       :obsfam:          Obs family (e.g. 'CH'
    !       :codtypelist:     Code type list 
    !       :filename:        Output file name
    !       :save_diagn:      Logical indicating gridded diagnostics are to be save
    !       :date:            YYYYMMDDHH
    !       :deltaLat:        Size of latitude bins for diagnostics (degrees)
    !       :deltaLon:        Size of longitude bins for diagnostics (degrees)
    !       :deltaPressure:   Size of vertical bins for diagnostics (Pascal)
    !       :pressmin:        bottom of top layer for diagnostics (in Pa).
    !       :anlm_mode:       Logical indicating if OmA diagnostics are to be generated. 
    !
    ! :Output: Content of ascii file with obs space diagnostics
    !
    ! :Comments:
    !
    !   - Although Jo_analysis is already calculated in OBS_JOBS obsSpaceData and can be passed
    !     to osd_obsspace_diagn_add, it is recalculated in osd_obsspace_diagn_add since OBS_JOBS
    !     will be set to zero for diagnostic-only observations.
    !

    implicit none

    !Arguments:
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columnTrlOnAnlIncLev
    character(len=*)        :: obsfam
    character(len=*)        :: filename
    integer, intent(in)     :: codtyplist(:)
    real(8), intent(in)     :: deltaLat
    real(8), intent(in)     :: deltaLon
    real(8), intent(in)     :: deltaPressure
    real(8), intent(in)     :: pressmin
    logical, intent(in)     :: anlm_mode
    logical, intent(in)     :: save_diagn

    !Locals:
    type(struct_osd_diagn) :: obs_diagn
    integer :: headerIndex,bodyIndex,vco,nlev_obs,ilev_obs,nlev_mod,ilev_mod

    integer, parameter :: nmax=100
    integer :: varno,varno_elemID(nmax)
    integer :: elemID,num_elemID,nset,iass
    character(len=9) :: stnid_elemID(nmax)
    logical :: unilev_elemID(nmax),unilevel,diagn_only,assim_obs,status_hpht
    character(len=256) :: label
    real(8) :: lat,lon
    character(len=12) :: stnid

    real(8), allocatable :: lev(:), omp(:), oma(:), obs(:)
    real(8), allocatable :: pres_mod(:), sigma_obs(:), sqrtHPHT(:)
    logical, allocatable :: success(:)
    integer, allocatable :: status(:)
    real(8), pointer :: height_mod(:)
    
    ! Get combination lists to group diagnostics by
    call oss_get_comboIdlist(obsSpaceData,stnid_elemID,varno_elemID,unilev_elemID,num_elemID,nset)
    
    if (num_elemID == 0) return

    if (mmpi_myid == 0) then
       write(*,*)
       write(*,*) "osd_obsDiagnostics: Observation-space diagnostics for chemical constituents"
       write(*,*)
    end if

    ! Read forecast error std. dev. at obs locations from obs files (saved in OBS_HPHT in obsSpaceData)
    !status_hpht=.false.
    if (anlm_mode) call osd_ReadSqrtHPHT(obsSpaceData,obsfam,codtyplist,status_hpht)
     
    ! Allocate memory for diagnostic arrays in obs_diagn
    call osd_obsspace_diagn_alloc(obs_diagn,deltaLat,deltaLon,deltaPressure,pressmin)

    ! Loop over all pairs in *_elemID lists
   
    do elemID=1,num_elemID

       ! Initialize the diagnostic arrays
       call osd_obsspace_diagn_init(obs_diagn)
       
       call obs_set_current_header_list(obsSpaceData,obsfam)
       HEADER: do

          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER
  
          ! Body info that we only need for first point in the profile
          bodyIndex = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)     
          vco = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)

          if (vco /= 1 .and. vco /= 2 .and. vco /= 4 .and. vco /= 5) then
             ! Vertical coordinate not handled
             write(*,*) 'osd_obsDiagnostics: Currently unaccounted VCO = ',vco
             cycle HEADER
          end if

          varno = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          nlev_obs = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)
          stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)

          ! Identify max number of profile points in the profile (exclude BUFR_SCALE_EXPONENT elements)
          call obs_set_current_body_list(obsSpaceData,headerIndex)
          do
             bodyIndex = obs_getBodyIndex(obsSpaceData)
             if (bodyIndex < 0) exit
             if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) then
                nlev_obs = nlev_obs-1
             else
                varno=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
             end if
          end do

          ! Determine if this observation should be added to this group (as specified by nset)
          if (.not.utl_stnid_equal(stnid_elemID(elemID),stnid)) cycle HEADER
          if (nset >= 2.and.varno /= varno_elemID(elemID)) cycle HEADER
          if (nset >= 3.and..not.(( nlev_obs == 1 .and. vco == 4 ).eqv.unilev_elemID(elemID))) cycle HEADER

          ! Accumulate for this combo
          
          lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)*MPC_DEGREES_PER_RADIAN_R8
          lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)*MPC_DEGREES_PER_RADIAN_R8
              
          allocate(lev(nlev_obs), omp(nlev_obs), oma(nlev_obs), obs(nlev_obs))
          allocate(sigma_obs(nlev_obs), success(nlev_obs), status(nlev_obs))          
          if (anlm_mode) allocate(sqrtHPHT(nlev_obs))

          lev(:) = 0.0d0
          omp(:) = 0.0d0
          oma(:) = 0.0d0
          obs(:) = 0.0d0
          sigma_obs(:) = 0.0d0
          if (anlm_mode) sqrtHPHT(:) = 0.0d0
          status(:) = 0
          assim_obs = .false.
          ilev_obs = 0   
          
          call obs_set_current_body_list(obsSpaceData,headerIndex)
          BODY: do
             
             bodyIndex = obs_getBodyIndex(obsSpaceData)
             if (bodyIndex < 0) exit BODY
             
             if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= varno) cycle BODY
             
             ilev_obs = ilev_obs+1

             iass = obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex)

             ! Indicates if diagnostics are to be calculated but observation not assimilated
             diagn_only = oopc_diagnOnly(obsfam,stnid,varno,nlev_obs,obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex))

             assim_obs = ((.not.diagn_only).and.anlm_mode) .or. assim_obs

             if (iass == obs_assimilated) then
                status(ilev_obs) = 1
             else if (diagn_only) then
                status(ilev_obs) = 2
             end if
             
             lev(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

             ! Include in the sum assimilated data and diagnostic only data
             if (status(ilev_obs) > 0) then
                
                omp(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
                obs(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
                sigma_obs(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
                if (anlm_mode) then
                   oma(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_OMA,bodyIndex)
                   if (status_hpht) then
                      sqrtHPHT(ilev_obs) = obs_bodyElem_r(obsSpaceData,OBS_HPHT,bodyIndex)
                      if (sqrtHPHT(ilev_obs) < 1.D-30) then
                         write(*,*) 'osd_obsDiagnostics: WARNING. sqrtHPHT not found for all obs'
                         write(*,*) 'Will not be used in Desroziers-based diagnostics.'       
                         status_hpht=.false.
                      end if
                   end if                      
                end if                   
                   
             end if

          end do BODY
             
          ! Convert to pressure if needed and identify unilevel observations
          unilevel = .false.
          select case(vco)
          case(1)
             ! Height coordinate
             
             nlev_mod = col_getNumLev(columnTrlOnAnlIncLev,'TH')  ! number of model levels     
             height_mod => col_getColumn(columnTrlOnAnlIncLev,headerIndex,'Z_T') ! geopotential
                
             allocate(pres_mod(nlev_mod))
               
             do ilev_mod=1,nlev_mod
                pres_mod(ilev_mod) = col_getPressure(columnTrlOnAnlIncLev,ilev_mod,headerIndex,'TH') ! model pressure
             end do

             ! Convert altidudes to pressure
             success = status.gt.0
             lev = phf_convert_z_to_pressure(lev,height_mod,pres_mod,nlev_obs,nlev_mod,lat/MPC_DEGREES_PER_RADIAN_R8,success)

             deallocate(pres_mod)
                
          case(4,5)
             ! Uni-level observations
             unilevel = .true.
          end select
             
          ! Add observation to diagnostic arrays
          if (anlm_mode) then 
             if (status_hpht) then
                call osd_obsspace_diagn_add(obs_diagn, lat, lon, lev, &
                          pressmin, omp, obs, sigma_obs, &
                          nlev_obs, unilevel, assim_obs, status, &
                          oma_opt=oma, sqrtHPHT_opt=sqrtHPHT)
             else 
                call osd_obsspace_diagn_add(obs_diagn, lat, lon, lev, &
                          pressmin, omp, obs, sigma_obs, &
                          nlev_obs, unilevel, assim_obs, status,oma_opt=oma)
              end if
          else
             call osd_obsspace_diagn_add(obs_diagn, lat, lon, lev, &
                          pressmin, omp, obs, sigma_obs, &
                          nlev_obs, unilevel, assim_obs, status)
          end if
          
          deallocate(lev,omp,oma,obs,sigma_obs,success,status)
          if (anlm_mode) deallocate(sqrtHPHT)
       
       end do HEADER
       
       ! Prepare output identification
       select case(nset)
       case(1)
          write(label,*) ">>> Statistics for STNID ",stnid_elemID(elemID)
       case(2)
          write(label,*) ">>> Statistics for STNID ",stnid_elemID(elemID)," and BUFR # ",varno_elemID(elemID)
       case(3)
          write(label,*) ">>> Statistics for STNID ",stnid_elemID(elemID)," and BUFR # ",varno_elemID(elemID)," and UNILEV = ",unilev_elemID(elemID)
       end select

       ! Sum over different processors
       call osd_obsspace_diagn_MPIreduce(obs_diagn)
       
       ! Output, and deallocate diagnostic arrays
       if (mmpi_myid == 0) call osd_obsspace_diagn_print(obs_diagn,filename, save_diagn, &
                              'stats', pressmin, status_hpht, label_opt=label, openfile_opt=.true.)
 
    end do
    
    ! Output diagnostics summary (over all CH observations)
    if (mmpi_myid == 0) call osd_obsspace_diagn_print(obs_diagn,filename, save_diagn, &
                            'summary', pressmin, status_hpht, openfile_opt=.true.)
 
    ! Deallocate arrays in obs_diagn
    call osd_obsspace_diagn_dealloc(obs_diagn)
   
    if (mmpi_myid == 0) then
       write(*,*)
       write(*,*) "End osd_obsDiagnostics"
       write(*,*)
    end if

  end subroutine osd_obsDiagnostics

  !--------------------------------------------------------------------------
  ! osd_ReadSqrtHPHT
  !--------------------------------------------------------------------------
  subroutine osd_ReadSqrtHPHT( obsSpaceData, obsfam, codtyplist, status_hpht )
    !
    ! :Purpose: Read background error std. dev. at obs locations from the obs files and store
    !           under OBS_HPHT in obsSpaceData
    !
    ! :Arguments:
    !   :obsSpaceData:     Observation space data
    !   :obsfam:           Obs family. e.g. 'CH'
    !   :codtyplist:       Code type list associated to obsfam
    !   :status_hpht:      logical indicating if successfully retrieved sqrtHPHT from obs file
    !

    implicit none

    !Arguments:
    logical, intent(out) :: status_hpht
    type(struct_obs), intent(inout) :: obsSpaceData
    character(len=*), intent(in) :: obsfam
    integer :: codtyplist(:)

    !Locals:
    integer :: bodyIndex,headerIndex,rln,nlv,kk
    integer :: stat,varno,icodtyp
    integer, parameter :: max_nlev=500
    integer, parameter :: ndim=1
    real(8) :: array(max_nlev)
    character(len=12), parameter :: stnid='************'
    
    type(struct_oss_obsdata) :: SqrtHPHT_struct
    
    write(*,*) 'osd_ReadSqrtHPHT: STARTING'

    ! Retrieve data from FGE blocks, i.e. sqrt(HPHT) 
    ! (with ndim=1, bkstp=15 and block_type='DATA')

    status_hpht=.true.
    SqrtHPHT_struct = obsf_obsSub_read(obsfam,stnid,-1,max_nlev,ndim,bkstp_opt=15,block_opt='DATA', &
                                    match_nlev_opt=.false.,codtyp_opt=codtyplist)

    if (SqrtHPHT_struct%nrep == 0) then
       write(*,*) 'osd_ReadSqrtHPHT: WARNING. sqrtHPHT not found in obs file(s).'
       write(*,*) 'Will not be used in Desroziers-based diagnostics.'       
       status_hpht=.false.
       call oss_obsdata_dealloc(SqrtHPHT_struct)
       return
    end if
    
    ! Save in OBS_HPHT of obsSpaceData
    
    call obs_set_current_header_list(obsSpaceData,obsfam)
    HEADER: do

       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER
  
       icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
       if (all(icodtyp /= codtyplist)) cycle HEADER
  
       ! Search for corresponding HPHT profile/element
       
       array(:)=0.0D0
       array=oss_obsdata_get_data1d(SqrtHPHT_struct, &
             obs_headElem_r(obsSpaceData,OBS_LON,headerIndex),obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex), &
             obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex), &
             obs_elem_c(obsSpaceData,'STID',headerIndex),stat_opt=stat) 
 
       if (stat == 0) then

          ! Store OBS_HPHT profile

          rln=obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
          nlv=obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)
          kk=0
          do bodyIndex = rln, nlv + rln -1
             varno=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
             if ( varno > 10000 .and. varno < 16000 ) then 
                kk=kk+1
                call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,array(kk))
             else
                call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)          
             end if
          end do
       end if
       
    enddo HEADER

    call oss_obsdata_dealloc(SqrtHPHT_struct)
    
    write(*,*) 'osd_ReadSqrtHPHT: DONE'
    
  end subroutine osd_ReadSqrtHPHT

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_alloc
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_alloc( obs_diagn, deltaLat, deltaLon, deltaPressure, pressmin )
    !
    ! :Purpose: Allocates diagnostic arrays in obs_diagn.
    !          
    ! :Arguments:
    !     :deltaLat:       latitude bin size in degrees
    !     :deltaLon:       longitutde in degrees
    !     :deltaPressure:  pressures bin size in Pa (approximate)
    !     :pressmin:       bottom of top layer for diagnostics (in Pa).
    !

    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout) :: obs_diagn
    real(8)               , intent(in)    :: deltaLat
    real(8)               , intent(in)    :: deltaLon
    real(8)               , intent(in)    :: deltaPressure
    real(8)               , intent(in)    :: pressmin
    
    !Locals:
    integer :: nlev,nlat,nlon,nbin,nstat

    obs_diagn%deltaLat = deltaLat
    obs_diagn%deltaLon = deltaLon
    obs_diagn%deltaLogPressure = deltaPressure/1.0d5  ! set constant delta ln(P) bin

    nlat = floor(180.0d0/deltaLat)
    nlon = floor(360.0d0/deltaLon)

    ! Add a last unequal size bin if remainder is larger than one degree
    if (180.0d0-nlat*deltaLat > 1.) nlat = nlat+1
    if (360.0d0-nlon*deltaLon > 1.) nlon = nlon+1

    ! Set number of levels for a pressure coordinate in hPa to cover the range
    ! of 0.01*pressmin (in hPa) to 1000 hPa for layers 2 to nlev-1. 
    ! Layer 1 covers all pressure levels <= 0.01*pressmin hPa = top of layer 2,
    ! The bottom of layer nlev-1 is provided by press_bins(nlev).
    nlev = 2 + nint(log(1.0d5/pressmin)/obs_diagn%deltaLogPressure) ! last index is for unilevel observations
    nbin = nlev*nlat*nlon

    ! Two different statistics held, stat=1 for RMS and stat=2 for mean
    nstat = 2

    obs_diagn%nlat = nlat
    obs_diagn%nlon = nlon
    obs_diagn%nlev = nlev
    obs_diagn%nbin = nbin
    obs_diagn%nstat = nstat

    if (.not. allocated(obs_diagn%OmP_stats)) allocate(obs_diagn%OmP_stats(nlat,nlon,nlev,nstat))
    if (.not. allocated(obs_diagn%OmA_stats)) allocate(obs_diagn%OmA_stats(nlat,nlon,nlev,nstat))
    if (.not. allocated(obs_diagn%obs_stats)) allocate(obs_diagn%obs_stats(nlat,nlon,nlev,nstat))
    if (.not. allocated(obs_diagn%Jo_stats))  allocate(obs_diagn%Jo_stats(nlat,nlon,nlev,nstat*2+1))
    if (.not. allocated(obs_diagn%Jpa_stats))  allocate(obs_diagn%Jpa_stats(nlat,nlon,nlev))
    if (.not. allocated(obs_diagn%counts))    allocate(obs_diagn%counts(nlat,nlon,nlev))
    if (.not. allocated(obs_diagn%nstatus))   allocate(obs_diagn%nstatus(nlat,nlon,nlev,0:2))
    if (.not. allocated(obs_diagn%diagR_stats))    allocate(obs_diagn%diagR_stats(nlat,nlon,nlev,3))
    if (.not. allocated(obs_diagn%diagHPHT_stats)) allocate(obs_diagn%diagHPHT_stats(nlat,nlon,nlev,3))

  end subroutine osd_obsspace_diagn_alloc

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_init
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_init(obs_diagn)
    !
    ! :Purpose: Initializes diagnostic arrays in obs_diagn.
    !
    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout) :: obs_diagn

    obs_diagn%OmP_stats(:,:,:,:) = 0.0d0
    obs_diagn%OmA_stats(:,:,:,:) = 0.0d0
    obs_diagn%obs_stats(:,:,:,:) = 0.0d0
    obs_diagn%Jo_stats(:,:,:,:)  = 0.0d0
    obs_diagn%Jpa_stats(:,:,:)  = 0.0d0
    obs_diagn%counts(:,:,:) = 0
    obs_diagn%nstatus(:,:,:,:) = 0
    obs_diagn%diagR_stats(:,:,:,:)    = 0.0d0
    obs_diagn%diagHPHT_stats(:,:,:,:) = 0.0d0

    obs_diagn%allow_print_summary = .false.
    obs_diagn%assim_mode = .false.

  end subroutine osd_obsspace_diagn_init

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_dealloc
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_dealloc(obs_diagn)
    !          
    ! :Purpose: Deallocates diagnostic arrays in obs_diagn.
    !
    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout) :: obs_diagn

    deallocate(obs_diagn%OmP_stats,obs_diagn%OmA_stats,obs_diagn%obs_stats)
    deallocate(obs_diagn%Jo_stats,obs_diagn%Jpa_stats,obs_diagn%counts,obs_diagn%nstatus)
    deallocate(obs_diagn%diagR_stats,obs_diagn%diagHPHT_stats)

  end subroutine osd_obsspace_diagn_dealloc

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_add
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_add( obs_diagn, lat, lon, pressure, pressmin, OmP, obs, sigma_obs, nlev_obs, &
                                     unilevel, assim_obs, status, OmA_opt, sqrtHPHT_opt)
    !        
    ! :Purpose: Adds an observation to the diagnostic arrays in obs_diagn.
    !
    ! :Arguments:
    !     :lat:            latitude in degrees
    !     :lon:            longitutde in degrees
    !     :pressure:       pressures of the profile (Pa)
    !     :pressmin:       bottom of top layer for diagnostics (in Pa).
    !     :OmP:            obs - background
    !     :OmA_opt:        obs - analysis
    !     :obs:            observations
    !     :Jo:             cost function
    !     :sigma_obs:      observation error standard deviation
    !     :sqrtHPHT_opt:   forecast error standard deviation in obs space
    !     :assim_obs:      indicates if the profile belongs to an assimilated data set
    !     :status:         indicates status of the observations, with values denoting:
    !
    !                      - 0 - observation has been rejected and not included in diagnostics
    !                      - 1 - observation has been assimilated
    !                      - 2 - observation has been used for diagnostics only (not assimilated)
    !                      only observations with status=1,2 will be added to the statistic arrays
    !     :nlev_obs:       number of observations in the profile
    !     :unilevel:       if the observation does not have a defined height coordinate
    ! 

    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout)        :: obs_diagn
    real(8)               , intent(in)           :: lat
    real(8)               , intent(in)           :: lon
    integer               , intent(in)           :: nlev_obs
    integer               , intent(in)           :: status(nlev_obs)
    real(8)               , intent(in)           :: pressure(nlev_obs)
    real(8)               , intent(in)           :: OmP(nlev_obs)
    real(8)               , intent(in)           :: obs(nlev_obs)
    real(8)               , intent(in)           :: sigma_obs(nlev_obs)
    real(8)               , intent(in)           :: pressmin
    logical               , intent(in)           :: unilevel
    logical               , intent(in)           :: assim_obs
    real(8)               , intent(in), optional :: OmA_opt(nlev_obs)
    real(8)               , intent(in), optional :: sqrtHPHT_opt(nlev_obs)

    !Locals:
    integer :: ilat,ilon,ilev,ilev_obs

    if (assim_obs) obs_diagn%assim_mode = .true.

    ! Put in first/list bin if lat,lon,level lower/higher than diagnostic range
    
    ilat = max(min(1 + floor((90.0d0+lat)/obs_diagn%deltaLat), obs_diagn%nlat), 1)
    ilon = max(min(1 + floor(lon/obs_diagn%deltaLon), obs_diagn%nlon), 1)

    LEVELS: do ilev_obs=1,nlev_obs
          
       if (unilevel) then
          ilev = obs_diagn%nlev  ! put unilevel data in last level index
       else
          ilev = max(min(2 + floor(log(pressure(ilev_obs)/pressmin)/obs_diagn%deltaLogPressure), obs_diagn%nlev-1), 1)
       end if
       
       obs_diagn%nstatus(ilat,ilon,ilev,status(ilev_obs)) = obs_diagn%nstatus(ilat,ilon,ilev,status(ilev_obs)) + 1

       if (status(ilev_obs) == 0) cycle LEVELS  ! skip adding of stats if the observation was rejected
   
       obs_diagn%counts(ilat,ilon,ilev) = obs_diagn%counts(ilat,ilon,ilev) + 1

       obs_diagn%OmP_stats(ilat,ilon,ilev,1) = obs_diagn%OmP_stats(ilat,ilon,ilev,1) + OmP(ilev_obs)**2
       obs_diagn%OmP_stats(ilat,ilon,ilev,2) = obs_diagn%OmP_stats(ilat,ilon,ilev,2) + OmP(ilev_obs)

       obs_diagn%obs_stats(ilat,ilon,ilev,1) = obs_diagn%obs_stats(ilat,ilon,ilev,1) + obs(ilev_obs)**2
       obs_diagn%obs_stats(ilat,ilon,ilev,2) = obs_diagn%obs_stats(ilat,ilon,ilev,2) + obs(ilev_obs)

       obs_diagn%Jo_stats(ilat,ilon,ilev,2)  = obs_diagn%Jo_stats(ilat,ilon,ilev,2)  + 0.5 * OmP(ilev_obs)**2 / sigma_obs(ilev_obs)**2
       if (status(ilev_obs) == 1) obs_diagn%Jo_stats(ilat,ilon,ilev,4)  = obs_diagn%Jo_stats(ilat,ilon,ilev,4)  + 0.5 * OmP(ilev_obs)**2 / sigma_obs(ilev_obs)**2

       if (present(OmA_opt)) then
          
          obs_diagn%OmA_stats(ilat,ilon,ilev,1) = obs_diagn%OmA_stats(ilat,ilon,ilev,1) + OmA_opt(ilev_obs)**2
          obs_diagn%OmA_stats(ilat,ilon,ilev,2) = obs_diagn%OmA_stats(ilat,ilon,ilev,2) + OmA_opt(ilev_obs)
          obs_diagn%Jo_stats(ilat,ilon,ilev,1)  = obs_diagn%Jo_stats(ilat,ilon,ilev,1)  + 0.5 * OmA_opt(ilev_obs)**2 / sigma_obs(ilev_obs)**2

          if (status(ilev_obs) == 1) then

             obs_diagn%diagR_stats(ilat,ilon,ilev,1) = obs_diagn%diagR_stats(ilat,ilon,ilev,1) &
                + OmP(ilev_obs)*OmA_opt(ilev_obs)/sigma_obs(ilev_obs)**2
             obs_diagn%diagR_stats(ilat,ilon,ilev,2) = obs_diagn%diagR_stats(ilat,ilon,ilev,2) &
                + OmP(ilev_obs)/sigma_obs(ilev_obs)
             obs_diagn%diagR_stats(ilat,ilon,ilev,3) = obs_diagn%diagR_stats(ilat,ilon,ilev,3) &
                + OmA_opt(ilev_obs)/sigma_obs(ilev_obs)

             if (present(sqrtHPHT_opt)) then
                if (sqrtHPHT_opt(ilev_obs) > 0.0) then
                   obs_diagn%diagHPHT_stats(ilat,ilon,ilev,1) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,1) &
                        + OmP(ilev_obs)*(OmP(ilev_obs)-OmA_opt(ilev_obs))/sqrtHPHT_opt(ilev_obs)**2
                   obs_diagn%diagHPHT_stats(ilat,ilon,ilev,2) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,2) &
                        + OmP(ilev_obs)/sqrtHPHT_opt(ilev_obs)
                   obs_diagn%diagHPHT_stats(ilat,ilon,ilev,3) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,3) &
                        + (OmP(ilev_obs)-OmA_opt(ilev_obs))/sqrtHPHT_opt(ilev_obs)
                   ! Following two diagnostics do not account for spatial correlations (spatial correlations of HPHT)! 
                   ! As a consequence, Jb likely overestimated.
                   obs_diagn%Jo_stats(ilat,ilon,ilev,5)  = obs_diagn%Jo_stats(ilat,ilon,ilev,5)  + 0.5 * (OmP(ilev_obs))**2 &
                                                  / (sigma_obs(ilev_obs)**2 + sqrtHPHT_opt(ilev_obs)**2) 
                   obs_diagn%Jpa_stats(ilat,ilon,ilev)  = obs_diagn%Jpa_stats(ilat,ilon,ilev)  + 0.5 * (OmA_opt(ilev_obs)-OmP(ilev_obs))**2 &
                                                  / sqrtHPHT_opt(ilev_obs)**2
                end if
             else

                ! No division by sqrtHPHT
                
                obs_diagn%diagHPHT_stats(ilat,ilon,ilev,1) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,1) + OmP(ilev_obs)*(OmP(ilev_obs)-OmA_opt(ilev_obs))
                obs_diagn%diagHPHT_stats(ilat,ilon,ilev,2) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,2) + OmP(ilev_obs)
                obs_diagn%diagHPHT_stats(ilat,ilon,ilev,3) = obs_diagn%diagHPHT_stats(ilat,ilon,ilev,3) + OmP(ilev_obs)-OmA_opt(ilev_obs)
             end if
             
          end if

       end if

    end do LEVELS

  end subroutine osd_obsspace_diagn_add

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_MPIreduce
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_MPIreduce(obs_diagn)
    !         
    ! :Purpose: Performs a MPI allreduce on diagnostic arrays in obs_diagn.
    !

    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout) :: obs_diagn

    !Locals:
    
    ! MPI global arrays
    real(8), allocatable :: OmP_global(:,:,:,:), OmA_global(:,:,:,:), obs_global(:,:,:,:), Jo_global(:,:,:,:)
    real(8), allocatable :: Jpa_global(:,:,:), diagR_global(:,:,:,:), diagHPHT_global(:,:,:,:)
    integer, allocatable :: counts_global(:,:,:),nstatus_global(:,:,:,:)
    logical :: assim_global

    integer :: nlat,nlon,nlev,nstat,nbin,ierr

    nlat = obs_diagn%nlat
    nlon = obs_diagn%nlon
    nlev = obs_diagn%nlev
    nstat = obs_diagn%nstat
    nbin = obs_diagn%nbin

    ! Allocate memory for mpi global arrays
    allocate(OmP_global(nlat,nlon,nlev,nstat))
    allocate(OmA_global(nlat,nlon,nlev,nstat))
    allocate(obs_global(nlat,nlon,nlev,nstat))
    allocate(Jo_global(nlat,nlon,nlev,nstat*2+1))
    allocate(Jpa_global(nlat,nlon,nlev))
    allocate(diagR_global(nlat,nlon,nlev,3))
    allocate(diagHPHT_global(nlat,nlon,nlev,3))
    allocate(counts_global(nlat,nlon,nlev))
    allocate(nstatus_global(nlat,nlon,nlev,0:2))
    
    ! Reduce from all mpi processes
    call rpn_comm_allreduce(obs_diagn%OmP_stats,OmP_global,nbin*nstat,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%OmA_stats,OmA_global,nbin*nstat,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%obs_stats,obs_global,nbin*nstat,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%Jo_stats,Jo_global,nbin*(nstat*2+1),"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%Jpa_stats,Jpa_global,nbin,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%diagR_stats,diagR_global,nbin*3,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%diagHPHT_stats,diagHPHT_global,nbin*3,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%counts,counts_global,nbin,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%nstatus,nstatus_global,nbin*3,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(obs_diagn%assim_mode,assim_global,1,"MPI_LOGICAL","MPI_LOR","GRID",ierr)

    ! save in struct_osd_diagn
    obs_diagn%OmP_stats = OmP_global
    obs_diagn%OmA_stats = OmA_global
    obs_diagn%obs_stats = obs_global
    obs_diagn%Jo_stats = Jo_global
    obs_diagn%Jpa_stats = Jpa_global
    obs_diagn%diagR_stats = diagR_global
    obs_diagn%diagHPHT_stats = diagHPHT_global
    obs_diagn%counts = counts_global
    obs_diagn%nstatus = nstatus_global
    obs_diagn%assim_mode = assim_global

    deallocate(OmP_global,OmA_global,obs_global,Jo_global,Jpa_global,diagR_global,diagHPHT_global, &
               counts_global,nstatus_global)
    
  end subroutine osd_obsspace_diagn_MPIreduce

  !--------------------------------------------------------------------------
  ! osd_obsspace_diagn_print
  !--------------------------------------------------------------------------
  subroutine osd_obsspace_diagn_print(obs_diagn, filename, save_diagn, print_type, pressmin, status_hpht, label_opt, openfile_opt )
    !        
    ! :Purpose: Prints observation space diagnostics. If called with print_type = 'stats', the
    !           printed statistics will be added to the total diagnostic arrays.
    !        
    ! :Arguments:
    !     :filename:       output file name
    !     :save_diagn:     Logical indicating gridded diagnostics are to be save
    !     :print_type:     Specifies which statistics to print, with possible values:
    !
    !                      - 'stats'  - prints statistics for the the arrays within obs_diagn
    !                      - 'summary'- prints total statistics held in the saved variables
    !                                 within this subrouine
    !     :pressmin:       min pressure level for output
    !     :label_opt:      label to print (only relevant if print_type = 'stats')
    !     :openfile_opt:   logical indicating if file filename is to be opened.
    !     :status_hpht:    logical indicating if sqrtHPHT were available.
    !
    
    implicit none

    !Arguments:
    type(struct_osd_diagn), intent(inout)        :: obs_diagn
    character(len=*)                             :: print_type
    character(len=*)                             :: filename
    real(8)                                      :: pressmin
    logical               , intent(in)           :: status_hpht
    logical               , intent(in)           :: save_diagn
    logical               , intent(in), optional :: openfile_opt
    character(len=256)    , intent(in), optional :: label_opt

    !Locals:
    integer, external :: fnom, fclos
    real(8) :: Jo_a,Jo_b,Jpa_assim, Jo_a_assim,Jo_b_assim, Jo_p_assim
    real(8), save :: Jo_a_total=0.0d0, Jo_b_total=0.0d0, Jpa_total_assim=0.0d0
    real(8), save :: Jo_a_total_assim=0.0d0, Jo_b_total_assim=0.0d0, Jo_p_total_assim=0.0d0
    integer, save :: counts_total=0, counts_total_assim=0

    integer :: ierr,unit,icount,nlat,nlon,nlev,ilev,ilat,ilon,icount_assim
    integer, allocatable :: ncounts(:), ncounts_assim(:)
    real(8), allocatable :: press_bins(:)
    logical :: fileout_exist,multilevel,unilevel
    
    select case(trim(print_type))
    case('stats','STATS')
       
       ! Print observation-space statistics to listing file or to output file obspace_diag_filename.
             
       ! Open and append to output file if requested
             
       if (present(openfile_opt)) then
          if (openfile_opt) then 
             inquire(file=filename, exist=fileout_exist)
             call utl_open_asciifile(filename,unit) 
          else
             unit=6
          end if
       else
          unit=6             
       end if
             
       if (present(label_opt)) then
          write(unit,*)
          write(unit,*) trim(label_opt)
       end if
                      
       if (any(obs_diagn%counts.gt.0)) then

          nlat = obs_diagn%nlat
          nlon = obs_diagn%nlon
          nlev = obs_diagn%nlev
          
          allocate(ncounts(nlev), ncounts_assim(nlev), press_bins(nlev+1))
          
          ! Pressure boundaries for each bin starting from a top layer with lower boundary of 0.01*pressmin (i=2) in hPa 
          ! and extending to the surface.
          
          press_bins(2:nlev-1) = (/ (0.01*pressmin*exp((ilev-2)*obs_diagn%deltaLogPressure), ilev=2,nlev-1) /)
          press_bins(1) = 0.0d0
          press_bins(nlev) = 1200.0d0 ! set to pressure larger than the largest expected surface pressure.
          press_bins(nlev+1) = 0.0d0  ! set to zero for unilevel

          ! Total counts for each level 
          ncounts = sum(sum(obs_diagn%counts,dim=1),dim=1)
          counts_total = counts_total + sum(ncounts)
          ncounts_assim = sum(sum(obs_diagn%nstatus(:,:,:,1),dim=1),dim=1)
          counts_total_assim = counts_total_assim + sum(ncounts_assim)

          ! Indicates if any multilevel or unilevel observations exist
          multilevel = any(obs_diagn%nstatus(:,:,1:nlev-1,:).gt.0)
          unilevel = any(obs_diagn%nstatus(:,:,nlev,:).gt.0)
          
          if (obs_diagn%assim_mode) then
              write(unit,*)
              write(unit,'(A)') "  Elements for calc of obs and background error standard deviation scaling factors via"
              write(unit,'(A)') "  the Desroziers approach can be found in the third and fourth blocks of statistics."
              write(unit,*)
          end if

          write(unit,*)
          write(unit,'(A)') " Global statistics"
          write(unit,'(A)') " -----------------"

          ! Multi-level data

          if (multilevel) then
            
             icount = sum(ncounts(1:nlev-1))
             icount_assim = sum(ncounts_assim(1:nlev-1))
             
             if (icount > 0) then
                Jo_a = sum(obs_diagn%Jo_stats(:,:,1:nlev-1,1))
                Jo_b = sum(obs_diagn%Jo_stats(:,:,1:nlev-1,2))
                Jo_a_assim = sum(obs_diagn%Jo_stats(:,:,1:nlev-1,3))
                Jo_b_assim = sum(obs_diagn%Jo_stats(:,:,1:nlev-1,4))
                Jo_p_assim = sum(obs_diagn%Jo_stats(:,:,1:nlev-1,5))
                Jpa_assim = sum(obs_diagn%Jpa_stats(:,:,1:nlev-1))
                Jo_a_total = Jo_a_total + Jo_a
                Jo_b_total = Jo_b_total + Jo_b
                Jo_a_total_assim = Jo_a_total_assim + Jo_a_assim
                Jo_b_total_assim = Jo_b_total_assim + Jo_b_assim
                Jo_p_total_assim = Jo_p_total_assim + Jo_p_assim
                Jpa_total_assim = Jpa_total_assim + Jpa_assim
             else
                Jo_a = 0.0d0
                Jo_a_assim = 0.0d0
                Jo_b_assim = 0.0d0
                Jo_p_assim = 0.0d0
                Jpa_assim = 0.0d0
             end if
                    
             obs_diagn%allow_print_summary = .true.

             call print_J(unit,"Multi-level data:",Jo_a,Jo_b,icount,Jo_a_assim,Jo_b_assim,Jpa_assim,Jo_p_assim,icount_assim)

             call print_stats(unit,obs_diagn,press_bins,1,nlat,1,nlon,1,nlev-1) 
             if (obs_diagn%assim_mode) call print_Desroziers(unit,obs_diagn,press_bins,1,nlat,1,nlon,1,nlev-1,status_hpht)
             
          end if
             
          ! Uni-level data
          
          if (unilevel) then
             
             if (ncounts(nlev) > 0) then
                Jo_a = sum(obs_diagn%Jo_stats(:,:,nlev,1))
                Jo_b = sum(obs_diagn%Jo_stats(:,:,nlev,2))
                Jo_a_assim = sum(obs_diagn%Jo_stats(:,:,nlev,3))
                Jo_b_assim = sum(obs_diagn%Jo_stats(:,:,nlev,4))
                Jo_p_assim = sum(obs_diagn%Jo_stats(:,:,nlev,5))
                Jpa_assim = sum(obs_diagn%Jpa_stats(:,:,nlev))
                Jo_a_total = Jo_a_total + Jo_a
                Jo_b_total = Jo_b_total + Jo_b
                Jo_a_total_assim = Jo_a_total_assim + Jo_a_assim
                Jo_b_total_assim = Jo_b_total_assim + Jo_b_assim
                Jo_p_total_assim = Jo_p_total_assim + Jo_p_assim
                Jpa_total_assim = Jpa_total_assim + Jpa_assim
             else
                Jo_a = 0.0d0
                Jo_b = 0.0d0
                Jo_a_assim = 0.0d0
                Jo_b_assim = 0.0d0
                Jo_p_assim = 0.0d0
                Jpa_assim = 0.0d0
             end if

             obs_diagn%allow_print_summary = .true.

             call print_J(unit,"Uni-level data:",Jo_a,Jo_b,ncounts(nlev),Jo_a_assim,Jo_b_assim,Jpa_assim,Jo_p_assim,ncounts_assim(nlev))
            
             call print_stats(unit,obs_diagn,press_bins,1,nlat,1,nlon,nlev,nlev) 
             if (obs_diagn%assim_mode) call print_Desroziers(unit,obs_diagn,press_bins,1,nlat,1,nlon,nlev,nlev,status_hpht)
             
          end if
                          
          deallocate(ncounts,ncounts_assim)
             
          ! Output lat,lon dependent averages to file if obsspace_diagn_filename is provided
          if (present(openfile_opt)) then
             if (openfile_opt .and. save_diagn .and. (nlat > 1 .or. nlon > 1)) then
                
                write(unit,*)
                write(unit,'(A)') " Lat-lon gridded statistics"
                write(unit,'(A)') " --------------------------"
                write(unit,*)
                write(unit,'(2X,3(A,I4))') "nlat = ",nlat," , nlon = ",nlon," , nlev = ",nlev
                write(unit,*)

                do ilat=1,nlat
                   do ilon=1,nlon
                      write(unit,'(2X,2(A,I6),3X,2(F8.1,A,F8.1))') "ilat = ",ilat," , ilon = ",ilon, &
                           (ilat-1.)*obs_diagn%deltaLat-90.," < lat < ",ilat*obs_diagn%deltaLat-90., &
                           (ilon-1.)*obs_diagn%deltaLon," < lon < ",ilon*obs_diagn%deltaLon
                      if (any(obs_diagn%nstatus(ilat,ilon,1:nlev-1,:) > 0)) then
                         write(unit,*)
                         write(unit,'(A)') " Multi-level data:"
                         write(unit,*)
                         call print_stats(unit,obs_diagn,press_bins,ilat,ilat,ilon,ilon,1,nlev-1) 
                         if (obs_diagn%assim_mode) call print_Desroziers(unit,obs_diagn,press_bins,ilat,ilat,ilon,ilon,1,nlev-1,status_hpht)
                      else if (multilevel) then
                         write(unit,*)
                         write(unit,'(A)') " No multi-level data."
                         write(unit,*)
                      end if
                      if (any(obs_diagn%nstatus(ilat,ilon,nlev,:) > 0)) then
                         write(unit,*)
                         write(unit,'(A)') " Uni-level data:"
                         write(unit,*)
                         call print_stats(unit,obs_diagn,press_bins,ilat,ilat,ilon,ilon,nlev,nlev) 
                         if (obs_diagn%assim_mode) call print_Desroziers(unit,obs_diagn,press_bins,ilat,ilat,ilon,ilon,nlev,nlev,status_hpht)
                      else if (unilevel) then
                         write(unit,*)
                         write(unit,'(A)') " No uni-level data."
                         write(unit,*)
                      end if
                   end do
                end do
                
             end if
          end if

          deallocate(press_bins)

       else
          write(unit,*)
          write(unit,*) "No data found for this combination."
          write(unit,*)
       end if
       
       if (present(openfile_opt)) then
          if (openfile_opt) ierr=fclos(unit)     
       end if
             
    case('summary','SUMMARY')

       if (.not.obs_diagn%allow_print_summary) then
          write(*,*) "osd_obsspace_diagn_print: allow_print_summary is set to false, no summary will be printed."
          return
       end if
          
       if (present(openfile_opt)) then
          if (openfile_opt) then 
             inquire(file=filename, exist=fileout_exist)
             call utl_open_asciifile(filename,unit) 
          else
             unit=6
          end if
       else
          unit=6
       end if
            
       call print_J(unit,"Total cost function for CH observations:",Jo_a_total,Jo_b_total,counts_total, &
                         Jo_a_total_assim,Jo_b_total_assim,Jpa_total_assim,Jo_a_total_assim,counts_total_assim)

       if (present(openfile_opt)) then
          if (openfile_opt) ierr=fclos(unit) 
       end if
    
    case default
       call utl_abort("osd_obsspace_diagn_print: Invalid print_type select of " // trim(print_type) )
    end select

  
  contains

    !--------------------------------------------------------------------------
    ! print_J
    !--------------------------------------------------------------------------
    subroutine print_J( unit, title, Jo_analysis, Jo_backgrnd, nobs, Jo_anal_assim, Jo_bgck_assim, Jpa_assim, Jo_p_assim, nobs_assim )
      !
      implicit none

      !Arguments:
      integer         , intent(in) :: unit
      integer         , intent(in) :: nobs
      integer         , intent(in) :: nobs_assim
      character(len=*), intent(in) :: title
      real(8)         , intent(in) :: Jo_analysis
      real(8)         , intent(in) :: Jo_backgrnd
      real(8)         , intent(in) :: Jo_anal_assim
      real(8)         , intent(in) :: Jo_bgck_assim
      real(8)         , intent(in) :: Jpa_assim
      real(8)         , intent(in) :: Jo_p_assim

      !Locals:
      real(8) :: Jo_analysis_norm,Jo_backgrnd_norm,Jt_assim
      real(8) :: Jpa_norm_assim,Jt_norm_assim,Jo_p_norm_assim
      character(len=100) :: fmt

      if (Jo_analysis < 1.0d6 .and. Jo_backgrnd < 1.0d6) then
         fmt = '(A,F24.8,A,F24.8)'
      else
         fmt = '(A,ES24.8,A,ES24.8)'
      end if

      if (nobs > 0) then
         Jo_analysis_norm = 2.*Jo_analysis/nobs
         Jo_backgrnd_norm = 2.*Jo_backgrnd/nobs
      else
         Jo_analysis_norm = 0.0d0
         Jo_backgrnd_norm = 0.0d0
      end if

      if ((nobs > 0 .and. nobs_assim > 0) .or. nobs_assim == 0) then
         write(unit,*)
         write(unit,'(A)') " " // title // " totals "
         write(unit,*)
         write(unit,fmt)       "   Jo(x_analysis) = ",Jo_analysis," ,   2*Jo(x_analysis)/N = ",Jo_analysis_norm
         write(unit,fmt)       "   Jo(x_backgrnd) = ",Jo_backgrnd," ,   2*Jo(x_backgrnd)/N = ",Jo_backgrnd_norm
         write(unit,'(A,I24)') "                N = ",nobs
         write(unit,*)
      end if

      if (nobs_assim > 0) then
         Jo_analysis_norm = 2.*Jo_anal_assim/nobs_assim
         Jo_backgrnd_norm = 2.*Jo_bgck_assim/nobs_assim
         Jpa_norm_assim = 2.*Jpa_assim/nobs_assim 
         Jt_norm_assim = Jpa_norm_assim+Jo_analysis_norm
         Jo_p_norm_assim = 2.*Jo_p_assim/nobs_assim

         write(unit,*)
         write(unit,'(A)') " " // title // " assimilated data "
         write(unit,*)
         write(unit,fmt)       "   J(A-P)+Jo(x_a) = ",Jt_assim,     " ,   2*[J(A-P)+Jo(x_a)]/N = ",Jt_norm_assim
         write(unit,fmt)       "   J(A-P)         = ",Jpa_assim,    " ,   2*J(A-P)/N           = ",Jpa_norm_assim
         write(unit,fmt)       "   Jo(x_analysis) = ",Jo_anal_assim," ,   2*Jo(x_analysis)/N   = ",Jo_analysis_norm
         write(unit,fmt)       "   Jo(x_backgrnd) = ",Jo_bgck_assim," ,   2*Jo(x_backgrnd)/N   = ",Jo_backgrnd_norm
         write(unit,fmt)       "   J(O-P)         = ",Jo_p_assim,   " ,   2*J(O-P)/N           = ",Jo_p_norm_assim
         write(unit,'(A,I24)') "                N = ",nobs_assim
         write(unit,*)
      end if
       
    end subroutine print_J

    !--------------------------------------------------------------------------
    ! print_stats
    !--------------------------------------------------------------------------
    subroutine print_stats( unit, obs_diagn, pressure, ilat_start, ilat_end, ilon_start, ilon_end, ilev_start, ilev_end )

      implicit none

      !Arguments:
      type(struct_osd_diagn), intent(in) :: obs_diagn
      real(8)               , intent(in) :: pressure(obs_diagn%nlev+1)
      integer               , intent(in) :: unit
      integer               , intent(in) :: ilat_start
      integer               , intent(in) :: ilat_end
      integer               , intent(in) :: ilon_start
      integer               , intent(in) :: ilon_end
      integer               , intent(in) :: ilev_start
      integer               , intent(in) :: ilev_end

      !Locals:
      integer :: ilev,level,counts(obs_diagn%nlev),N_assim,N_diagn,N_rej
      real(8) :: pres1,pres2,jo_a,jo_b,jo_a_norm,jo_b_norm,obs_sum,obs_mean,obs_std,OmP_mean,OmP_rms,OmA_mean,OmA_rms
      logical :: skip(obs_diagn%nlev)

      skip(:) = .false.

      write(unit,'(A)') "  Layer     Pressure (hPa)      Counts (N)    N_assim    N_diagn      N_rej    Jo(O-A)     2*Jo(O-A)/N   Jo(O-P)     2*Jo(O-P)/N"
      write(unit,'(A)') "  -----     --------------      ----------    -------    -------      -----    -------     -----------   -------     -----------"

      do ilev = ilev_start, ilev_end

         counts(ilev) = sum(obs_diagn%counts(ilat_start:ilat_end,ilon_start:ilon_end,ilev))

         N_rej   = sum(obs_diagn%nstatus(ilat_start:ilat_end,ilon_start:ilon_end,ilev,0))
         N_assim = sum(obs_diagn%nstatus(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))
         N_diagn = sum(obs_diagn%nstatus(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2))

         skip(ilev) = counts(ilev) == 0 .and. N_rej == 0

         if (skip(ilev)) cycle

         if (ilev < nlev) then
            pres1 = pressure(ilev)
            pres2 = pressure(ilev+1)
            level = ilev
         else
            pres1 = 0.0d0
            pres2 = 0.0d0
            level = 0
         end if
         
         jo_a = sum(obs_diagn%Jo_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))
         jo_b = sum(obs_diagn%Jo_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2))

         if ( counts(ilev) == 0 ) then
            jo_a_norm = 0.0d0
            jo_b_norm = 0.0d0
         else
            jo_a_norm = 2.*jo_a/counts(ilev)
            jo_b_norm = 2.*jo_b/counts(ilev)
         end if

         write(unit,'(2X, I3, 2(2X, F11.4), 4(2X,I9), 4(2X, ES11.4))') &
              level,pres1,pres2,counts(ilev),N_assim,N_diagn,N_rej,jo_a,jo_a_norm,jo_b,jo_b_norm

      end do
    
      write(unit,*)
      ! Division by <O> removed to prevent unlikely division by zero (or <= zero)
      ! write(unit,'(A)') "  Layer     Pressure (hPa)       Counts   obs (mean)   obs (std)  rms(O-P)/<O>   <O-P>/<O>   rms(O-A)/<O>  <O-A>/<O>"
      ! write(unit,'(A)') "  -----     --------------       ------   ----------   ---------  ------------   ---------   ------------  ---------"
      write(unit,'(A)') "  Layer     Pressure (hPa)       Counts   obs (mean)   obs (std)    rms(O-P)       <O-P>       rms(O-A)      <O-A>"
      write(unit,'(A)') "  -----     --------------       ------   ----------   ---------    --------       -----       --------      -----"


      do ilev=ilev_start,ilev_end

         if (skip(ilev)) cycle

         if (ilev < nlev) then
            pres1 = pressure(ilev)
            pres2 = pressure(ilev+1)
            level = ilev
         else
            pres1 = 0.0d0
            pres2 = 0.0d0
            level = 0
         end if
         
         if (counts(ilev) == 0) then
            obs_mean = 0.0d0
            obs_std = 0.0d0
            OmP_rms = 0.0d0
            OmP_mean = 0.0d0
            OmA_rms = 0.0d0
            OmA_mean = 0.0d0
    
            write(unit,'(2X, I3, 2(2X, F11.4), 2X, I6, 6(2X, ES11.4))') &
                 level,pres1,pres2,counts(ilev),obs_mean,obs_std,OmP_rms,OmP_mean,OmA_rms,OmA_mean
         else
            obs_sum = sum(obs_diagn%obs_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2))
            obs_mean = obs_sum / counts(ilev)
            obs_std = sqrt(max(0.0D0, sum(obs_diagn%obs_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))/counts(ilev) - (obs_sum/counts(ilev))**2 ))

            OmP_rms = sqrt(max(0.0D0, sum(obs_diagn%OmP_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1)) / counts(ilev) ))
            OmP_mean = sum(obs_diagn%OmP_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2)) / obs_sum
            
            OmA_rms = sqrt(max(0.0D0, sum(obs_diagn%OmA_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1)) / counts(ilev) ))
            OmA_mean = sum(obs_diagn%OmA_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2)) / counts(ilev)
    
            ! write(unit,'(2X, I3, 2(2X, F11.4), 6(2X, ES11.4))') &
            !     level,pres1,pres2,obs_mean,obs_std,OmP_rms/obs_mean,OmP_mean/obs_sum,OmA_rms/obs_mean,OmA_mean/obs_sum
    
            write(unit,'(2X, I3, 2(2X, F11.4), 2X, I6, 6(2X, ES11.4))') &
                  level,pres1,pres2,counts(ilev),obs_mean,obs_std,OmP_rms,OmP_mean,OmA_rms,OmA_mean

         end if
 
      end do

      write(unit,*)

    end subroutine print_stats

    !--------------------------------------------------------------------------
    ! print_Desroziers
    !--------------------------------------------------------------------------
    subroutine print_Desroziers( unit, obs_diagn, pressure, ilat_start, ilat_end, ilon_start, ilon_end, ilev_start, ilev_end, status_hpht )
      !
      ! :Purpose: Prints elements contributing to the calc of scaling factors for observation and background
      !           error std. dev. based on the Desroziers approach.
      !
      implicit none

      !Arguments:
      type(struct_osd_diagn), intent(in) :: obs_diagn
      real(8)               , intent(in) :: pressure(obs_diagn%nlev+1)
      integer               , intent(in) :: unit
      integer               , intent(in) :: ilat_start
      integer               , intent(in) :: ilat_end
      integer               , intent(in) :: ilon_start
      integer               , intent(in) :: ilon_end
      integer               , intent(in) :: ilev_start
      integer               , intent(in) :: ilev_end
      logical               , intent(in) :: status_hpht

      !Locals:
      integer :: ilev,level,N_assim(obs_diagn%nlev)
      real(8) :: pres1,pres2,sum_prod,sum_OmP,sum_OmA,sum_AmP,scaling

      N_assim(:) = 0

      write(unit,'(A)') "  Layer     Pressure (hPa)         N_assim   sum[(O-P)(O-A)/var(O)]  sum[(O-P)/sig(O)]  sum[(O-A)/sig(O)]  mean[(O-P)(O-A)/var(O)]-mean[(O-P)/sig(O)]*mean[(O-A)/sig(O)]"
      write(unit,'(A)') "  -----     --------------         -------   ----------------------  -----------------  -----------------  -------------------------------------------------------------"

      do ilev=ilev_start,ilev_end
         
         N_assim(ilev) = sum(obs_diagn%nstatus(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))

         if (N_assim(ilev) == 0) cycle

         if (ilev < nlev) then
            pres1 = pressure(ilev)
            pres2 = pressure(ilev+1)
            level = ilev
         else
            pres1 = 0.0d0
            pres2 = 0.0d0
            level = 0
         end if
         
         sum_prod = sum(obs_diagn%diagR_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))
         sum_OmP  = sum(obs_diagn%diagR_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2))
         sum_OmA  = sum(obs_diagn%diagR_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,3))

         scaling = sum_prod/N_assim(ilev) - sum_OmP*sum_OmA/N_assim(ilev)**2

         write(unit,'(2X, I3, 2(2X, F11.4), 2X, I9, 11X, ES11.4, 2(8X,ES11.4), 10X, ES11.4)') level,pres1,pres2,N_assim(ilev),sum_prod,sum_OmP,sum_OmA,scaling
            
      end do

      write(unit,*)
      if (status_hpht) then
         write(unit,'(A)') "  Layer     Pressure (hPa)         N_assim   sum[(O-P)(A-P)/var(P)]  sum[(O-P)/sig(P)]  sum[(A-P)/sig(P)]  mean[(O-P)(A-P)/var(P)]-mean[(O-P)/sig(P)]*mean[(A-P)/sig(P)]"
      else
         write(unit,'(A)') "  Layer     Pressure (hPa)         N_assim       sum[(O-P)(A-P)]         sum[O-P]           sum[A-P]            mean[(O-P)(A-P)]-mean[O-P]*mean[A-P]"
      end if
      write(unit,'(A)') "  -----     --------------         -------   ----------------------  -----------------  -----------------  -------------------------------------------------------------"

      do ilev=ilev_start,ilev_end
         
         if (N_assim(ilev) == 0) cycle
         
         if (ilev < nlev) then
            pres1 = pressure(ilev)
            pres2 = pressure(ilev+1)
            level = ilev
         else
            pres1 = 0.0d0
            pres2 = 0.0d0
            level = 0
         end if
         
         sum_prod = sum(obs_diagn%diagHPHT_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,1))
         sum_OmP  = sum(obs_diagn%diagHPHT_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,2))
         sum_AmP  = sum(obs_diagn%diagHPHT_stats(ilat_start:ilat_end,ilon_start:ilon_end,ilev,3))
         
         scaling = sum_prod/N_assim(ilev) - sum_OmP*sum_AmP/N_assim(ilev)**2
         
         write(unit,'(2X, I3, 2(2X, F11.4), 2X, I9, 11X, ES11.4, 2(8X,ES11.4), 10X, ES11.4)') level,pres1,pres2,N_assim(ilev),sum_prod,sum_OmP,sum_AmP,scaling
         
      end do

      write(unit,*)

    end subroutine print_Desroziers

  end subroutine osd_obsspace_diagn_print

end module obsSpaceDiag_mod
