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

MODULE localizationSpectral_mod
  ! MODULE localizationSpectral_mod (prefix='lsp' category='2. B and R matrices')
  !
  ! :Purpose: To compute localized 3D gridpoint amplitude fields for each
  !           ensemble member from a given (1D) control vector of
  !           SPECTRAL ELEMENTS
  !
  use midasMpi_mod
  use utilities_mod
  use globalSpectralTransform_mod
  use lamSpectralTransform_mod
  use localizationFunction_mod
  use horizontalCoord_mod
  use earthConstants_mod
  use ensembleStatevector_mod
  implicit none
  save
  private

  ! public derived type
  public :: struct_lsp
  ! public procedures
  public :: lsp_setup, lsp_Lsqrt, lsp_LsqrtAd, lsp_finalize
  public :: lsp_reducetompilocal, lsp_reducetompilocal_r4
  public :: lsp_expandtompiglobal, lsp_expandtompiglobal_r4

  type :: struct_lsp
     logical             :: initialized = .false.
     real(8),allocatable :: LhorizSqrt(:,:)
     real(8),allocatable :: LvertSqrt(:,:)
     type(struct_lst)    :: lst    ! Spectral transform Parameters
     real(8)             :: dlon
     integer             :: gstID
     integer, pointer    :: ilaList_mpiglobal(:)
     integer, pointer    :: ilaList_mpilocal(:)
     integer             :: cvDim_mpilocal
     integer             :: cvDim_mpiglobal
     integer             :: nla_mpilocal
     integer             :: nla_mpiglobal
     integer             :: ntrunc
     integer             :: nphase
     integer             :: ni, nj
     integer             :: myLatBeg, myLatEnd
     integer             :: myLonBeg, myLonEnd
     integer             :: nEnsOverDimension, nEns
     integer             :: nLev
     integer             :: mymBeg, mymEnd, mymSkip
     integer             :: mynBeg, mynEnd, mynSkip
     logical             :: global
  end type struct_lsp

  real(8), parameter :: rsq2 = sqrt(2.0d0)

  logical, parameter :: verbose = .false.

CONTAINS

!--------------------------------------------------------------------------
! lsp_setup
!--------------------------------------------------------------------------
  SUBROUTINE lsp_setup(hco_loc, nEns, nLev, pressureProfile, ntrunc, locType, &
                       locMode, horizLengthScale1, horizLengthScale2, vertLengthScale, &
                       cvDim_out, lsp, nEnsOverDimension_out)
    implicit none
  
    type(struct_lsp), pointer ::  lsp

    type(struct_hco), pointer :: hco_loc

    integer, intent(in) :: nEns
    integer, intent(in) :: nLev
    integer, intent(in) :: nTrunc

    real(8), intent(in) :: pressureProfile(nLev)
    real(8), intent(in) :: horizLengthScale1
    real(8), intent(in) :: horizLengthScale2
    real(8), intent(in) :: vertLengthScale

    character(len=*), intent(in) :: locType, locMode

    integer, intent(out) :: cvDim_out
    integer, intent(out) :: nEnsOverDimension_out

    integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax, mymCount, mynCount, mIndex, nIndex, maxMyNla
    integer :: myMemberBeg, myMemberEnd, myMemberCount, maxMyMemberCount, ierr

    if (verbose) write(*,*) 'Entering lsp_Setup'

    !
    !- 1.  Allocation
    !
    if (.not. associated(lsp)) then
      allocate(lsp)
    else
      call utl_abort('lsp_setup: supplied lsp must be null')
    endif

    !
    !- 2.  Settings
    !
    
    !- 2.1 Mode
    if ( trim(locType) == 'spectral') then
       if (mmpi_myid == 0) write(*,*)
       if (mmpi_myid == 0) write(*,*) 'lsp_setup: LocType = ', trim(locType)
    else
       write(*,*)
       write(*,*) 'locType = ', trim(locType)
       call utl_abort('lsp_setup: unknown locType')
    end if

    !- 2.2 Ensemble members and Levels
    lsp%nEns=nEns
    lsp%nLev=nLev
 
    !- 2.3 Horizontal grid
    lsp%ni   = hco_loc%ni
    lsp%nj   = hco_loc%nj

    call mmpi_setup_latbands(lsp%nj, latPerPE, latPerPEmax, lsp%myLatBeg, lsp%myLatEnd)
    call mmpi_setup_lonbands(lsp%ni, lonPerPE, lonPerPEmax, lsp%myLonBeg, lsp%myLonEnd)

    lsp%global = hco_loc%global
    if (lsp%global) then
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'lsp_setup: GLOBAL mode activated'
    else
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'lsp_setup: LAM mode activated'
    end if

    !- 2.4 Spectral Transform
    lsp%nTrunc=nTrunc
    lsp%dlon = hco_loc%dlon

    call mmpi_setup_levels(lsp%nEns,myMemberBeg,myMemberEnd,myMemberCount)
    call rpn_comm_allreduce(myMemberCount, maxMyMemberCount, &
                              1,"MPI_INTEGER","mpi_max","GRID",ierr)
    nEnsOverDimension_out     = mmpi_npex * maxMyMemberCount
    lsp%nEnsOverDimension = nEnsOverDimension_out

    if (lsp%global) then
       ! Global Mode
       lsp%nphase = 2
       lsp%nla_mpiglobal = (lsp%ntrunc+1)*(lsp%ntrunc+2)/2
       
       lsp%gstID = gst_setup(lsp%ni,lsp%nj,lsp%ntrunc,lsp%nEnsOverDimension)
       if (mmpi_myid == 0) write(*,*) 'lsp_setup: returned value of gstID = ',lsp%gstID

    else
       ! LAM mode
       call lst_Setup(lsp%lst,                                          & ! OUT
                      lsp%ni, lsp%nj, lsp%dlon, lsp%ntrunc, & ! IN
                      'LatLonMN', maxlevels_opt=lsp%nEnsOverDimension,  & ! IN
                       gridDataOrder_opt='kij')                               ! IN

       lsp%nphase       = lsp%lst%nphase
       lsp%nla_mpilocal = lsp%lst%nla
       
    end if

    !- 2.5 Distribute control vector over mpi processes according to member index and m
    call mmpi_setup_m(lsp%ntrunc,lsp%mymBeg,lsp%mymEnd,lsp%mymSkip,mymCount)
    call mmpi_setup_n(lsp%ntrunc,lsp%mynBeg,lsp%mynEnd,lsp%mynSkip,mynCount)

    if (lsp%global) then
       ! compute arrays to facilitate conversions between ila_mpilocal and ila_mpiglobal
       call gst_ilaList_mpiglobal(lsp%ilaList_mpiglobal,lsp%nla_mpilocal,maxMyNla,lsp%gstID, &
                                  lsp%mymBeg,lsp%mymEnd,lsp%mymSkip,lsp%mynBeg, &
                                  lsp%mynEnd,lsp%mynSkip)
       call gst_ilaList_mpilocal(lsp%ilaList_mpilocal,lsp%gstID,lsp%mymBeg,lsp%mymEnd, &
                                 lsp%mymSkip,lsp%mynBeg,lsp%mynEnd,lsp%mynSkip)
       write(*,*) 'lsp_setup: nla_mpiglobal, nla_mpilocal, maxMyNla = ', lsp%nla_mpiglobal, &
            lsp%nla_mpilocal, maxMyNla
    end if

    !- 2.6 Localization
    call lfn_setup('FifthOrder') ! IN

    call setupLocalizationMatrices(lsp, horizLengthScale1, horizLengthScale2, & ! IN
                                   vertLengthScale, pressureProfile, locMode)  ! IN

    if (lsp%global) then
       lsp%cvDim_mpiglobal = (lsp%ntrunc+1)*(lsp%ntrunc+1)*lsp%nLev*lsp%nEns
       lsp%cvDim_mpilocal  = 0
       
       do mIndex = lsp%mymBeg, lsp%mymEnd, lsp%mymSkip
          do nIndex = lsp%mynBeg, lsp%mynEnd, lsp%mynSkip
             if (mIndex.le.nIndex) then
                if (mIndex == 0) then
                   ! controlVector only contains real part for mIndex=0
                   lsp%cvDim_mpilocal = lsp%cvDim_mpilocal + 1*lsp%nLev*lsp%nEns
                else
                   ! controlVector contains real and imag parts for mIndex>0
                   lsp%cvDim_mpilocal = lsp%cvDim_mpilocal + 2*lsp%nLev*lsp%nEns
                end if
             end if
          end do
       end do
    else
       lsp%cvDim_mpiglobal = lsp%lst%nlaGlobal * lsp%nphase * lsp%nLev * lsp%nEns
       lsp%cvDim_mpilocal  = lsp%nla_mpilocal      * lsp%nphase * lsp%nLev * lsp%nEns
       write(*,*) 'cvDim_mpiglobal ', lsp%cvDim_mpiglobal, lsp%lst%nlaGlobal, &
            lsp%nphase, lsp%nLev, lsp%nEns
       write(*,*) 'cvDim_mpilocal  ', lsp%cvDim_mpilocal, lsp%nla_mpilocal, &
            lsp%nphase, lsp%nLev, lsp%nEns
    end if
    cvDim_out = lsp%cvDim_mpilocal

    !
    !- 3.  Ending
    !
    lsp%initialized = .true.

  END SUBROUTINE lsp_setup

!--------------------------------------------------------------------------
! setupLocalizationMatrices
!--------------------------------------------------------------------------
  SUBROUTINE setupLocalizationMatrices(lsp,horizLengthScale1,horizLengthScale2 ,vertLengthScale,&
                                       pressureProfile,localizationMode)
    implicit none

    type(struct_lsp), pointer    :: lsp

    real(8), intent(in) :: horizLengthScale1,horizLengthScale2 
    real(8), intent(in) :: vertLengthScale
    real(8), intent(in) :: pressureProfile(lsp%nLev)
    character(len=*), intent(in) :: localizationMode

    real(8)  :: zr,zcorr

    integer :: levIndex,levIndex1,levIndex2,ierr

    real(8) :: horizLengthScaleAll(lsp%nLev)

    if (verbose) write(*,*) 'Entering setupLocalizationMatrices'

    !
    !- 1.  Allocation
    !
    allocate(lsp%LhorizSqrt(0:lsp%nTrunc,lsp%nLev),stat=ierr)
    if (ierr.ne.0 ) then
       write(*,*) 'lsp_setup: Problem allocating memory! id=9',ierr
       call utl_abort('setupLocalizationMatrices')
    end if
    
    allocate(lsp%LvertSqrt(lsp%nLev,lsp%nLev),stat=ierr)
    if (ierr.ne.0 ) then
       write(*,*) 'bmatrixEnsemble: Problem allocating memory! id=10',ierr
       call utl_abort('setupLocalizationMatrices')
    end if
    
    !
    !- 2.  Compute HORIZONTAL localization correlation matrix
    !

    !  2.1 Determine localization length scale for each vertical level
    if ( trim(localizationMode) == 'LevelDependent' .and. horizLengthScale2 > 0.0d0 ) then
      ! vertically varying horizontal localization (linear in log P)
      do levIndex = 1, lsp%nLev
        horizLengthScaleAll(levIndex) = ( horizLengthScale1*( log(pressureProfile(levIndex)) - &
                                                              log(pressureProfile(1       )) ) +    &
                                          horizLengthScale2*( log(pressureProfile(lsp%nLev   )) - &
                                                              log(pressureProfile(levIndex)) ) ) /  &
                                         ( log(pressureProfile(lsp%nLev))-log(pressureProfile(1)) )
        if (mmpi_myid == 0) then
          write(*,*) 'loc: localization length scale (',levIndex,') = ',horizLengthScaleAll(levIndex)
        end if
      end do
    else
      ! vertically constant horizontal localization
       horizLengthScaleAll(:) = horizLengthScale1
    end if
    
    !- 2.2 Compute the matrix
    if (lsp%global) then
       call setupGlobalSpectralHLoc(lsp,horizLengthScaleAll) ! IN
    else
       call setupLamSpectralHLoc(lsp,horizLengthScaleAll) ! IN
    end if
    
    !
    !- 3.  Compute VERTICAL localization correlation matrix
    !
    if (vertLengthScale > 0.0d0 .and. lsp%nLev > 1) then

      !-  3.1 Calculate 5'th order function
      do levIndex1 = 1, lsp%nLev
        do levIndex2 = 1, lsp%nLev
          ZR = abs(log(pressureProfile(levIndex2)) - log(pressureProfile(levIndex1)))
          zcorr = lfn_response(zr,vertLengthScale)
          lsp%LvertSqrt(levIndex1,levIndex2) = zcorr
        end do
      end do

      !- 3.2 Compute sqrt of the matrix
      if (vertLengthScale > 0.0d0) then
        call utl_matSqrt(lsp%LvertSqrt(1,1),lsp%nLev,1.0d0,.false.)
      end if

    else
      lsp%LvertSqrt(:,:) = 1.d0 ! no vertical localization
    end if

  END SUBROUTINE setupLocalizationMatrices

!--------------------------------------------------------------------------
! setupGlobalSpectralHLoc
!--------------------------------------------------------------------------
  SUBROUTINE setupGlobalSpectralHLoc(lsp, local_length)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(8), intent(in)  :: local_length(lsp%nLev)

    real(8) ::   zlc,zr,zpole,zcorr

    ! NOTE: arrays passed to spectral transform are dimensioned as follows
    !       gd: lat/lon tiles and sp: member index
    real(8) :: sp_mpilocal(lsp%nla_mpilocal,lsp%nphase,lsp%nLev)
    real(8) :: sp_mpiglobal(lsp%nla_mpiglobal,lsp%nphase,lsp%nLev)
    real(8) :: sp_mympiglobal(lsp%nla_mpiglobal,lsp%nphase,lsp%nLev)
    real(8) :: zgd_gst(lsp%myLonBeg:lsp%myLonEnd,lsp%myLatBeg:lsp%myLatEnd,lsp%nLev)

    integer :: nIndex,latIndex,lonIndex,levIndex,nsize,ierr
    integer :: ila_mpiglobal,jla_mpilocal,gstID2

    if (local_length(1).gt.0.0d0) then

      gstID2 = gst_setup(lsp%ni,lsp%nj,lsp%ntrunc,lsp%nLev)

      zgd_gst(:,:,:) = 0.0d0
      do levIndex = 1, lsp%nLev
        ! Calculate 5th Order Correlation Functions in Physical Space
        zlc = 1000.0d0*local_length(levIndex)
        do latIndex = lsp%myLatBeg, lsp%myLatEnd
          zr = ec_ra * acos(gst_getrmu(latIndex,gstID2))
          zcorr = lfn_response(zr,zlc)
          do lonIndex = lsp%myLonBeg, lsp%myLonEnd
            zgd_gst(lonIndex,latIndex,levIndex) = zcorr
          end do
        end do
      end do

      ! Transform to spectral space
      sp_mpilocal(:,:,:) = 0.0d0
      call gst_setID(gstID2)
      call gst_reespe(sp_mpilocal,zgd_gst)

      ! Make mpiglobal in spectral space
      sp_mympiglobal(:,:,:) = 0.0d0
      do jla_mpilocal = 1, lsp%nla_mpilocal
        ila_mpiglobal = lsp%ilaList_mpiglobal(jla_mpilocal)
        sp_mympiglobal(ila_mpiglobal,:,:) = sp_mpilocal(jla_mpilocal,:,:)
      end do
      nsize = lsp%nla_mpiglobal*lsp%nphase*lsp%nLev
      sp_mpiglobal(:,:,:) = 0.0d0
      call rpn_comm_allreduce(sp_mympiglobal,sp_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
      
      do levIndex = 1, lsp%nLev
        do nIndex = 0, lsp%ntrunc
          lsp%LhorizSqrt(nIndex,levIndex) = sp_mpiglobal(nIndex+1,1,levIndex)
        end do
      end do

      ! Make sure it's one at the pole
      do levIndex = 1, lsp%nLev
        do  nIndex = 0, lsp%ntrunc
          lsp%LhorizSqrt(nIndex,levIndex) = abs(lsp%LhorizSqrt(nIndex,levIndex))
        end do
      end do
      do levIndex = 1, lsp%nLev
        zpole = 0.d0
        do  nIndex = 0, lsp%ntrunc
          zpole = zpole + lsp%LhorizSqrt(nIndex,levIndex)*sqrt((2.d0*nIndex+1.d0)/2.d0)
        end do
        if (zpole.le.0.d0) then
          write(*,*)'POLE VALUE NEGATIVE IN setupGlobalSpectralHLoc levIndex=',levIndex
          call utl_abort('setupGlobalSpectralHLoc')
        end if
        do nIndex = 0, lsp%ntrunc
          lsp%LhorizSqrt(nIndex,levIndex) = lsp%LhorizSqrt(nIndex,levIndex)/zpole
        end do
      end do

      ! Convert back to correlations and take sqrt
      do levIndex = 1, lsp%nLev
        do nIndex = 0, lsp%ntrunc
          lsp%LhorizSqrt(nIndex,levIndex) = sqrt( 0.5d0*lsp%LhorizSqrt(nIndex,levIndex) * &
                                                       ((2.0d0/(2.0d0*nIndex+1.0d0))**0.5d0) )
        end do
      end do

    else

       ! NO HORIZONTAL LOCALIZATION, set lsp%LhorizSqrt to 1.0 for wavenumber 0
       do levIndex = 1, lsp%nLev
          lsp%LhorizSqrt(:,levIndex) = 0.0d0
          lsp%LhorizSqrt(0,levIndex) = 1.0d0
       end do
    end if

  END SUBROUTINE setupGlobalSpectralHLoc

!--------------------------------------------------------------------------
! setupLamSpectralHLoc
!--------------------------------------------------------------------------
  SUBROUTINE setupLamSpectralHLoc(lsp, local_length)
    implicit none

    type(struct_lsp), pointer     :: lsp
    
    real(8), intent(in)  :: local_length(lsp%nLev)

    real(8), allocatable :: sp(:,:,:)
    real(8), allocatable :: gd(:,:,:)
    real(8), allocatable :: SumWeight(:)

    real(8) :: sum

    type(struct_lst)     :: lst_hloc    ! Spectral transform Parameters

    integer :: k, p, e, ila, totwvnb

    character(len=19)   :: kind

    if ( local_length(1) > 0.d0 ) then
      !
      !- 1. Enforce HORIZONTAL LOCALIZATION
      !

      !- 1.1 Setup a non-MPI spectral transform
      call lst_Setup(lst_hloc,                             & ! OUT
                     lsp%ni, lsp%nj, lsp%dlon, & ! IN
                     lsp%ntrunc, 'NoMpi')                ! IN

      !- 1.2 Create a correlation function in physical space
      allocate (gd(lsp%ni,lsp%nj,lsp%nLev))

      call lfn_CreateBiPerFunction( gd,                                  & ! OUT
                                    local_length, lsp%dlon,          & ! IN
                                    lsp%ni, lsp%nj, lsp%nLev ) ! IN

      !- 1.3 Transform to spectral space
      allocate (sp(lst_hloc%nla, lsp%nphase, lsp%nLev))

      kind = 'GridPointToSpectral'
      call lst_VarTransform(lst_hloc,           & ! IN
                            sp,                 & ! OUT
                            gd,                 & ! IN
                            kind, lsp%nLev)   ! IN
 
      !- 1.4 Compute band mean
      allocate(SumWeight(0:lsp%nTrunc))
      SumWeight  (:)  = 0.d0

      lsp%LhorizSqrt(:,:) = 0.d0
      do totwvnb = 0, lsp%ntrunc
         do e = 1, lst_hloc%nePerK(totwvnb)
            ila = lst_hloc%ilaFromEK(e,totwvnb)
            do p = 1, lst_hloc%nphase
               SumWeight(totwvnb) = SumWeight(totwvnb) + lst_hloc%Weight(ila)
               do k = 1, lsp%nLev
                  lsp%LhorizSqrt(totwvnb,k) = lsp%LhorizSqrt(totwvnb,k) + &
                                                             lst_hloc%Weight(ila) * abs(sp(ila,p,k))
                end do
             end do
         end do
      end do

      do totwvnb = 0, lsp%ntrunc
         if (SumWeight(totwvnb) /= 0.d0) then
            lsp%LhorizSqrt(totwvnb,:) = lsp%LhorizSqrt(totwvnb,:) / SumWeight(totwvnb)
         else
            lsp%LhorizSqrt(totwvnb,:) = 0.d0
         end if
      end do

      deallocate(SumWeight)

      !- 1.5 Normalization to one of correlation function from spectral densities: Part 1
!$OMP PARALLEL DO PRIVATE (totwvnb,k,sum)
      do k = 1, lsp%nLev
         sum = 0.0d0
         do totwvnb = 0, lsp%ntrunc
            sum = sum + real(totwvnb,8) * lsp%LhorizSqrt(totwvnb,k)
         end do
         do totwvnb = 0, lsp%ntrunc
            if ( sum /= 0.0d0 ) then
               lsp%LhorizSqrt(totwvnb,k) = lsp%LhorizSqrt(totwvnb,k) / sum
            else
               lsp%LhorizSqrt(totwvnb,k) = 0.d0
            end if
         end do
      end do
!$OMP END PARALLEL DO

      !- 1.6 Normalization to one of correlation function from spectral densities: Part 2

      !- 1.6.1 Spectral transform of a delta function (at the center of the domain)
      gd(:,:,:) = 0.d0
      gd(lsp%ni/2,lsp%nj/2,:) = 1.d0

      kind = 'GridPointToSpectral'
      call lst_VarTransform(lst_hloc,           & ! IN
                            sp,                 & ! OUT
                            gd,                 & ! IN
                            kind, lsp%nLev )  ! IN

      !- 1.6.2 Apply the correlation function
!$OMP PARALLEL DO PRIVATE (totwvnb,e,ila,p,k)
      do totwvnb = 0, lsp%ntrunc
         do e = 1, lst_hloc%nePerK(totwvnb)
            ila = lst_hloc%ilaFromEK(e,totwvnb)
            do p = 1, lsp%nphase
               do k = 1, lsp%nLev
                  sp(ila,p,k) = sp(ila,p,k) * lsp%LhorizSqrt(totwvnb,k) * &
                                lst_hloc%NormFactor(ila,p) * lst_hloc%NormFactorAd(ila,p)
               end do
            end do
         end do
      end do
!$OMP END PARALLEL DO

      !- 1.6.3 Move back to physical space
      kind = 'SpectralToGridPoint'
      call lst_VarTransform(lst_hloc,           & ! IN
                            sp,                 & ! IN
                            gd,                 & ! OUT
                            kind, lsp%nLev )  ! IN

      !- 1.6.4 Normalize to 1
      do k = 1, lsp%nLev
         if ( gd(lsp%ni/2,lsp%nj/2,k) <= 0.d0 ) then
            write(*,*) 'setupLamSpectralHLoc: Problem in normalization ',gd(lsp%ni/2,lsp%nj/2,k)
            call utl_abort('setupLamSpectralHLoc')
         end if
         if ( mmpi_myid == 0 ) then
           write(*,*) 'setupLamSpectralHLoc: Normalization factor = ', k, gd(lsp%ni/2,lsp%nj/2,k), 1.d0 / gd(lsp%ni/2,lsp%nj/2,k)
         end if
         lsp%LhorizSqrt(:,k) = lsp%LhorizSqrt(:,k) / gd(lsp%ni/2,lsp%nj/2,k)
      end do

      !- 1.7 Take sqrt
      lsp%LhorizSqrt(:,:) = sqrt(lsp%LhorizSqrt(:,:))

      deallocate(sp)
      deallocate(gd)

    else
      !
      !- 2. NO HORIZONTAL LOCALIZATION: set ensLocal%LhorizSqrt to 1.0 for wavenumber 0
      !
      lsp%LhorizSqrt(:,:) = 0.0d0
      lsp%LhorizSqrt(0,:) = 1.0d0
    end if

  END SUBROUTINE setupLamSpectralHLoc

!--------------------------------------------------------------------------
! lsp_Lsqrt
!--------------------------------------------------------------------------
  SUBROUTINE lsp_Lsqrt(lsp, controlVector, ensAmplitude, stepIndex)
    implicit none

    type(struct_lsp), pointer     :: lsp

    integer, intent(in)  :: stepIndex
    real(8), intent(in)  :: controlVector(lsp%cvDim_mpilocal)
    type(struct_ens)     :: ensAmplitude

    integer :: levIndex1,levIndex2,jla,p,memberIndex
    integer :: levIndex, latIndex
    character(len=19) :: kind

    real(8) ,allocatable :: sp_hLoc(:,:,:,:),sp_vhLoc(:,:,:,:) 

    real(8), pointer :: ensAmplitude_oneLev(:,:,:,:)

    allocate(sp_vhLoc(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEnsOverDimension))
    allocate(sp_hLoc (lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns             ))

    if (verbose) write(*,*) 'Entering lsp_Lsqrt'
    call lsp_check(lsp)

    !
    !- 1.  Horizontal Localization
    !
    if (lsp%global) then
       call globalSpectralHLoc( lsp,           & ! IN
                                sp_hLoc,      & ! OUT
                                controlVector ) ! IN
    else
       call lamSpectralHLoc( lsp,           & ! IN
                             sp_hLoc,      & ! OUT
                             controlVector ) ! IN
    end if

    !
    !- 2.  Vertical localization
    !
    if (lsp%nEnsOverDimension > lsp%nEns) then
       sp_vhLoc(:,:,:,lsp%nEns+1:lsp%nEnsOverDimension) = 0.0d0
    end if

!$OMP PARALLEL DO PRIVATE (memberIndex,levIndex1,levIndex2,p,jla)
    do memberIndex = 1, lsp%nEns
      sp_vhLoc(:,:,:,memberIndex) = 0.0d0
      do levIndex1 = 1, lsp%nLev
        do levIndex2 = 1, lsp%nLev
          do p = 1, lsp%nphase
            do jla = 1, lsp%nla_mpilocal
              sp_vhLoc(jla,p,levIndex1,memberIndex) = sp_vhLoc(jla,p,levIndex1,memberIndex) +  &
                          lsp%LvertSqrt(levIndex1,levIndex2)*sp_hLoc(jla,p,levIndex2,memberIndex)
            end do
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

    !
    !- 3.  Transform to gridpoint space all ensemble amplitudes
    !
    do levIndex = 1, lsp%nLev ! loop over all levels for the amplitudes

      ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,levIndex)

!$OMP PARALLEL DO PRIVATE (latIndex)
       do latIndex = lsp%myLatBeg, lsp%myLatEnd
          ensAmplitude_oneLev(:,stepIndex,:,latIndex) = 0.0d0
       end do
!$OMP END PARALLEL DO

       if (lsp%global) then
          call gst_setID(lsp%gstID)
          call gst_speree_kij(sp_vhLoc(:,:,levIndex,:),ensAmplitude_oneLev(:,stepIndex,:,:))
       else
          kind = 'SpectralToGridPoint'
          call lst_VarTransform(lsp%lst,                          & ! IN
                                sp_vhLoc(:,:,levIndex,:),             & ! IN
                                ensAmplitude_oneLev(:,stepIndex,:,:), & ! OUT
                                kind, lsp%nEnsOverDimension)        ! IN
       end if

    end do ! Loop on levels

    deallocate(sp_hLoc )
    deallocate(sp_vhLoc)
    
  END SUBROUTINE lsp_Lsqrt

!--------------------------------------------------------------------------
! globalSpectralHLoc
!--------------------------------------------------------------------------
  SUBROUTINE globalSpectralHLoc(lsp, sp_all, controlVector)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(8), intent(in)  :: controlVector(lsp%cvDim_mpilocal)
    real(8), intent(out) :: sp_all(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns)

    integer :: levIndex, mIndex, nIndex, ila_mpilocal, ila_mpiglobal, dimIndex, memberIndex 

    if (verbose) write(*,*) 'Entering globalSpectralHloc'
    call lsp_check(lsp)

    dimIndex = 0

    do memberIndex = 1, lsp%nEns

      do levIndex = 1, lsp%nLev
        do mIndex = lsp%mymBeg, lsp%mymEnd, lsp%mymSkip
          do nIndex = lsp%mynBeg, lsp%mynEnd, lsp%mynSkip
            if (mIndex .le. nIndex) then

              ila_mpiglobal = gst_getnind(mIndex,lsp%gstID) + nIndex - mIndex
              ila_mpilocal  = lsp%ilaList_mpilocal(ila_mpiglobal)
              if (mIndex == 0) then
                ! controlVector only contain real part for mIndex=0
                dimIndex = dimIndex + 1
                sp_all(ila_mpilocal,1,levIndex,memberIndex) = controlVector(dimIndex)*lsp%LhorizSqrt(nIndex,levIndex)*rsq2
                sp_all(ila_mpilocal,2,levIndex,memberIndex) = 0.0d0
              else
                ! controlVector contains real and imag parts for mIndex>0
                dimIndex = dimIndex + 1
                sp_all(ila_mpilocal,1,levIndex,memberIndex) = controlVector(dimIndex)*lsp%LhorizSqrt(nIndex,levIndex)
                dimIndex = dimIndex + 1
                sp_all(ila_mpilocal,2,levIndex,memberIndex) = controlVector(dimIndex)*lsp%LhorizSqrt(nIndex,levIndex)
              end if

            end if
          end do
        end do
      end do

      if (dimIndex.gt.lsp%cvDim_mpilocal) then
        write(*,*) 'loc globalSpectralHLoc: dimIndex > cvDim_mpilocal! ',dimIndex,memberIndex,lsp%cvDim_mpilocal
        call utl_abort('globalSpectralHLoc')
      end if

    end do

  END SUBROUTINE globalSpectralHLoc

!--------------------------------------------------------------------------
! lamSpectralHLoc
!--------------------------------------------------------------------------
  SUBROUTINE lamSpectralHLoc(lsp, sp_all, controlVector)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(8), intent(in)  :: controlVector(lsp%cvDim_mpilocal)
    real(8), intent(out) :: sp_all(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns)

    integer :: levIndex,jla, dimIndex, memberIndex, p 

    if (verbose) write(*,*) 'Entering lamSpectralHloc'
    call lsp_check(lsp)

    !
    !- Reshape + Horizontal localization + Scaling (parseval)
    !
    dimIndex = 0

    do memberIndex = 1, lsp%nEns

       do levIndex = 1, lsp%nLev
         do jla = 1, lsp%nla_mpilocal
            do p = 1, lsp%nphase
              dimIndex = dimIndex + 1
              sp_all(jla,p,levIndex,memberIndex) = controlVector(dimIndex)           * &
                                                lsp%LhorizSqrt(lsp%lst%k(jla),levIndex) * &
                                                lsp%lst%NormFactor(jla,p)
            end do
         end do
       end do
       if (dimIndex > lsp%cvDim_mpilocal ) then
          write(*,*) 'loc lamSpectralHLoc: dimIndex > cvDim! ',dimIndex,memberIndex,lsp%cvDim_mpilocal
          call utl_abort('lamSpectralHLoc')
       end if

    end do

  END SUBROUTINE lamSpectralHLoc
  
!--------------------------------------------------------------------------
! spectralLocalizationSqrtAd
!--------------------------------------------------------------------------
  SUBROUTINE lsp_LsqrtAd(lsp, ensAmplitude, controlVector, stepIndex)
    implicit none

    type(struct_lsp), pointer      :: lsp

    integer, intent(in)   :: stepIndex
    real(8), intent(out)  :: controlVector(lsp%cvDim_mpilocal)
    type(struct_ens)      :: ensAmplitude

    integer :: levIndex1,levIndex2,jla,memberIndex,p
    integer :: levIndex
    character(len=19) :: kind

    real(8) ,allocatable :: sp_hLoc(:,:,:,:),sp_vhLoc(:,:,:,:)

    real(8), pointer :: ensAmplitude_oneLev(:,:,:,:)

    if (verbose) write(*,*) 'Entering lsp_LsqrtAd'
    call lsp_check(lsp)

    allocate(sp_vhLoc(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEnsOverDimension))
    allocate(sp_hLoc (lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns             ))

    !
    !- 3.  Transform to gridpoint space all ensemble amplitudes
    !
    do levIndex = 1, lsp%nLev ! loop over all levels for the amplitudes

      ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,levIndex)

      sp_vhLoc(:,:,levIndex,:) = 0.d0 ! needed, not everything is set

      if (lsp%global) then
        call gst_setID(lsp%gstID)
        call gst_speree_kij_ad(sp_vhLoc(:,:,levIndex,:),ensAmplitude_oneLev(:,stepIndex,:,:))
      else
        kind = 'GridPointToSpectral'
        call lst_VarTransform(lsp%lst,                          & ! IN
                              sp_vhLoc(:,:,levIndex,:),             & ! OUT
                              ensAmplitude_oneLev(:,stepIndex,:,:), & ! IN
                              kind, lsp%nEnsOverDimension )       ! IN
        end if

     end do

     !
     !- 2.  Vertical Localization
     !
!$OMP PARALLEL DO PRIVATE (memberIndex,levIndex1,levIndex2,p,jla)
     do memberIndex = 1, lsp%nEns
        sp_hLoc(:,:,:,memberIndex) = 0.0d0
        do levIndex1 = 1, lsp%nLev
           do levIndex2 = 1, lsp%nLev
              do p = 1, lsp%nphase
                 do jla = 1, lsp%nla_mpilocal
                    sp_hLoc(jla,p,levIndex2,memberIndex) = sp_hLoc(jla,p,levIndex2,memberIndex) +  & 
                         lsp%LvertSqrt(levIndex2,levIndex1)*sp_vhLoc(jla,p,levIndex1,memberIndex)
                 end do
              end do
           end do
        end do
     end do
!$OMP END PARALLEL DO

    !
    !- 1.  Horizontal Localization
    !
    if (lsp%global) then
      call globalSpectralHLocAd( lsp, sp_hLoc,  & ! IN
                                 controlVector ) ! OUT
    else
      call lamSpectralHLocAd( lsp, sp_hLoc,  & ! IN
                              controlVector ) ! OUT
    end if

    deallocate(sp_hLoc )
    deallocate(sp_vhLoc)

  END SUBROUTINE lsp_LsqrtAd

!--------------------------------------------------------------------------
! globalSpectralHLocAd
!--------------------------------------------------------------------------
  SUBROUTINE globalSpectralHLocAd(lsp, sp_all, controlVector)
    implicit none

    type(struct_lsp), pointer       :: lsp

    real(8), intent(out)   :: controlVector(lsp%cvDim_mpilocal)
    real(8), intent(in)    :: sp_all(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns)

    integer :: levIndex, mIndex, nIndex, ila_mpilocal, ila_mpiglobal, dimIndex, memberIndex 

    if (verbose) write(*,*) 'Entering globalSpectralHLocAd'
    call lsp_check(lsp)

    dimIndex = 0

    do memberIndex = 1, lsp%nEns

       do levIndex = 1, lsp%nLev
          do mIndex = lsp%mymBeg, lsp%mymEnd, lsp%mymSkip
            do nIndex = lsp%mynBeg, lsp%mynEnd, lsp%mynSkip
              if (mIndex .le. nIndex) then

                ila_mpiglobal = gst_getnind(mIndex,lsp%gstID) + nIndex - mIndex
                ila_mpilocal  = lsp%ilaList_mpilocal(ila_mpiglobal)
                if (mIndex == 0) then
                  ! controlVector only contain real part for mIndex=0
                  dimIndex = dimIndex + 1
                  controlVector(dimIndex) = controlVector(dimIndex) +  &
                                            sp_all(ila_mpilocal,1,levIndex,memberIndex)*lsp%LhorizSqrt(nIndex,levIndex)*rsq2
                else
                  ! controlVector contains real and imag parts for mIndex>0
                  dimIndex = dimIndex + 1
                  controlVector(dimIndex) = controlVector(dimIndex) +  &
                                            sp_all(ila_mpilocal,1,levIndex,memberIndex)*lsp%LhorizSqrt(nIndex,levIndex)*2.0d0
                  dimIndex = dimIndex + 1
                  controlVector(dimIndex) = controlVector(dimIndex) +  &
                                            sp_all(ila_mpilocal,2,levIndex,memberIndex)*lsp%LhorizSqrt(nIndex,levIndex)*2.0d0
                end if

             end if
           end do
         end do
       end do

       if (dimIndex.gt.lsp%cvDim_mpilocal) then
          write(*,*) 'loc globalSpectralHLocAd: dimIndex > cvDim_mpilocal! ',dimIndex,memberIndex,lsp%cvDim_mpilocal
          call utl_abort('globalSpectralHLocAd')
       end if
    
    end do

  END SUBROUTINE globalSpectralHLocAd

!--------------------------------------------------------------------------
! lamSpectralHLocAd
!--------------------------------------------------------------------------
  SUBROUTINE lamSpectralHLocAd(lsp, sp_all, controlVector)
    implicit none

    type(struct_lsp), pointer       :: lsp

    real(8), intent(out)   :: controlVector(lsp%cvDim_mpilocal)
    real(8), intent(in)    :: sp_all(lsp%nla_mpilocal,lsp%nphase,lsp%nLev,lsp%nEns)

    integer :: jla, levIndex, dimIndex, memberIndex, p

    if (verbose) write(*,*) 'Entering lamSpectralHLocAd'
    call lsp_check(lsp)

    !
    !- Reshape + Horizontal localization + Scaling (parseval)
    !
    dimIndex = 0

    do memberIndex = 1, lsp%nEns

       do levIndex = 1, lsp%nLev
         do jla = 1, lsp%nla_mpilocal
           do p = 1, lsp%nphase
             dimIndex = dimIndex + 1
             controlVector(dimIndex) = controlVector(dimIndex) +           &
                                       ( sp_all(jla,p,levIndex,memberIndex)  * &
                                         lsp%LhorizSqrt(lsp%lst%k(jla),levIndex) * &
                                         lsp%lst%NormFactorAd(jla,p)    )
           end do
         end do
       end do
       if (dimIndex > lsp%cvDim_mpilocal ) then
          write(*,*) 'BEN: lamSpectralHLocAD: dimIndex > cvDim! ',dimIndex, memberIndex, lsp%cvDim_mpilocal
          call utl_abort('lamSpectralHLocAd')
       end if

    end do

  END SUBROUTINE lamSpectralHLocAd

!--------------------------------------------------------------------------
! lsp_finalize
!--------------------------------------------------------------------------
  SUBROUTINE lsp_finalize(lsp)
    implicit none

    type(struct_lsp), pointer      :: lsp

    if (verbose) write(*,*) 'Entering lsp_finalize'
    call lsp_check(lsp)

    deallocate(lsp%LhorizSqrt)
    deallocate(lsp%LvertSqrt )

  end SUBROUTINE lsp_finalize

!--------------------------------------------------------------------------
! lsp_reduceToMPILocal
!--------------------------------------------------------------------------
  SUBROUTINE lsp_reduceToMPILocal(lsp,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(8), intent(out) :: cv_mpilocal(lsp%cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    real(8), allocatable :: cv_allmaxmpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc, cvDim_maxmpilocal
    integer :: dimIndex_mpilocal, dimIndex_mpiglobal, ila_mpilocal, ila_mpiglobal
    integer :: mIndex, nIndex, memberIndex, levIndex, ierr, p, nlaMax

    if (verbose) write(*,*) 'Entering lsp_reduceToMPILocal'
    call lsp_check(lsp)

    call rpn_comm_allreduce(lsp%cvDim_mpilocal, cvDim_maxmpilocal, &
         1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    allocate(cvDim_allMpiLocal(mmpi_nprocs))
    call rpn_comm_allgather(lsp%cvDim_mpiLocal   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ierr)

    ! assign part of mpiglobal vector from current mpi process

    if (lsp%global) then

       ! Global

       allocate(allnBeg(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynBeg,1,"mpi_integer",       &
                               allnBeg,1,"mpi_integer","GRID",ierr)
       allocate(allnEnd(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynEnd,1,"mpi_integer",       &
                               allnEnd,1,"mpi_integer","GRID",ierr)
       allocate(allnSkip(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynSkip,1,"mpi_integer",       &
                               allnSkip,1,"mpi_integer","GRID",ierr)

       allocate(allmBeg(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymBeg,1,"mpi_integer",       &
                               allmBeg,1,"mpi_integer","GRID",ierr)
       allocate(allmEnd(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymEnd,1,"mpi_integer",       &
                               allmEnd,1,"mpi_integer","GRID",ierr)
       allocate(allmSkip(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymSkip,1,"mpi_integer",       &
                               allmSkip,1,"mpi_integer","GRID",ierr)


       if (mmpi_myid == 0) then

          allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,dimIndex_mpilocal,memberIndex,levIndex,mIndex,nIndex,ila_mpiglobal,dimIndex_mpiglobal)
          do jproc = 0, (mmpi_nprocs-1)
             cv_allmaxmpilocal(:,jproc+1) = 0.d0

             dimIndex_mpilocal = 0
             do memberIndex = 1, lsp%nEns

                do levIndex = 1, lsp%nLev
                   do mIndex = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                      do nIndex = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                         if (mIndex.le.nIndex) then

                            ! figure out index into global control vector
                            ila_mpiglobal = gst_getNIND(mIndex,lsp%gstID) + nIndex - mIndex
                            if (mIndex == 0) then
                               ! for mIndex=0 only real part
                               dimIndex_mpiglobal = ila_mpiglobal
                            else
                               ! for mIndex>0 both real and imaginary part
                               dimIndex_mpiglobal = 2*ila_mpiglobal-1 - (lsp%ntrunc+1)
                            end if
                            ! add offset for level
                            dimIndex_mpiglobal = dimIndex_mpiglobal + (levIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)
                            ! add offset for member index
                            dimIndex_mpiglobal = dimIndex_mpiglobal + (memberIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)*lsp%nLev
                            
                            if (mIndex == 0) then
                               ! controlVector only contain real part for mIndex=0
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)
                            else
                               ! controlVector contains real and imag parts for mIndex>0
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal+1)
                            end if
                           
                            if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                               write(*,*)
                               write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                               write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                               call utl_abort('lsp_reduceToMPILocal')
                            end if
                            if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                               write(*,*)
                               write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                               write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                               call utl_abort('lsp_reduceToMPILocal')
                            end if
 
                         end if
                      end do
                   end do
                end do
                
             end do

          end do ! procs
!$OMP END PARALLEL DO

       else
          allocate(cv_allmaxmpilocal(1,1))
       end if

       deallocate(allnBeg)
       deallocate(allnEnd)
       deallocate(allnSkip)
       deallocate(allmBeg)
       deallocate(allmEnd)
       deallocate(allmSkip)

    else
       
      ! LAM
      call rpn_comm_allreduce(lsp%lst%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ierr)

      if (mmpi_myid == 0) then
         allocate(allnlaLocal(mmpi_nprocs))
         allocate(allilaGlobal(nlaMax,mmpi_nprocs))
      else
         allocate(allnlaLocal(1))
         allocate(allilaGlobal(1,1))
      end if
      
      allocate(ilaGlobal(nlaMax))
      ilaGlobal(:)             = -1
      ilaGlobal(1:lsp%lst%nla) = lsp%lst%ilaGlobal(:)
      
      call rpn_comm_gather(lsp%lst%nla, 1, "mpi_integer",       &
                           allnlaLocal, 1, "mpi_integer", 0, "GRID", ierr)
      call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                           allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ierr)

      deallocate(ilaGlobal)

      if (mmpi_myid == 0) then

         allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

         do jproc = 0, (mmpi_nprocs-1)
            cv_allmaxmpilocal(:,jproc+1) = 0.d0
            do memberIndex = 1, lsp%nEns
               do levIndex = 1, lsp%nLev
                  do ila_mpilocal = 1, allnlaLocal(jproc+1)
                     do p = 1, lsp%lst%nphase

                        dimIndex_mpilocal = ( (levIndex-1) * lsp%nEns * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                        ( (memberIndex-1) * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                                               ( (ila_mpilocal-1) * lsp%lst%nphase ) + p

                        ila_mpiglobal = allilaGlobal(ila_mpilocal,jproc+1)
                        if ( ila_mpiglobal <= 0 ) then 
                           write(*,*) 'lsp_reduceToMPILocal: invalid ila_mpiglobal index ', ila_mpiglobal
                           call utl_abort('lsp_reduceToMPILocal')
                        end if
                        dimIndex_mpiglobal = ( (levIndex-1) * lsp%nEns * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                         ( (memberIndex-1) * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                                           ( (ila_mpiglobal-1) * lsp%lst%nphase ) + p
  
                        if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_reduceToMPILocal')
                         end if
                         if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_reduceToMPILocal')
                         end if
                  
                         cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)

                     end do
                  end do
               end do
            end do
         end do
      else
         allocate(cv_allmaxmpilocal(1,1))
      end if

     deallocate(allnlaLocal)
     deallocate(allilaGlobal) 
      
   end if

   !- Distribute
   allocate(displs(mmpi_nprocs))
   do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
   end do

   call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_double_precision", &
                          cv_mpiLocal, lsp%cvDim_mpiLocal, "mpi_double_precision", &
                          0, "GRID", ierr)

   deallocate(displs) 
   deallocate(cv_allMaxMpiLocal)
   deallocate(cvDim_allMpiLocal)
   
 END SUBROUTINE Lsp_reduceToMPILocal

!--------------------------------------------------------------------------
! lsp_reduceToMPILocal_r4
!--------------------------------------------------------------------------
  SUBROUTINE lsp_reduceToMPILocal_r4(lsp,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(4), intent(out) :: cv_mpilocal(lsp%cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    real(4), allocatable :: cv_allmaxmpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc, cvDim_maxmpilocal
    integer :: dimIndex_mpilocal, dimIndex_mpiglobal, ila_mpilocal, ila_mpiglobal
    integer :: mIndex, nIndex, memberIndex, levIndex, ierr, p, nlaMax

   if (verbose) write(*,*) 'Entering lsp_reduceToMPILocal_r4'
    call lsp_check(lsp)

    call rpn_comm_allreduce(lsp%cvDim_mpilocal, cvDim_maxmpilocal, &
         1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    allocate(cvDim_allMpiLocal(mmpi_nprocs))
    call rpn_comm_allgather(lsp%cvDim_mpiLocal   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ierr)

    ! assign part of mpiglobal vector from current mpi process

    if (lsp%global) then

       ! Global

       allocate(allnBeg(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynBeg,1,"mpi_integer",       &
                               allnBeg,1,"mpi_integer","GRID",ierr)
       allocate(allnEnd(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynEnd,1,"mpi_integer",       &
                               allnEnd,1,"mpi_integer","GRID",ierr)
       allocate(allnSkip(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mynSkip,1,"mpi_integer",       &
                               allnSkip,1,"mpi_integer","GRID",ierr)

       allocate(allmBeg(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymBeg,1,"mpi_integer",       &
                               allmBeg,1,"mpi_integer","GRID",ierr)
       allocate(allmEnd(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymEnd,1,"mpi_integer",       &
                               allmEnd,1,"mpi_integer","GRID",ierr)
       allocate(allmSkip(mmpi_nprocs))
       call rpn_comm_allgather(lsp%mymSkip,1,"mpi_integer",       &
                               allmSkip,1,"mpi_integer","GRID",ierr)


       if (mmpi_myid == 0) then

          allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,dimIndex_mpilocal,memberIndex,levIndex,mIndex,nIndex,ila_mpiglobal,dimIndex_mpiglobal)
          do jproc = 0, (mmpi_nprocs-1)
             cv_allmaxmpilocal(:,jproc+1) = 0.d0

             dimIndex_mpilocal = 0
             do memberIndex = 1, lsp%nEns

                do levIndex = 1, lsp%nLev
                   do mIndex = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                      do nIndex = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                         if (mIndex.le.nIndex) then

                            ! figure out index into global control vector
                            ila_mpiglobal = gst_getNIND(mIndex,lsp%gstID) + nIndex - mIndex
                            if (mIndex == 0) then
                               ! for mIndex=0 only real part
                               dimIndex_mpiglobal = ila_mpiglobal
                            else
                               ! for mIndex>0 both real and imaginary part
                               dimIndex_mpiglobal = 2*ila_mpiglobal-1 - (lsp%ntrunc+1)
                            end if
                            ! add offset for level
                            dimIndex_mpiglobal = dimIndex_mpiglobal + (levIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)
                            ! add offset for member index
                            dimIndex_mpiglobal = dimIndex_mpiglobal + (memberIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)*lsp%nLev
                            
                            if (mIndex == 0) then
                               ! controlVector only contain real part for mIndex=0
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)
                            else
                               ! controlVector contains real and imag parts for mIndex>0
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)
                               dimIndex_mpilocal = dimIndex_mpilocal + 1
                               cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal+1)
                            end if
                           
                            if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                               write(*,*)
                               write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                               write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                               call utl_abort('lsp_reduceToMPILocal')
                            end if
                            if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                               write(*,*)
                               write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                               write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                               call utl_abort('lsp_reduceToMPILocal')
                            end if
 
                         end if
                      end do
                   end do
                end do
                
             end do

          end do ! procs
!$OMP END PARALLEL DO

       else
          allocate(cv_allmaxmpilocal(1,1))
       end if

       deallocate(allnBeg)
       deallocate(allnEnd)
       deallocate(allnSkip)
       deallocate(allmBeg)
       deallocate(allmEnd)
       deallocate(allmSkip)

    else
       
      ! LAM
      call rpn_comm_allreduce(lsp%lst%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ierr)

      if (mmpi_myid == 0) then
         allocate(allnlaLocal(mmpi_nprocs))
         allocate(allilaGlobal(nlaMax,mmpi_nprocs))
      else
         allocate(allnlaLocal(1))
         allocate(allilaGlobal(1,1))
      end if
      
      allocate(ilaGlobal(nlaMax))
      ilaGlobal(:)             = -1
      ilaGlobal(1:lsp%lst%nla) = lsp%lst%ilaGlobal(:)
      
      call rpn_comm_gather(lsp%lst%nla, 1, "mpi_integer",       &
                           allnlaLocal, 1, "mpi_integer", 0, "GRID", ierr)
      call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                           allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ierr)

      deallocate(ilaGlobal)

      if (mmpi_myid == 0) then

         allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

         do jproc = 0, (mmpi_nprocs-1)
            cv_allmaxmpilocal(:,jproc+1) = 0.d0
            do memberIndex = 1, lsp%nEns
               do levIndex = 1, lsp%nLev
                  do ila_mpilocal = 1, allnlaLocal(jproc+1)
                     do p = 1, lsp%lst%nphase

                        dimIndex_mpilocal = ( (levIndex-1) * lsp%nEns * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                        ( (memberIndex-1) * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                                               ( (ila_mpilocal-1) * lsp%lst%nphase ) + p

                        ila_mpiglobal = allilaGlobal(ila_mpilocal,jproc+1)
                        if ( ila_mpiglobal <= 0 ) then 
                           write(*,*) 'lsp_reduceToMPILocal: invalid ila_mpiglobal index ', ila_mpiglobal
                           call utl_abort('lsp_reduceToMPILocal')
                        end if
                        dimIndex_mpiglobal = ( (levIndex-1) * lsp%nEns * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                         ( (memberIndex-1) * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                                           ( (ila_mpiglobal-1) * lsp%lst%nphase ) + p
  
                        if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_reduceToMPILocal')
                         end if
                         if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_reduceToMPILocal')
                         end if
                  
                         cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1) = cv_mpiglobal(dimIndex_mpiglobal)

                     end do
                  end do
               end do
            end do
         end do
      else
         allocate(cv_allmaxmpilocal(1,1))
      end if

     deallocate(allnlaLocal)
     deallocate(allilaGlobal) 
      
   end if

   !- Distribute
   allocate(displs(mmpi_nprocs))
   do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
   end do

   call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_real4", &
                          cv_mpiLocal, lsp%cvDim_mpiLocal, "mpi_real4", &
                          0, "GRID", ierr)

   deallocate(displs) 
   deallocate(cv_allMaxMpiLocal)
   deallocate(cvDim_allMpiLocal)
   
 END SUBROUTINE Lsp_reduceToMPILocal_r4

!--------------------------------------------------------------------------
! lsp_expandToMPIGlobal
!--------------------------------------------------------------------------
  SUBROUTINE lsp_expandToMPIGlobal(lsp,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(8), intent(in)  :: cv_mpilocal(lsp%cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer     :: cv_allmaxmpilocal(:,:) => null()

    integer, allocatable :: cvDim_allMpilocal(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: dimIndex_mpilocal, dimIndex_mpiglobal, ila_mpiglobal, ila_mpilocal, cvDim_maxmpilocal
    integer :: mIndex, nIndex, jproc, memberIndex, levIndex, ierr, p, nlaMax

    if (verbose) write(*,*) 'Entering lsp_expandToMPIGlobal'
    call lsp_check(lsp)

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    allocate(cvDim_allMpiLocal(mmpi_nprocs))
    call rpn_comm_allgather(lsp%cvDim_mpiLocal   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ierr)

    call rpn_comm_allreduce(lsp%cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:lsp%cvDim_mpilocal) = cv_mpilocal(1:lsp%cvDim_mpilocal)

    nullify(cv_allmaxmpilocal)
    if (mmpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if
    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if (lsp%global) then

       ! Global
       if (mmpi_myid == 0) then
          allocate(allnBeg(mmpi_nprocs))
          allocate(allnEnd(mmpi_nprocs))
          allocate(allnSkip(mmpi_nprocs))
          allocate(allmBeg(mmpi_nprocs))
          allocate(allmEnd(mmpi_nprocs))
          allocate(allmSkip(mmpi_nprocs))
       else
          allocate(allnBeg(1))
          allocate(allnEnd(1))
          allocate(allnSkip(1))
          allocate(allmBeg(1))
          allocate(allmEnd(1))
          allocate(allmSkip(1))
       end if

       call rpn_comm_gather(lsp%mynBeg  ,1,"mpi_integer",       &
                            allnBeg ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mynEnd  ,1,"mpi_integer",       &
                            allnEnd ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mynSkip ,1,"mpi_integer",       &
                            allnSkip,1,"mpi_integer",0,"GRID",ierr)

       call rpn_comm_gather(lsp%mymBeg  ,1,"mpi_integer",       &
                            allmBeg ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mymEnd  ,1,"mpi_integer",       &
                            allmEnd ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mymSkip ,1,"mpi_integer",       &
                            allmSkip,1,"mpi_integer",0,"GRID",ierr)

       ! reorganize gathered mpilocal control vectors into the mpiglobal control vector
       if (mmpi_myid == 0) then
         cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,dimIndex_mpilocal,memberIndex,levIndex,mIndex,nIndex,ila_mpiglobal,dimIndex_mpiglobal)
         do jproc = 0, (mmpi_nprocs-1)
           dimIndex_mpilocal = 0
           do memberIndex = 1, lsp%nEns

             do levIndex = 1, lsp%nLev
               do mIndex = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                 do nIndex = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
                   if (mIndex.le.nIndex) then

                     ! figure out index into global control vector
                     ila_mpiglobal = gst_getNIND(mIndex,lsp%gstID) + nIndex - mIndex
                     if (mIndex == 0) then
                       ! for mIndex=0 only real part
                       dimIndex_mpiglobal = ila_mpiglobal
                     else
                       ! for mIndex>0 both real and imaginary part
                       dimIndex_mpiglobal = 2*ila_mpiglobal-1 - (lsp%ntrunc+1)
                     end if
                     ! add offset for level
                     dimIndex_mpiglobal = dimIndex_mpiglobal + (levIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)
                     ! add offset for member index
                     dimIndex_mpiglobal = dimIndex_mpiglobal + (memberIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)*lsp%nLev

                     ! index into local control vector
                     if (mIndex == 0) then
                       ! only real component for mIndex=0
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                     else
                       ! both real and imaginary components for mIndex>0
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal+1) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                     end if

                     if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                        write(*,*)
                        write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                        write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                        call utl_abort('lsp_expandToMPIGlobal')
                     end if
                     if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                        write(*,*)
                        write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                        write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                        call utl_abort('lsp_expandToMPIGlobal')
                     end if

                   end if
                 end do
               end do
             end do
           end do
         end do ! jproc
!$OMP END PARALLEL DO

      end if ! myid == 0 

      deallocate(allnBeg)
      deallocate(allnEnd)
      deallocate(allnSkip)
      deallocate(allmBeg)
      deallocate(allmEnd)
      deallocate(allmSkip)

    else

      ! LAM
       call rpn_comm_allreduce(lsp%lst%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ierr)

       if (mmpi_myid == 0) then
          allocate(allnlaLocal(mmpi_nprocs))
          allocate(allilaGlobal(nlaMax,mmpi_nprocs))
       else
          allocate(allnlaLocal(1))
          allocate(allilaGlobal(1,1))
       end if

       allocate(ilaGlobal(nlaMax))
       ilaGlobal(:)             = -1
       ilaGlobal(1:lsp%lst%nla) = lsp%lst%ilaGlobal(:)

       call rpn_comm_gather(lsp%lst%nla, 1, "mpi_integer",       &
                            allnlaLocal, 1, "mpi_integer", 0, "GRID", ierr)
       call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                            allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ierr)

       deallocate(ilaGlobal)

       if (mmpi_myid == 0) then
          cv_mpiglobal(:) = 0.0d0

          do jproc = 0, (mmpi_nprocs-1)
             do memberIndex = 1, lsp%nEns
                do levIndex = 1, lsp%nLev
                   do ila_mpilocal = 1, allnlaLocal(jproc+1)
                      do p = 1, lsp%lst%nphase

                         dimIndex_mpilocal = ( (levIndex-1) * lsp%nEns * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                         ( (memberIndex-1) * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                                      ( (ila_mpilocal-1) * lsp%lst%nphase ) + p

                         ila_mpiglobal = allilaGlobal(ila_mpilocal,jproc+1)
                         if ( ila_mpiglobal <= 0 ) then 
                            write(*,*) 'lsp_expandToMPIGlobal: invalid ila_mpiglobal index ', ila_mpiglobal
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if

                         dimIndex_mpiglobal = ( (levIndex-1) * lsp%nEns * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                          ( (memberIndex-1) * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                                            ( (ila_mpiglobal-1) * lsp%lst%nphase ) + p

                         if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if
                         if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if

                         cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                         
                      end do
                   end do
                end do
             end do
          end do

       end if

       deallocate(allnlaLocal)
       deallocate(allilaGlobal)

    end if

    deallocate(cv_allmaxmpilocal)
    deallocate(cvDim_allMpiLocal)

  end SUBROUTINE Lsp_expandToMPIGlobal

!--------------------------------------------------------------------------
! lsp_expandToMPIGlobal_r4
!--------------------------------------------------------------------------
  SUBROUTINE lsp_expandToMPIGlobal_r4(lsp,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_lsp), pointer     :: lsp

    real(4), intent(in)  :: cv_mpilocal(lsp%cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer     :: cv_allmaxmpilocal(:,:) => null()

    integer, allocatable :: cvDim_allMpilocal(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: dimIndex_mpilocal, dimIndex_mpiglobal, ila_mpiglobal, ila_mpilocal, cvDim_maxmpilocal
    integer :: mIndex, nIndex, jproc, memberIndex, levIndex, ierr, p, nlaMax

    if (verbose) write(*,*) 'Entering lsp_expandToMPIGlobal_r4'
    call lsp_check(lsp)

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    allocate(cvDim_allMpiLocal(mmpi_nprocs))
    call rpn_comm_allgather(lsp%cvDim_mpiLocal   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ierr)

    call rpn_comm_allreduce(lsp%cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:lsp%cvDim_mpilocal) = cv_mpilocal(1:lsp%cvDim_mpilocal)

    nullify(cv_allmaxmpilocal)
    if (mmpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if
    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if (lsp%global) then

       ! Global
       if (mmpi_myid == 0) then
          allocate(allnBeg(mmpi_nprocs))
          allocate(allnEnd(mmpi_nprocs))
          allocate(allnSkip(mmpi_nprocs))
          allocate(allmBeg(mmpi_nprocs))
          allocate(allmEnd(mmpi_nprocs))
          allocate(allmSkip(mmpi_nprocs))
       else
          allocate(allnBeg(1))
          allocate(allnEnd(1))
          allocate(allnSkip(1))
          allocate(allmBeg(1))
          allocate(allmEnd(1))
          allocate(allmSkip(1))
       end if

       call rpn_comm_gather(lsp%mynBeg  ,1,"mpi_integer",       &
                            allnBeg ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mynEnd  ,1,"mpi_integer",       &
                            allnEnd ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mynSkip ,1,"mpi_integer",       &
                            allnSkip,1,"mpi_integer",0,"GRID",ierr)

       call rpn_comm_gather(lsp%mymBeg  ,1,"mpi_integer",       &
                            allmBeg ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mymEnd  ,1,"mpi_integer",       &
                            allmEnd ,1,"mpi_integer",0,"GRID",ierr)
       call rpn_comm_gather(lsp%mymSkip ,1,"mpi_integer",       &
                            allmSkip,1,"mpi_integer",0,"GRID",ierr)

       ! reorganize gathered mpilocal control vectors into the mpiglobal control vector
       if (mmpi_myid == 0) then
         cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,dimIndex_mpilocal,memberIndex,levIndex,mIndex,nIndex,ila_mpiglobal,dimIndex_mpiglobal)
         do jproc = 0, (mmpi_nprocs-1)
           dimIndex_mpilocal = 0
           do memberIndex = 1, lsp%nEns

             do levIndex = 1, lsp%nLev
               do mIndex = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                 do nIndex = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
                   if (mIndex.le.nIndex) then

                     ! figure out index into global control vector
                     ila_mpiglobal = gst_getNIND(mIndex,lsp%gstID) + nIndex - mIndex
                     if (mIndex == 0) then
                       ! for mIndex=0 only real part
                       dimIndex_mpiglobal = ila_mpiglobal
                     else
                       ! for mIndex>0 both real and imaginary part
                       dimIndex_mpiglobal = 2*ila_mpiglobal-1 - (lsp%ntrunc+1)
                     end if
                     ! add offset for level
                     dimIndex_mpiglobal = dimIndex_mpiglobal + (levIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)
                     ! add offset for member index
                     dimIndex_mpiglobal = dimIndex_mpiglobal + (memberIndex-1) * (lsp%ntrunc+1)*(lsp%ntrunc+1)*lsp%nLev

                     ! index into local control vector
                     if (mIndex == 0) then
                       ! only real component for mIndex=0
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                     else
                       ! both real and imaginary components for mIndex>0
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                       dimIndex_mpilocal = dimIndex_mpilocal + 1
                       cv_mpiglobal(dimIndex_mpiglobal+1) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                     end if

                     if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                        write(*,*)
                        write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                        write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                        call utl_abort('lsp_expandToMPIGlobal')
                     end if
                     if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                        write(*,*)
                        write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                        write(*,*) '       proc, levIndex, nIndex, mIndex = ',jproc, levIndex, nIndex, mIndex
                        call utl_abort('lsp_expandToMPIGlobal')
                     end if

                   end if
                 end do
               end do
             end do
           end do
         end do ! jproc
!$OMP END PARALLEL DO

      end if ! myid == 0 

      deallocate(allnBeg)
      deallocate(allnEnd)
      deallocate(allnSkip)
      deallocate(allmBeg)
      deallocate(allmEnd)
      deallocate(allmSkip)

    else

      ! LAM
       call rpn_comm_allreduce(lsp%lst%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ierr)

       if (mmpi_myid == 0) then
          allocate(allnlaLocal(mmpi_nprocs))
          allocate(allilaGlobal(nlaMax,mmpi_nprocs))
       else
          allocate(allnlaLocal(1))
          allocate(allilaGlobal(1,1))
       end if

       allocate(ilaGlobal(nlaMax))
       ilaGlobal(:)             = -1
       ilaGlobal(1:lsp%lst%nla) = lsp%lst%ilaGlobal(:)

       call rpn_comm_gather(lsp%lst%nla, 1, "mpi_integer",       &
                            allnlaLocal, 1, "mpi_integer", 0, "GRID", ierr)
       call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                            allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ierr)

       deallocate(ilaGlobal)

       if (mmpi_myid == 0) then
          cv_mpiglobal(:) = 0.0d0

          do jproc = 0, (mmpi_nprocs-1)
             do memberIndex = 1, lsp%nEns
                do levIndex = 1, lsp%nLev
                   do ila_mpilocal = 1, allnlaLocal(jproc+1)
                      do p = 1, lsp%lst%nphase

                         dimIndex_mpilocal = ( (levIndex-1) * lsp%nEns * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                         ( (memberIndex-1) * allnlaLocal(jproc+1) * lsp%lst%nphase ) + &
                                                      ( (ila_mpilocal-1) * lsp%lst%nphase ) + p

                         ila_mpiglobal = allilaGlobal(ila_mpilocal,jproc+1)
                         if ( ila_mpiglobal <= 0 ) then 
                            write(*,*) 'lsp_expandToMPIGlobal: invalid ila_mpiglobal index ', ila_mpiglobal
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if

                         dimIndex_mpiglobal = ( (levIndex-1) * lsp%nEns * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                          ( (memberIndex-1) * lsp%lst%nlaGlobal * lsp%lst%nphase ) + &
                                                            ( (ila_mpiglobal-1) * lsp%lst%nphase ) + p

                         if (dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpilocal > cvDim_allMpiLocal(jproc+1)', dimIndex_mpilocal, cvDim_allMpiLocal(jproc+1)
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if
                         if (dimIndex_mpiglobal > lsp%cvDim_mpiglobal) then
                            write(*,*)
                            write(*,*) 'ERROR: dimIndex_mpiglobal > cvDim_mpiglobal', dimIndex_mpiglobal, lsp%cvDim_mpiglobal
                            write(*,*) '       proc, memberIndex, levIndex, ila, p = ',jproc,memberIndex,levIndex,ila_mpilocal,p
                            call utl_abort('lsp_expandToMPIGlobal')
                         end if

                         cv_mpiglobal(dimIndex_mpiglobal) = cv_allmaxmpilocal(dimIndex_mpilocal,jproc+1)
                         
                      end do
                   end do
                end do
             end do
          end do

       end if

       deallocate(allnlaLocal)
       deallocate(allilaGlobal)

    end if

    deallocate(cv_allmaxmpilocal)
    deallocate(cvDim_allMpiLocal)

  end SUBROUTINE lsp_expandToMPIGlobal_r4

!--------------------------------------------------------------------------
!   LSP_CHECK
!--------------------------------------------------------------------------
  subroutine lsp_check(lsp)
    implicit none

    type(struct_lsp), pointer :: lsp

    if ( .not. lsp%initialized) then
       call utl_abort('lsp_check: structure not initialized')
    end if

  end subroutine lsp_check

end module localizationSpectral_mod
