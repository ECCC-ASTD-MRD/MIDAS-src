
module BmatrixChem_mod 
  ! MODULE BmatrixChem_mod (prefix='bchm' category='2. B and R matrices')
  !
  ! :Purpose: Contains routines involving the application of
  !           background-error covariance matrix(ces). Matrix based on
  !           horizontally homogeneous/isotropic correlations. This module
  !           includes the transform from control vector (spectral space) to 
  !           analysis increments, related utilites, and the transform's adjoint.
  !
  !           Based on elements of bmatrixHI_mod.ftn90
  !
  ! :Comments:
  !
  !   1. Covariances uncoupled from weather variable.
  !
  !   2. One could potentially make public the functions/routines which are
  !      identical to those in bmatrixhi_mod.ftn90 (except possibly in name) so
  !      that one copy is present in the code.
  !
  !
  ! Public Subroutines:
  !
  !    bchm_setupCH:  Must be called first. 
  !                   Acquire constituents backgound error standard 
  !                   deviations and spectral space correlations which are 
  !                   read and prepared by bcsc_setupCH. 
  !    bchm_finalize: Deallocate internal module arrays.
  !    bchm_BSqrt:    Transformations from control vector to analysis
  !                   increments in the minimization process.
  !    bchm_BSqrtAd:  Adjoint of bchm_BSqrt.
  !    bchm_expand*   MPI manipulations of contol vector(s)
  !    bchm_reduce*   MPI manipulations related to contol vector(s)
  !
  
  use midasMpi_mod
  use gridStateVector_mod
  use gridVariableTransforms_mod
  use globalSpectralTransform_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use bCovarSetupChem_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: bchm_setupCH,bchm_finalize,bchm_BSqrt,bchm_BSqrtAd
  public :: bchm_expandToMPIglobal,bchm_expandToMPIglobal_r4,bchm_reduceToMPIlocal,bchm_reduceToMPIlocal_r4
  
  logical             :: initialized = .false.                    
  integer             :: nla_mpiglobal,nla_mpilocal           
  integer             :: cvDim_mpilocal,cvDim_mpiglobal           
  integer             :: gstID, gstID2          

  integer             :: mymBeg,mymEnd,mymSkip,mymCount
  integer             :: mynBeg,mynEnd,mynSkip,mynCount
  integer             :: maxMyNla
  integer             :: myLatBeg,myLatEnd
  integer             :: myLonBeg,myLonEnd
  integer, pointer    :: ilaList_mpiglobal(:)
  integer, pointer    :: ilaList_mpilocal(:)
                            
  integer, external   :: get_max_rss

  ! Bacgkround error covariance matrix elements.
  ! One could add an additional dimension to corns  
  ! for separate block-univariate correlation matrices.
  ! This would also permit merging of bmatrixhi_mod and bmatrixchem_mod
  ! into one module.

  character(len=20)   :: TransformVarKindCH


  type(struct_bcsc_bgStats) :: bgStats ! Background covariances
                                       ! and related elements

  !*************************************************************************
    
  contains

  !--------------------------------------------------------------------------
  ! bchm_setupCH
  !--------------------------------------------------------------------------
  subroutine bchm_setupCH(hco_in,vco_in,cvDim_out)
   
    !:Purpose: Acquire constituents backgound error standard deviations (stddev),
    !          spectral space correlations (corns) and related elements
    !          which are read and prepared by lower level module routine
    !          bcsc_setupCH.
    !

    implicit none

    !Arguments
    type(struct_hco), intent(in), pointer :: hco_in
    type(struct_vco), intent(in) ,pointer :: vco_in
    integer, intent(out) :: cvDim_out

    !Locals
    integer :: jn,jm
    integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax, nlev_T_even, ierr
    logical :: covarNeeded

    ! Read and prepare covariances and related elements
    
    call bcsc_setupCH(hco_in,vco_in,covarNeeded,'Analysis')
    if (.not.covarNeeded) then
      ! Assumes CH covariances not needed.
      cvDim_out=0
      return
    end if
    
    ! Get covariances and related elements required for assimilation
        
    call bcsc_getCovarCH(bgStats,transformVarKind_opt=transformVarKindCH)
        
    ! Set vertical dimension
    ! Need an even number of levels for spectral transform
    
    if (mod(bgStats%nlev,2) /= 0) then
      nLev_T_even = bgStats%nlev+1
    else
      nLev_T_even = bgStats%nlev
    end if
     
    ! Spectral transform and MPI setup parameters.

    nla_mpiglobal = (bgStats%ntrunc+1)*(bgStats%ntrunc+2)/2    
    gstID  = gst_setup(bgStats%ni,bgStats%nj,bgStats%ntrunc,bgStats%nkgdim)
    gstID2 = gst_setup(bgStats%ni,bgStats%nj,bgStats%ntrunc,nlev_T_even)
    if (mmpi_myid == 0) write(*,*) 'bchm_setupCH: returned value of gstID =',gstID
    if (mmpi_myid == 0) write(*,*) 'bchm_setupCH: returned value of gstID2=',gstID2

    call mmpi_setup_latbands(bgStats%nj, latPerPE, latPerPEmax, myLatBeg, &
      myLatEnd)
    call mmpi_setup_lonbands(bgStats%ni, lonPerPE, lonPerPEmax, myLonBeg, &
      myLonEnd)

    call mmpi_setup_m(bgStats%ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
    call mmpi_setup_n(bgStats%ntrunc,mynBeg,mynEnd,mynSkip,mynCount)

    call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID, &
      mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
    call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip, &
      mynBeg,mynEnd,mynSkip)

    ! Compute mpilocal control vector size
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if (jm <= jn) then
          if (jm == 0) then
            ! only real component for jm=0
            cvDim_mpilocal = cvDim_mpilocal + 1*bgStats%nkgdim
          else
            ! both real and imaginary components for jm>0
            cvDim_mpilocal = cvDim_mpilocal + 2*bgStats%nkgdim
          end if
        end if
      end do
    end do
    cvDim_out = cvDim_mpilocal

    ! Also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_mpiglobal, 1, &
                            "mpi_integer", "mpi_sum", "GRID", ierr)

    initialized = .true.
   
  end subroutine bchm_setupCH

  !--------------------------------------------------------------------------
  ! bchm_bSqrt
  !--------------------------------------------------------------------------
  subroutine bchm_bSqrt(controlvector_in,statevector, stateVectorRef_opt)
    !
    !:Purpose: To apply B_CHM^1/2 to a control vector.
    !
    ! Based on bhi_bsqrt
    
    implicit none

    !Arguments
    real(8), intent(inout) :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv), intent(inout) :: statevector
    type(struct_gsv), intent(in) , optional :: stateVectorRef_opt
    
    !Locals
    real(8) ,allocatable :: gd_out(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,bgStats%nkgdim)
    character(len=30) :: transform
    integer :: varIndex
    
    if (.not.initialized) return

    if (mmpi_myid == 0) write(*,*) 'bchm_bsqrt: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim))

    call bchm_cain(controlVector_in,hiControlVector)

    call bchm_spa2gd(hiControlVector,gd_out)
    
    call bchm_copyToStatevector(statevector,gd_out)

    if ( trim(transformVarKindCH) /= '' ) then  
     
      transform = trim(transformVarKindCH)//'CH_tlm'
      do varIndex= 1, bgStats%numvar3d+bgStats%numvar2d
   
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
    
        if ( present(stateVectorRef_opt) ) then
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              stateVectorRef_opt=stateVectorRef_opt, & ! IN
			      varName_opt=bgStats%varNameList(varIndex) ) ! IN
        else
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              varName_opt=bgStats%varNameList(varIndex) ) ! IN
        end if

      end do
    end if
    
    deallocate(gd_out)

    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mmpi_myid == 0) write(*,*) 'bchm_bsqrt: done'

  end subroutine bchm_bSqrt

  !--------------------------------------------------------------------------
  ! bchm_bSqrtAd
  !--------------------------------------------------------------------------
  subroutine bchm_bSqrtAd(statevector,controlVector_out, stateVectorRef_opt)
    !
    !:Purpose: To apply adjoint of B_CHM^1/2 to a statevector.
    !
    !
    ! Based on bhi_bSqrtAd.
    
    implicit none

    !Arguments
    real(8), intent(inout) :: controlVector_out(cvDim_mpilocal)
    type(struct_gsv), intent(inout) :: statevector
    type(struct_gsv), intent(in), optional :: stateVectorRef_opt
    
    !Locals
    real(8), allocatable :: gd_in(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,bgStats%nkgdim)
    character(len=30) :: transform
    integer :: varIndex  

    if ( .not.initialized ) then
      if (mmpi_myid == 0) write(*,*) 'bMatrixChem not initialized'
      return
    end if

    if (mmpi_myid == 0) write(*,*) 'bchm_bsqrtad: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim))

    if ( trim(transformVarKindCH) /= '' ) then  
          
      transform = trim(transformVarKindCH)//'CH_ad'            
      do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
   
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
    
        if ( present(stateVectorRef_opt) ) then
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              stateVectorRef_opt=stateVectorRef_opt, &  ! IN
			      varName_opt=bgStats%varNameList(varIndex) ) ! IN
        else
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              varName_opt=bgStats%varNameList(varIndex) ) ! IN
        end if

      end do
    end if

    call bchm_copyFromStatevector(statevector,gd_in)

    call bchm_spa2gdad(gd_in,hiControlVector)

    call bchm_cainad(hiControlVector,controlVector_out)

    deallocate(gd_in)

    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mmpi_myid == 0) write(*,*) 'bchm_bsqrtad: done'

  end subroutine bchm_bSqrtAd

  !--------------------------------------------------------------------------
  ! bchm_cain
  !--------------------------------------------------------------------------
  subroutine bchm_cain(controlVector_in,hiControlVector_out)
    !
    implicit none

    !Arguments
    real(8), intent(inout) :: controlVector_in(cvDim_mpilocal)
    real(8), intent(inout) :: hiControlVector_out(nla_mpilocal,2,bgStats%nkgdim)

    !Locals
    integer :: jdim, levelIndex, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    hiControlVector_out(:,:,:) = 0.0d0
    do levelIndex = 1, bgStats%nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if (jm <= jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if (jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,levelIndex) = controlVector_in(jdim)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,levelIndex) = controlVector_in(jdim)
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,2,levelIndex) = controlVector_in(jdim)
            end if
          end if
        end do
      end do
    end do

  end subroutine bchm_cain

  !--------------------------------------------------------------------------
  ! bchm_cainAd
  !--------------------------------------------------------------------------
  subroutine bchm_cainAd(hiControlVector_in,controlVector_out)
    implicit none

    !Arguments
    real(8), intent(inout) :: controlVector_out(cvDim_mpilocal)
    real(8), intent(inout) :: hiControlVector_in(nla_mpilocal,2,bgStats%nkgdim)

    !Locals
    integer :: jdim, levelIndex, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    do levelIndex = 1, bgStats%nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if (jm <= jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if (jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,levelIndex)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,levelIndex)*2.0d0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,2,levelIndex)*2.0d0
            end if
          end if
        end do
      end do
    end do

  end subroutine bchm_cainAd

  !--------------------------------------------------------------------------
  ! bchm_spa2gd
  !--------------------------------------------------------------------------
  subroutine bchm_spa2gd(hiControlVector_in,gd_out)
    implicit none

    !Arguments
    real(8), intent(inout) :: hiControlVector_in(nla_mpilocal,2,bgStats%nkgdim)
    real(8), intent(inout) :: gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)

    !Locals
    real(8) :: sp(nla_mpilocal,2,bgStats%nkgdim)

    integer :: jn,jm,ila_mpilocal,ila_mpiglobal,icount
    real(8) :: sq2
    real(8) , allocatable :: zsp(:,:,:), zsp2(:,:,:)
    integer :: levelIndex, lonIndex, latIndex
    real(8), target  :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)

    ! maybe not needed:
    sp(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)

    allocate(zsp(bgStats%nkgdim,2,mymCount))
    allocate(zsp2(bgStats%nkgdim,2,mymCount))

    !$OMP PARALLEL DO PRIVATE(jn,jm,levelIndex,ila_mpiglobal,ila_mpilocal, &
      zsp2,zsp,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if (jm <= jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, bgStats%nkgdim
            zsp(levelIndex,1,icount) = hiControlVector_in(ila_mpilocal,1, &
	      levelIndex)
            zsp(levelIndex,2,icount) = hiControlVector_in(ila_mpilocal,2, &
	      levelIndex)
          end do
        end if
      end do
      if (icount > 0) then

        CALL DGEMM('N','N',bgStats%nkgdim,2*icount,bgStats%nkgdim,1.0d0, &
	  bgStats%corns(1,1,jn),bgStats%nkgdim,zsp(1,1,1), &
	  bgStats%nkgdim,0.0d0,zsp2(1,1,1),bgStats%nkgdim)

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if (jm <= jn) then
            icount = icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do levelIndex = 1, bgStats%nkgdim
              sp(ila_mpilocal,1,levelIndex) = zsp2(levelIndex,1,icount)
              sp(ila_mpilocal,2,levelIndex) = zsp2(levelIndex,2,icount)
            end do
          end if
        end do

      end if

      ! make adjustments for jm=0
      if (mymBeg == 0) then

        ila_mpiglobal = gst_getNind(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do levelIndex = 1, bgStats%nkgdim
          sp(ila_mpilocal,1,levelIndex) = sp(ila_mpilocal,1,levelIndex)*sq2
          sp(ila_mpilocal,2,levelIndex) = 0.0d0
        end do

      end if

    end do
    !$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)


    !$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, bgStats%nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          gd(lonIndex,latIndex,levelIndex) = 0.0d0
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree(sp,gd)
    call gst_setID(gstID2)

    !$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, bgStats%nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          gd(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex) &
	    *bgStats%stddev(lonIndex,latIndex,levelIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, bgStats%nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          gd_out(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine bchm_spa2gd

  !--------------------------------------------------------------------------
  ! bchm_spa2gdad
  !--------------------------------------------------------------------------
  subroutine bchm_spa2gdad(gd_in,hiControlVector_out)
    implicit none

    !Arguments
    real(8), intent(inout) :: hiControlVector_out(nla_mpilocal,2,bgStats%nkgdim)
    real(8), intent(inout) :: gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)

    !Locals
    real(8) :: sp(nla_mpilocal,2,bgStats%nkgdim)

    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, icount
    real(8) :: sq2
    real(8), allocatable :: zsp(:,:,:), zsp2(:,:,:)

    integer :: levelIndex, lonIndex, latIndex
    real(8), target :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)

    !$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, bgStats%nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd                                                      
          gd(lonIndex,latIndex,levelIndex) = gd_in(lonIndex,latIndex,levelIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, bgStats%nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          gd(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex)* &
 	     bgStats%stddev(lonIndex,latIndex,levelIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree_ad(sp,gd)

    hiControlVector_out(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)
    allocate(zsp(bgStats%nkgdim,2,mymCount))
    allocate(zsp2(bgStats%nkgdim,2,mymCount))

    !$OMP PARALLEL DO PRIVATE(JN,JM,levelIndex,ILA_MPILOCAL,ILA_MPIGLOBAL,zsp, &
      zsp2,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if (jm <= jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNind(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, bgStats%nkgdim
            zsp2(levelIndex,1,icount) = sp(ila_mpilocal,1,levelIndex)
            zsp2(levelIndex,2,icount) = sp(ila_mpilocal,2,levelIndex)
          end do
        end if
      end do

      if (icount > 0) then

        CALL DGEMM('T','N',bgStats%nkgdim,2*icount,bgStats%nkgdim,1.0d0, &
	  bgStats%corns(1,1,jn),bgStats%nkgdim,zsp2(1,1,1), &
	  bgStats%nkgdim,0.0d0,zsp(1,1,1),bgStats%nkgdim)

        icount = 0
        do jm = mymBeg, jn, mymSkip
          icount=icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, bgStats%nkgdim
            hiControlVector_out(ila_mpilocal,1,levelIndex) = zsp(levelIndex,1,icount)
            hiControlVector_out(ila_mpilocal,2,levelIndex) = zsp(levelIndex,2,icount)
          end do
        end do

      end if

      ! make adjustments for jm=0
      if (mymBeg == 0) then

        ila_mpiglobal = gst_getNIND(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do levelIndex = 1, bgStats%nkgdim
          hiControlVector_out(ila_mpilocal,1,levelIndex) = &
	    hiControlVector_out(ila_mpilocal,1,levelIndex)*sq2
          hiControlVector_out(ila_mpilocal,2,levelIndex) = &
	    hiControlVector_out(ila_mpilocal,2,levelIndex)*sq2
        end do

      end if

    end do
    !$OMP END PARALLEL DO

    deallocate(zsp)
    deallocate(zsp2)

  end subroutine bchm_spa2gdad

  !--------------------------------------------------------------------------
  ! bchm_copyToStatevector
  !--------------------------------------------------------------------------
  subroutine bchm_copyToStatevector(statevector,gd)
    implicit none
    
    !Arguments
    type(struct_gsv), intent(inout) :: statevector
    real(8), intent(inout) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)
    
    !Locals
    integer :: lonIndex, levelIndex, levelIndex2, latIndex, varIndex, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      call gsv_getField(statevector,field,bgStats%varNameList(varIndex))
      ilev1 = bgStats%nsposit(varIndex)
      ilev2 = bgStats%nsposit(varIndex+1)-1 
        
      !!!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,levelIndex2,lonIndex)
      do levelIndex = ilev1, ilev2
        levelIndex2 = levelIndex-ilev1+1
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            field(lonIndex,latIndex,levelIndex2) = gd(lonIndex,latIndex,levelIndex)
          end do
        end do
      end do
      !!!$OMP END PARALLEL DO
    end do
  end subroutine bchm_copyToStatevector

  !--------------------------------------------------------------------------
  ! bchm_copyFromStatevector
  !--------------------------------------------------------------------------
  subroutine bchm_copyFromStatevector(statevector,gd)
    implicit none

    !Arguments
    type(struct_gsv), intent(inout) :: statevector
    real(8), intent(inout) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,bgStats%nkgdim)
    
    !Locals
    integer :: lonIndex, levelIndex, levelIndex2, latIndex, varIndex, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      call gsv_getField(statevector,field,bgStats%varNameList(varIndex))

      ilev1 = bgStats%nsposit(varIndex)
      ilev2 = bgStats%nsposit(varIndex+1)-1 

      !!!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,levelIndex2,lonIndex)
      do levelIndex = ilev1, ilev2
        levelIndex2 = levelIndex-ilev1+1
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            gd(lonIndex,latIndex,levelIndex) = field(lonIndex,latIndex,levelIndex2)
          end do
        end do
      end do
      !!!$OMP END PARALLEL DO
     end do

  end subroutine bchm_copyFromStatevector

  !--------------------------------------------------------------------------
  ! bchm_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine bchm_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    !Locals
    real(8), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: levelIndex,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    if (.not.initialized) return

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if (mmpi_myid == 0) then
      allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
      allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

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

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare to data to be distributed
    if (mmpi_myid == 0) then

      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        cv_allmaxmpilocal(:,jproc+1) = 0.d0
        jdim_mpilocal = 0

        do levelIndex = 1, bgStats%nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

              if (jm <= jn) then
                      
                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                ! figure out index into global control vector
                if (jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (bgStats%ntrunc+1)
                end if
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * &
		                 (bgStats%ntrunc+1)*(bgStats%ntrunc+1)
                      
                ! index into local control vector computer as in cain
                if (jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                end if
                      
                if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                  write(*,*)
                  write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                  write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                  call utl_abort('bchm_reduceToMPILocal')
                end if
                if (jdim_mpiglobal > cvDim_mpiglobal) then
                  write(*,*)
                  write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                  write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                  call utl_abort('bchm_reduceToMPILocal')
                end if

              end if
            end do
          end do
        end do
 
      end do ! jproc
      !$OMP END PARALLEL DO

    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mmpi_nprocs))
    do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                ! to take the outgoing data to process jproc
    end do

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_double_precision", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_double_precision", &
                           0, "GRID", ierr)

    !- End
    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine bchm_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! bchm_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  subroutine bchm_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    !Locals
    real(4), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: levelIndex,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    if (.not.initialized) return

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if (mmpi_myid == 0) then
      allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
      allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

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

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare to data to be distributed
    if (mmpi_myid == 0) then

      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        cv_allmaxmpilocal(:,jproc+1) = 0.d0
        jdim_mpilocal = 0

        do levelIndex = 1, bgStats%nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

              if (jm <= jn) then
                      
                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                ! figure out index into global control vector
                if  (jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (bgStats%ntrunc+1)
                end if
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * &
		                 (bgStats%ntrunc+1)*(bgStats%ntrunc+1)
                      
                ! index into local control vector computer as in cain
                if (jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                end if
                      
                if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                  write(*,*)
                  write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                  write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                  call utl_abort('bchm_reduceToMPILocal')
                end if
                if (jdim_mpiglobal > cvDim_mpiglobal) then
                  write(*,*)
                  write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                  write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                  call utl_abort('bchm_reduceToMPILocal')
                end if

              end if
            end do
          end do
        end do
 
      end do ! jproc
      !$OMP END PARALLEL DO

    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mmpi_nprocs))
    do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                ! to take the outgoing data to process jproc
    end do

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_real4", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_real4", &
                           0, "GRID", ierr)

    !- End
    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine bchm_reduceToMPILocal_r4

  !--------------------------------------------------------------------------
  ! bchm_expandToMPIGlobal
  !--------------------------------------------------------------------------
  subroutine bchm_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    !Locals
    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: levelIndex, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    if (.not.initialized) return

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if (mmpi_myid == 0) then
      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
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

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if (mmpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do levelIndex = 1, bgStats%nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if (jm <= jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if (jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (bgStats%ntrunc+1)
                end if
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * &
		                 (bgStats%ntrunc+1)*(bgStats%ntrunc+1)

                ! index into local control vector
                if (jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                end if

                if (jdim_mpiglobal > cvDim_mpiglobal) then
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=', &
		             jdim_mpiglobal,cvDim_mpiglobal,levelIndex,jn,jm
		end if

              end if
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
    deallocate(cv_allmaxmpilocal)

  end subroutine bchm_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! bchm_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  subroutine bchm_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    !Locals
    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: levelIndex, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    if (.not.initialized) return

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if (mmpi_myid == 0) then
      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
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

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if (mmpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do levelIndex = 1, bgStats%nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if (jm <= jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if (jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (bgStats%ntrunc+1)
                end if
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * &
		                 (bgStats%ntrunc+1)*(bgStats%ntrunc+1)

                ! index into local control vector
                if (jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                end if

                if (jdim_mpiglobal > cvDim_mpiglobal) then
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=', &
		             jdim_mpiglobal,cvDim_mpiglobal,levelIndex,jn,jm
                end if 
              end if
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
    deallocate(cv_allmaxmpilocal)

  end subroutine bchm_expandtompiglobal_r4

  !--------------------------------------------------------------------------
  ! bchm_finalize
  !--------------------------------------------------------------------------
  subroutine bchm_finalize()
    implicit none

    if (initialized) call bcsc_finalize()

  end subroutine bchm_finalize

end module BmatrixChem_mod
