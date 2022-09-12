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

module minimization_mod
  ! MODULE minimization_mod (prefix='min' category='1. High-level functionality')
  !
  ! :Purpose: Minimization for variational assimilation, including the
  !           subroutine that evaluates the cost function and its gradient.
  !
  use codePrecision_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use verticalCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use midasMpi_mod
  use horizontalCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bmatrix_mod
  use bMatrix1DVar_mod
  use bmatrixhi_mod
  use bmatrixchem_mod
  use bmatrixEnsemble_mod
  use stateToColumn_mod
  use varNameList_mod
  use varqc_mod
  use randomNumber_mod
  use rmatrix_mod
  use costFunction_mod
  use residual_mod
  use obsOperators_mod
  use innovation_mod
  use quasinewton_mod
  use utilities_mod
  use biasCorrectionSat_mod
  use columnVariableTransforms_mod
  implicit none
  save
  private

  ! public variables
  public              :: min_niter, min_nsim

  ! public procedures
  public              :: min_Setup, min_minimize, min_writeHessian

  type(struct_obs)       , pointer :: obsSpaceData_ptr         => null()
  type(struct_columnData), pointer :: columnAnlInc_ptr         => null()
  type(struct_columnData), pointer :: columnTrlOnAnlIncLev_ptr => null()
  type(struct_hco)       , pointer :: hco_anl                  => null()

  logical             :: initialized = .false.

  integer             :: nmtra,nwork,min_nsim
  integer             :: nvadim_mpilocal ! for mpi
  integer             :: min_niter
  integer,external    :: get_max_rss
  logical             :: preconFileExists
  character(len=20)   :: preconFileName    = './preconin'  
  character(len=20)   :: preconFileNameOut = './pm1q'  
  character(len=20)   :: preconFileNameOut_pert = './pm1q_pert'  
  integer             :: n1gc = 3

  ! variables stored for later call to min_writeHessian
  real(8), allocatable :: vatra(:)
  real(8), pointer     :: controlVectorIncrSum_ptr(:)
  real(8), allocatable :: controlVectorIncrSumZero(:)
  real(8) :: zeps0, zdf1
  integer :: itertot, isimtot, iztrl(5), imode
  integer :: outerLoopIndex
  logical :: llvazx
  logical :: initializeForOuterLoop
  logical :: deallocHessian
  logical :: isMinimizationFinalCall

  ! namelist variables
  real(8) :: REPSG, rdf1fac
  integer :: numIterMaxInnerLoopUsed
  integer :: NVAMAJ, NITERMAX, NSIMMAX, nwoqcv
  integer :: numAnalyses
  logical :: lxbar, lwrthess, lgrtest, lvazx
  logical :: lvarqc, writeAnalysis
  logical :: oneDVarMode
  character(len=256) :: ensPathName

  NAMELIST /NAMMIN/ NVAMAJ, NITERMAX, NSIMMAX
  NAMELIST /NAMMIN/ LGRTEST
  NAMELIST /NAMMIN/ lxbar, lwrthess, lvazx
  NAMELIST /NAMMIN/ REPSG, rdf1fac
  NAMELIST /NAMMIN/ LVARQC, NWOQCV
  NAMELIST /NAMMIN/ numAnalyses, ensPathName

CONTAINS

  subroutine min_setup( nvadim_mpilocal_in, hco_anl_in, oneDVarMode_opt, &
                        varqc_opt, nwoqcv_opt )
    !
    ! :Purpose: Reading NAMMIN namelist to setup variables for minimization.
    !
    implicit none

    ! Arguments:
    integer, intent(in)                   :: nvadim_mpilocal_in
    type(struct_hco), pointer, intent(in) :: hco_anl_in
    logical, intent(in), optional         :: oneDVarMode_opt
    logical, intent(out), optional        :: varqc_opt
    integer, intent(out), optional        :: nwoqcv_opt

    ! Locals:
    integer :: ierr,nulnam
    integer,external :: fnom,fclos

    call utl_tmg_start(90,'--Minimization')

    if ( nvadim_mpilocal_in /= cvm_nvadim ) then
      write(*,*) 'nvadim_mpilocal_in,cvm_nvadim=',nvadim_mpilocal_in,cvm_nvadim
      call utl_abort('min_setup: control vector dimension not consistent')
    endif

    nvadim_mpilocal = nvadim_mpilocal_in

    if ( present(oneDVarMode_opt) ) then
      oneDVarMode = oneDVarMode_opt
    else
      oneDVarMode = .false.
    end if

    hco_anl => hco_anl_in

    ! set default values for namelist variables
    nvamaj = 6
    nitermax = 0
    rdf1fac  = 0.25d0
    nsimmax  = 500
    lgrtest  = .false.
    lwrthess = .true.
    lxbar    = .false.
    lvazx    = .false.
    repsg    = 1.0d-5
    lvarqc   = .false.
    nwoqcv   = 5
    numAnalyses = 20
    ensPathName = './ensemble'

    ! read in the namelist NAMMIN
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nammin,iostat=ierr)
    if(ierr.ne.0) call utl_abort('min_setup: Error reading namelist')
    write(*,nml=nammin)
    ierr=fclos(nulnam)

    IF(N1GC == 3)THEN
      NMTRA = (4 + 2*NVAMAJ)*nvadim_mpilocal
    ELSE
      call utl_abort('min_setup: only N1GC=3 currently supported!') 
    END IF
    WRITE(*,9401)N1GC,NVAMAJ,NMTRA
 9401 FORMAT(4X,'N1GC = ',I2,4X,'NVAMAJ = ',I3,/5X,"NMTRA =",1X,I14)

    if ( present(varqc_opt) ) varqc_opt = lvarqc
    if ( present(nwoqcv_opt) ) nwoqcv_opt = nwoqcv

    if(LVARQC .and. mmpi_myid == 0) write(*,*) 'VARIATIONAL QUALITY CONTROL ACTIVATED.'

    initialized=.true.

    call utl_tmg_stop(90)

  end subroutine min_setup


  subroutine min_minimize( outerLoopIndex_in, columnTrlOnAnlIncLev, obsSpaceData, controlVectorIncrSum, &
                           vazx, numIterMaxInnerLoop, deallocHessian_opt, &
                           isMinimizationFinalCall_opt, numIterMaxInnerLoopUsed_opt )
    !
    ! :Purpose: Minimizing cost function to get the increments.
    !           The maximum number of inner-loop iterations is set to nitermax if the
    !           namelist variable nitermax is provided. Otherwise, it is set to the
    !           numIterMaxInnerLoop supplied by the calling subroutine/program.
    !           numIterMaxInnerLoopUsed_opt is passing the maximum number of inner-loop
    !           iterations to the calling subroutine/program.
    !
    implicit none

    ! Arguments:
    integer, intent(in)                    :: outerLoopIndex_in
    type(struct_columnData), intent(inout) :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    real(8), intent(inout)        , target :: controlVectorIncrSum(:)
    real(8), intent(inout)                 :: vazx(:)
    integer, intent(in)                    :: numIterMaxInnerLoop
    integer, intent(out),         optional :: numIterMaxInnerLoopUsed_opt
    logical, intent(in),          optional :: deallocHessian_opt
    logical, intent(in),          optional :: isMinimizationFinalCall_opt

    ! Locals:
    type(struct_columnData) :: columnAnlInc
    integer                 :: get_max_rss

    call utl_tmg_start(90,'--Minimization')

    write(*,*) '--------------------------------'
    write(*,*) '--Starting subroutine minimize--'
    write(*,*) '--------------------------------'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    initializeForOuterLoop = .true.
    outerLoopIndex = outerLoopIndex_in

    if ( present(deallocHessian_opt) ) then
      deallocHessian = deallocHessian_opt
    else
      deallocHessian = .true.
    end if
    if ( present(isMinimizationFinalCall_opt) ) then
      isMinimizationFinalCall = isMinimizationFinalCall_opt
    else
      isMinimizationFinalCall = .true.
    end if

    if ( (nitermax > 0 .and. numIterMaxInnerLoop > 0) .or. &
         (nitermax == 0 .and. numIterMaxInnerLoop == 0) ) then
      call utl_abort('min_minimize: one of nitermax or numIterMaxInnerLoop should be zero and the other positive')
    end if

    if ( nitermax > 0 ) then
      numIterMaxInnerLoopUsed = nitermax
    else if ( numIterMaxInnerLoop > 0 ) then
      numIterMaxInnerLoopUsed = numIterMaxInnerLoop
    else
      call utl_abort('min_minimize: one of the variables nIterMax and numIterMaxInnerLoop must be positive')
    end if
    if ( present(numIterMaxInnerLoopUsed_opt) ) numIterMaxInnerLoopUsed_opt = numIterMaxInnerLoopUsed

    controlVectorIncrSum_ptr => controlVectorIncrSum

    ! zero array for writing to hessian
    allocate(controlVectorIncrSumZero(cvm_nvadim))
    controlVectorIncrSumZero(:) = 0.0d0

    call col_setVco(columnAnlInc,col_getVco(columnTrlOnAnlIncLev))
    call col_allocate(columnAnlInc,col_getNumCol(columnTrlOnAnlIncLev),mpiLocal_opt=.true.)

    write(*,*) 'oti_timeBinning: For 4D increment'
    call oti_timeBinning(obsSpaceData,tim_nstepobsinc)

    call quasiNewtonMinimization( columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData, vazx )

    call col_deallocate(columnAnlInc)

    deallocate(controlVectorIncrSumZero)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) '--Done subroutine minimize--'

    call utl_tmg_stop(90)

  end subroutine min_minimize


  subroutine quasiNewtonMinimization( columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData, vazx )
      !
      !:Purpose: 3D/En-VAR minimization
      !
      implicit none

      ! Arguments:
      type(struct_columnData), target :: columnAnlInc
      type(struct_columnData), target :: columnTrlOnAnlIncLev
      type(struct_obs),        target :: obsSpaceData
      real(8)                         :: vazx(:)

      ! Locals:
      integer              :: nulout = 6
      integer              :: impres
      INTEGER              :: NGRANGE = 10 ! range of powers of 10 used for gradient test

      real    :: zzsunused(1)
      integer :: intUnused(1)

      real(8),allocatable :: vazg(:)

      real(8) :: dlds(1)
      logical :: llvarqc, lrdvatra, llxbar

      integer :: itermax, iterdone, itermaxtodo, isimmax, indic, iitnovqc
      integer :: ierr, isimdone, jdata, isimnovqc
      integer :: ibrpstamp, isim3d
      real(8) :: zjsp, zxmin, zeps1
      real(8) :: dlgnorm, dlxnorm

      real(8) :: zeps0_000,zdf1_000
      integer :: iterdone_000,isimdone_000

      if (lvarqc .and. outerLoopIndex==1) call vqc_setup(obsSpaceData)

      min_nsim=0 

      if(mmpi_myid == 0) then
        impres=5
      else 
        impres=0
      endif

      ! Check for preconditioning file
      inquire(file=preconFileName,exist=preconFileExists)
      if(preconFileExists) then
        if(mmpi_myid == 0) write(*,*) 'PRECONDITIONING FILE FOUND:',preconFileName
      else
        if(mmpi_myid == 0) write(*,*) 'NO PRECONDITIONING FILE FOUND:',preconFileName
      endif

      ! Initialize Hessian only at first outerLoop (mpilocal)
      if ( outerLoopIndex == 1 ) then
        allocate(vatra(nmtra))
        vatra(:) = 0.0d0
      end if

      allocate(vazg(nvadim_mpilocal),stat=ierr)
      if(ierr.ne.0) then
        write(*,*) 'minimization: Problem allocating memory! id=2',ierr
        call utl_abort('min quasiNewtonMinimization')
      endif

      ! set module variable pointers for obsspacedata and the two column objects
      obsSpaceData_ptr         => obsSpaceData
      columnAnlInc_ptr         => columnAnlInc
      columnTrlOnAnlIncLev_ptr => columnTrlOnAnlIncLev

      ! Set-up the minimization

      ! initialize iteration/simulation counters to zero
      ITERTOT = 0
      isimtot = 0

      ! initialize control vector related arrays to zero
      vazx(:) = 0.0d0
      vazg(:) = 0.0d0


      ! If minimization start without qcvar : turn off varqc to compute
      ! innovations and test the gradients

      if (outerLoopIndex == 1) zeps0 = repsg

      if (preconFileExists) then
        if ( mmpi_myid == 0 ) write(*,*) 'quasiNewtonMinimization : Preconditioning mode'
        lrdvatra = .true.
        imode = 2
        llvazx = lvazx ! from namelist (default is .false.)
        llxbar = lxbar ! from namelist (default is .false.)
      else
        lrdvatra = .false.
        imode = 0
      endif

      ! read the hessian from preconin file at first outer-loop iteration
      if ( lrdvatra .and. outerLoopIndex == 1 ) then
        ibrpstamp = tim_getDatestamp() ! ibrpstamp is a I/O argument of hessianIO

        call hessianIO (preconFileName,0,                    &
          isim3d,ibrpstamp,zeps0_000,zdf1_000,iterdone_000,  &
          isimdone_000,iztrl,vatra,controlVectorIncrSumZero, &
          vazx,llxbar,llvazx,n1gc,imode)
      endif

      iterdone = 0
      isimdone = 0
      itermax = numIterMaxInnerLoopUsed
      itermaxtodo = itermax
      isimmax = nsimmax

      ! do the gradient test for the starting point of minimization
      if ( lgrtest .and. outerLoopIndex == 1 ) then
        ! save user-requested varqc switch
        llvarqc = lvarqc

        if ( nwoqcv > 0 ) lvarqc = .false.
        call utl_tmg_start(91,'----QuasiNewton')
        call grtest2(simvar,nvadim_mpilocal,vazx,ngrange)
        call utl_tmg_stop(91)

        lvarqc = llvarqc
      endif

      itertot = iterdone
      isimtot = isimdone

      llvarqc = lvarqc
      if ( nwoqcv > 0 .and. outerLoopIndex == 1 ) lvarqc = .false.
      INDIC =2
      call utl_tmg_start(91,'----QuasiNewton')
      call simvar(indic,nvadim_mpilocal,vazx,zjsp,vazg)
      call utl_tmg_stop(91)
      lvarqc = llvarqc
      
      if ( outerLoopIndex == 1 ) zdf1 = rdf1fac * ABS(zjsp)

      CALL PRSCAL(nvadim_mpilocal,VAZG,VAZG,DLGNORM)
      DLGNORM = DSQRT(DLGNORM)
      CALL PRSCAL(nvadim_mpilocal,VAZX,VAZX,DLXNORM)
      DLXNORM = DSQRT(DLXNORM)
      WRITE(*,*)' |X| = ', DLXNORM
      WRITE(*,FMT=9220) ZJSP, DLGNORM
 9220 FORMAT(/4X,'J(X) = ',G23.16,4X,'|Grad J(X)| = ',G23.16)

      ! Iterations of the minimization algorithm

      ZXMIN = epsilon(ZXMIN)
      WRITE(*,FMT=9320)ZXMIN,ZDF1,ZEPS0,IMPRES,ITERMAX,NSIMMAX

 9320 FORMAT(//,10X,' Minimization QNA_N1QN3 starts ...',/  &
             10x,'DXMIN =',G23.16,2X,'DF1 =',G23.16,2X,'EPSG =',G23.16  &
             /,10X,'IMPRES =',I3,2X,'NITER = ',I3,2X,'NSIM = ',I3)

      ! Begin the minimization
      if ( numIterMaxInnerLoopUsed > 0 ) then

        ! First do iterations without var-QC only at the beginning of first outer-loop.
        if (lvarqc .and. nwoqcv > 0 .and. iterdone < nwoqcv .and. outerLoopIndex == 1 ) then
          iitnovqc = min(nwoqcv - iterdone,itermax)
          isimnovqc = isimmax
          lvarqc = .false.
          call utl_tmg_start(91,'----QuasiNewton')

          zeps1 = zeps0

          call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, vazx,  &
              zjsp,vazg, zxmin, zdf1, zeps1, impres, nulout, imode,       &
              iitnovqc, isimnovqc ,iztrl, vatra, nmtra, intUnused,   &
              zzsunused, dlds)
          call utl_tmg_stop(91)
          call fool_optimizer(obsSpaceData)

          isimnovqc = isimnovqc - 1
          itermaxtodo = itermaxtodo - iitnovqc + 1
          isimmax = isimmax - isimnovqc + 1

          itertot = itertot + iitnovqc
          isimtot = isimtot + isimnovqc

          zeps0 = zeps0/zeps1
          lvarqc = .true.

          if ((imode == 4 .or. imode == 1) .and. itertot < itermax) then
            imode = 2
            INDIC = 2
            call utl_tmg_start(91,'----QuasiNewton')
            call simvar(indic,nvadim_mpilocal,vazx,zjsp,vazg)
            call utl_tmg_stop(91)
          else
            write(*,*) 'minimization_mod: qna_n1qn3 imode = ', imode
            call utl_abort('minimization_mod: qna_n1qn3 mode not equal to 1 or 4')
          endif
        endif

        ! Now do main minimization with var-QC
        call utl_tmg_start(91,'----QuasiNewton')

        if ( outerLoopIndex > 1 ) imode = 2

        zeps1 = zeps0

        call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, vazx,  &
            zjsp,vazg, zxmin, zdf1, zeps1, impres, nulout, imode,   &
            itermaxtodo,isimmax, iztrl, vatra, nmtra, intUnused, zzsunused,   &
            dlds)
        call utl_tmg_stop(91)
        call fool_optimizer(obsSpaceData)

        itertot = itertot + itermaxtodo
        isimtot = isimtot + isimmax

        zeps0 = zeps0/zeps1


        WRITE(*,FMT=9500) imode,iterdone,itertot-iterdone,itertot,isimdone,isimtot-isimdone,isimtot
 9500   FORMAT(//,20X,20('*'),2X    &
        ,/,20X,'              Minimization ended with MODE:',I4  &
        ,/,20X,' Number of iterations done in previous job:',I4  &
        ,/,20X,'          Number of iterations in this job:',I4  &
        ,/,20X,'                Total number of iterations:',I4  &
        ,/,20X,'Number of simulations done in previous job:',I4  &
        ,/,20X,'         Number of simulations in this job:',I4  &
        ,/,20X,'               Total number of simulations:',I4)

        min_niter = itertot

        ! Test the gradient at the final point
        if ( lgrtest .and. isMinimizationFinalCall ) then
          WRITE(*,FMT=9400)
 9400     FORMAT(//,12X,40('**'),/,12X,'TESTING THE GRADIENT AT THE FINAL POINT',/,40('**'))
          call utl_tmg_start(91,'----QuasiNewton')
          call grtest2(simvar,nvadim_mpilocal,vazx,ngrange)
          call utl_tmg_stop(91)
        end if

        ! Print some contents of obsSpaceData to the listing
        if(mmpi_myid == 0) then
          do jdata = 1, min(1,obs_numHeader(obsSpaceData))
            call obs_prnthdr(obsSpaceData,jdata)
            call obs_prntbdy(obsSpaceData,jdata)
          end do
        endif

      else

        ! no analysis performed for ensemble mean, ensure mean increment is zero
        vazx(:) = 0.0d0
        min_niter = 0

      endif ! if numIterMaxInnerLoopUsed > 0

      ! deallocate the gradient
      deallocate(vazg)
      if ( deallocHessian .and. .not. lwrthess ) deallocate(vatra)

  end subroutine quasiNewtonMinimization


  subroutine min_writeHessian(vazx)
    implicit none

    real(8) :: vazx(:)

    call utl_tmg_start(90,'--Minimization')

    if ( lwrthess ) then
      ! Write out the Hessian to file
      if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! zero array for writing to hessian
      if ( .not. allocated(controlVectorIncrSumZero) ) then
        allocate(controlVectorIncrSumZero(cvm_nvadim))
      end if
      controlVectorIncrSumZero(:) = 0.0d0

      call hessianIO (preconFileNameOut,1,  &
        min_nsim,tim_getDatestamp(),zeps0,zdf1,itertot,isimtot,  &
        iztrl,vatra,controlVectorIncrSumZero,vazx,.true.,llvazx,n1gc,imode)

      deallocate(controlVectorIncrSumZero)

      if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    endif

    if ( lwrthess ) then
      deallocate(vatra)
    end if

    call utl_tmg_stop(90)

  end subroutine min_writeHessian


  subroutine simvar(na_indic,na_dim,da_v,da_J,da_gradJ)
    !
    !:Purpose: To implement the Variational solver as described in
    !          Courtier, 1997, Dual Formulation of Four-Dimentional Variational
    !          Assimilation, Q.J.R., pp2449-2461.
    !
    !:Arguments:
    !    :na_indic:
    !               =1 No action taken
    !
    !               =2 Same as 4 (compute J and gradJ) but do not interrupt
    !                  timer of the minimizer.
    !
    !               =3 Compute Jo and gradJo only.
    !
    !               =4 Both J(u) and its gradient are computed.
    !
    !               .. Note:: 1 and 4 are reserved values for call back from
    !                         m1qn3. For direct calls, use a value other than 1
    !                         and 4.
    implicit none

    ! Arguments:
    integer :: na_indic
    integer :: na_dim! Control-vector dimension, forecast-error covariance space
    real(8) :: da_v(na_dim) ! Control variable, forecast-error covariance space
    real*8  :: da_J ! Cost function of the Variational algorithm
    real(8) :: da_gradJ(na_dim) ! Gradient of the Variational Cost funtion

    ! Locals:
    real*8, dimension(na_dim) :: dl_v
    real*8 :: dl_Jb, dl_Jo
    type(struct_gsv), save :: statevector
    type(struct_vco), pointer :: vco_anl

    call utl_tmg_stop(91)
    call utl_tmg_stop(90)

    if (na_indic .ne. 1) then ! No action taken if na_indic == 1
       min_nsim = min_nsim + 1

       if(mmpi_myid == 0) then
         write(*,*) 'Entering simvar for simulation ',min_nsim
         write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
       endif

       ! note: controlVectorIncrSum_ptr is sum of previous outer-loops
       dl_v(1:nvadim_mpilocal) = da_v(1:nvadim_mpilocal) + controlVectorIncrSum_ptr(1:nvadim_mpilocal)     

       ! Computation of background term of cost function:
       dl_Jb = dot_product(dl_v(1:nvadim_mpilocal),dl_v(1:nvadim_mpilocal))/2.d0  
       call mmpi_allreduce_sumreal8scalar(dl_Jb,"GRID")

       if (oneDVarMode) then
         call bmat1D_sqrtB(da_v, nvadim_mpilocal, columnAnlInc_ptr, obsSpaceData_ptr)
         call cvt_transform(columnAnlInc_ptr, 'ZandP_tl', columnTrlOnAnlIncLev_ptr)
       else
         if (.not.gsv_isAllocated(statevector)) then
           write(*,*) 'min-simvar: allocating increment stateVector'
           vco_anl => col_getVco(columnTrlOnAnlIncLev_ptr)
           call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                dataKind_opt=pre_incrReal, mpi_local_opt=.true.)
           call gio_readMaskFromFile(statevector,'./analysisgrid')
         end if

         call bmat_sqrtB(da_v,nvadim_mpilocal,statevector)

         ! put in columnAnlInc H_horiz dx
         call s2c_tl(statevector,columnAnlInc_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)
       end if

       ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
       call utl_tmg_start(10,'--Observations')
       call utl_tmg_start(18,'----ObsOper_TL')
       call oop_Htl(columnAnlInc_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr,min_nsim, &
                    initializeLinearization_opt=initializeForOuterLoop) 
       call utl_tmg_stop(18)
       call utl_tmg_stop(10)

       ! Calculate OBS_OMA from OBS_WORK : d-Hdx
       call res_compute(obsSpaceData_ptr)  

       call bcs_calcbias_tl(da_v,OBS_OMA,obsSpaceData_ptr,columnTrlOnAnlIncLev_ptr)

       ! Save as OBS_WORK : R**-1/2 (d-Hdx)
       call utl_tmg_start(10,'--Observations')
       call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_OMA)  
       call utl_tmg_stop(10)

       ! Store J-obs in OBS_JOBS : 1/2 * R**-1 (d-Hdx)**2
       call cfn_calcJo(obsSpaceData_ptr) 

       ! Store modified J_obs in OBS_JOBS : -ln((gamma-exp(J))/(gamma+1)) 
       IF ( LVARQC ) THEN
         call vqc_NlTl(obsSpaceData_ptr)
       endif

       dl_Jo = 0.d0
       call utl_tmg_start(90,'--Minimization')
       call utl_tmg_start(92,'----SumCostFunction')
       call cfn_sumJo(obsSpaceData_ptr,dl_Jo)
       da_J = dl_Jb + dl_Jo
       if (na_indic  ==  3) then
          da_J = dl_Jo
          IF(mmpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  JO = ",G23.16,6X)') dl_Jo
       else
          da_J = dl_Jb + dl_Jo
          IF(mmpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  Jb = ",G23.16,6X,"JO = ",G23.16,6X,"Jt = ",G23.16)') dl_Jb,dl_Jo,da_J
       endif
       call utl_tmg_stop(92)
       call utl_tmg_stop(90)

       ! Modify OBS_WORK : R**-1 (d-Hdx)
       call utl_tmg_start(10,'--Observations')
       call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  
       call utl_tmg_stop(10)

       IF ( LVARQC ) THEN
         call vqc_ad(obsSpaceData_ptr)
       endif

       ! Calculate adjoint of d-Hdx (mult OBS_WORK by -1)
       call res_computeAd(obsSpaceData_ptr)  

       call col_zero(columnAnlInc_ptr)

       ! Put in column : -H_vert**T R**-1 (d-Hdx)
       call utl_tmg_start(10,'--Observations')
       call utl_tmg_start(19,'----ObsOper_AD')
       call oop_Had(columnAnlInc_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr, &
                    initializeLinearization_opt=initializeForOuterLoop)
       call utl_tmg_stop(19)
       call utl_tmg_stop(10)

       ! Put in statevector -H_horiz**T H_vert**T R**-1 (d-Hdx)
       if (oneDVarMode) then
         ! no interpolation needed for 1Dvar case
       else
         call s2c_ad(statevector,columnAnlInc_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)
       end if

       da_gradJ(:) = 0.d0
       call bcs_calcbias_ad(da_gradJ,OBS_WORK,obsSpaceData_ptr)

       if (oneDVarMOde) then
         call cvt_transform( columnAnlInc_ptr, 'ZandP_ad', columnTrlOnAnlIncLev_ptr)      ! IN
         call bmat1D_sqrtBT(da_gradJ, nvadim_mpilocal, columnAnlInc_ptr, obsSpaceData_ptr)
       else
         call bmat_sqrtBT(da_gradJ,nvadim_mpilocal,statevector)
       end if

       if (na_indic .ne. 3) then
         da_gradJ(1:nvadim_mpilocal) = dl_v(1:nvadim_mpilocal) + da_gradJ(1:nvadim_mpilocal)
       endif
    endif

    call utl_tmg_start(90,'--Minimization')
    call utl_tmg_start(91,'----QuasiNewton')

    initializeForOuterLoop = .false.

    if(mmpi_myid == 0) write(*,*) 'end of simvar'

  end subroutine simvar


  subroutine dscalqn(kdim,px,py,ddsc)
    !:Purpose: Interface for the inner product to be used by the minimization
    !          subroutines QNA_N1QN3.
    !
    !:Arguments:
    !   i : kdim
    !   i : px,py
    !   o : ddsc
    !   i : kzs
    !   i : pzs
    !   i : ddzs
    implicit none

    ! Arguments:
    integer:: kdim     ! dimension of the vectors
    real*8 :: px(kdim) ! vector for which <PX,PY> is being calculated
    real*8 :: py(kdim) ! vector for which <PX,PY> is being calculated
    real*8 :: ddsc     ! result of the inner product

    CALL PRSCAL(KDIM,PX,PY,DDSC)
    RETURN
  end subroutine dscalqn


  subroutine prscal(kdim,px,py,ddsc)
    !
    !:Purpose: To evaluate the inner product used in the minimization
    !
    !:Arguments:
    !    i : KDIM
    !    i : PX, PY
    !    o : DDSC
    !
    implicit none

    ! Arguments:
    integer :: kdim     ! dimension of the vectors
    real*8  :: px(kdim) ! vector for which <PX,PY> is being calculated
    real*8  :: py(kdim) ! vector for which <PX,PY> is being calculated
    real*8  :: ddsc     ! result of the inner product

    ! Locals:
    INTEGER J

    DDSC = 0.D0

    do j=1,nvadim_mpilocal
      DDSC = DDSC + PX(J)*PY(J)
    ENDDO

    call mmpi_allreduce_sumreal8scalar(ddsc,"GRID")

  end subroutine prscal


  subroutine dcanab(kdim,py,px)
    !
    !:Purpose: Change of variable associated with the canonical inner product:
    !
    !    * PX = L^-1 * Py with L related to the inner product
    !    * <PX,PY> = PX^t  L^t  L PY
    !    * (see the modulopt documentation aboutn DTCAB)
    !    * NOTE: L is assumed to be the identity!
    !
    implicit none

    ! Arguments:
    integer kdim
    real*8 px(kdim), py(kdim)

    ! Locals
    INTEGER JDIM

    DO JDIM = 1, KDIM
      PX(JDIM) = PY(JDIM)
    ENDDO

    RETURN
  end subroutine DCANAB


  subroutine dcanonb(kdim,px,py)
    !
    !:Purpose: Change of variable associated with the canonical inner product:
    !
    !    * PY = L * PX with L related to the inner product
    !    * <PX,PY> = PX^t  L^t  L PY
    !    * (see the modulopt documentation about DTONB)
    !
    implicit none

    ! Arguments:
    integer kdim
    real*8 px(kdim), py(kdim)

    ! Locals:
    INTEGER JDIM

    DO JDIM = 1, KDIM
      PY(JDIM) = PX(JDIM)
    ENDDO

    RETURN
  end subroutine DCANONB


  subroutine hessianIO (cfname,status,                             &
                        nsim,kbrpstamp,zeps1,zdf1,itertot,isimtot, &
                        iztrl,vatra,vazxbar,vazx,llxbar,llvazx,k1gc,imode)
    !
    !:Purpose: Read-Write Hessian and increment (possibly for outer loop) on a
    !          file
    !
    !:Arguments:
    !   * i   cfname
    !   * i   status
    !   * i   nsim
    !   * io  kbrpstamp
    !   * i   zeps1
    !   * i   zdf1
    !   * i   itertot
    !   * i   isimtot
    !   * i   iztrl
    !   * i   vatra
    !   * i   vazxbar
    !   * i   vazx
    !   * i   llxbar
    !   * i   llvazx
    !   * i   k1gc
    !   * o   imode
    !
    implicit none

    ! Arguments:
    character(len=*) :: cfname ! precon file
    integer status  ! = 0 if READ, = 1 if WRITE
    integer nsim    ! Number of simulations in QNA_N1QN3
    integer kbrpstamp ! Date
    real*8 :: zeps1    ! Parameter in QNA_N1QN3
    real*8 :: zdf1     ! Parameter in QNA_N1QN3
    integer itertot ! Parameter in QNA_N1QN3
    integer isimtot ! Parameter in QNA_N1QN3
    integer, target:: iztrl(5)     ! Localisation parameters for Hessian
    real*8, target :: vatra(nmtra) ! Hessian
    real*8, target :: vazxbar(nvadim_mpilocal) ! Vazx of previous loop
    real*8, target :: vazx(nvadim_mpilocal) ! Current state of the minimization
    logical llxbar  ! read in vaxzbar if dates are compatible
    logical llvazx  ! Logical to read vazx
    integer k1gc    ! Minimizer ID (2: m1qn2, 3: m1qn3)
    integer imode   ! If status=0, set imode=0 (no prec) or 2 (prec)

    ! Locals:
    real*4, allocatable :: vatravec_r4_mpiglobal(:)
    real*4, allocatable :: vatra_r4(:)
    real*8, allocatable :: vazxbar_mpiglobal(:),vazx_mpiglobal(:)

    integer :: ibrpstamp,ireslun, ierr, fnom, fclos
    integer :: nvadim_mpiglobal,nmtra_mpiglobal
    integer :: ivadim, itrunc
    integer :: ivamaj
    integer :: jvec, i1gc,ictrlvec,ii
    integer, dimension(10), target, save :: iztrl_io

    character(len=3) :: cl_version

    if (status == 0) then
      if (mmpi_myid == 0) write(*,*) 'Read  Hessian in min_hessianIO from file ', cfname
      call utl_tmg_start(93,'----ReadHess')
    elseif (status == 1) then
      if (mmpi_myid == 0) write(*,*) 'Write Hessian in min_hessianIO to file ', cfname
      call utl_tmg_start(94,'----WriteHess')
    else
      call utl_abort('min_hessianIO: status not valid ')
    endif

    call rpn_comm_allreduce(nvadim_mpilocal,nvadim_mpiglobal,1,"mpi_integer","mpi_sum","GRID",ierr)
    call rpn_comm_allreduce(nmtra          ,nmtra_mpiglobal, 1,"mpi_integer","mpi_sum","GRID",ierr)

    ireslun=0

    if (status == 0) then
      !
      !- 1.  Read Hessian on processor 0 and distribute the data to the other processors
      !

      if (mmpi_myid == 0) then
         !- Open the Hessian matrix file
         ierr = fnom(ireslun,cfname,'FTN+SEQ+UNF+OLD+R/O',0)

         !- Checking version number
         read (ireslun) cl_version, i1gc
         if (trim(cl_version) == 'V5') then
            write(*,*) 'min_hessianIO: using single precision V5 Hessian'
         else if (trim(cl_version) == 'V4') then
            write(*,*) 'min_hessianIO: using DEPRECIATED single precision V5 Hessian'
         else
            write(*,*)
            write(*,*) 'min_hessianIO : Only V5 Hessian are supported, found ', trim(cl_version)
            call utl_abort('min_hessianIO')
         endif

         if (i1gc == 3 .and. i1gc == k1gc) then
            write(*,*) trim(cl_version),' M1QN3'
         else
            write(*,*) 'Version, n1gc =',trim(cl_version),i1gc
            call utl_abort('min_hessianIO: Inconsistant input hessian')
         endif

         rewind (ireslun)

         read(ireslun) cl_version,i1gc,nsim,ibrpstamp,zeps1,zdf1,itertot,isimtot,ivadim,itrunc
         read(ireslun) ivamaj,iztrl_io
         if ((ivamaj.ne.nvamaj).or.(nvadim_mpiglobal.ne.ivadim)) then
            write(*,*) nvamaj,ivamaj,nvadim_mpiglobal,ivadim
            call utl_abort('min_hessianIO : ERROR, size of V5 Hessian not consistent')
         endif

      end if

      ! ibrpstamp and iztrl_io must be broadcasted
      call rpn_comm_bcast(ibrpstamp,  1, "MPI_INTEGER", 0, "GRID", ierr)
      call rpn_comm_bcast(iztrl_io , 10, "MPI_INTEGER", 0, "GRID", ierr)

      !- Read the Hessian
      if(mmpi_myid == 0) then 
         write(*,*) 'min_hessianIO : reading Hessian'
         allocate(vatravec_r4_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vatravec_r4_mpiglobal(1))
      end if
      allocate(vatra_r4(nvadim_mpilocal))

      if (k1gc == 3) ictrlvec = 2*nvamaj+1
      do jvec = 1, ictrlvec
 
         if (mmpi_myid == 0) then
            read(ireslun) vatravec_r4_mpiglobal
         end if

         call bmat_reduceToMPILocal_r4( vatra_r4,            & ! OUT
                                        vatravec_r4_mpiglobal) ! IN (contains data only on proc 0)
         !$OMP PARALLEL DO PRIVATE(ii)
         do ii = 1, nvadim_mpilocal
           vatra((jvec-1)*nvadim_mpilocal+ii) = real(vatra_r4(ii),8)
         enddo
         !$OMP END PARALLEL DO

      end do

      deallocate(vatravec_r4_mpiglobal)       
      deallocate(vatra_r4)       

      imode = 2
      iztrl(:) = iztrl_io(1:5)
      if(k1gc == 3) then
         iztrl(1) = nvadim_mpilocal
         iztrl(2) = 0
         iztrl(3) = nvamaj
         iztrl(4) = iztrl_io(4)
         iztrl(5) = iztrl_io(5)
      end if

      !- Read VAZXBAR 
      if (ibrpstamp == kbrpstamp .and. llxbar) then
         if (mmpi_myid == 0) then 
            write(*,*) 'min_hessianIO : reading vazxbar'
            allocate(vazxbar_mpiglobal(nvadim_mpiglobal))
            read(ireslun) vazxbar_mpiglobal
         else
            allocate(vazxbar_mpiglobal(1))
         end if
         call bmat_reduceToMPILocal( vazxbar,          & ! OUT
                                     vazxbar_mpiglobal ) ! IN (contains data only on proc 0)
         deallocate(vazxbar_mpiglobal)
      end if

      !- Read VAZX
      if (ibrpstamp == kbrpstamp .and. llvazx) then
         if (mmpi_myid == 0) then 
            write(*,*) 'min_hessianIO : reading vazx'
            allocate(vazx_mpiglobal(nvadim_mpiglobal))
            read(ireslun) vazx_mpiglobal
         else
            allocate(vazx_mpiglobal(1))
         end if
         call bmat_reduceToMPILocal( vazx,           & ! OUT
                                     vazx_mpiglobal )  ! IN (contains data only on proc 0)
         deallocate(vazx_mpiglobal)
      end if

      if (ibrpstamp /= kbrpstamp) then
         kbrpstamp = ibrpstamp
      end if

      if (mmpi_myid == 0) ierr = fclos(ireslun)

    elseif (status == 1) then
      !
      !- 2.  Write Hessian
      !

      !- Open the Hessian matrix file on processor 0 and write some metadata 
      if (mmpi_myid == 0) ierr = fnom(ireslun,cfname, 'FTN+SEQ+UNF' , 0)
      cl_version = 'V5'
      itrunc=0
      iztrl_io(1: 5) = iztrl(:)
      iztrl_io(6:10) = 0 ! dummy values for hessian size legacy purposes
      if(mmpi_myid == 0) write(ireslun) cl_version,k1gc,nsim,kbrpstamp,zeps1,zdf1,itertot,isimtot,nvadim_mpiglobal,itrunc
      if(mmpi_myid == 0) write(ireslun) nvamaj,iztrl_io

      !- Write the Hessian
      if (mmpi_myid == 0) then
         allocate(vatravec_r4_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vatravec_r4_mpiglobal(1))
      endif
      allocate(vatra_r4(nvadim_mpilocal))

      if (k1gc == 3) ictrlvec = 2*nvamaj+1
      do jvec = 1, ictrlvec
        !$OMP PARALLEL DO PRIVATE(ii)
        do ii = 1, nvadim_mpilocal
          vatra_r4(ii) = vatra((jvec-1)*nvadim_mpilocal+ii)
        enddo
        !$OMP END PARALLEL DO
        call bmat_expandToMPIGlobal_r4( vatra_r4,              & ! IN
                                        vatravec_r4_mpiglobal )  ! OUT
        if (mmpi_myid == 0) write(ireslun) vatravec_r4_mpiglobal
      enddo
      deallocate(vatravec_r4_mpiglobal)
      deallocate(vatra_r4)

      !- Write VAZXBAR
      if (mmpi_myid == 0) then
         allocate(vazxbar_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vazxbar_mpiglobal(1)) ! dummy
      end if
      call bmat_expandToMPIGlobal( vazxbar,          & ! IN
                                   vazxbar_mpiglobal ) ! OUT
      if (mmpi_myid == 0) write(ireslun) vazxbar_mpiglobal(1:nvadim_mpiglobal)
      deallocate(vazxbar_mpiglobal)

      !- Write VAZX
      if (mmpi_myid == 0) then
         allocate(vazx_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vazx_mpiglobal(1))
      end if
      call bmat_expandToMPIGlobal( vazx,          & ! IN
                                   vazx_mpiglobal ) ! OUT
      if (mmpi_myid == 0) write(ireslun) vazx_mpiglobal(1:nvadim_mpiglobal)
      deallocate(vazx_mpiglobal)

      !- Close the Hessian matrix file
      if (mmpi_myid == 0) ierr = fclos(ireslun)
    else
      call utl_abort('min_hessianIO: status not valid')
    endif

    if (status == 0) then
      call utl_tmg_stop(93)
    elseif (status == 1) then
      call utl_tmg_stop(94)
    endif

  end subroutine hessianIO


  subroutine grtest2(simul,na_dim,da_x0,na_range)
  !
  !:Purpose: To compare the variation of the functional against what the
  !          gradient gives for small changes in the control variable. This test
  !          should be accurate for values as small as DLALPHA =  SQRT(machine
  !          precision). (see Courtier, 1987)
  !
  !:Arguments:
  !    :na_range:  the test will be carried over values of ALPHA ranging between
  !                10**(-NA_RANGE) < ALPHA < 0.1
  !
  implicit none

  ! Arguments:
  external simul ! simulator: return cost function estimate and its gradient
  integer, intent(in) :: na_dim ! Size of the control vector
  real*8,  intent(in) :: da_x0(na_dim) ! Control vector
  integer, intent(in) :: na_range

  ! Locals:
  integer :: nl_indic, nl_j
  real*8  :: dl_wrk(na_dim),dl_gradj0(na_dim), dl_x(na_dim)
  real*8  :: dl_J0, dl_J, dl_test, dl_start,dl_end
  real*8  :: dl_alpha, dl_gnorm0


  ! 1. Initialize dl_gradj0 at da_x0
  !    ------------------------------------

  nl_indic = 2
  call simul(nl_indic,na_dim,da_x0,dl_j0,dl_gradj0)
  dl_gnorm0 = dot_product(dl_gradj0,dl_gradj0)
  call mmpi_allreduce_sumreal8scalar(dl_gnorm0,"GRID")
  dl_start = 1.d0
  dl_end   = 10.0d0**(-na_range)
  write(*,FMT=9100) dl_start,dl_end, dl_j0, dl_gnorm0
  ! 2. Perform the test
  !    ----------------

  if(dl_gnorm0 == 0.d0)then
     write(*,FMT=9101)
     return
  end if
  write(*,FMT=9200)
  do  nl_j = 1, na_range
     dl_alpha = 10.0d0**(- nl_j)
     dl_x(:) = da_x0(:) - dl_alpha*dl_gradJ0(:)
     call simul(nl_indic,na_dim,dl_x,dl_j,dl_wrk)
     dl_test = (dl_j-dl_j0)/(-dl_alpha * dl_gnorm0)
     write(*,FMT=9201)nl_j, dl_alpha, dl_j, dl_test
  end do

9100 format(//,4X,&
          'GRTEST- The gradient is being tested for',&
          G23.16,' <= ALPHA <= ',G23.16,/,12X,&
          'Initial value of J(X):',1x,G23.16,4x,&
          'Norm of GRAD J(X)**2: ',G23.16)
9101 format(/,4X,'-In GRTEST: gradient vanishes exactly',&
          '. Gradient test cannot be performed at this point')
9200 format(/,4X,'J',8X,'ALPHA',11X,'J(X)',12X,'TEST')

9201 format(2X,'GRTEST: step',2X,I3,4X,G23.16,4X,G23.16,4X,&
          G23.16)
  return
  end subroutine grtest2

end module minimization_mod
