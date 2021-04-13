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
  use MathPhysConstants_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use verticalCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use mpi_mod
  use mpivar_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use gridStateVector_mod
  use bmatrix_mod
  use var1D_mod
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
  public              :: min_niter,min_nsim
  ! public procedures
  public              :: min_Setup, min_minimize, min_writeHessian

  type(struct_obs)       , pointer :: obsSpaceData_ptr => null()
  type(struct_columnData), pointer :: column_ptr       => null()
  type(struct_columnData), pointer :: columng_ptr      => null()

  type(struct_gsv), pointer :: stateVectorRefHU_ptr => null()

  logical             :: initialized = .false.

  integer             :: envar_loop   ! environment variable

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
  real(8), pointer :: dg_vbar(:)
  real(8) :: zeps1, zdf1
  integer :: itertot, isimtot, iztrl(5), imode
  logical :: llvazx

  ! namelist variables
  INTEGER NVAMAJ,NITERMAX,NSIMMAX
  logical :: lxbar,lwrthess,lgrtest,lvazx
  REAL*8 REPSG,rdf1fac
  logical :: lvarqc,pertBhiOnly, writeAnalysis
  integer :: nwoqcv
  integer :: numIterMax_pert, numAnalyses, ntrunc_pert
  character(len=256) :: ensPathName
  real(8) :: e1_scaleFactor, e2_scaleFactor
  integer,parameter   :: maxNumLevels=200
  real(8) :: pertScaleFactor_UV(maxNumLevels)
  logical :: oneDVarMode=.false.

  NAMELIST /NAMMIN/NVAMAJ,NITERMAX,NSIMMAX
  NAMELIST /NAMMIN/LGRTEST
  NAMELIST /NAMMIN/lxbar,lwrthess,lvazx
  NAMELIST /NAMMIN/REPSG,rdf1fac
  NAMELIST /NAMMIN/LVARQC,NWOQCV
  NAMELIST /NAMMIN/numIterMax_pert,numAnalyses,ensPathName
  NAMELIST /NAMMIN/e1_scaleFactor,e2_scaleFactor,pertBhiOnly
  NAMELIST /NAMMIN/pertScaleFactor_UV,ntrunc_pert

CONTAINS

  subroutine min_setup(nvadim_mpilocal_in, oneDVarMode_opt)
    implicit none
    integer,intent(in) :: nvadim_mpilocal_in
    logical,intent(in),optional :: oneDVarMode_opt

    integer :: ierr,nulnam
    integer,external :: fnom,fclos
    character(len=32) :: envVariable
    integer :: length_envVariable, status

    if(nvadim_mpilocal_in.ne.cvm_nvadim) then
      write(*,*) 'nvadim_mpilocal,cvm_nvadim=',nvadim_mpilocal,cvm_nvadim
      call utl_abort('min_setup: control vector dimension not consistent')
    endif

    nvadim_mpilocal=nvadim_mpilocal_in

    if ( present(oneDVarMode_opt) ) oneDVarMode = oneDVarMode_opt

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
    numIterMax_pert = 0
    numAnalyses = 20
    ensPathName = './ensemble'
    e1_scaleFactor = 0.66d0
    e2_scaleFactor = 0.33d0
    pertBhiOnly = .true.
    pertScaleFactor_UV(:) = 1.0d0
    ntrunc_pert = 0

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

    if(LVARQC .and. mpi_myid == 0) write(*,*) 'VARIATIONAL QUALITY CONTROL ACTIVATED.'

    ! Retrieve environment variables related to doing an ensemble of perturbed analyses
    status = 0
    call get_environment_variable('envar_loop',envVariable,length_envVariable,status,.true.)
    if (status.gt.1) then
      write(*,*) 'min_analysisPert: Problem when getting the environment variable envar_loop'
      envar_loop = 1
    elseif (status == 1) then
      envar_loop = 1
    else
      write(*,*) 'min_analysisPert: The environment variable envar_loop has been detected: ',envVariable
      read(envVariable,'(i8)') envar_loop
      write(*,*) 'envar_loop = ',envar_loop
    endif

    initialized=.true.

  end subroutine min_setup

  subroutine min_minimize(columng,obsSpaceData,vazx,stateVectorRef_opt)
    implicit none

    real*8 :: vazx(:)
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: columng
    type(struct_gsv), target, optional :: stateVectorRef_opt

    type(struct_columnData)   :: column
    integer :: get_max_rss

    write(*,*) '--------------------------------'
    write(*,*) '--Starting subroutine minimize--'
    write(*,*) '--------------------------------'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call tmg_start(3,'MIN')

    if ( present(statevectorRef_opt) ) then
      if ( statevectorRef_opt%allocated ) stateVectorRefHU_ptr => stateVectorRef_opt
    end if

    call col_setVco(column,col_getVco(columng))
    call col_allocate(column,col_getNumCol(columng),mpiLocal_opt=.true.)

    write(*,*) 'oti_timeBinning: For 4D increment'
    call oti_timeBinning(obsSpaceData,tim_nstepobsinc)

    call quasiNewtonMinimization(column,columng,obsSpaceData,vazx)

    call col_deallocate(column)
    call tmg_stop(3)

    ! Memory deallocations for non diagonal R matrices for radiances
    call rmat_cleanup()
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) '--Done subroutine minimize--'

  end subroutine min_minimize



  subroutine quasiNewtonMinimization(column,columng,obsSpaceData,vazx)
      !
      !:Purpose: 3D/En-VAR minimization
      implicit none

      ! Arguments:
      type(struct_columnData),target :: column,columng
      type(struct_obs),target :: obsSpaceData
      real*8 :: vazx(:)

      ! Locals:
      integer              :: nulout = 6
      integer              :: impres
      INTEGER              :: NGRANGE = 10 ! range of powers of 10 used for gradient test

      real    :: zzsunused(1)
      integer :: intUnused(1)

      real*8,allocatable :: vazg(:)

      real*8 :: dlds(1)
      logical :: llvarqc, lldf1, lrdvatra, llxbar

      integer :: itermax, iterdone, itermaxtodo, isimmax, indic, iitnovqc
      integer :: ierr, isimdone, jdata, isimnovqc
      integer :: ibrpstamp, isim3d
      real*8 :: zjsp, zxmin, zeps0
      real*8 :: dlgnorm, dlxnorm

      real*8 :: zeps0_000,zdf1_000
      integer :: iterdone_000,isimdone_000

      if (lvarqc) call vqc_setup(obsSpaceData)

      min_nsim=0 

      if(mpi_myid == 0) then
        impres=5
      else 
        impres=0
      endif

      ! Check for preconditioning file
      inquire(file=preconFileName,exist=preconFileExists)
      if(preconFileExists) then
        if(mpi_myid == 0) write(*,*) 'PRECONDITIONING FILE FOUND:',preconFileName
      else
        if(mpi_myid == 0) write(*,*) 'NO PRECONDITIONING FILE FOUND:',preconFileName
      endif

      ! allocate control vector related arrays (these are all mpilocal)

      allocate(dg_vbar(nvadim_mpilocal),stat=ierr)
      if(ierr.ne.0) then
        write(*,*) 'minimization: Problem allocating memory! id=1',ierr
        call utl_abort('min quasiNewtonMinimization')
      endif

      allocate(vazg(nvadim_mpilocal),stat=ierr)
      if(ierr.ne.0) then
        write(*,*) 'minimization: Problem allocating memory! id=2',ierr
        call utl_abort('min quasiNewtonMinimization')
      endif

      allocate(vatra(nmtra),stat=ierr)
      if(ierr.ne.0) then
        write(*,*) 'minimization: Problem allocating memory! id=4',ierr
        call utl_abort('min quasiNewtonMinimization')
      endif

      ! set module variable pointers for obsspacedata and the two column objects
      obsSpaceData_ptr => obsSpaceData
      column_ptr       => column
      columng_ptr      => columng

      ! Set-up the minimization

      ! initialize iteration/simulation counters to zero
      ITERTOT  = 0
      isimtot = 0

      ! initialize control vector related arrays to zero
      vazx(:)=0.0d0
      dg_vbar(:)=0.0d0
      vazg(:)=0.0d0
      vatra(:)=0.0d0

      ! save user-requested varqc switch
      llvarqc = lvarqc

      ! If minimization start without qcvar : turn off varqc to compute
      ! innovations and test the gradients

      lldf1 = .true.
      if (preconFileExists) then
        if(mpi_myid == 0) write(*,*) 'quasiNewtonMinimization : Preconditioning mode'
        lrdvatra = .true.
        imode = 2
        llvazx = lvazx ! from namelist (default is .false.)
        llxbar = lxbar ! from namelist (default is .false.)
      else
        lrdvatra = .false.
        imode = 0
        zeps0 = repsg
      endif

      ! read the hessian from preconin file
      if (lrdvatra) then
        ibrpstamp = tim_getDatestamp() ! ibrpstamp is a I/O argument of hessianIO

        call hessianIO (preconFileName,0,                     &
          isim3d,ibrpstamp,zeps0_000,zdf1_000,iterdone_000,  &
          isimdone_000,iztrl,vatra,dg_vbar,                  &
          vazx,llxbar,llvazx,n1gc,imode)

        if (ibrpstamp == tim_getDatestamp() .and. lxbar) then
          ! use vbar for true outer loop
          zeps0  = zeps0_000
          zdf1   = zdf1_000
          lldf1 = .false.     ! don't re-compute df1 base on Cost function
        else
          zeps0 = repsg
          lldf1 = .true.      ! Compute df1 base on Cost function
        endif
      endif

      iterdone = 0
      isimdone = 0
      itermax = nitermax
      itermaxtodo = itermax
      isimmax = nsimmax

      if (nwoqcv > 0) lvarqc = .false.

      ! do the gradient test for the starting point of minimization
      if(lgrtest) then
        call grtest2(simvar,nvadim_mpilocal,vazx,ngrange)
      endif

      zeps1 = zeps0

      itertot = iterdone
      isimtot = isimdone

      INDIC =2
      call simvar(indic,nvadim_mpilocal,vazx,zjsp,vazg)

      if (lldf1) ZDF1     =  rdf1fac * ABS(ZJSP)

!     Put QCVAR logical to its original values
      lvarqc=llvarqc

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
      if(nitermax.gt.0) then

        ! First do iterations without var-QC
        if (lvarqc .and. nwoqcv > 0 .and. iterdone < nwoqcv) then
          iitnovqc = min(nwoqcv - iterdone,itermax)
          isimnovqc = isimmax
          lvarqc = .false.
          call tmg_start(70,'QN')
          call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, vazx,  &
              zjsp,vazg, zxmin, zdf1, zeps1, impres, nulout, imode,       &
              iitnovqc, isimnovqc ,iztrl, vatra, nmtra, intUnused,   &
              zzsunused, dlds)
          call tmg_stop(70)
          call fool_optimizer(obsSpaceData)

          isimnovqc = isimnovqc - 1
          itermaxtodo = itermaxtodo - iitnovqc + 1
          isimmax = isimmax - isimnovqc + 1

          itertot = itertot + iitnovqc
          isimtot = isimtot + isimnovqc

          zeps1 = zeps0/zeps1
          zeps0 = zeps1
          lvarqc = .true.

          if ((imode == 4 .or. imode == 1) .and. itertot < itermax) then
            imode = 2
            INDIC = 2
            call simvar(indic,nvadim_mpilocal,vazx,zjsp,vazg)
          else
            write(*,*) 'minimization_mod: qna_n1qn3 imode = ', imode
            call utl_abort('minimization_mod: qna_n1qn3 mode not equal to 1 or 4')
          endif
        endif

        ! Now do main minimization with var-QC
        call tmg_start(70,'QN')
        call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, vazx,  &
            zjsp,vazg, zxmin, zdf1, zeps1, impres, nulout, imode,   &
            itermaxtodo,isimmax, iztrl, vatra, nmtra, intUnused, zzsunused,   &
            dlds)
        call tmg_stop(70)
        call fool_optimizer(obsSpaceData)

        itertot = itertot + itermaxtodo
        isimtot = isimtot + isimmax

        zeps1 = zeps0/zeps1


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
        if (lgrtest) then
          WRITE(*,FMT=9400)
 9400     FORMAT(//,12X,40('**'),/,12X,'TESTING THE GRADIENT AT THE FINAL POINT',/,40('**'))
          call grtest2(simvar,nvadim_mpilocal,vazx,ngrange)
        end if

        do jdata = 1, nvadim_mpilocal
          dg_vbar(jdata) = vazx(jdata) + dg_vbar(jdata)
        enddo

        ! Print some contents of obsSpaceData to the listing
        if(mpi_myid == 0) then
          do jdata = 1, min(1,obs_numHeader(obsSpaceData))
            call obs_prnthdr(obsSpaceData,jdata)
            call obs_prntbdy(obsSpaceData,jdata)
          end do
        endif

      else

        ! no analysis performed for ensemble mean, ensure mean increment is zero
        vazx(:) = 0.0d0
        min_niter = 0

      endif ! if nitermax .gt. 0


      ! If requested, compute analysis increment for ensemble perturbations
      if(numIterMax_pert.gt.0) then
        dg_vbar(:) = 0.0d0
        call tmg_start(4,'MINPERT')
        call min_analysisPert(vatra,iztrl,zdf1,column,columng,obsSpaceData)
        call tmg_stop(4)
      endif

      ! Set the QC flags to be consistent with VAR-QC if control analysis
      if(lvarqc) call vqc_listrej(obsSpaceData)

      ! deallocate the gradient
      deallocate(vazg)
      if ( .not. lwrthess ) then
        deallocate(vatra)
        deallocate(dg_vbar)
      end if

  end subroutine quasiNewtonMinimization


  subroutine min_writeHessian(vazx)
    implicit none

    real(8) :: vazx(:)

    if ( nitermax > 0 .and. lwrthess) then
      ! Write out the Hessian to file
      if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      call hessianIO (preconFileNameOut,1,  &
        min_nsim,tim_getDatestamp(),zeps1,zdf1,itertot,isimtot,  &
        iztrl,vatra,dg_vbar,vazx,.true.,llvazx,n1gc,imode)
      if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    endif

    if ( lwrthess ) then
      deallocate(vatra)
      deallocate(dg_vbar)
    end if

  end subroutine min_writeHessian


  subroutine min_analysisPert(vatra,iztrl,zdf1,column,columng, &
                              obsSpaceData)
    !
    !:Purpose: To use QNA_N1QN3 minimization to perform analysis step on
    !          ensemble perturbations
    implicit none

    ! Arguments
    real(8)                        :: vatra(:)
    integer                        :: iztrl(:)
    real(8)                        :: zdf1
    type(struct_columnData),target :: column,columng
    type(struct_obs),target        :: obsSpaceData

    ! Locals
    type(struct_gsv) :: statevector_ens(numAnalyses)
    type(struct_gsv) :: statevector_mean, statevector_incr, statevector_incr_perturbed, statevector_randpert
    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl
    real(8), allocatable :: incr_cv(:)
    real(8)           :: scalefactor
    integer :: indexAnalysis,stepIndex
    character(len=80) :: fileName
    character(len=4)  :: censnumber
    character(len=8)   :: datestr_last
    character(len=2)   :: hourstr_last
    character(len=2)   :: trialTimeIndex_str
    integer :: stamp_last, ndate, ntime, ierr, newdate
    real(8),allocatable :: vazg(:)
    real(8) :: zjsp, zxmin
    real(8) :: dlds(1)
    real :: zzsunused(1)
    integer :: intUnused(1)
    integer :: nulout, impres, simtot, indic
    logical :: llvazx = .false.

    ! check if user wants to compute any perturbed analyses
    if(numAnalyses.le.0) then
      write(*,*) 'min_analysisPert: numAnalyses not positive, do nothing'
      return
    endif

    ! initialization
    hco_anl => agd_getHco('ComputationalGrid')
    vco_anl => col_getVco(columng)

    call gsv_allocate(statevector_mean, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    call gsv_allocate(statevector_incr, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    call gsv_allocate(statevector_incr_perturbed, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    call gsv_allocate(statevector_randpert, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! allocate local arrays
    allocate(incr_cv(nvadim_mpilocal))
    allocate(vazg(nvadim_mpilocal))
    vazg(:)=0.0d0

    ! figure out date strings for origin time of trials
    call incdatr(stamp_last,tim_getDatestamp(),-6.0d0)
    ierr = newdate(stamp_last,ndate,ntime,-3)
    write(datestr_last,'(i8.8)') ndate
    write(hourstr_last,'(i2.2)') ntime/1000000

    ! get all background ensemble members
    do indexAnalysis = 1, numAnalyses
      if(mpi_myid == 0) write(*,*) ' '
      if(mpi_myid == 0) write(*,*) 'min_analysisPert: reading member #', indexAnalysis
      if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      call gsv_allocate(statevector_ens(indexAnalysis), tim_nstepobsinc, hco_anl, vco_anl, &
                        datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                        allocHeight_opt=.false., allocPressure_opt=.false.)

      write(censnumber,'(i4.4)') indexAnalysis + (envar_loop-1)*numAnalyses
      fileName = trim(ensPathName) // '/' // datestr_last // hourstr_last // '_006_' // trim(censnumber)
      if(mpi_myid == 0) write(*,*) 'Reading from file: ', fileName
      do stepIndex = 1, statevector_ens(indexAnalysis)%numStep
        call gsv_readFromFile( statevector_ens(indexAnalysis),fileName,' ','P',stepIndex_opt=stepIndex, &
                               containsFullField_opt=.true. )
      enddo

    enddo

    ! read the background ensemble mean
    ! NOTE: assume it is supplied to VAR task as the trial file
    if(mpi_myid == 0) write(*,*) 'min_analysisPert: reading ensemble mean'
    do stepIndex = 1, statevector_mean%numStep
      write(trialTimeIndex_str,'(i2.2)') stepIndex
      fileName = './trlm_' // trim(trialTimeIndex_str)
      call gsv_readFromFile( statevector_mean,fileName,' ','P',stepIndex_opt=stepIndex, &
                             containsFullField_opt=.true. )
    enddo

    ! remove mean
    do indexAnalysis = 1, numAnalyses
      scaleFactor = -1.0d0
      call gsv_add(statevector_mean,statevector_ens(indexAnalysis),scaleFactor)
    enddo
    
    ! do perturbation minimizations
    lvarqc = .false.
    do indexAnalysis = 1, numAnalyses

      ! write out the original background ensemble perturbation
      call writeToFile4D(statevector_ens(indexAnalysis),'./bgpert','BGPERT',indexAnalysis)

      ! multiply by -1 (-xb')
      call gsv_scale(statevector_ens(indexAnalysis),-1.0d0)

      ! compute -H*xb', put in OBS_WORK
      call s2c_tl(statevector_ens(indexAnalysis),column,columng,obsSpaceData)  ! put in column H_horiz
      call oop_Htl(column,columng,obsSpaceData,min_nsim)

      ! undo the multiply by -1 (-xb')
      call gsv_scale(statevector_ens(indexAnalysis),-1.0d0)

      ! obs perturbation added to OBS_WORK and copy into OBS_OMP
      call inn_perturbObs(obsSpaceData,numAnalyses,indexAnalysis,envar_loop,OBS_WORK,OBS_OMP)

      ! compute initial gradient
      incr_cv(:) = 0.0d0
      indic = 2
      call simvar(indic, nvadim_mpilocal, incr_cv, zjsp, vazg)

      ! Initialization for call to QNA_N1QN3
      itertot = numIterMax_pert
      simtot = 2*itertot
      if(indexAnalysis == 1) then
        if(preconFileExists) then
          ! warm start with Hessian from file
          imode = 2
        else
          ! cold start and initially double number of iterations
          imode = 0
          vatra(:) = 0.0d0
          itertot = 2*itertot
          simtot = 2*itertot
        endif
      else
        ! warm start with Hessian from previous member
        imode=2
      endif
      zeps1 = 1.0d-5
      zxmin = epsilon(zxmin)
      nulout = 6
      if(mpi_myid == 0) then
        impres=5
      else 
        impres=0
      endif
      ! call QNA_N1QN3 minimization for perturbation
      call tmg_start(70,'QN')
      call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, incr_cv,  &
          zjsp, vazg, zxmin, zdf1, zeps1, impres, nulout, imode,       &
          itertot, simtot ,iztrl, vatra, nmtra, intUnused,   &
          zzsunused, dlds)
      call tmg_stop(70)
      call fool_optimizer(obsSpaceData)

      ! multiply by B^1/2
      call bmat_sqrtB(incr_cv,nvadim_mpilocal,statevector_incr)

      ! write perturbation analysis increment to file before adding random model-error
      !call writeToFile4D(statevector_incr,'./pert_inc0','PERT_INC0',indexAnalysis)

      ! compute random model-error
      call calcRandomPert(statevector_randpert,numAnalyses,indexAnalysis)

      ! add E1 random model-error to increment and write to file
      call gsv_copy(statevector_randpert,statevector_incr_perturbed)
      call gsv_scale(statevector_incr_perturbed,e1_scaleFactor)
      call gsv_add(statevector_incr,statevector_incr_perturbed)
      call writeToFile4D(statevector_incr_perturbed,'./pert_inc1','PERT_INC1',indexAnalysis)

      ! add E2 random model-error to increment and write to file
      call gsv_copy(statevector_randpert,statevector_incr_perturbed)
      call gsv_scale(statevector_incr_perturbed,e2_scaleFactor)
      call gsv_add(statevector_incr,statevector_incr_perturbed)
      call writeToFile4D(statevector_incr_perturbed,'./pert_inc2','PERT_INC2',indexAnalysis)

      ! deallocate statevector_ens(indexAnalysis), no longer needed
      call gsv_deallocate(statevector_ens(indexAnalysis))

    enddo

    ! Write out the final Hessian to file
    if(lwrthess) then
      call hessianIO (preconFileNameOut_pert,1,  &
        min_nsim,tim_getDatestamp(),zeps1,zdf1,itertot,simtot,  &
        iztrl,vatra,dg_vbar,incr_cv,.true.,llvazx,n1gc,imode)
    endif

    ! deallocate all local arrays
    deallocate(incr_cv)
    deallocate(vazg)
    call gsv_deallocate(statevector_incr)
    call gsv_deallocate(statevector_incr_perturbed)
    call gsv_deallocate(statevector_randpert)
    call gsv_deallocate(statevector_mean)

  end subroutine min_analysisPert


  subroutine writeToFile4D(statevector,fileName,cetiket,indexAnalysis)
    implicit none
    ! arguments
    type(struct_gsv)   :: statevector
    character(len=*)   :: fileName
    character(len=*)   :: cetiket
    integer            :: indexAnalysis
    ! locals
    character(len=100) :: fileNameFull
    integer            :: stepIndex

    do stepIndex = 1, statevector%numStep
      fileNameFull = trim(fileName) // trim(fileNameExt(statevector,stepIndex,indexAnalysis))
      call gsv_writeToFile(statevector,fileNameFull,cetiket,scaleFactor_opt=1.0d0, &
                           ip3_opt=0, stepIndex_opt=stepIndex)
    enddo

  end subroutine writeToFile4D


  function fileNameExt(statevector,stepIndex,indexAnalysis) result(fileNameExtStr)
    implicit none
    ! arguments
    type(struct_gsv)  :: statevector
    integer           :: stepIndex, indexAnalysis
    character(len=20) :: fileNameExtStr
    ! locals
    real(8) :: deltaHours
    character(len=4) :: coffset, cmember

    call difdatr(gsv_getDateStamp(statevector,stepIndex),tim_getDatestamp(),deltaHours)
    if(nint(deltaHours*60.0d0).lt.0) then
      write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
    else
      write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
    endif

    write(cmember,'(I4.4)') indexAnalysis + (envar_loop-1)*numAnalyses

    fileNameExtStr = '_' // trim(coffset) // 'm' // '_' // trim(cmember)

  end function fileNameExt


  subroutine calcRandomPert(statevector_randpert,numAnalyses,indexAnalysis)
    !
    !:Purpose: To compute additive inflation random perturbations for VarEnKF
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector_randpert
    integer :: numAnalyses, indexAnalysis

    ! Locals:
    integer :: iseed,jj,nlev_T,nlev_M,jvar,jlev,indexAnalysis2
    real(pre_incrReal), pointer :: field(:,:,:,:)
    real(8), allocatable :: cv_pert_mpiglobal(:), cv_pert_mpilocal(:)
    real(8), pointer :: cv_pert_bens_mpilocal(:), cv_pert_bhi_mpilocal(:)
    real(8), pointer :: cv_pert_bchm_mpilocal(:)
    real(8), allocatable :: scaleFactorBhi(:),scaleFactorBchm(:,:)
    logical, save :: firstTime = .true.
    real(8), allocatable, save :: cv_pert_mean_mpilocal(:)
    integer :: lon1, lon2, lat1, lat2, icount, nlev, nlevmax

    allocate(cv_pert_mpiglobal(cvm_nvadim_mpiglobal))
    allocate(cv_pert_mpilocal(cvm_nvadim))

    if(firstTime) then
      firstTime = .false.
      ! compute mean perturbation for this batch
      cv_pert_mpiglobal(:) = 0.0d0
      do indexAnalysis2 = 1, numAnalyses
        iseed = abs(tim_getDatestamp()) + indexAnalysis2 + (envar_loop-1)*numAnalyses
        call rng_setup(iseed) ! JFC: should be called only once, no???
        do jj = 1, cvm_nvadim_mpiglobal
          cv_pert_mpiglobal(jj) = cv_pert_mpiglobal(jj) + rng_gaussian()
        enddo
      enddo
      cv_pert_mpiglobal(:) = cv_pert_mpiglobal(:)/real(numAnalyses,8)

      allocate(cv_pert_mean_mpilocal(cvm_nvadim))
      call bmat_reduceToMPILocal( cv_pert_mean_mpilocal, & ! OUT
                                  cv_pert_mpiglobal )      ! IN

    endif ! firstTime

    ! compute perturbation and make mpilocal
    iseed = abs(tim_getDatestamp()) + indexAnalysis + (envar_loop-1)*numAnalyses
    write(*,*) 'min_calcRandomPert: indexAnalysis, iseed=', indexAnalysis, iseed
    call rng_setup(iseed) ! JFC : why re-initializing the seed??? 
    do jj = 1, cvm_nvadim_mpiglobal
      cv_pert_mpiglobal(jj) = rng_gaussian()
    enddo

    call bmat_reduceToMPILocal( cv_pert_mpilocal,  & ! OUT
                                cv_pert_mpiglobal )  ! IN
    deallocate(cv_pert_mpiglobal)

    ! remove the ensemble mean
    cv_pert_mpilocal(:) = cv_pert_mpilocal(:) - cv_pert_mean_mpilocal(:)

    if(pertBhiOnly) then
      ! set Bensemble component of control vector to zero
      if(cvm_subVectorExists('B_ENS')) then
        cv_pert_bens_mpilocal => cvm_getSubVector(cv_pert_mpilocal,'B_ENS')
        if(associated(cv_pert_bens_mpilocal)) cv_pert_bens_mpilocal(:) = 0.0d0
      endif
    endif

    ! do spectral truncation of control vector
    if(ntrunc_pert.gt.0) then
      ! Check for weather field static covariances
      if(cvm_subVectorExists('B_HI')) then
        cv_pert_bhi_mpilocal => cvm_getSubVector(cv_pert_mpilocal,'B_HI')
        if (associated(cv_pert_bhi_mpilocal)) call bhi_truncateCV(cv_pert_bhi_mpilocal,ntrunc_pert)
      endif

      ! Check for constituent field static covariances
      if(cvm_subVectorExists('B_CHM')) then
        cv_pert_bchm_mpilocal => cvm_getSubVector(cv_pert_mpilocal,'B_CHM')
        if (associated(cv_pert_bchm_mpilocal)) call bchm_truncateCV(cv_pert_bchm_mpilocal,ntrunc_pert)
      endif

    endif

    call bmat_sqrtB(cv_pert_mpilocal,cvm_nvadim,statevector_randpert)

    if(ntrunc_pert.gt.0) then
      write(*,*) 'WARNING: No scaleFactor applied to truncated perturbation!!!'
    endif

    deallocate(cv_pert_mpilocal)

    lon1=statevector_randpert%myLonBeg
    lon2=statevector_randpert%myLonEnd
    lat1=statevector_randpert%myLatBeg
    lat2=statevector_randpert%myLatEnd

    ! undo the Bhi (and Bchm) scaleFactor(s)
    if(pertBhiOnly) then
      nlev_T = gsv_getNumLev(statevector_randpert,'TH')
      nlev_M = gsv_getNumLev(statevector_randpert,'MM')
      nlevmax=max(nlev_T,nlev_M)
      allocate(scaleFactorBhi(nlevmax))
      call bhi_getScaleFactor(scaleFactorBhi)
      icount=0
      do jvar=1,vnl_numvarmax 
        if(gsv_varExist(statevector_randpert,vnl_varNameList(jvar))) then
           call gsv_getField(statevector_randpert,field,vnl_varNameList(jvar))
           nlev=gsv_getNumLev(statevector_randpert,vnl_varLevelFromVarname(vnl_varNameList(jvar))) 
           if (vnl_varKindFromVarname(vnl_varNameList(jvar)).eq.'MT') then
             write(*,*) 'min_calcRandomPert: undo Bhi scaleFactor varname= ',vnl_varNameList(jvar)
             if (nlev.gt.1) then 
                do jlev = 1,nlev 
                  if(scaleFactorBhi(jlev).gt.0.0d0) then
                    field(lon1:lon2,lat1:lat2,jlev,:)=field(lon1:lon2,lat1:lat2,jlev,:)/scaleFactorBhi(jlev)
                  endif
                enddo
             else
                if(scaleFactorBhi(nlevmax).gt.0.0d0) then
                  field(lon1:lon2,lat1:lat2,1,:)=field(lon1:lon2,lat1:lat2,1,:)/scaleFactorBhi(nlevmax)
                endif
             end if
           else if (vnl_varKindFromVarname(vnl_varNameList(jvar)).eq.'CH') then
             if (icount.eq.0) then
                allocate(scaleFactorBchm(nlevmax,100))
                call bchm_getScaleFactor(scaleFactorBchm)
             end if
             icount=icount+1
             write(*,*) 'min_calcRandomPert: undo Bchm scaleFactor varname= ',vnl_varNameList(jvar)
             if (nlev.gt.1) then 
                do jlev = 1,nlev  
                  if(scaleFactorBchm(jlev,icount).gt.0.0d0) then
                    field(lon1:lon2,lat1:lat2,jlev,:)=field(lon1:lon2,lat1:lat2,jlev,:)/scaleFactorBChm(jlev,icount)
                  endif
                enddo
             else 
                if(scaleFactorBchm(1,icount).gt.0.0d0) then
                   field(lon1:lon2,lat1:lat2,1,:)=field(lon1:lon2,lat1:lat2,1,:)/scaleFactorBChm(1,icount)
                endif
             end if
           endif
        endif
      enddo
      deallocate(scaleFactorBhi)
      if (icount.gt.0) deallocate(scaleFactorBchm)
    endif

    ! apply additional scaling of random perturbation (initially only for UU, VV)
    nlev_T = gsv_getNumLev(statevector_randpert,'TH')
    nlev_M = gsv_getNumLev(statevector_randpert,'MM')
    ! for 3D variables
    do jvar=1,vnl_numvarmax3D 
      if(gsv_varExist(statevector_randpert,vnl_varNameList3D(jvar)).and.  &
            (trim(vnl_varNameList3D(jvar)) == 'UU' .or.  &
             trim(vnl_varNameList3D(jvar)) == 'VV') ) then
        write(*,*) 'min_calcRandomPert: pertScaleFactor_UV varname= ',vnl_varNameList3D(jvar)
        call gsv_getField(statevector_randpert,field,vnl_varNameList3D(jvar))
        do jlev = 1, gsv_getNumLev(statevector_randpert,vnl_varLevelFromVarname(vnl_varNameList3D(jvar)))   
          write(*,*) 'min_calcRandomPert: pertScaleFactor_UV= ',jlev,pertScaleFactor_UV(jlev)
          field(lon1:lon2,lat1:lat2,jlev,:)=field(lon1:lon2,lat1:lat2,jlev,:)*pertScaleFactor_UV(jlev)
        enddo
      endif
    enddo

  end subroutine calcRandomPert


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
    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl

    if (na_indic  ==  1 .or. na_indic  ==  4) call tmg_stop(70)
    call tmg_start(80,'MIN_SIMVAR')
    if (na_indic .ne. 1) then ! No action taken if na_indic == 1
       min_nsim = min_nsim + 1

       if(mpi_myid == 0) then
         write(*,*) 'Entering simvar for simulation ',min_nsim
         write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
       endif

       ! note: dg_vbar = sum(v) of previous outer-loops
       dl_v(1:nvadim_mpilocal) = da_v(1:nvadim_mpilocal) + dg_vbar(1:nvadim_mpilocal)
     
       ! Computation of background term of cost function:
       dl_Jb = dot_product(dl_v(1:nvadim_mpilocal),dl_v(1:nvadim_mpilocal))/2.d0  
       call mpi_allreduce_sumreal8scalar(dl_Jb,"GRID")

       if (oneDVarMode) then
         call var1D_sqrtB(da_v, nvadim_mpilocal, column_ptr, obsSpaceData_ptr)
         call cvt_transform(column_ptr, columng_ptr, 'PsfcToP_tl')
       else
         if (.not.statevector%allocated) then
           write(*,*) 'min-simvar: allocating increment stateVector'
           hco_anl => agd_getHco('ComputationalGrid')
           vco_anl => col_getVco(columng_ptr)
           call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                dataKind_opt=pre_incrReal, mpi_local_opt=.true.)
           call gsv_readMaskFromFile(statevector,'./analysisgrid')
       
           if ( associated(stateVectorRefHU_ptr) ) then
             call bmat_sqrtB(da_v,nvadim_mpilocal,statevector, &
                  stateVectorRef_opt=stateVectorRefHU_ptr)
           else
             call bmat_sqrtB(da_v,nvadim_mpilocal,statevector)
           end if

           call tmg_start(30,'OBS_INTERP')
           call s2c_tl(statevector,column_ptr,columng_ptr,obsSpaceData_ptr)  ! put in column H_horiz dx
           call tmg_stop(30)
         end if
       end if
      
       call tmg_start(40,'OBS_TL')
       call oop_Htl(column_ptr,columng_ptr,obsSpaceData_ptr,min_nsim)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
       call tmg_stop(40)

       call res_compute(obsSpaceData_ptr)  ! Calculate OBS_OMA from OBS_WORK : d-Hdx

       call bcs_calcbias_tl(da_v,OBS_OMA,obsSpaceData_ptr,columng_ptr)

       call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_OMA)  ! Save as OBS_WORK : R**-1/2 (d-Hdx)

       call cfn_calcJo(obsSpaceData_ptr)  ! Store J-obs in OBS_JOBS : 1/2 * R**-1 (d-Hdx)**2

       IF (LVARQC) THEN
          call vqc_tl(obsSpaceData_ptr)  ! Store modified J_obs in OBS_JOBS : -ln((gamma-exp(J))/(gamma+1)) 
       endif

       dl_Jo = 0.d0
       call cfn_sumJo(obsSpaceData_ptr,dl_Jo)
       da_J = dl_Jb + dl_Jo
       if (na_indic  ==  3) then
          da_J = dl_Jo
          IF(mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  JO = ",G23.16,6X)') dl_Jo
       else
          da_J = dl_Jb + dl_Jo
          IF(mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  Jb = ",G23.16,6X,"JO = ",G23.16,6X,"Jt = ",G23.16)') dl_Jb,dl_Jo,da_J
       endif

       call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  ! Modify OBS_WORK : R**-1 (d-Hdx)

       IF (LVARQC) THEN
          call vqc_ad(obsSpaceData_ptr)
       endif

       call res_computeAd(obsSpaceData_ptr)  ! Calculate adjoint of d-Hdx (mult OBS_WORK by -1)

       call col_zero(column_ptr)

       call tmg_start(41,'OBS_AD')
       call oop_Had(column_ptr,columng_ptr,obsSpaceData_ptr)   ! Put in column : -H_vert**T R**-1 (d-Hdx)
       call tmg_stop(41)

       if (oneDVarMode) then
         ! no interpolation needed for 1Dvar case
       else
         call tmg_start(31,'OBS_INTERPAD')
         call s2c_ad(statevector,column_ptr,columng_ptr,obsSpaceData_ptr)  ! Put in statevector -H_horiz**T H_vert**T R**-1 (d-Hdx)
         call tmg_stop(31)
       end if

       da_gradJ(:) = 0.d0
       call bcs_calcbias_ad(da_gradJ,OBS_WORK,obsSpaceData_ptr)

       if (oneDVarMOde) then
         call cvt_transform( column, columng, 'PsfcToP_ad')      ! IN
         call var1D_sqrtBT(da_gradJ, nvadim_mpilocal, column_ptr, obsSpaceData_ptr)
       else
         if ( associated(stateVectorRefHU_ptr) ) then
           call bmat_sqrtBT(da_gradJ,nvadim_mpilocal,statevector, &
                            stateVectorRef_opt=stateVectorRefHU_ptr)
         else
           call bmat_sqrtBT(da_gradJ,nvadim_mpilocal,statevector)
           !call gsv_deallocate(statevector)
         end if
       end if

       if (na_indic .ne. 3) then
         da_gradJ(1:nvadim_mpilocal) = dl_v(1:nvadim_mpilocal) + da_gradJ(1:nvadim_mpilocal)
       endif
    endif
    call tmg_stop(80)
    if (na_indic  ==  1 .or. na_indic  ==  4) call tmg_start(70,'QN')

    if(mpi_myid == 0) write(*,*) 'end of simvar'

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

    call mpi_allreduce_sumreal8scalar(ddsc,"GRID")

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
      if (mpi_myid == 0) write(*,*) 'Read  Hessian in min_hessianIO from file ', cfname
      call tmg_start(88,'MIN_READHESS')
    elseif (status == 1) then
      if (mpi_myid == 0) write(*,*) 'Write Hessian in min_hessianIO to file ', cfname
      call tmg_start(94,'MIN_WRITEHESS')
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

      if (mpi_myid == 0) then
         !- Open the Hessian matrix file
         ierr = fnom(ireslun,cfname,'FTN+SEQ+UNF+OLD+R/O',0)

         !- Checking version number
         call tmg_start(98,'MIN_READHESS_IO')
         read (ireslun) cl_version, i1gc
         call tmg_stop(98)
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

         call tmg_start(98,'MIN_READHESS_IO')
         rewind (ireslun)

         read(ireslun) cl_version,i1gc,nsim,ibrpstamp,zeps1,zdf1,itertot,isimtot,ivadim,itrunc
         read(ireslun) ivamaj,iztrl_io
         call tmg_stop(98)
         if ((ivamaj.ne.nvamaj).or.(nvadim_mpiglobal.ne.ivadim)) then
            write(*,*) nvamaj,ivamaj,nvadim_mpiglobal,ivadim
            call utl_abort('min_hessianIO : ERROR, size of V5 Hessian not consistent')
         endif

      end if

      ! ibrpstamp and iztrl_io must be broadcasted
      call rpn_comm_bcast(ibrpstamp,  1, "MPI_INTEGER", 0, "GRID", ierr)
      call rpn_comm_bcast(iztrl_io , 10, "MPI_INTEGER", 0, "GRID", ierr)

      !- Read the Hessian
      if(mpi_myid == 0) then 
         write(*,*) 'min_hessianIO : reading Hessian'
         allocate(vatravec_r4_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vatravec_r4_mpiglobal(1))
      end if
      allocate(vatra_r4(nvadim_mpilocal))

      if (k1gc == 3) ictrlvec = 2*nvamaj+1
      do jvec = 1, ictrlvec
 
         if (mpi_myid == 0) then
            call tmg_start(98,'MIN_READHESS_IO')
            read(ireslun) vatravec_r4_mpiglobal
            call tmg_stop(98)
         end if

         call tmg_start(119,'MIN_READHESS_REDUCE')
         call bmat_reduceToMPILocal_r4( vatra_r4,            & ! OUT
                                        vatravec_r4_mpiglobal) ! IN (contains data only on proc 0)
         call tmg_stop(119)
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
         if (mpi_myid == 0) then 
            write(*,*) 'min_hessianIO : reading vazxbar'
            allocate(vazxbar_mpiglobal(nvadim_mpiglobal))
            call tmg_start(98,'MIN_READHESS_IO')
            read(ireslun) vazxbar_mpiglobal
            call tmg_stop(98)
         else
            allocate(vazxbar_mpiglobal(1))
         end if
         call bmat_reduceToMPILocal( vazxbar,          & ! OUT
                                     vazxbar_mpiglobal ) ! IN (contains data only on proc 0)
         deallocate(vazxbar_mpiglobal)
      end if

      !- Read VAZX
      if (ibrpstamp == kbrpstamp .and. llvazx) then
         if (mpi_myid == 0) then 
            write(*,*) 'min_hessianIO : reading vazx'
            allocate(vazx_mpiglobal(nvadim_mpiglobal))
            call tmg_start(98,'MIN_READHESS_IO')
            read(ireslun) vazx_mpiglobal
            call tmg_stop(98)
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

      if (mpi_myid == 0) ierr = fclos(ireslun)

    elseif (status == 1) then
      !
      !- 2.  Write Hessian
      !

      !- Open the Hessian matrix file on processor 0 and write some metadata 
      if (mpi_myid == 0) ierr = fnom(ireslun,cfname, 'FTN+SEQ+UNF' , 0)
      cl_version = 'V5'
      itrunc=0
      iztrl_io(1: 5) = iztrl(:)
      iztrl_io(6:10) = 0 ! dummy values for hessian size legacy purposes
      call tmg_start(76,'MIN_WRITEHESS_IO')
      if(mpi_myid == 0) write(ireslun) cl_version,k1gc,nsim,kbrpstamp,zeps1,zdf1,itertot,isimtot,nvadim_mpiglobal,itrunc
      if(mpi_myid == 0) write(ireslun) nvamaj,iztrl_io
      call tmg_stop(76)

      !- Write the Hessian
      if (mpi_myid == 0) then
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
        call tmg_start(76,'MIN_WRITEHESS_IO')
        if (mpi_myid == 0) write(ireslun) vatravec_r4_mpiglobal
        call tmg_stop(76)
      enddo
      deallocate(vatravec_r4_mpiglobal)
      deallocate(vatra_r4)

      !- Write VAZXBAR
      if (mpi_myid == 0) then
         allocate(vazxbar_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vazxbar_mpiglobal(1)) ! dummy
      end if
      call bmat_expandToMPIGlobal( vazxbar,          & ! IN
                                   vazxbar_mpiglobal ) ! OUT
      call tmg_start(76,'MIN_WRITEHESS_IO')
      if (mpi_myid == 0) write(ireslun) vazxbar_mpiglobal(1:nvadim_mpiglobal)
      call tmg_stop(76)
      deallocate(vazxbar_mpiglobal)

      !- Write VAZX
      if (mpi_myid == 0) then
         allocate(vazx_mpiglobal(nvadim_mpiglobal))
      else
         allocate(vazx_mpiglobal(1))
      end if
      call bmat_expandToMPIGlobal( vazx,          & ! IN
                                   vazx_mpiglobal ) ! OUT
      call tmg_start(76,'MIN_WRITEHESS_IO')
      if (mpi_myid == 0) write(ireslun) vazx_mpiglobal(1:nvadim_mpiglobal)
      call tmg_stop(76)
      deallocate(vazx_mpiglobal)

      !- Close the Hessian matrix file
      if (mpi_myid == 0) ierr = fclos(ireslun)
    else
      call utl_abort('min_hessianIO: status not valid')
    endif

    if (status == 0) then
      call tmg_stop(88)
    elseif (status == 1) then
      call tmg_stop(94)
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
  integer :: nl_indic, nl_j !, i, ierr
  real*8  :: dl_wrk(na_dim),dl_gradj0(na_dim), dl_x(na_dim)
  real*8  :: dl_J0, dl_J, dl_test, dl_start,dl_end
  real*8  :: dl_alpha, dl_gnorm0 !, xsave,xpert


  ! 1. Initialize dl_gradj0 at da_x0
  !    ------------------------------------

  nl_indic = 2
<<<<<<< HEAD
  call simul(nl_indic,na_dim,da_x0,dl_j0,dl_gradj0)
=======
  call simul(nl_indic, na_dim, da_x0, dl_j0,dl_gradj0,dataptr(1))
!  dl_x(:) = da_x0(:)
!  xpert = 1.d-6
!  do i=1, 320
!    if (mpi_myId ==0 ) then
!       xsave = dl_x(i)
!       dl_x(i) = xpert + dl_x(i)
!    end if 
!    call simul(nl_indic,na_dim,dl_x,dl_j,dl_wrk,dataptr(1))
!    if (mpi_myId ==0 ) then
!      write(200,'(5e26.18)') da_x0(i), dl_x(i), dl_gradj0(i), dl_J, dl_j0
!      dl_x(i) = xsave
!   end if
!  end do
!  if (mpi_myId ==0 ) flush(200)
!  call rpn_comm_barrier('GRID',ierr)
!  call utl_abort('test Sylvain')

>>>>>>> Issue #309: introduction of code for 1Dvar program
  dl_gnorm0 = dot_product(dl_gradj0,dl_gradj0)
  call mpi_allreduce_sumreal8scalar(dl_gnorm0,"GRID")
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
