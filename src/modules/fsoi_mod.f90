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

module fsoi_mod
  ! MODULE fsoi_mod (prefix='fso' category='1. High-level functionality')
  !
  ! :Purpose: Observation impact (FSOI) library
  !
  use midasMpi_mod
  use codePrecision_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use obsSpaceData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use columnData_mod
  use controlVector_mod
  use bmatrix_mod
  use bmatrixensemble_mod
  use stateToColumn_mod
  use obsOperators_mod
  use rMatrix_mod
  use MathPhysConstants_mod
  use quasinewton_mod
  use costFunction_mod
  use tovs_nl_mod
  use timeCoord_mod
  use utilities_mod
  use rttov_const, only: inst_name, platform_name
  implicit none
  save
  private

  ! public subroutines and functions
  public :: fso_setup, fso_ensemble

  ! module private variables
  type(struct_obs),        pointer :: obsSpaceData_ptr
  type(struct_columnData), pointer :: columnTrlOnAnlIncLev_ptr
  type(struct_columnData), pointer :: column_ptr
  real(8),allocatable :: vhat(:)
  integer,external    :: get_max_rss
  integer             :: fso_nsim, nvadim_mpilocal

  type(struct_hco), pointer :: hco_anl => null()

  ! namelist variables
  integer             :: nvamaj, nitermax, nsimmax
  real(8)             :: leadTime, repsg, rdf1fac,latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm
  logical             :: includeUVnorm, includeTTnorm, includeP0norm, includeHUnorm, includeTGnorm
  character(len=256)  :: forecastPath
  character(len=4)    :: fsoMode

  contains

  !--------------------------------------------------------------------------
  ! fso_setup
  !--------------------------------------------------------------------------
  subroutine fso_setup(hco_anl_in)
    !
    ! :Purpose: Initialise the FSOI module: read the namelist and initialise
    !           global variables and structure
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco_anl_in

    ! Locals:
    integer :: ierr,nulnam
    integer :: fnom,fclos

    NAMELIST /NAMFSO/leadTime, nvamaj, nitermax, nsimmax
    NAMELIST /NAMFSO/repsg, rdf1fac, forecastPath, fsoMode
    NAMELIST /NAMFSO/latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm
    NAMELIST /NAMFSO/includeUVnorm, includeTTnorm, includeP0norm, includeHUnorm, includeTGnorm

    ! set default values for namelist variables
    leadtime = 12.0d0
    nvamaj = 6
    nitermax = 100
    nsimmax  = 120
    repsg    = 1d-5
    rdf1fac  = 0.25d0
    latMinNorm = -95.0d0
    latMaxNorm = 95.0d0
    lonMinNorm = -185.0d0
    lonMaxNorm = 365.0d0
    includeUVnorm=.true.
    includeTTnorm=.true.
    includeP0norm=.true.
    includeHUnorm=.false.
    includeTGnorm=.false.
    forecastPath = './forecasts'
    fsoMode  = 'HFSO'

    ! read in the namelist NAMFSO
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namfso,iostat=ierr)
    if(ierr /= 0) call utl_abort('fso_setup: Error reading namelist')
    write(*,nml=namfso)
    ierr = fclos(nulnam)
    ! convert Latmin,max and Lonmin,max from degres into RAD
    latMinNorm = latMinNorm*MPC_RADIANS_PER_DEGREE_R8
    latMaxNorm = latMaxNorm*MPC_RADIANS_PER_DEGREE_R8
    lonMinNorm = lonMinNorm*MPC_RADIANS_PER_DEGREE_R8
    lonMaxNorm = lonMaxNorm*MPC_RADIANS_PER_DEGREE_R8

    call ben_setFsoLeadTime(leadTime)
    fso_nsim = 0

    hco_anl => hco_anl_in

  end subroutine fso_setup

  !--------------------------------------------------------------------------
  ! fso_ensemble
  !--------------------------------------------------------------------------
  subroutine fso_ensemble(columnTrlOnAnlIncLev,obsSpaceData)
    !
    ! :Purpose: Perform forecast sensitivity to observation calculation using
    !           ensemble approach
    !
    implicit none

    ! Arguments:
    type(struct_columnData), target, intent(in)     :: columnTrlOnAnlIncLev
    type(struct_obs),        target, intent(inout)  :: obsSpaceData

    ! Locals:
    type(struct_columnData),target  :: column
    type(struct_gsv)                :: statevector_FcstErr, statevector_fso, statevector_HUreference
    type(struct_vco), pointer       :: vco_anl
    real(8),allocatable             :: ahat(:), zhat(:)
    integer                         :: dateStamp_fcst, dateStamp

    !for Observation space
    integer                         :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex
    real(8)                         :: fso_ori, fso_fin

    if (mmpi_myid == 0) write(*,*) 'fso_ensemble: starting'

    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    nvadim_mpilocal = cvm_nvadim

    ! initialize column object for storing "increment"
    call col_setVco(column,col_getVco(columnTrlOnAnlIncLev))
    call col_allocate(column,col_getNumCol(columnTrlOnAnlIncLev),mpiLocal_opt=.true.)

    ! compute dateStamp_fcst
    call incdatr(dateStamp_fcst, tim_getDatestamp(), leadTime)
    dateStamp = tim_getDatestamp()
    write(*,*) 'fso_ensemble: analysis datestamp = ',dateStamp
    write(*,*) 'fso_ensemble: forecast datestamp = ',dateStamp_fcst

    ! allocate control vector related arrays (these are all mpilocal)
    allocate(ahat(nvadim_mpilocal))
    allocate(vhat(nvadim_mpilocal))
    allocate(zhat(nvadim_mpilocal))

    ! initialize control vector related arrays to zero
    ahat(:) = 0.0d0
    vhat(:) = 0.0d0

    ! for statevector_FcstErr
    call gsv_allocate(statevector_FcstErr, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! for statevector_fso
    call gsv_allocate(statevector_fso, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.)

    ! for statevector_HUreference (verifying analysis)
    call gsv_allocate(statevector_HUreference, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! compute forecast error = C * (error_t^fa + error_t^fb)
    call calcFcstError(columnTrlOnAnlIncLev,statevector_FcstErr,statevector_HUreference)


    ! compute vhat = B_t^T/2 * C * (error_t^fa + error_t^fb)
    call bmat_sqrtBT(vhat, nvadim_mpilocal, statevector_FcstErr, useFSOFcst_opt = .true., &
    stateVectorRef_opt=statevector_HUreference)

    if (mmpi_myid == 0) write(*,*) 'fso: B_t^T/2 * C * (error_t^fa + error_t^fb) max,min:', &
        maxval(vhat),minval(vhat)

    if( trim(fsoMode) == 'HFSO' ) then
      call minimize(nvadim_mpilocal, zhat, column, columnTrlOnAnlIncLev, obsSpaceData)
      ahat = zhat + vhat
      call bmat_sqrtB(ahat, nvadim_mpilocal, statevector_fso)
    elseif( trim(fsoMode) == 'EFSO' ) then
      call bmat_sqrtB(vhat, nvadim_mpilocal, statevector_fso)
    end if

    ! Compute yhat = [R^-1 H B^1/2 ahat], and put in OBS_FSO
    call s2c_tl(statevector_fso,column,columnTrlOnAnlIncLev,obsSpaceData)  ! put in column H_horiz B^1/2 ahat
    call oop_Htl(column,columnTrlOnAnlIncLev,obsSpaceData,1)          ! Save as OBS_WORK: H_vert H_horiz B^1/2 vhat = H B^1/2 ahat
    call rmat_RsqrtInverseAllObs(obsSpaceData,OBS_FSO,OBS_WORK) ! Save as OBS_FSO : R**-1/2 H B^1/2 ahat
    call rmat_RsqrtInverseAllObs(obsSpaceData,OBS_FSO,OBS_FSO)  ! Save as OBS_FSO : R**-1 H B^1/2 ahat\

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
          fso_ori = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
          fso_fin = fso_ori * obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
          call obs_bodySet_r(obsSpaceData,OBS_FSO,bodyIndex, fso_fin)
        end if
      end do

    end do

    ! print out the information of total FSO for each family
    call sumFSO(obsSpaceData)

    ! deallocate the control vector related arrays
    deallocate(ahat)
    deallocate(vhat)
    deallocate(zhat)
    call col_deallocate(column)

    if (mmpi_myid == 0) write(*,*) 'fso_ensemble: Finished'

  end subroutine fso_ensemble

  !--------------------------------------------------------------------------
  ! calcFcstError
  !--------------------------------------------------------------------------
  subroutine calcFcstError(columnTrlOnAnlIncLev,statevector_out,statevector_verifAnalysis)
    !
    ! :Purpose: Reads the forecast from background and analysis, the verifying
    !           analysis based on these inputs, calculates the Forecast error
    !
    implicit none

    ! Arguments:
    type(struct_columnData), target, intent(in)     :: columnTrlOnAnlIncLev
    type(struct_gsv)       , target, intent(inout)  :: statevector_out, statevector_verifAnalysis

    ! Locals:
    type(struct_gsv)                :: statevector_fa, statevector_fb, statevector_a
    character(len=256)              :: fileName_fa, fileName_fb, fileName_a
    logical                         :: faExists
    type(struct_gsv)                :: statevector_tempfa, statevector_tempfb
    type(struct_vco), pointer       :: vco_anl
    integer                         :: dateStamp_fcst, dateStamp

    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    ! compute dateStamp_fcst
    call incdatr(dateStamp_fcst, tim_getDatestamp(), leadTime)
    dateStamp = tim_getDatestamp()
    write(*,*) 'fso_ensemble: analysis datestamp = ',dateStamp
    write(*,*) 'fso_ensemble: forecast datestamp = ',dateStamp_fcst

    ! read forecasts from the analysis and background state
    fileName_fa = trim(forecastPath) // '/forecast_a'
    inquire(file=trim(fileName_fa),exist=faExists)
    write(*,*) 'faExists', faExists
    call gsv_allocate(statevector_fa, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gio_readFromFile(statevector_fa, fileName_fa, ' ', 'P', containsFullField_opt=.true.)

    !for statevecotr_tempfa
    call gsv_allocate(statevector_tempfa, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    !for statevector_fb
    fileName_fb = trim(forecastPath) // '/forecast_b'
    call gsv_allocate(statevector_fb, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gio_readFromFile(statevector_fb, fileName_fb, ' ', 'P', containsFullField_opt=.true.)

    !for statevecotr_tempfb
    call gsv_allocate(statevector_tempfb, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! read verifying analysis
    fileName_a = trim(forecastPath) // '/analysis'
    call gsv_allocate(statevector_a, 1,hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gio_readFromFile(statevector_a, fileName_a, ' ', 'A', containsFullField_opt=.true.)

    ! compute error of both forecasts (overwrite forecasts with error)
    call gsv_add(statevector_a, statevector_fa, -1.0d0)
    call gsv_add(statevector_a, statevector_fb, -1.0d0)

    call gsv_copy(statevector_fa,statevector_tempfa)
    call gsv_copy(statevector_fb,statevector_tempfb)
    call multEnergyNorm(statevector_tempfa, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, &
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state
    call multEnergyNorm(statevector_tempfb, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, &
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state

    ! compute error Norm =  C * (error_t^fa + error_t^fb)
    call gsv_add(statevector_fa, statevector_fb, 1.0d0)
    call multEnergyNorm(statevector_fb, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, &
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state
    call gsv_copy(statevector_fb,statevector_out)
    call gsv_copy(statevector_a,statevector_verifAnalysis)

  end subroutine calcFcstError

  !--------------------------------------------------------------------------
  ! minimize
  !--------------------------------------------------------------------------
  subroutine minimize(nvadim,zhat,column,columnTrlOnAnlIncLev,obsSpaceData)
    !
    ! :Purpose: Performs HFSO quasi-Newton minimization
    !
    implicit none

    ! Arguments:
    integer,                         intent(in)   :: nvadim
    real(8), dimension(nvadim),      intent(out)  :: zhat
    type(struct_columnData), target, intent(in)   :: column, columnTrlOnAnlIncLev
    type(struct_obs),        target, intent(in)   :: obsSpaceData

    ! Locals:
    integer                         :: nulout = 6
    integer                         :: imode, itermax, isimmax, indic, nmtra
    integer                         :: impres, iztrl(10), intunused(1)
    real                            :: rspunused(1)
    real(8)                         :: zjsp, zxmin, zdf1, zeps, dlgnorm, dlxnorm,zspunused(1)
    real(8),allocatable             :: gradJ(:), vatra(:)

    call utl_tmg_start(90,'--Minimization')

    if (mmpi_myid == 0) write(*,*) 'minimize: starting'

    nmtra = (4 + 2*nvamaj)*nvadim
    write(*,'(4X,"NVAMAJ = ",I3,/5X,"NMTRA =",I14)') nvamaj,nmtra

    columnTrlOnAnlIncLev_ptr => columnTrlOnAnlIncLev
    column_ptr  => column
    obsSpaceData_ptr => obsSpaceData

    allocate(gradJ(nvadim))
    allocate(vatra(nmtra))

    gradJ(:) = 0.0d0
    vatra(:) = 0.0d0
    zhat(:)  = 0.0d0

    ! Compute zhat by performing variational minimization
    ! Set-up for the minimization
    if (mmpi_myid == 0) then
      impres = 5
    else
      impres = 0
    end if

    imode = 0
    zeps = repsg
    itermax = nitermax
    isimmax = nsimmax
    zxmin = epsilon(zxmin)
    ! initial gradient calculation
    indic = 2
    call utl_tmg_start(91,'----QuasiNewton')
    call simvar(indic,nvadim,zhat,zjsp,gradJ)
    call utl_tmg_stop(91)
    zdf1 =  rdf1fac * abs(zjsp)

    ! print amplitude of initial gradient and cost function value
    call prscal(nvadim,gradJ,gradJ,dlgnorm)
    dlgnorm = dsqrt(dlgnorm)
    call prscal(nvadim,zhat,zhat,dlxnorm)
    dlxnorm = dsqrt(dlxnorm)
    write(*,*)' |X| = ', dlxnorm
    write(*,'(/4X,"J(X) = ",G23.16,4X,"|Grad J(X)| = ",G23.16)') zjsp, dlgnorm
    write(*,'(//,10X," Minimization QNA_N1QN3 starts ...",/  &
             10x,"DXMIN =",G23.16,2X,"DF1 =",G23.16,2X,"EPSG =",G23.16  &
             /,10X,"IMPRES =",I3,2X,"NITER = ",I3,2X,"NSIM = ",I3)') zxmin,zdf1,zeps,impres,itermax,isimmax

    ! Do the minimization
    call utl_tmg_start(91,'----QuasiNewton')
    call qna_n1qn3(simvar, prscal, dcanonb, dcanab, nvadim, zhat,  &
                   zjsp, gradJ, zxmin, zdf1, zeps, impres, nulout, imode,   &
                   itermax,isimmax, iztrl, vatra, nmtra, intunused, rspunused,  &
                   zspunused)
    call utl_tmg_stop(91)
    call fool_optimizer(obsSpaceData)

    write(*,'(//,20X,20("*"),2X    &
        ,/,20X,"              Minimization ended with MODE:",I4  &
        ,/,20X,"                Total number of iterations:",I4  &
        ,/,20X,"               Total number of simulations:",I4)' ) imode,itermax,isimmax

    if (mmpi_myid == 0) write(*,*) 'minimize: Finished'

    deallocate(vatra)
    deallocate(gradJ)

    call utl_tmg_stop(90)

  end subroutine minimize

  !--------------------------------------------------------------------------
  ! sumFSO
  !--------------------------------------------------------------------------
  subroutine sumFSO(obsSpaceData)
    !
    ! :Purpose: Print out the information of total FSO for each family
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData

    ! Locals:
    real(8)            :: pfso_1
    integer            :: bodyIndex,itvs,isens,headerIndex
    integer            :: bodyIndexBeg, bodyIndexEnd

    integer, parameter :: numFamily = 10
    character(len=2), parameter :: familyList(numFamily) = (/'UA','AI','SF','SC','TO','SW','PR','RO','GP','CH'/)
    real(8)            :: tfso(numFamily), tfsotov_sensors(tvs_nsensors),totFSO
    integer            :: numAss_local(numFamily), numAss_global(numFamily)
    integer            :: numAss_sensors_loc(tvs_nsensors), numAss_sensors_glb(tvs_nsensors)
    integer            :: ierr, familyIndex

    if (mmpi_myid == 0) write(*,*) 'sumFSO: Starting'

    ! initialize
    do familyIndex = 1, numFamily
      tfso(familyIndex) = 0.d0
      numAss_local(familyIndex) = 0
      numAss_global(familyIndex) = 0
    end do

    tfsotov_sensors(:) = 0.d0
    numAss_sensors_loc(:) = 0
    numAss_sensors_glb(:) = 0
    totFSO = 0.d0

    do familyIndex = 1, numFamily
      do bodyIndex = 1, obs_numbody(obsSpaceData)

        pfso_1 = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 ) then
          ! FSO for each family
          if (obs_getFamily(obsSpaceData,bodyIndex=bodyIndex) == familyList(familyIndex) ) then
            tfso(familyIndex) = tfso(familyIndex) + pfso_1
            numAss_local(familyIndex) = numAss_local(familyIndex) + 1
          end if
        end if

      end do
    end do

    do itvs = 1, tvs_nobtov
      headerIndex  = tvs_headerIndex(itvs)
      if (headerIndex > 0 ) then
        bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          pfso_1 = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 ) then
            isens = tvs_lsensor (itvs)
            tfsotov_sensors(isens) =  tfsotov_sensors(isens) + pfso_1
            numAss_sensors_loc(isens) = numAss_sensors_loc(isens) + 1
          end if
        end do
      end if
    end do

    do familyIndex = 1, numFamily
      call mmpi_allreduce_sumreal8scalar(tfso(familyIndex),'GRID')
      totFSO = totFSO + tfso(familyIndex)
      call rpn_comm_allreduce(numAss_local(familyIndex), numAss_global(familyIndex) ,1,'MPI_INTEGER','MPI_SUM','GRID',ierr)
    end do

    do isens = 1, tvs_nsensors
      call mmpi_allreduce_sumreal8scalar(tfsotov_sensors(isens),'GRID')
      call rpn_comm_allreduce(numAss_sensors_loc(isens), numAss_sensors_glb(isens) ,1,'MPI_INTEGER','MPI_SUM','GRID',ierr)
    end do

    if (mmpi_myid == 0) then

      write(*,*) ' '
      write(*,'(a15,f15.8)') 'Total FSO=', totFSO
      write(*,*) ' '

      do familyIndex = 1, numFamily
        write(*,'(a4,a2,a2,f15.8,a16,i10)') 'FSO-', familyList(familyIndex), '=', tfso(familyIndex),'  Count Number=', numAss_global(familyIndex)
      end do
      write(*,*) ' '

      if (tvs_nsensors > 0) then
        write(*,'(1x,a)') 'For TOVS decomposition by sensor:'
        write(*,'(1x,a)') '#  plt sat ins    FSO'
        do isens = 1, tvs_nsensors
          write(*,'(i2,1x,a,1x,a,1x,i2,1x,f15.8,i10)') isens,inst_name(tvs_instruments(isens)), &
                platform_name(tvs_platforms(isens)),tvs_satellites(isens),tfsotov_sensors(isens), numAss_sensors_glb(isens)
        end do
        write(*,*) ' '
      end if

    end if

  end subroutine sumFSO

  !--------------------------------------------------------------------------
  ! simvar
  !--------------------------------------------------------------------------
  subroutine simvar(indic,nvadim,zhat,Jtotal,gradJ)
    !
    ! :Purpose: Implement the Variational solver as described in
    !           Courtier, 1997, Dual formulation of four-dimentional variational
    !           assimilation, Q.J.R., pp2449-2461.
    !
    ! :Arguments:
    !   :indic:   Value of indic
    !             Note: 1 and 4 are reserved values for call back from m1qn3.
    !             For direct calls use other value than 1 and 4.
    !             =1 No action taken; =4 Both J(u) and its gradient are computed.
    !             =2 Same as 4 (compute J and gradJ) but do not interrupt timer
    !             of the minimizer.
    !             =3 Compute Jo and gradJo only.
    !
    !   :nvadim:  Dimension of the control vector in forecast error covariances space
    !
    !   :zhat:    Control variable in forecast error covariances space
    !
    !   :Jtotal:  Cost function of the Variational algorithm
    !
    !   :gradJ:   Gradient of the Variational Cost funtion
    !
    implicit none

    ! Arguments:
    integer,                    intent(in)    :: indic
    integer,                    intent(in)    :: nvadim
    real(8), dimension(nvadim), intent(inout) :: zhat
    real(8),                    intent(out)   :: Jtotal
    real(8), dimension(nvadim), intent(out)   :: gradJ

    ! Locals:
    real(8) :: ahat_vhat(nvadim)
    real(8) :: Jb, Jobs
    type(struct_gsv)          :: statevector
    type(struct_vco), pointer :: vco_anl

    call utl_tmg_stop(91)
    call utl_tmg_stop(90)

    if (indic /= 1) then ! No action taken if indic == 1
      fso_nsim = fso_nsim + 1

      if (mmpi_myid == 0) then
        write(*,*) 'simvar: entering for simulation ',fso_nsim
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
        call flush(6)
      end if

      ! note: vhat = B_t^T/2 hat(del x_t)
      ahat_vhat(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + vhat(1:nvadim_mpilocal)

      ! Computation of background term of cost function:
      Jb = dot_product(zhat(1:nvadim_mpilocal),zhat(1:nvadim_mpilocal))/2.d0
      call mmpi_allreduce_sumreal8scalar(Jb,'GRID')

      vco_anl => col_getVco(columnTrlOnAnlIncLev_ptr)
      call gsv_allocate(statevector,tim_nstepobsinc, hco_anl, vco_anl, &
                        dataKind_opt=pre_incrReal, mpi_local_opt=.true.)

      call bmat_sqrtB(ahat_vhat,nvadim_mpilocal,statevector)

      call s2c_tl(statevector,column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)  ! put in column H_horiz dx
      call utl_tmg_start(10,'--Observations')
      call utl_tmg_start(18,'----ObsOper_TL')
      call oop_Htl(column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr,fso_nsim)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx

      call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  ! Save as OBS_WORK : R**-1/2 (Hdx)
      call utl_tmg_stop(18)
      call utl_tmg_stop(10)

      call cfn_calcJo(obsSpaceData_ptr)  ! Store J-obs in OBS_JOBS : 1/2 * R**-1 (Hdx)**2

      Jobs = 0.d0
      call utl_tmg_start(90,'--Minimization')
      call utl_tmg_start(92,'----SumCostFunction')
      call cfn_sumJo(obsSpaceData_ptr,Jobs)
      Jtotal = Jb + Jobs
      if (indic == 3) then
        Jtotal = Jobs
        if (mmpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  JO = ",G23.16,6X)') Jobs
      else
        Jtotal = Jb + Jobs
        if (mmpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  Jb = ",G23.16,6X,"JO = ",G23.16,6X,"Jt = ",G23.16)') Jb,Jobs,Jtotal
      end if
      call utl_tmg_stop(92)
      call utl_tmg_stop(90)

      call utl_tmg_start(10,'--Observations')
      call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  ! Modify OBS_WORK : R**-1 (Hdx)
      call utl_tmg_stop(10)

      call col_zero(column_ptr)

      call utl_tmg_start(10,'--Observations')
      call utl_tmg_start(19,'----ObsOper_AD')
      call oop_Had(column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)   ! Put in column : H_vert**T R**-1 (Hdx)
      call utl_tmg_stop(19)
      call utl_tmg_stop(10)

      call s2c_ad(statevector,column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)  ! Put in statevector H_horiz**T H_vert**T R**-1 (Hdx)

      gradJ(:) = 0.d0
      call bmat_sqrtBT(gradJ,nvadim_mpilocal,statevector)
      call gsv_deallocate(statevector)

      if (indic /= 3) then
        gradJ(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + gradJ(1:nvadim_mpilocal)
      end if
    end if

    call utl_tmg_start(90,'--Minimization')
    call utl_tmg_start(91,'----QuasiNewton')

    if (mmpi_myid == 0) write(*,*) 'simvar: Finished'

  end subroutine simvar

  !--------------------------------------------------------------------------
  ! prscal
  !--------------------------------------------------------------------------
  subroutine prscal(kdim,px,py,ddsc)
    !
    ! :Purpose: evaluation of the inner product in canonical space used in the minimization
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: kdim               ! dimension of the vectors
    real(8), intent(in)  :: px(kdim), py(kdim) ! vector components for which <px,py> is being calculated
    real(8), intent(out) :: ddsc               ! result of the inner product

    ! Locals:
    integer ::  cvIndex

    ddsc = 0.d0

    do cvIndex=1,nvadim_mpilocal
      ddsc = ddsc + px(cvIndex)*py(cvIndex)
    end do

    call mmpi_allreduce_sumreal8scalar(ddsc,'GRID')

  end subroutine prscal

  !--------------------------------------------------------------------------
  ! dcanab
  !--------------------------------------------------------------------------
  subroutine dcanab(kdim,py,px)
    !
    ! :Purpose: Change of variable associated with the canonical inner product
    !           to compute PX = L^-1 * Py with L related to the inner product
    !           <PX,PY> = PX^t  L^t  L PY
    !           (see the modulopt documentation aboutn DTCAB)
    !           Double precision version based on single precision CTCAB.
    !           Refered to  as dummy argument DTCAB by N1QN3 minimization
    !           package.
    !
    implicit none

    ! Arguments:
    integer, intent(in)     :: kdim     ! dimension of the vectors
    real(8), intent(inout)  :: px(kdim)
    real(8), intent(in)     :: py(kdim)

    ! Locals:
    integer :: cvIndex

    do cvIndex = 1, kdim
      px(cvIndex) = py(cvIndex)
    end do

  end subroutine dcanab

  !--------------------------------------------------------------------------
  ! dcanonb
  !--------------------------------------------------------------------------
  subroutine dcanonb(kdim,px,py)
    !
    ! :Purpose: Change of variable associated with the canonical inner product
    !           to compute PY = L * PX with L related to the inner product
    !           <PX,PY> = PX^t  L^t  L PY
    !           (see the modulopt documentation about DTONB)
    !           Double precision version based on single precision CTCAB.
    !           Refered to  as dummy argument DTCAB by N1QN3 minimization
    !           package.
    !
    implicit none

    ! Arguments:
    integer, intent(in)     :: kdim     ! dimension of the vectors
    real(8), intent(in)     :: px(kdim)
    real(8), intent(inout)  :: py(kdim)

    ! Locals:
    integer :: cvIndex

    do cvIndex = 1, kdim
      py(cvIndex) = px(cvIndex)
    end do

  end subroutine dcanonb

  !--------------------------------------------------------------------------
  ! multEnergyNorm
  !--------------------------------------------------------------------------
  subroutine multEnergyNorm(statevector_inout, statevector_ref,  &
                                latMin, latMax, lonMin, lonMax,      &
                                uvNorm,ttNorm,p0Norm,huNorm,tgNorm)
    !
    ! :Purpose: Computes energy norms
    !           For some positive definite symmetric matrix defining the energy,
    !             total energy = x^T * C * x
    !             statevector_inout = C * statevector_inout
    !           (Buehner, Du and Bedard, 2018)
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)  :: statevector_inout
    type(struct_gsv), intent(in)     :: statevector_ref
    real(8),          intent(in)     :: latMin, latMax, lonMin, lonMax
    logical,          intent(in)     :: uvNorm, ttNorm, p0Norm, huNorm, tgNorm

    ! Locals:
    integer              :: stepIndex, lonIndex, levIndex, latIndex, lonIndex2, latIndex2, status, nLev_M, nLev_T
    real(8)              :: scaleFactor, scaleFactorConst, scaleFactorLat, scaleFactorLon, scaleFactorLev
    real(8)              :: pfac, tfac, qfac
    real(8)              :: sumScale , sumeu, sumev, sumep, sumet, sumeq
    real(8), pointer     :: field_UU(:,:,:,:), field_VV(:,:,:,:), field_T(:,:,:,:), field_LQ(:,:,:,:)
    real(8), pointer     :: field_Psfc(:,:,:,:), field_TG(:,:,:,:),Psfc_ptr(:,:,:)
    real(8), pointer     :: Press_T(:,:,:)
    real(8), pointer     :: Press_M(:,:,:)
    real(8), allocatable :: Psfc_ref(:,:)
    real(8), parameter   :: T_r = 280.0D0
    real(8), parameter   :: Psfc_r = 100000.0D0 ! unit Pa
    real(8), parameter   :: sigma = 0.3 ! weight factor for humidity

    if (mmpi_myid == 0) write(*,*) 'multEnergyNorm: START'
    nullify(Press_T,Press_M)

    ! the factors for TT, HU and Ps (for wind is 1)
    tfac = mpc_cp_dry_air_r8/T_r                                 ! temperature factor (c_p/T_r)
    qfac = sigma*mpc_heat_condens_water_r8**2/(mpc_cp_dry_air_r8*T_r)  ! humidity factor ( (l_p*l_p)/(c_p*T_r) )
    pfac = mpc_rgas_dry_air_r8*T_r/(Psfc_r**2)                   ! surface pressure factor (R*T_r/Psfc_r^2)

    if (.not. gsv_isAllocated(statevector_inout)) then
      call utl_abort('multEnergyNorm: gridStateVector_inout not yet allocated')
    end if

    nLev_M = gsv_getNumLev(statevector_inout,'MM')
    nLev_T = gsv_getNumLev(statevector_inout,'TH')

    ! compute 3D log pressure fields
    call gsv_getField(statevector_ref,Psfc_ptr,'P0')
    allocate(Psfc_ref(statevector_inout%lonPerPEmax,statevector_inout%latPerPEmax))
    Psfc_ref(:,:) =  &
                  Psfc_ptr(statevector_inout%myLonBeg:statevector_inout%myLonEnd,  &
                  statevector_inout%myLatBeg:statevector_inout%myLatEnd, 1)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_T,  &
                        levels=Press_T,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_M,  &
                        levels=Press_M,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    ! dlat * dlon
    scaleFactorConst = statevector_inout%hco%dlat*statevector_inout%hco%dlon

    ! for wind components if to include in Norm calculation
    call gsv_getField(statevector_inout,field_UU,'UU')
    call gsv_getField(statevector_inout,field_VV,'VV')
    sumeu = 0.0D0
    sumev = 0.0D0
    sumScale = 0.0D0
    if (uvNorm) then
      do levIndex = 1, nLev_M
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if ( levIndex == nLev_M) then
                scaleFactorLev = Press_M(lonIndex2, latIndex2, nLev_M)-Press_T(lonIndex2, latIndex2, nLev_T-1)
              else if ( Press_T(lonIndex2, latIndex2, levIndex) < 10000.0D0) then
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_T(lonIndex2, latIndex2, levIndex+1) -  Press_T(lonIndex2, latIndex2, levIndex)
              end if

              scaleFactor = scaleFactorConst * scaleFactorLat* scaleFactorLon * scaleFactorLev
              sumScale = sumScale + scaleFactor

              sumeu = sumeu + &
                      0.5 * field_UU(lonIndex,latIndex,levIndex,stepIndex) * field_UU(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor
              sumev = sumev + &
                      0.5 * field_VV(lonIndex,latIndex,levIndex,stepIndex) * field_VV(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor

              field_UU(lonIndex,latIndex,levIndex,stepIndex) = &
                   field_UU(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor
              field_VV(lonIndex,latIndex,levIndex,stepIndex) = &
                   field_VV(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor
            end do !lonIndex
          end do !latIndex
        end do ! stepIndex
      end do ! levIndex
      call mmpi_allreduce_sumreal8scalar(sumeu,'grid')
      call mmpi_allreduce_sumreal8scalar(sumev,'grid')
      call mmpi_allreduce_sumreal8scalar(sumScale,'grid')

      sumeu = sumeu/sumScale
      sumev = sumev/sumScale

      field_UU(:,:,:,:) = field_UU(:,:,:,:)/sumScale
      field_VV(:,:,:,:) = field_VV(:,:,:,:)/sumScale
    else
      field_UU(:,:,:,:) = field_UU(:,:,:,:)*0.0D0
      field_VV(:,:,:,:) = field_VV(:,:,:,:)*0.0D0
    end if ! if uvNorm

    if (mmpi_myid == 0)  write(*,*) 'energy for UU=', sumeu
    if (mmpi_myid == 0)  write(*,*) 'energy for VV=', sumev

    ! for Temperature
    call gsv_getField(statevector_inout,field_T,'TT')
    sumScale = 0.0D0
    sumet = 0.0D0
    if (ttNorm) then
      do levIndex = 1, nLev_T
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if (levIndex == nLev_T) then  !surface
                scaleFactorLev =  Press_T(lonIndex2, latIndex2, nLev_T)-Press_T(lonIndex2, latIndex2, nLev_T-1)
              else if (levIndex == 1)  then  ! top
                scaleFactorLev = 0.0D0
              else if ( Press_M(lonIndex2, latIndex2, levIndex-1) < 10000.0D0) then
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_M(lonIndex2, latIndex2, levIndex ) - Press_M(lonIndex2, latIndex2, levIndex-1)
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon * scaleFactorLev
              sumet = sumet + &
                   0.5 * tfac * field_T(lonIndex,latIndex,levIndex,stepIndex) * field_T(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor
              sumScale = sumScale + scaleFactor
              field_T(lonIndex,latIndex,levIndex,stepIndex) = &
                           field_T(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * tfac * scaleFactor
            end do
          end do
        end do ! stepIndex
      end do ! levIndex
      call mmpi_allreduce_sumreal8scalar(sumet,'grid')
      call mmpi_allreduce_sumreal8scalar(sumScale,'grid')
      sumet = sumet/sumScale

      field_T(:,:,:,:) = field_T(:,:,:,:)/sumScale
    else
      field_T(:,:,:,:) = field_T(:,:,:,:)*0.0D0
    end if ! if ttNorm

    if (mmpi_myid == 0)  write(*,*) 'energy for TT=', sumet

    ! humidity (set to zero, for now)
    call gsv_getField(statevector_inout,field_LQ,'HU')
    sumScale = 0.0D0
    sumeq = 0.0D0
    if (huNorm) then
      do levIndex = 1, nLev_T
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if ( levIndex == nLev_T) then !surface
                scaleFactorLev =  Press_T(lonIndex2, latIndex2, nLev_T) - Press_T(lonIndex2, latIndex2, nLev_T-1)
              else if (levIndex == 1)  then  ! top
                scaleFactorLev = 0.0D0
              else if ( Press_M(lonIndex2, latIndex2, levIndex-1) < 10000.0D0) then
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_M(lonIndex2, latIndex2, levIndex ) - Press_M(lonIndex2, latIndex2, levIndex-1)
              end if

              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon * scaleFactorLev
              sumScale = sumScale + scaleFactor

              sumeq = sumeq + 0.5 * qfac * &
                    field_LQ(lonIndex,latIndex,levIndex,stepIndex) * field_LQ(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor

              field_LQ(lonIndex,latIndex,levIndex,stepIndex) = &
                       field_LQ(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor * qfac * 0.0

            end do
          end do
        end do ! stepIndex
      end do ! latIndex
      call mmpi_allreduce_sumreal8scalar(sumScale,'grid')
      call mmpi_allreduce_sumreal8scalar(sumeq,'grid')
      sumeq = sumeq/sumScale
      field_LQ(:,:,:,:) = field_LQ(:,:,:,:)/sumScale
    else
      field_LQ(:,:,:,:) = field_LQ(:,:,:,:)*0.0D0
    end if ! if huNorm

    if (mmpi_myid == 0)  write(*,*) 'energy for HU=', sumeq

    ! surface pressure
    call gsv_getField(statevector_inout,field_Psfc,'P0')
    sumScale = 0.0D0
    sumep = 0.0
    if (p0Norm) then
      do stepIndex = 1, statevector_inout%numStep
        do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon
              sumScale = sumScale + scaleFactor
              sumep = sumep + 0.5 * pfac * &
                  field_Psfc(lonIndex,latIndex,1,stepIndex) * field_Psfc(lonIndex,latIndex,1,stepIndex) * scaleFactor
              field_Psfc(lonIndex,latIndex,1,stepIndex) = &
              field_Psfc(lonIndex,latIndex,1,stepIndex) * 0.5 * scaleFactor * pfac
          end do
        end do ! latIndex
      end do ! stepIndex

      call mmpi_allreduce_sumreal8scalar(sumep,'grid')
      call mmpi_allreduce_sumreal8scalar(sumScale,'grid')
      sumep = sumep/sumScale
      field_Psfc(:,:,:,:) =  field_Psfc(:,:,:,:)/sumScale
    else
      field_Psfc(:,:,:,:) =  field_Psfc(:,:,:,:)*0.0D0
    end if ! if p0Norm

    if (mmpi_myid == 0)  write(*,*) 'energy for Ps=', sumep

    ! skin temperature (set to zero for now)
    call gsv_getField(statevector_inout,field_TG,'TG')
    sumScale = 0.0D0
    if (tgNorm) then
      do stepIndex = 1, statevector_inout%numStep
        do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon
              sumScale = sumScale + scaleFactor
              field_TG(lonIndex,latIndex,1,stepIndex) = &
              field_TG(lonIndex,latIndex,1,stepIndex) * 0.5 * scaleFactor * 0.0
          end do
        end do ! latIndex
      end do ! stepIndex
      call mmpi_allreduce_sumreal8scalar(sumScale,'grid')
      field_TG(:,:,:,:) = field_TG(:,:,:,:)/sumScale
    else
      field_TG(:,:,:,:) = field_TG(:,:,:,:)*0.0D0
    end if ! if tgNorm

    if (mmpi_myid == 0) write(*,*) 'energy for total=', sumeu + sumev + sumet + sumep + sumeq
    deallocate(Press_T,Press_M)
    deallocate(Psfc_ref)

    if (mmpi_myid == 0) write(*,*) 'multEnergyNorm: END'

  end subroutine multEnergyNorm

end module fsoi_mod
