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
!!  MODULE fso (Forecast sensitivity to observations.  prefix="fso")
!!
!!  *Purpose*:  FSO calculatio using a forecast ensemble
!!
!!  Subroutines (public):
!!    fso_setup
!!    fso_ensemble
!!
!!  Dependencies: 
!!
!--------------------------------------------------------------------------
module fso_mod
  use MathPhysConstants_mod
  use timeCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use controlVector_mod
  use mpivar_mod
  use horizontalCoord_mod
  use gridStateVector_mod
  use bmatrix_mod
  use bmatrixensemble_mod
  use tovs_nl_mod
  use stateToColumn_mod
  use analysisGrid_mod
  use utilities_mod
  use obsOperators_mod
  use costFunction_mod
  use quasinewton_mod
  implicit none
  save
  private

  ! public variables

  ! public procedures
  public              :: fso_setup, fso_ensemble

  type struct_dataptr
    type(struct_obs),pointer        :: obsSpaceData
    type(struct_columnData),pointer :: column
    type(struct_columnData),pointer :: columng
  end type struct_dataptr

  logical             :: initialized = .false.
  integer             :: nmtra, fso_nsim, nvadim_mpilocal
  integer             :: dataptr_int_size=0
  integer,external    :: get_max_rss
  real*8,allocatable  :: vhat(:)


  ! namelist variables
  integer             :: nvamaj, nitermax, nsimmax
  real(8)             :: leadTime, repsg, rdf1fac
  character(len=256)  :: forecastPath

  NAMELIST /NAMFSO/leadTime, nvamaj, nitermax, nsimmax
  NAMELIST /NAMFSO/repsg, rdf1fac, forecastPath

CONTAINS

!--------------------------------------------------------------------------
! FSO_setup
!--------------------------------------------------------------------------
  subroutine fso_setup

    implicit none

    integer :: ierr,nulnam
    integer :: fnom,fclos

    ! set default values for namelist variables
    leadtime = 12.0d0
    nvamaj = 6
    nitermax = 100
    nsimmax  = 120
    repsg    = 1d-5
    rdf1fac  = 0.25d0
    forecastPath = './forecasts'

    ! read in the namelist NAMFSO
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namfso,iostat=ierr)
    if(ierr.ne.0) call utl_abort('fso_setup: Error reading namelist')
    write(*,nml=namfso)
    ierr=fclos(nulnam)

    call ben_setFsoLeadTime(leadTime)
    fso_nsim=0 
    initialized=.true.

  end subroutine fso_setup

!--------------------------------------------------------------------------
! FSO_ensemble
!--------------------------------------------------------------------------
  subroutine fso_ensemble(columng,obsSpaceData)


    implicit none

    type(struct_columnData),target  :: columng
    type(struct_obs),target         :: obsSpaceData
    type(struct_columnData),target  :: column
    type(struct_dataptr)            :: dataptr 
    integer,allocatable             :: dataptr_int(:) ! obs array used to transmit pointer
    type(struct_gsv)                :: statevector_fa, statevector_fb, statevector_a, statevector_fso
    type(struct_gsv)                :: statevector_tempfa, statevector_tempfb !DPP
    type(struct_hco), pointer       :: hco_anl
    type(struct_vco), pointer       :: vco_anl
    integer                         :: nulout = 6
    integer                         :: impres, izs(1), iztrl(10)
    real                            :: zzsunused(1)
    real*8,allocatable              :: gradJ(:), ahat(:), vatra(:), zhat(:)
    integer                         :: imode, itermax, isimmax, indic
    real*8                          :: zjsp, zxmin, zdf1, zeps, dlgnorm, dlxnorm,  dlds(1)
    integer                         :: fnom,fclos,ierr
    integer                         :: dateStamp_fcst
    character(len=256)              :: fileName_fa, fileName_fb, fileName_a 
    character(len=8)                :: unitconversion

    !for Observation space 
    integer                         :: index_header, IDATA, IDATEND,  index_body 
    real*8                          :: fso_ori
  
    logical                         :: faExists

    !********************* some words should be here

   
      !hco_anl => agd_getHco('Analysis')
      hco_anl => agd_getHco('ComputationalGrid')
      vco_anl => col_getVco(columng)

      nvadim_mpilocal = cvm_nvadim
      nmtra = (4 + 2*nvamaj)*nvadim_mpilocal
      write(*,9401) nvamaj,nmtra
  9401 format(4X,'NVAMAJ = ',I3,/5X,"NMTRA =",1X,I14)

     
      ! initialize column object for storing "increment"
      call col_setVco(column,col_getVco(columng))
      call col_allocate(column,col_getNumCol(columng),mpi_local=.true.)
      call col_copyLatLon(columng,column)

      write(*,*) 'PRDATABIN: For 4D increment'
      call tim_sutimeinterp(obsSpaceData)

      ! compute dateStamp_fcst
      call incdatr(dateStamp_fcst, tim_getDatestamp(), leadTime)
      write(*,*) 'fso_ensemble: analysis datestamp = ',tim_getDatestamp()
      write(*,*) 'fso_ensemble: forecast datestamp = ',dateStamp_fcst

      ! allocate control vector related arrays (these are all mpilocal)
      allocate(ahat(nvadim_mpilocal))
      allocate(vhat(nvadim_mpilocal))
      allocate(zhat(nvadim_mpilocal)) 
      allocate(gradJ(nvadim_mpilocal))
      allocate(vatra(nmtra))

      ! initialize control vector related arrays to zero
      ahat(:)=0.0d0
      vhat(:)=0.0d0
      zhat(:)=0.0d0
      gradJ(:)=0.0d0
      vatra(:)=0.0d0

     
      ! ------------------------------------------------------
      ! Compute vhat
      ! ------------------------------------------------------

      ! read forecasts from the analysis and background state
      fileName_fa = trim(forecastPath) // '/forecast_a'
      inquire(file=trim(fileName_fa),exist=faExists)
      write(*,*) 'faExists', faExists
      write(*,*) 'DPP', fileName_fa
      call gsv_allocate(statevector_fa, 1, hco_anl, vco_anl, &
                        datestamp=datestamp_fcst, mpi_local=.true.)
      call gsv_readFromFile(statevector_fa, fileName_fa, ' ', 'P')

      !DPP for statevector_tempfa and statevector_tempfb
      !for statevecotr_tempfa
      call gsv_allocate(statevector_tempfa, 1, hco_anl, vco_anl, &
                        datestamp=datestamp_fcst, mpi_local=.true.)

      !for statevecotr_tempfb
      call gsv_allocate(statevector_tempfb, 1, hco_anl, vco_anl, &
                        datestamp=datestamp_fcst, mpi_local=.true.)

      fileName_fb = trim(forecastPath) // '/forecast_b'
      call gsv_allocate(statevector_fb, 1, hco_anl, vco_anl, &
                        datestamp=datestamp_fcst, mpi_local=.true.)
      call gsv_readFromFile(statevector_fb, fileName_fb, ' ', 'P')

      ! read verifying analysis
      fileName_a = trim(forecastPath) // '/analysis'
      call gsv_allocate(statevector_a, 1,hco_anl, vco_anl, &
                        datestamp=datestamp_fcst, mpi_local=.true.)
      call gsv_readFromFile(statevector_a, fileName_a, ' ', 'A')


      ! compute error of both forecasts (overwrite forecasts with error)
      call gsv_add(statevector_a, statevector_fa, -1.0d0)
      call gsv_add(statevector_a, statevector_fb, -1.0d0)


  !DPP
      call gsv_copy(statevector_fa,statevector_tempfa)
      call gsv_copy(statevector_fb,statevector_tempfb)
      call gsv_multEnergyNorm(statevector_tempfa, statevector_a) ! use analysis as reference state
      call gsv_multEnergyNorm(statevector_tempfb, statevector_a) ! use analysis as reference state


      ! compute error Norm =  C * (error_t^fa + error_t^fb)
      call gsv_add(statevector_fa, statevector_fb, 1.0d0)
      call gsv_multEnergyNorm(statevector_fb, statevector_a) ! use analysis as reference state

      ! compute vhat = B_t^T/2 * C * (error_t^fa + error_t^fb)  
      call bmat_sqrtBT(vhat, nvadim_mpilocal, statevector_fb, useForecast = .true.) 

     
      write(*,*) 'DPP-vhat',maxval(vhat),minval(vhat) 
   
 
      ! ------------------------------------------------------
      ! Compute zhat by performing variational minimization
      ! ------------------------------------------------------    

      ! Set-up for the minimization
      if(mpi_myid.eq.0) then
       impres=5
      else 
       impres=0
      endif

      ! recast pointer to obsSpaceData as an integer array, so it can be passed through n1qn3 to simvar
      dataptr%obsSpaceData => obsSpaceData
      dataptr%column       => column
      dataptr%columng      => columng
      dataptr_int_size = size(transfer(dataptr,dataptr_int))
      allocate(dataptr_int(dataptr_int_size))
      dataptr_int(1:dataptr_int_size)=transfer(dataptr,dataptr_int)

      imode = 0
      zeps = repsg
      itermax = nitermax
      isimmax = nsimmax
      zxmin = epsilon(zxmin)

      ! initial gradient calculation
      indic = 2
      call simvar(indic,nvadim_mpilocal,zhat,zjsp,gradJ,dataptr_int(1))
      zdf1 =  rdf1fac * abs(zjsp)

      ! print amplitude of initial gradient and cost function value
      call prscal(nvadim_mpilocal,gradJ,gradJ,dlgnorm)
      dlgnorm = dsqrt(dlgnorm)
      call prscal(nvadim_mpilocal,zhat,zhat,dlxnorm)
      dlxnorm = dsqrt(dlxnorm)
      write(*,*)' |X| = ', dlxnorm
      write(*,fmt=9220) zjsp, dlgnorm
 9220 format(/4X,'J(X) = ',G23.16,4X,'|Grad J(X)| = ',G23.16)

      write(*,fmt=9320)zxmin,zdf1,zeps,impres,itermax,isimmax
 9320 format(//,10X,' Minimization QNA_N1QN3 starts ...',/  &
             10x,'DXMIN =',G23.16,2X,'DF1 =',G23.16,2X,'EPSG =',G23.16  &
             /,10X,'IMPRES =',I3,2X,'NITER = ',I3,2X,'NSIM = ',I3)

      ! Do the minimization
      call tmg_start(70,'QN')
      call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim_mpilocal, zhat,  &
          zjsp, gradJ, zxmin, zdf1, zeps, impres, nulout, imode,   &
          itermax,isimmax, iztrl, vatra, nmtra, dataptr_int(1), zzsunused,   &
          dlds)
      call tmg_stop(70)
      call fool_optimizer(obsSpaceData)

      write(*,FMT=9500) imode,itermax,isimmax
 9500 FORMAT(//,20X,20('*'),2X    &
        ,/,20X,'              Minimization ended with MODE:',I4  &
        ,/,20X,'                Total number of iterations:',I4  &
        ,/,20X,'               Total number of simulations:',I4)

      ! Compute ahat = zhat + vhat
      ahat = zhat + vhat

      ! ------------------------------------------------------
      ! Compute yhat = [R^-1 H B^1/2 ahat], and put in OBS_FSO
      ! ------------------------------------------------------    

      call gsv_allocate(statevector_fso, tim_nstepobsinc, hco_anl, vco_anl, &
                        datestamp=tim_getDatestamp(), mpi_local=.true.)

      call bmat_sqrtB(ahat, nvadim_mpilocal, statevector_fso) 
      call s2c_tl(statevector_fso,column,columng,obsSpaceData)  ! put in column H_horiz B^1/2 ahat
      call oop_Htl(column,columng,obsSpaceData,1) !DPP  ! Save as OBS_WORK: H_vert H_horiz B^1/2 vhat = H B^1/2 ahat
      call cfn_RsqrtInverse(obsSpaceData,OBS_FSO,OBS_WORK) ! Save as OBS_FSO : R**-1/2 H B^1/2 ahat
      call cfn_RsqrtInverse(obsSpaceData,OBS_FSO,OBS_FSO)  ! Save as OBS_FSO : R**-1 H B^1/2 ahat\

      ! Due to the very small value of FSO, here it is enlarged by 1e6 therefore in the script file to extract FSO it should be divided by 1e6
      do index_header =1, obs_numHeader(obsSpaceData)
        IDATA   = obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
        IDATEND = obs_headElem_i(obsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
        do index_body=idata,idatend
           fso_ori= obs_bodyElem_r(obsSpaceData,OBS_FSO,index_body)
           call obs_bodySet_r(obsSpaceData,OBS_FSO,index_body, fso_ori*1e6)
        enddo
      enddo


      ! deallocate the control vector related arrays
      deallocate(ahat)
      deallocate(gradJ)
      deallocate(vatra)
      deallocate(vhat)
      deallocate(zhat)
      deallocate(dataptr_int)
      call col_deallocate(column)

  end subroutine fso_ensemble

  subroutine simvar(indic,nvadim,zhat,Jtotal,gradJ,dataptr_int)
    implicit none
    ! Argument declarations
    integer :: nvadim ! Dimension of the control vector in forecast error coraviances space
    ! Value of indic
    ! Note: 1 and 4 are reserved values for call back from m1qn3.
    !       For direct calls use other value than 1 and 4.
    ! =1 No action taken; =4 Both J(u) and its gradient are computed.
    ! =2 Same as 4 (compute J and gradJ) but do not interrupt timer of the
    !    minimizer.
    ! =3 Compute Jo and gradJo only.
    integer :: indic 
    real*8  :: Jtotal ! Cost function of the Variational algorithm
    real*8, dimension(nvadim) :: gradJ ! Gradient of the Variational Cost funtion
    !!real*8, dimension(nvadim) :: ahat ! Control variable in forecast error covariances space
    real*8, dimension(nvadim) :: zhat ! Control variable in forecast error covariances space
    integer :: dataptr_int(dataptr_int_size)  ! integer work area used to transmit a pointer to the obsSpaceData
    !
    ! Purpose: Implement the Variational solver as described in
    ! Courtier, 1997, Dual formulation of four-dimentional variational assimilation,
    ! Q.J.R., pp2449-2461.
    !
    ! Author : Simon Pellerin *ARMA/MSC October 2005
    !          (Based on previous versions of evaljo.ftn, evaljg.ftn and evaljgns.ftn).
    !
    ! Local declaration
    real*8 :: ahat_vhat(nvadim)
    real*8 :: Jb, Jobs
    type(struct_gsv) :: statevector
    type(struct_dataptr) :: dataptr
    type(struct_obs),pointer :: obsSpaceData
    type(struct_columnData),pointer :: column,columng
    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl


    ! Convert the integer array dataptr_int back into a pointer to the obsSpaceData
    dataptr=transfer(dataptr_int(1:dataptr_int_size),dataptr)
    obsSpaceData => dataptr%obsSpaceData
    column       => dataptr%column
    columng      => dataptr%columng

    if (indic .eq. 1 .or. indic .eq. 4) call tmg_stop(70)

    call tmg_start(80,'MIN_SIMVAR')
    if (indic .ne. 1) then ! No action taken if indic == 1
       fso_nsim = fso_nsim + 1

       if(mpi_myid == 0) then
         write(*,*) 'Entering simvar for simulation ',fso_nsim
         write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
         call flush(6)
       endif

       ! note: vhat = B_t^T/2 hat(del x_t)
        ahat_vhat(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + vhat(1:nvadim_mpilocal)     

       ! Computation of background term of cost function:
       Jb = dot_product(zhat(1:nvadim_mpilocal),zhat(1:nvadim_mpilocal))/2.d0  
       call tmg_start(89,'MIN_COMM')       
       call mpi_allreduce_sumreal8scalar(Jb,"GRID")
       call tmg_stop(89)

             
       hco_anl => agd_getHco('ComputationalGrid')
       vco_anl => col_getVco(columng)
       call gsv_allocate(statevector,tim_nstepobsinc, hco_anl, vco_anl, &
                         mpi_local=.true.)

       call bmat_sqrtB(ahat_vhat,nvadim_mpilocal,statevector)

       call tmg_start(30,'OBS_INTERP')
       call s2c_tl(statevector,column,columng,obsSpaceData)  ! put in column H_horiz dx
       call tmg_stop(30)

       call tmg_start(40,'OBS_TL')
       call oop_Htl(column,columng,obsSpaceData,fso_nsim)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
       call tmg_stop(40)

       call cfn_RsqrtInverse(obsSpaceData,OBS_WORK,OBS_WORK)  ! Save as OBS_WORK : R**-1/2 (Hdx)
     
       call cfn_calcJo(obsSpaceData)  ! Store J-obs in OBS_JOBS : 1/2 * R**-1 (Hdx)**2
     
       Jobs = 0.d0
       call cfn_sumJo(obsSpaceData,Jobs)
       Jtotal = Jb + Jobs
       if (indic .eq. 3) then
          Jtotal = Jobs
          IF(mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  JO = ",G23.16,6X)') Jobs
       else
          Jtotal = Jb + Jobs
          IF(mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  Jb = ",G23.16,6X,"JO = ",G23.16,6X,"Jt = ",G23.16)') Jb,Jobs,Jtotal
       endif

       call cfn_RsqrtInverse(obsSpaceData,OBS_WORK,OBS_WORK)  ! Modify OBS_WORK : R**-1 (Hdx)

       call col_zero(column)

       call tmg_start(41,'OBS_AD')
       call oop_Had(column,columng,obsSpaceData)   ! Put in column : H_vert**T R**-1 (Hdx)
       call tmg_stop(41)

       call tmg_start(31,'OBS_INTERPAD')
       call s2c_ad(statevector,column,columng,obsSpaceData)  ! Put in statevector H_horiz**T H_vert**T R**-1 (Hdx)
       call tmg_stop(31)

       gradJ(:) = 0.d0
       call bmat_sqrtBT(gradJ,nvadim_mpilocal,statevector)
       call gsv_deallocate(statevector)

       if (indic .ne. 3) then
          gradJ(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + gradJ(1:nvadim_mpilocal)
       endif
    endif
    call tmg_stop(80)
    if (indic .eq. 1 .or. indic .eq. 4) call tmg_start(70,'QN')

    if(mpi_myid.eq.0) write(*,*) 'end of simvar'

  end subroutine simvar


  SUBROUTINE DSCALQN(KDIM,PX,PY,DDSC,KZS, PZS, DDZS)
    !***s/r DSCALQN: inner product in canonical space
    !*    ------------------- 
    !**    Purpose: interface for the inner product to be used
    !*     .        by the minimization subroutines N1QN3.
    !*
    !*Arguments
    !*     i : KDIM      : dimension of the vectors
    !*     i : PX, PY    : vector for which <PX,PY> is being calculated
    !*     o : DDSC      : result of the inner product
    !*     --------------
    !*     i :  KZS(1)   : unused working space for INTEGER  (not used)
    !*     i :  PZS(1)   : unused working space for REAL     (not used)
    !*     i : PDZS(1)   : unused working space for REAL*8   (not used)
    IMPLICIT NONE

    REAL PZS(1)
    INTEGER KZS(1)
    REAL*8  DDZS(1)

    INTEGER KDIM
    REAL*8 PX(KDIM), PY(KDIM)
    REAL*8 DDSC

    CALL PRSCAL(KDIM,PX,PY,DDSC)
    RETURN
  END SUBROUTINE DSCALQN


  SUBROUTINE PRSCAL(KDIM,PX,PY,DDSC)
    !***s/r PRSCAL: inner product in canonical space
    !*
    !*Author  : P. Gauthier *ARMA/AES  January 27, 1993
    !**    Purpose: evaluation of the inner product used in the
    !*     .        minimization
    !*
    !*Arguments
    !*     i : KDIM     : dimension of the vectors
    !*     i : PX, PY   : vector for which <PX,PY> is being calculated
    !*     o : DDSC     : result of the inner product
    !*
    !* Implicit argument: SCALP(KDIM) assumed to be unity

    IMPLICIT NONE

    INTEGER KDIM, J, RR
    REAL*8 PX(KDIM), PY(KDIM)
    REAL*8 DDSC
    REAL*8 partialsum(128)
    INTEGER mythread,numthreads,jstart,jend
    INTEGER omp_get_thread_num,omp_get_num_threads

    call tmg_start(71,'QN_PRSCAL')
    DDSC = 0.D0

    do j=1,nvadim_mpilocal
      DDSC = DDSC + PX(J)*PY(J)
    ENDDO

    call tmg_start(79,'QN_COMM')
    call mpi_allreduce_sumreal8scalar(ddsc,"GRID")
    call tmg_stop(79)

    call tmg_stop(71)

    RETURN
  END SUBROUTINE PRSCAL


  SUBROUTINE DCANAB(KDIM,PY,PX,KZS,PZS,PDZS)
    !***s/r DCANAB  - Change of variable associated with the canonical
    !*     .         inner product
    !*
    !*Author    JM Belanger CMDA/SMC   May 2001
    !*     .    Double precision version based on single precision CTCAB.
    !*          Refered to  as dummy argument DTCAB by N1QN3 minimization
    !*          package.
    !*    -------------------
    !**    Purpose: to compute PX = L^-1 * Py with L related to the inner product
    !*     .        <PX,PY> = PX^t  L^t  L PY
    !*     .        (see the modulopt documentation aboutn DTCAB)
    !*     NOTE: L is assumed to be the identity!
    IMPLICIT NONE

    INTEGER KDIM, KZS(1)
    REAL PZS(1)
    REAL*8 PX(KDIM), PY(KDIM)
    REAL*8 PDZS(1)

    INTEGER JDIM

    DO JDIM = 1, KDIM
      PX(JDIM) = PY(JDIM)
    ENDDO

    RETURN
  END SUBROUTINE DCANAB


  SUBROUTINE DCANONB(KDIM,PX,PY,KZS,PZS,PDZS)
    !***s/r DCANONB  - Change of variable associated with the canonical
    !*     .          inner product
    !*
    !*Author    JM Belanger CMDA/SMC  May 2001
    !*     .    Double precision version based on single precision CANONB.
    !*          Refered to as dummy argument DTONB by N1QN3 minimization
    !*          package.
    !*    -------------------
    !**    Purpose: to compute PY = L * PX with L related to the inner product
    !*     .        <PX,PY> = PX^t  L^t  L PY
    !*     .        (see the modulopt documentation about DTONB)
    !*     .

    IMPLICIT NONE
    INTEGER KDIM, KZS(1)
    REAL PZS(1)
    REAL*8 PX(KDIM), PY(KDIM)
    REAL*8 PDZS(1)

    INTEGER JDIM

    DO JDIM = 1, KDIM
      PY(JDIM) = PX(JDIM)
    ENDDO

    RETURN
  END SUBROUTINE DCANONB

end module fso_mod
