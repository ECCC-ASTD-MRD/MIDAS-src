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
!! MODULE bMatrix (prefix="bmat")
!!
!! *Purpose*: A higher-level module that takes care of calling subroutines 
!!            in the lower-level modules bmatrixHI/lambmatrixHI and bmatrixEnsemble.
!!
!! Revisions: Pind Du, CMDA, 2014
!!v            - Added 'use bMatrixChem_mod'
!!v            Y. Rochon, ARQI. July 2015
!!v            - Added 'public bchm_getScaleFactor' for constituents. 
!!
!! Comments:
!!v            - Considerations for ensemble-based and regional static covariances
!!v              for constituents not yet included.
!--------------------------------------------------------------------------
MODULE BMatrix_mod

  use mpivar_mod
  use bMatrixHI_mod
  use bMatrixEnsemble_mod
  use bMatrixChem_mod
  use controlVector_mod
  use gridStateVector_mod
  use LAMbMatrixHI_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use globalSpectralTransform_mod
  use utilities_mod
  implicit none
  save
  private
  
  ! public procedures
  public :: bmat_setup, bmat_finalize, bmat_sqrtB, bmat_sqrtBT
  public :: bmat_reduceToMPILocal, bmat_reduceToMPILocal_r4, bmat_expandToMPIGlobal, bmat_expandToMPIGlobal_r4
  ! public procedures through inheritance
  public :: bhi_getScaleFactor,bhi_truncateCV,ben_getScaleFactor,ben_getnEns,ben_getPerturbation
  public :: bchm_getScaleFactor,bchm_truncateCV
  public :: ben_setFsoLeadTime


  type(struct_hco), pointer :: hco_anl

contains
  
!--------------------------------------------------------------------------
! bmat_setup
!--------------------------------------------------------------------------
  SUBROUTINE bmat_setup(hco_anl_in, vco_anl_in)
    !
    !- bmat_setup - Initializes the analysis Background term for the 
    !               specific analysis configuration used.
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvdimchm and 
    !             section 2.3)
    !
    ! Comments:
    !
    IMPLICIT NONE

    type(struct_vco), pointer :: vco_anl_in
    type(struct_hco), pointer :: hco_anl_in

    integer :: cvdimens, cvdimhi
    integer :: get_max_rss
    integer :: cvdimchm
   
    !
    !- 1.  Get/Check the analysis grid info
    !

    !- 1.1 Horizontal Grid info
    hco_anl => hco_anl_in

    !
    !- 2.  Set the B matrices
    !
    cvdimhi  = 0
    cvdimens = 0
    cvdimchm = 0

    !- 2.1 Time-Mean Homogeneous and Isotropic...
    if ( hco_anl%global ) then
      write(*,*)
      write(*,*) 'Setting up the modular GLOBAL HI covariances...'
      call bhi_Setup( hco_anl, vco_anl_in, & ! IN
                      cvdimhi )              ! OUT
    else
      write(*,*)
      write(*,*) 'Setting up the modular LAM HI covariances...'
      call lbhi_Setup( hco_anl, vco_anl_in, & ! IN
                       cvdimhi )              ! OUT
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'Dimension of HI  control vector returned:',cvdimhi

    !- 2.2 Flow-dependent Ensemble-Based
    write(*,*)
    write(*,*) 'Setting up the modular ENSEMBLE covariances...'
    call ben_Setup( hco_anl,             & ! IN
                    vco_anl_in,          & ! IN
                    cvdimens )             ! OUT

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'Dimension of ENS control vector returned:',cvdimens

    !-2.3  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
    if ( hco_anl % global ) then
      write(*,*)
      write(*,*) 'Setting up the modular GLOBAL HI-chm covariances...'
      call bchm_Setup( hco_anl, vco_anl_in, & ! IN
                       cvdimchm )             ! OUT

    else
      ! Done in lbhi_Setup 
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'Dimension of CH static control vector returned:',cvdimchm

    !
    !- 3.  Setup the control vector
    !
    call cvm_Setup( cvdimhi, cvdimens, cvdimchm ) ! IN

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'Dimension of TOTAL control vector:',cvm_nvadim
    
  END SUBROUTINE bmat_setup

!--------------------------------------------------------------------------
! bmat_sqrtB
!-------------------------------------------------------------------------- 
  SUBROUTINE bmat_sqrtB(controlVector,cvdim,statevector,useForecast_opt)
    implicit none
    !
    !- Purpose: Transforms model state from error covariance space
    !           to grid point space.
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvBchm and 
    !             section 2.3)
    !
    ! Comments:
    !
    ! LAM and Ensemble cases not done for constituents. TBD 
    !
    integer         :: cvdim
    real(8)         :: controlVector(cvdim)
    real(8),pointer :: cvBhi(:), cvBen(:), field(:,:,:), field4d(:,:,:,:)
    real(8),pointer :: cvBchm(:)
    logical,optional :: useForecast_opt

    type(struct_gsv) :: statevector, statevector_temp

    !
    !- 1.  Set analysis increment to zero
    !
    call gsv_zero(statevector)

    !
    !- 2.  Compute the analysis increment
    !

    !- 2.1 Allocate and set to zero another temporary statevector
    call gsv_allocate(statevector_temp, statevector%numStep, hco_anl, gsv_getVco(statevector), &
                      mpi_local=.true.)
    call gsv_zero(statevector_temp)

    !- 2.2 Compute 3D contribution to increment from BmatrixHI
    call tmg_start(50,'B_HI')
    if ( cvm_subVectorExists(cvm_BHI) ) then
      cvBhi => cvm_getSubVector(controlVector,cvm_BHI)
      if ( statevector%hco%global ) then
        !- 2.2.1 Global mode
        call bhi_bsqrt( cvBhi,      & ! IN
                        statevector ) ! OUT
      else
        !- 2.2.2 LAM mode
        call lbhi_bSqrt( cvBhi,      & ! IN
                         statevector ) ! OUT
      end if
    end if
    call tmg_stop(50)

    !- 2.3  Compute 3D contribution to increment from BmatrixChem
    call tmg_start(123,'B_CHM')
    if ( cvm_subVectorExists(cvm_BCHM) ) then
      cvBchm => cvm_getSubVector(controlVector,cvm_BCHM)
      if ( statevector % hco % global ) then
        !- 2.3.1 Global mode
        call bchm_bsqrt( cvBchm,      & ! IN
                        statevector )   ! OUT
      else
        !- 2.3.2 LAM mode
        call utl_abort('bmat_sqrtB: local routine currently unavailable for chemical constituents, to be available via lbhi_bsqrt in the future')
      end if
    end if
    call tmg_stop(123)

    !- 2.4 copy 3D increment to other timesteps to create 4D increment
    if ( cvm_subVectorExists(cvm_BHI) .or. cvm_subVectorExists(cvm_BCHM) ) call gsv_3dto4d(statevector)

    !- 2.5 compute 4D contribution to increment from BmatrixEnsemble
    call tmg_start(60,'B_ENS')
    if ( cvm_subVectorExists(cvm_BEN) ) then
      cvBen => cvm_getSubVector(controlVector,cvm_BEN)
      if( present(useForecast_opt) ) then
        call ben_bsqrt(cvBen, statevector_temp, useForecast_opt)
      else
        call ben_bsqrt(cvBen, statevector_temp)
      end if
    end if
    call tmg_stop(60)

    !- 2.6 Add the two contributions together, result in statevector
    call gsv_add(statevector_temp,statevector)

    call gsv_deallocate(statevector_temp)

  END SUBROUTINE bmat_sqrtB

!--------------------------------------------------------------------------
! bmat_sqrtBT
!--------------------------------------------------------------------------
  SUBROUTINE bmat_sqrtBT(controlVector,cvdim,statevector,useForecast_opt)
    implicit none
    !
    !- Purpose: Transforms model state from grid point space 
    !           to error covariance space.
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvBchm and
    !             section 1.5)
    !
    ! Comments:
    !
    integer :: cvdim
    real(8) :: controlVector(cvdim)
    real(8),pointer :: cvBhi(:),cvBen(:),cvBchm(:)
    type(struct_gsv) :: statevector
    logical,optional :: useForecast_opt

    !- 1.1 set gradient to zero
    controlVector(:)=0.0d0
    
    !- 1.2 Add contribution to gradient from BmatrixEnsemble
    call tmg_start(61,'B_ENS_T')
    if ( cvm_subVectorExists(cvm_BEN) ) then
      cvBen=>cvm_getSubVector(controlVector,cvm_BEN)
      if ( present(useForecast_opt) ) then
        call ben_bsqrtad(statevector,cvBen,useForecast_opt)
      else
        call ben_bsqrtad(statevector,cvBen)
      end if
    end if
    call tmg_stop(61)

    !- 1.3 adjoint of copy 3D increment to 4D increment
    if ( cvm_subVectorExists(cvm_BHI) .or. cvm_subVectorExists(cvm_BCHM)) call gsv_3dto4dAdj(statevector)

    !- 1.4 add contribution to gradient from BmatrixChem
    call tmg_start(124,'B_CHM_T')
    if ( cvm_subVectorExists(cvm_BCHM) ) then
      cvBchm=>cvm_getSubVector(controlVector,cvm_BCHM)
      if ( statevector%hco%global ) then
        !- 1.4.1 add contribution to gradient from GLOBAL BmatrixChem
        call bchm_bsqrtad( statevector, & ! IN
                          cvBchm )        ! OUT
      else
        !- 1.4.2 add contribution to gradient from LAM BmatrixChem
        call utl_abort('bmat_sqrtBT: local routine currently unavailable for chemical constituents, to be available via lbhi_bSqrtAdj in the future')
      end if
    end if
    call tmg_stop(124)

    !- 1.5 add contribution to gradient from BmatrixHI
    call tmg_start(51,'B_HI_T')
    if ( cvm_subVectorExists(cvm_BHI) ) then
      cvBhi=>cvm_getSubVector(controlVector,cvm_BHI)
      if ( statevector%hco%global ) then
        !- 1.5.1 add contribution to gradient from GLOBAL BmatrixHI
        call bhi_bsqrtad( statevector, & ! IN
                          cvBhi )        ! OUT
      else
        !- 1.5.2 add contribution to gradient from LAM BmatrixHI
        call lbhi_bSqrtAdj( statevector, & ! IN
                            cvBhi )        ! OUT
      end if
    end if
    call tmg_stop(51)

  END SUBROUTINE bmat_sqrtBT

!--------------------------------------------------------------------------
! bmat_finalize
!--------------------------------------------------------------------------
  SUBROUTINE bmat_finalize()
    !
    !- Purpose: Releases memory used by B matrices.
    !
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvBchm)
    !
    ! Comments:
    !
    ! - LAM and ensemble components for constituents tbd.
    !
    implicit none    

    call bhi_finalize()
    call ben_finalize()
    call bchm_finalize()
    call lbhi_finalize()

  END SUBROUTINE bmat_finalize

!--------------------------------------------------------------------------
! BMAT_reduceToMPILocal
!--------------------------------------------------------------------------
  SUBROUTINE BMAT_reduceToMPILocal(cv_mpilocal,cv_mpiglobal,cvDim_mpilocal_out)
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, Dec 2014
    !           - Additions for chemical constituents (see cv*Bchm*)
    !
    ! Comments:
    !
    ! - LAM and ensemble components for constituents tbd.
    !
    implicit none
    real(8), intent(out) :: cv_mpilocal(:)
    real(8), intent(in)  :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpilocal_out

    integer :: cvDim_Bhi_mpilocal,cvDim_Ben_mpilocal,cvDim_Bchm_mpilocal

    real(8),pointer :: cvBhi_mpilocal(:) ,cvBen_mpilocal(:),cvBchm_mpilocal(:)
    real(8),pointer :: cvBhi_mpiglobal(:),cvBen_mpiglobal(:),cvBchm_mpiglobal(:)


    cvDim_Bhi_mpilocal = 0
    if(cvm_subVectorExists(cvm_BHI)) then
      cvBhi_mpilocal =>cvm_getSubVector(cv_mpilocal,cvm_BHI)
      if (mpi_myid == 0) then
         cvBhi_mpiglobal=>cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BHI)
      else
         cvBhi_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then 
         call bhi_reduceToMPILocal (cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpilocal)
      else
         call lbhi_reduceToMPILocal(cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpilocal)
      end if
    end if

    cvDim_Ben_mpilocal = 0
    if(cvm_subVectorExists(cvm_BEN)) then
       cvBen_mpilocal => cvm_getSubVector(cv_mpilocal,cvm_BEN)
       if (mpi_myid == 0) then
          cvBen_mpiglobal => cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BEN)
       else
          cvBen_mpiglobal => null()
       end if
       call ben_reduceToMPILocal(cvBen_mpilocal,cvBen_mpiglobal,cvDim_Ben_mpilocal)
    end if

    cvDim_Bchm_mpilocal = 0
    if(cvm_subVectorExists(cvm_BCHM)) then
      cvBchm_mpilocal => cvm_getSubVector(cv_mpilocal,cvm_BCHM)
      if (mpi_myid == 0) then
         cvBchm_mpiglobal => cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BCHM)
      else
         cvBchm_mpiglobal => null()
      end if
      if ( hco_anl%global ) then 
         call bchm_reduceToMPILocal(cvBchm_mpilocal,cvBchm_mpiglobal,cvDim_Bchm_mpilocal)
      else
!         Done in lbhi_reducetoMPILocal
      end if
    end if


    cvDim_mpilocal_out = cvDim_Bhi_mpilocal + cvDim_Ben_mpilocal + &
                         cvDim_Bchm_mpilocal

  END SUBROUTINE BMAT_reduceToMPILocal

!--------------------------------------------------------------------------
! BMAT_reduceToMPILocal_r4
!--------------------------------------------------------------------------
  SUBROUTINE BMAT_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal,cvDim_mpilocal_out)
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, Dec 2014
    !           - Additions for chemical constituents (see cv*Bchm*)
    !
    ! Comments:
    !
    ! - LAM and ensemble components for constituents tbd.
    !
    implicit none
    real(4), intent(out) :: cv_mpilocal(:)
    real(4), intent(in)  :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpilocal_out

    integer :: cvDim_Bhi_mpilocal,cvDim_Ben_mpilocal,cvDim_Bchm_mpilocal

    real(4),pointer :: cvBhi_mpilocal(:) ,cvBen_mpilocal(:),cvBchm_mpilocal(:)
    real(4),pointer :: cvBhi_mpiglobal(:),cvBen_mpiglobal(:),cvBchm_mpiglobal(:)

    cvDim_Bhi_mpilocal = 0
    if(cvm_subVectorExists(cvm_BHI)) then
      cvBhi_mpilocal =>cvm_getSubVector_r4(cv_mpilocal,cvm_BHI)
      if (mpi_myid == 0) then
         cvBhi_mpiglobal=>cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BHI)
      else
         cvBhi_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then 
         call bhi_reduceToMPILocal_r4 (cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpilocal)
      else
         call lbhi_reduceToMPILocal_r4(cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpilocal)
      end if
    end if

    cvDim_Ben_mpilocal = 0
    if(cvm_subVectorExists(cvm_BEN)) then
       cvBen_mpilocal =>cvm_getSubVector_r4(cv_mpilocal,cvm_BEN)
       if (mpi_myid == 0) then
          cvBen_mpiglobal=>cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BEN)
       else
          cvBen_mpiglobal=>null()
       end if
       call ben_reduceToMPILocal_r4(cvBen_mpilocal,cvBen_mpiglobal,cvDim_Ben_mpilocal)
    end if

    cvDim_Bchm_mpilocal = 0
    if(cvm_subVectorExists(cvm_BCHM)) then
      cvBchm_mpilocal =>cvm_getSubVector_r4(cv_mpilocal,cvm_BCHM)
      if (mpi_myid == 0) then
         cvBchm_mpiglobal=>cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BCHM)
      else
         cvBchm_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then 
         call bchm_reduceToMPILocal_r4(cvBchm_mpilocal,cvBchm_mpiglobal,cvDim_Bchm_mpilocal)
      else
!         Done in lbhi_reducetoMPILocal
      end if
    end if


    cvDim_mpilocal_out = cvDim_Bhi_mpilocal + cvDim_Ben_mpilocal + &
                         cvDim_Bchm_mpilocal

  END SUBROUTINE BMAT_reduceToMPILocal_r4

!--------------------------------------------------------------------------
! BMAT_expandToMPIGlobal
!--------------------------------------------------------------------------
  SUBROUTINE BMAT_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal,cvDim_mpiglobal_out)
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvBchm*)
    !
    ! Comments:
    !
    ! - LAM and ensemble components for constituents tbd.
    !
    implicit none

    real(8), intent(in)  :: cv_mpilocal(:)
    real(8), intent(out) :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpiglobal_out

    integer :: cvDim_Bhi_mpiglobal,cvDim_Ben_mpiglobal,cvDim_Bchm_mpiglobal

    real(8), pointer :: cvBhi_mpilocal(:) ,cvBen_mpilocal(:), cvBchm_mpilocal(:)
    real(8), pointer :: cvBhi_mpiglobal(:),cvBen_mpiglobal(:),cvBchm_mpiglobal(:)

    cvDim_Bhi_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BHI)) then
      cvBhi_mpilocal =>cvm_getSubVector(cv_mpilocal,cvm_BHI)
      if (mpi_myid == 0) then
         cvBhi_mpiglobal=>cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BHI)
      else
         cvBhi_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then
         call bhi_expandToMPIGlobal (cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpiglobal)
      else
         call lbhi_expandToMPIGlobal(cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpiglobal)
      end if
    end if

    cvDim_Ben_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BEN)) then
       cvBen_mpilocal =>cvm_getSubVector(cv_mpilocal,cvm_BEN)
       if (mpi_myid == 0) then
          cvBen_mpiglobal=>cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BEN)
       else
          cvBen_mpiglobal=>null()
       end if
       call ben_expandToMPIGlobal(cvBen_mpilocal,cvBen_mpiglobal,cvDim_Ben_mpiglobal)
    end if

    cvDim_Bchm_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BCHM)) then
      cvBchm_mpilocal =>cvm_getSubVector(cv_mpilocal,cvm_BCHM) 
      if (mpi_myid == 0) then
         cvBchm_mpiglobal=>cvm_getSubVector_mpiglobal(cv_mpiglobal,cvm_BCHM)
      else
         cvBchm_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then
         call bchm_expandToMPIGlobal(cvBchm_mpilocal,cvBchm_mpiglobal,cvDim_Bchm_mpiglobal)
      else
!         Done in lbhi_expandToMPIGlobal
      end if
    end if

    cvDim_mpiglobal_out = cvDim_Bhi_mpiglobal + cvDim_Ben_mpiglobal + cvDim_Bchm_mpiglobal

  end SUBROUTINE BMAT_expandToMPIGlobal

!--------------------------------------------------------------------------
! BMAT_expandToMPIGlobal_r4
!--------------------------------------------------------------------------
  SUBROUTINE BMAT_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal,cvDim_mpiglobal_out)
    !
    ! Revision:
    !           Ping Du, CMDA/MSC, July 2014
    !           - Additions for chemical constituents (see cvBchm*)
    !
    ! Comments:
    !
    ! - LAM and ensemble components for constituents tbd.
    !
    implicit none

    real(4), intent(in)  :: cv_mpilocal(:)
    real(4), intent(out) :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpiglobal_out

    integer :: cvDim_Bhi_mpiglobal,cvDim_Ben_mpiglobal,cvDim_Bchm_mpiglobal

    real(4), pointer :: cvBhi_mpilocal(:) ,cvBen_mpilocal(:), cvBchm_mpilocal(:)
    real(4), pointer :: cvBhi_mpiglobal(:),cvBen_mpiglobal(:),cvBchm_mpiglobal(:)

    cvDim_Bhi_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BHI)) then
      cvBhi_mpilocal =>cvm_getSubVector_r4(cv_mpilocal,cvm_BHI)
      if (mpi_myid == 0) then
         cvBhi_mpiglobal=>cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BHI)
      else
         cvBhi_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then
         call bhi_expandToMPIGlobal_r4 (cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpiglobal)
      else
         call lbhi_expandToMPIGlobal_r4(cvBhi_mpilocal,cvBhi_mpiglobal,cvDim_Bhi_mpiglobal)
      end if
    end if

    cvDim_Ben_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BEN)) then
       cvBen_mpilocal => cvm_getSubVector_r4(cv_mpilocal,cvm_BEN)
       if (mpi_myid == 0) then
          cvBen_mpiglobal => cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BEN)
       else
          cvBen_mpiglobal => null()
       end if
       call ben_expandToMPIGlobal_r4(cvBen_mpilocal,cvBen_mpiglobal,cvDim_Ben_mpiglobal)
    end if

    cvDim_Bchm_mpiglobal = 0
    if(cvm_subVectorExists(cvm_BCHM)) then
      cvBchm_mpilocal =>cvm_getSubVector_r4(cv_mpilocal,cvm_BCHM)
      if (mpi_myid == 0) then
         cvBchm_mpiglobal=>cvm_getSubVector_mpiglobal_r4(cv_mpiglobal,cvm_BCHM)
      else
         cvBchm_mpiglobal=>null()
      end if
      if ( hco_anl%global ) then
         call bchm_expandToMPIGlobal_r4(cvBchm_mpilocal,cvBchm_mpiglobal,cvDim_Bchm_mpiglobal)
      else
!         Done in lbhi_expandToMPIGlobal
      end if
    end if

    cvDim_mpiglobal_out = cvDim_Bhi_mpiglobal + cvDim_Ben_mpiglobal + cvDim_Bchm_mpiglobal

  end SUBROUTINE BMAT_expandToMPIGlobal_r4

END MODULE BMatrix_mod
