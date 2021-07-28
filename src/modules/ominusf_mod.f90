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

module oMinusF_mod
  ! MODULE oMinusF_mod (prefix='omf' category='1. High-level functionality')
  !
  ! :Purpose: Module for Observation minus Forecast (O-F) computation
  !
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpiVar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod
  use stateToColumn_mod
  use gridStateVector_mod
  use obsFiles_mod
  use obsFilter_mod
  use innovation_mod
  use tovs_nl_mod
  use obsErrors_mod
  use obsOperators_mod
  use biasCorrectionConv_mod
  use obsSpaceErrorStdDev_mod
  implicit none
  private

  ! public subroutines and functions
  public :: omf_oMinusF

  contains


    subroutine omf_oMinusF(trlColumnOnAnlLev, trlColumnOnTrlLev, obsSpaceData, &
                           varMode, addHBHT, addSigmaO)
      !
      ! :Purpose: compute Observation-minus-Forecast (OmF)
      !

      implicit none
      ! Arguments:
      type(struct_columnData),target, intent(inout)  :: trlColumnOnAnlLev
      type(struct_columnData),target, intent(inout)  :: trlColumnOnTrlLev
      type(struct_obs),       target, intent(inout)  :: obsSpaceData
      character(len=*), intent(in) :: varMode
      logical, intent(in) :: addHBHT
      logical, intent(in) :: addSigmaO

      ! locals
      type(struct_vco),       pointer :: vco_anl  => null()
      type(struct_hco),       pointer :: hco_anl  => null()
      type(struct_hco),       pointer :: hco_core => null()

      character(len=48) :: obsMpiStrategy
      character(len=3)  :: obsColumnMode
      character(len=10) :: trialFileName

      integer :: datestamp, get_max_rss

      write(*,*) " ---------------------------------------"
      write(*,*) " ---  START OF SUBROUTINE oMinusF    ---"
      write(*,*) " ---  Computation of the innovation  ---"
      write(*,*) " ---------------------------------------"

      !
      !- 1.  Settings and module initializations
      !
      write(*,*)
      write(*,*) '> omf_oMinusF: setup - START'

      obsMpiStrategy = 'LIKESPLITFILES'
      obsColumnMode  = 'VAR'
      trialFileName  = './trlm_01'

      !- 1.3 RAM disk usage
      call ram_setup

      !- 1.4 Temporal grid
      call tim_setup( fileNameForDate_opt=trim(trialFileName) )

      !- 1.5 Observation file names and get datestamp, but do not use it
      call obsf_setup(dateStamp, trim(varMode))

      !- 1.6 Constants
      if (mpi_myid == 0) call mpc_printConstants(6)

      !- 1.7 Variables of the model states
      call gsv_setup

      !- 1.8 Set the horizontal domain
      if ( addHBHT ) then
        call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
        if ( hco_anl % global ) then
          hco_core => hco_anl
        else
          !- Iniatilized the core (Non-Exteded) analysis grid
          call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
        end if
      else
        call hco_SetupFromFile(hco_anl, trim(trialFileName), ' ') ! IN
        hco_core => hco_anl
      end if

      ! 1.9 Setup a column vector following the analysis vertical grid
      if ( addHBHT ) then
        call vco_SetupFromFile(vco_anl,        & ! OUT
                               './analysisgrid') ! IN
        call col_setVco(trlColumnOnAnlLev,vco_anl)
      end if

      !- 1.10 Setup and read observations
      call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy,trim(varMode))!IN

      ! Apply optional bias corrections when namelist logicals aiBiasActive, gpBiasActive are TRUE
      ! (Only reverse existing corrections when namelist logicals aiRevOnly, gpRevOnly are TRUE)
      call bcc_applyAIBcor(obsSpaceData)    
      call bcc_applyGPBcor(obsSpaceData)      
      
      !- 1.11 Basic setup of columnData module
      call col_setup

      !- 1.12 Memory allocation for background column data
      if ( addHBHT ) then
        call col_allocate(trlColumnOnAnlLev, obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
      end if

      if ( addSigmaO ) then
        !- 1.13 Initialize the observation error covariances
        write(*,*)
        write(*,*) '> omf_oMinusF: Adding sigma_O'
        call oer_setObsErrors(obsSpaceData, trim(varMode))
      end if

      !- 1.14 Reading, horizontal interpolation and unit conversions of the 3D background fields
      call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData,hco_core)

      write(*,*)
      write(*,*) '> omf_oMinusF: setup - END'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      !
      !- 2.  O-F computation
      !

      !- 2.1 Compute observation innovations
      write(*,*)
      write(*,*) '> omf_oMinusF: compute innovation'
      call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)

      if ( addHBHT ) then
        write(*,*)
        write(*,*) '> omf_oMinusF: Adding HBH^T'
        !- 2.2 Interpolate background columns to analysis levels and setup for linearized H
        call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)
        !- 2.3 Compute the background errors in observation space
        call ose_computeStddev(trlColumnOnAnlLev,hco_anl,obsSpaceData)
      end if

    end subroutine omf_oMinusF

end module oMinusF_mod
