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

MODULE biasCorrectionConv_mod
  ! MODULE biasCorrectionConv_mod (prefix="bcc" category='1. High-level functionality')
  !
  ! :Purpose: Performs bias correction for conventional observations
  !
  use utilities_mod
  use obsSpaceData_mod
  

  implicit none
  save
  private
  integer, parameter :: nPhases=3, nLevels=5, nStationMax=100000 
  integer  :: nb_aircraft_bias
  real     :: corrects_TT(nStationMax,nPhases,nLevels)
  character(len=9) :: aircraft_ID(nStationMax)

  public               :: bcc_readConfig, bcc_readAIBcor
 

  
  integer, external            :: fnom, fclos
 
  namelist /nambiasconv/ 
  
CONTAINS
 
  !-----------------------------------------------------------------------
  ! bcc_readConfig
  !-----------------------------------------------------------------------
  subroutine bcc_readConfig()
    !
    ! :Purpose: Read nambiasconv namelist section
    !
    implicit none
    !Locals:
    integer  :: ierr,nulnam
  
    ! set default values for namelist variables
   
    ! read in the namelist NAMBIASCONV
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambiasconv,iostat=ierr)
    if ( ierr /= 0 .and. mpi_myid == 0 )  &
         write(*,*) 'bcs_readConfig: WARNING: Error reading namelist, assume it will not be used!'
    if ( mpi_myid == 0 ) write(*,nml=nambiasconv)
    ierr = fclos(nulnam)
    

  end subroutine bcc_readConfig

  !-----------------------------------------------------------------------
  ! bcc_readAIBcor
  !-----------------------------------------------------------------------
  subroutine bcc_readAIBcor(file_cor)
    !
    ! :Purpose: Read TT bias corrections 
    !     The first line of the file is the number of aircraft plus one.
    !     The rest of the file gives 15 values of Mean O-A for each aircraft, with each (AC,value) line written in format "a9,1x,f6.2".
    !     The order is the same as what is written by genbiascorr script genbc.aircraft_bcor.py.
    !     The first "aircraft" (AC name = BULKBCORS) values are the bulk corrections by layer for All-AC (first 5 values), AIREP/ADS 
    !     (second 5 values) and AMDAR/BUFR (last 5 values).
    !     Missing value = 99.0.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: file_cor
    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex, phaseIndex, levelIndex
    real    :: cor_ligne
    character(len=9) :: id_ligne

    nulcoeff = 0
    ierr = fnom(nulcoeff, file_cor, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBcor: unable to open airplanes bias correction file ' // file_cor )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nb_aircraft_bias
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBcor: error 1 while reading airplanes bias correction file ' // file_cor )
    end if
    do stationIndex=1,nb_aircraft_bias
      do phaseIndex=1,3
        do levelIndex=1,5
          read (nulcoeff, *, iostat=ierr) id_ligne,cor_ligne
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readAIBcor: error 2 while reading airplanes bias correction file ' // file_cor )
          end if
          corrects_TT(stationIndex,phaseIndex,levelIndex) = cor_ligne
          aircraft_ID(stationIndex)                       = id_ligne
          !print*, stationIndex, phaseIndex, levelIndex, aircraft_ID(stationIndex), corrects_TT(stationIndex,phaseIndex,levelIndex)
        end do
      end do
    end do
    ierr = fclos(nulcoeff)

    ! Check for bulk bias corrections at start of file
    if ( aircraft_ID(1) /= "BULKBCORS" ) then
      call utl_abort('bcc_readAIBcor: ERROR: Bulk bias corrections are missing in bias correction file!' )
    end if

  end subroutine bcc_readAIBcor


end MODULE biasCorrectionSat_mod

