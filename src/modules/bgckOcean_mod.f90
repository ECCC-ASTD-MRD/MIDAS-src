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
!-------------------------------------- LICENCE end --------------------------------------

module bgckOcean_mod
  ! MODULE bgckOcean_mod (prefix='ocebg' category='1. High-level functionality')
  !
  ! :Purpose: to perform ocean data background Check
  !
  use mpi_mod
  use utilities_mod
  use obsSpaceData_mod
  use columnData_mod
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use statetocolumn_mod 
  use bufr_mod
  use mathPhysConstants_mod

  implicit none

  
  save
  private

  ! Public functions/subroutines
  public :: ocebg_bgCheckSST
 
  ! External functions
  integer, external :: fnom, fclos  

  character(len=20) :: timeInterpType_nl       ! 'NEAREST' or 'LINEAR'
  integer           :: numObsBatches           ! number of batches for calling interp setup
  namelist /namOceanBGcheck/ timeInterpType_nl, numObsBatches

  contains

  !----------------------------------------------------------------------------------------
  ! ocebg_bgCheckSST
  !----------------------------------------------------------------------------------------
  subroutine ocebg_bgCheckSST(obsData, columnTrlOnTrlLev, hco)
    !
    !: Purpose: to compute SST data background Check  
    !           
    
    implicit none

    ! Arguments:
    type(struct_obs)       , intent(inout)       :: obsData           ! obsSpaceData object
    type(struct_columnData), intent(inout)       :: columnTrlOnTrlLev ! column data on trl levels
    type(struct_hco)       , intent(in), pointer :: hco               ! horizontal trl grid

    ! Locals:
    type(struct_gsv)            :: stateVector   ! state vector containing std B estimation field
    integer                     :: nulnam, ierr, headerIndex, bodyIndex, obsFlag, obsVarno
    integer                     :: numberObs, numberObsRejected  
    real(8)                     :: OER, OmP, FGE, bgCheck
    logical                     :: llok
    type(struct_columnData)     :: columnFGE
    
    write(*,*) 'ocebg_bgCheckSST: performing background check for the SST data...'
    
    ! Setting default namelist variable values
    timeInterpType_nl = 'NEAREST'
    numObsBatches = 20

    ! Read the namelist
    if (.not. utl_isNamelistPresent('namOceanBGcheck','./flnml')) then
      if (mpi_myid == 0) then
        write(*,*) 'ocebg_bgCheckSST: namOceanBGcheck is missing in the namelist.'
        write(*,*) 'ocebg_bgCheckSST: the default values will be taken.'
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam, nml = namOceanBGcheck, iostat = ierr)
      if (ierr /= 0) call utl_abort(myName//': Error reading namelist')
      ierr = fclos(nulnam)
    end if
    write(*,*) 'ocebg_bgCheckSST: interpolation type: ', timeInterpType_nl
    write(*,*) 'ocebg_bgCheckSST: number obs batches: ', numObsBatches

    ! Read First Guess Error (FGE) and put it into stateVector
    call gsv_allocate(stateVector, 1, hco, columnTrlOnTrlLev % vco, dataKind_opt = 4, &
                      hInterpolateDegree_opt = 'NEAREST', &
                      datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gsv_readFromFile(stateVector, './bgstddev', 'STDDEV', 'X', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    
    call col_setVco(columnFGE, col_getVco(columnTrlOnTrlLev))
    call col_allocate(columnFGE, col_getNumCol(columnTrlOnTrlLev))
   
    ! Convert stateVector to column object
    call s2c_nl(stateVector, obsData, columnFGE, hco, timeInterpType = timeInterpType_nl, &
                moveObsAtPole_opt = .true., numObsBatches_opt = numObsBatches, dealloc_opt = .true.)

    numberObs = 0
    numberObsRejected = 0
    do headerIndex = 1, obs_numheader(obsData)
      
      bodyIndex = obs_headElem_i(obsData, obs_rln, headerIndex)
      obsVarno  = obs_bodyElem_i(obsData, obs_vnm, bodyIndex)
      llok = (obs_bodyElem_i(obsData, obs_ass, bodyIndex) == obs_assimilated)
      if (llok) then
        if (obsVarno == bufr_sst) then
       
	  FGE = col_getElem(columnFGE, 1, headerIndex, 'TM')
	  OmP = obs_bodyElem_r(obsData, OBS_OMP , bodyIndex)
          OER = obs_bodyElem_r(obsData, OBS_OER , bodyIndex)
	    
	  if (FGE /= MPC_missingValue_R8 .and. OmP /= MPC_missingValue_R8) then 
	    
	    numberObs = numberObs + 1
	    call obs_bodySet_r(obsData, OBS_HPHT, bodyIndex, FGE)
	    bgCheck = (OmP)**2 / (FGE**2 + OER**2)
	    obsFlag = ocebg_setFlag(obsVarno, bgCheck)
	
            if (obsFlag >= 2) then
              numberObsRejected = numberObsRejected + 1
	      write(*,'(a,i7,a,i7)')'ocebg_bgCheckSST:*********** ', numberObsRejected, ', header index: ', headerIndex
	      write(*,'(a)') 'ocebg_bgCheckSST: rejected '//obs_elem_c(obsData, 'STID' , headerIndex)//' data:'
	      write(*,'(a,i5,a,4f10.4)') 'ocebg_bgCheckSST: codtype: ', obs_headElem_i(obsData, obs_ity, headerIndex), &
              ', lon/lat/obs.value/OmP: ', &
              obs_headElem_r(obsData, obs_lon, headerIndex) * MPC_DEGREES_PER_RADIAN_R8, &
              obs_headElem_r(obsData, obs_lat, headerIndex) * MPC_DEGREES_PER_RADIAN_R8, &
              obs_bodyElem_r(obsData, obs_var, bodyIndex), OmP
            end if
	    	      
	    ! update background check flags based on bgCheck
            ! (element flags + global header flags)  
	    if (obsFlag == 1) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 13))
            else if (obsFlag == 2) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 14))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 16))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 09))
              call obs_headSet_i(obsData, obs_st1, headerIndex, ibset(obs_headElem_i(obsData, obs_st1, headerIndex), 06))
            else if (obsFlag == 3) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 15))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 16))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 09))
              call obs_headSet_i(obsData, obs_st1, headerIndex, ibset(obs_headElem_i(obsData, obs_st1, headerIndex), 06))
            end if
	   
          end if
	end if
      end if
      
    end do 

    if (numberObs > 0) then
      write(*,*)' '
      write(*,*) 'ocebg_bgCheckSST: background check of TM data is computed'
      write(*,'(a, i7,a,i7,a)') 'ocebg_bgCheckSST:   ', numberObsRejected, ' observations out of ', numberObs,' rejected'
      write(*,*)' '
    end if
    
    call gsv_deallocate(stateVector)
    call col_deallocate(columnFGE)
    
  end subroutine ocebg_bgCheckSST

  !--------------------------------------------------------------------------
  ! ocebg_setFlag
  !--------------------------------------------------------------------------
  function ocebg_setFlag(obsVarno, bgCheck) result(obsFlag)
    !
    !:Purpose: Set background-check flags according to values set in a table.
    !          Original values in table come from ECMWF.
    !

    implicit none
    
    integer             :: obsFlag  ! obs flag 

    ! Arguments:
    integer, intent(in) :: obsVarno ! obsVarno, Universal Field-Identity Numbers defined in bufr_mod
    real(8), intent(in) :: bgCheck  ! normalized background departure

    ! Locals:      
    real(8), parameter :: multipleSST(3) = (/  5.d0, 25.d0, 30.d0 /)

    obsFlag = 0
 
    if (obsVarno == bufr_sst) then
      if (bgCheck >= multipleSST(1) .and. bgCheck < multipleSST(2)) then
        obsFlag = 1
      else if (bgCheck >= multipleSST(2) .and. bgCheck < multipleSST(3)) then
        obsFlag = 2
      else if (bgCheck >= multipleSST(3)) then
        obsFlag = 3
      end if
    end if

  end function ocebg_setFlag
  
end module bgckOcean_mod  
