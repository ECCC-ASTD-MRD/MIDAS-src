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

!--------------------------------------------------------------------------
!! *Purpose*: Do background check on all conventionnal observations
!!
!! @author P. Koclas *CMC/CMDA  Nov 1998
!!
!! Revision:  M. Sitwell (ARQI/AQRD) May 2015, March 2016
!!           - Added call to BGCDATA for chemical constituents
!!           - Added loop over bgfam list
!!           Y. Rochon (ARQI/AQRD) June 2016
!!           - Added call to osd_ObsSpaceDiag
!--------------------------------------------------------------------------
SUBROUTINE BGCHECK_CONV(columng,columnhr,obsSpaceData)

  use mpivar_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsSpaceDiag_mod
  use burpFiles_mod
  IMPLICIT NONE

  type(struct_obs)        :: obsSpaceData  ! Observation-related data
  type(struct_columnData) :: columng       !
  type(struct_columnData) :: columnhr      ! 

  INTEGER J,JDATA
  REAL*8 ZJO

  INTEGER*4      :: NULNAM,IER,FNOM,FCLOS
  CHARACTER *256 :: NAMFILE
  LOGICAL        :: NEW_BGCK_SW

  character(len=2), dimension(9) :: bgfam = (/ 'UA', 'AI', 'HU', 'SF', 'ST', 'SC', 'PR', 'GP', 'CH' /)
      
!
  call tmg_start(3,'BGCHECK_CONV')

  WRITE(*,FMT=9000)
9000 FORMAT(//,3(" **********"),/," BEGIN CONVENTIONNAL BACKGROUND CHECK",/,3(" **********"),/)

!stl
  NEW_BGCK_SW = .false.

  NAMELIST /NAMBGCKCONV/NEW_BGCK_SW
  NAMFILE=trim("flnml")
  nulnam=0
  IER=FNOM(NULNAM,NAMFILE,'R/O',0)

  READ(NULNAM,NML=NAMBGCKCONV,IOSTAT=IER)
  if(IER.ne.0) then
    write(*,*) 'selectb: No valid namelist NAMBGCKCONV found'
  endif

  iER=FCLOS(NULNAM)

  write(*,*) 'new_bgck_sw = ',new_bgck_sw
!stl


!     CALCULATE HBHT (sigma_B in observation space)
!     ----------------------------------------------
!
  call compute_HBHT( columng,   &   ! IN
                     columnhr,  &   ! IN
                     obsSpaceData ) ! INOUT

!
!     DO A BACKGROUND CHECK ON ALL THE OBSERVATIONS
!     ----------------------------------------------

  do j=1,size(bgfam)
     if (obs_famExist(obsSpaceData,bgfam(j))) CALL BGCDATA(ZJO,bgfam(j),obsSpaceData)
  end do

  if (obs_famExist(obsSpaceData,'RO')) CALL BGCGPSRO(columnhr,obsSpaceData)

!stl
  if (obs_famExist(obsSpaceData,'SW')) then
    if(new_bgck_sw) then
      CALL BGCDATA(ZJO,'SX',obsSpaceData)
    else
      CALL BGCDATA(ZJO,'SW',obsSpaceData)
    endif
  endif
!stl

! Conduct obs-space post-processing diagnostic tasks (some diagnostic 
! computations controlled by NAMOSD namelist in flnml)

  call osd_ObsSpaceDiag(obsSpaceData,columng,analysis_mode=.false.)

!
!     Write out contents of obsSpaceData into BURP files
!
  CALL burp_updateFiles(obsSpaceData)

  if (mpi_myid == 0 ) then
     DO j =1, min(1,obs_numHeader(obsSpaceData))
        call obs_prnthdr(obsSpaceData,j)
        call obs_prntbdy(obsSpaceData,j)
     END DO
  end if

  ! deallocate obsSpaceData
  call obs_finalize(obsSpaceData)

  call tmg_stop(3)

END SUBROUTINE BGCHECK_CONV
