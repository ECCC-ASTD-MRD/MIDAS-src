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
!!
!! *Purpose*: read observation files that are in bgckalt format
!!         (i.e. output from the background check) and
!!        transform them to the ObsSpaceData format. All observations will
!!        collected in a single ObSSpaceData structure. The structure is
!!        output in binary format.
!!
!--------------------------------------------------------------------------
program midas_prepcma
  !
  ! *Purpose*: read observation files that are in bgckalt format
  !         (i.e. output from the background check) and
  !        transform them to the ObsSpaceData format. All observations will
  !        collected in a single ObSSpaceData structure. The structure is
  !        output in binary format.
  !
  use obsSpaceData_mod
  use mathPhysConstants_mod
  use obsFiles_mod
  use utilities_mod
  use mpi_mod
  use ramDisk_mod
  use obsUtil_mod

  implicit none

  ! Namelist
  NAMELIST /NAMPREPCMA/ infl_sigo,thinning,dyn_sw_oer
  logical :: infl_sigo,thinning,dyn_sw_oer

  integer :: fnom, fclos, get_max_rss, nulnam, ierr, datestamp
  type(struct_obs),       target  :: obsSpaceData

  real(kind=8), dimension(1,1) :: hx_dummy
  integer, parameter :: zero_ens=0
  integer :: ncmahdr,ncmahx,ncmabdy,ncmadim,nobsout

  write(*,*) " --------------------------------------------"
  write(*,*) " ---  START OF MAIN PROGRAM midas-prepcma ---"
  write(*,*) " ---  Computation of the innovation       ---"
  write(*,*) " ---  Revision: GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE "
  write(*,*) " --------------------------------------------"

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_PREPCMA' )
  call tmg_start(1,'MAIN')

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('PREPCMA_BEG')
  endif

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namprepcma,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-prepcma: Error reading namelist')
  if (mpi_myid == 0) write(*,nml=namprepcma)
  ierr = fclos(nulnam)

  !- RAM disk usage
  call ram_setup

  !- Observation file names and get datestamp
  call obsf_setup(dateStamp, 'analysis' )

  write(*,*)
  write(*,*) '> midas-obsIO: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  call obs_class_initialize('ENKF')
  call obs_initialize( obsSpaceData, mpi_local=obsf_filesSplit() )

  call tmg_start(11,'READ_OBS')
  call obsf_readFiles( obsSpaceData )
  call tmg_stop(11)

  WRITE(*,*) 'midas_prepcma: obs_numheader(obsdat)', obs_numheader(obsdat)
  WRITE(*,*) 'midas_prepcma: obs_numbody(obsdat)  ', obs_numbody  (obsdat)

  nobsout = -1
  call openfile(nobsout,'cma.ascii.1','NEW','FORMATTED')
  call obs_print(obsSpaceData,nobsout)
  close(nobsout)

  ! 2.3 Write the results in CMA format
  write(*,*)
  write(*,*) '> midas-prepcma: writing to file'

  ncmahdr = 0
  call openfile(ncmahdr,'cmaheader','NEW','FORMATTED')
  ncmabdy = 0
  call openfile(ncmabdy,'cmabdy','NEW','FORMATTED')
  ncmadim = 0
  call openfile(ncmadim,'cmadim','NEW','FORMATTED')
  ncmahx  = 0

  call obs_write(obsSpaceData,hx_dummy,zero_ens,ncmahdr,ncmabdy,ncmahx,ncmadim)

  close(ncmahdr)
  close(ncmabdy)
  close(ncmadim)

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-prepcma: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_PREPCMA' )

  call rpn_comm_finalize(ierr)

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('PREPCMA_END')
  endif

contains

  subroutine openfile(unitnum,filename,filestat,fileform)
    !
    ! s/r openfile - open file filename and verify the exit
    !     status of the command.
    !
    ! author P. Houtekamer and H. Mitchell May 2005.
    !
    ! input:
    !     unitnum:  unit number for the file
    !     filename: file name
    !     filestat: status ('OLD' or 'NEW')
    !     fileform: format ('FORMATTED' or 'UNFORMATTED')
    !
    implicit none
    integer   :: unitnum,ierror
    character (len=*) :: filename
    character (len=*)  :: filestat
    character (len=*)  :: fileform

    write(*,*) 'open file: ',filename
    open(unit=unitnum,file=trim(filename),access='sequential', &
         form=trim(fileform),status=trim(filestat),IOSTAT=ierror)
    if (ierror.ne.0) then
      write(*,*) 'could not open unit ',unitnum, &
           ' for file: ',filename
      write(*,*) 'system error number: ',ierror
      call qqexit(1)
      stop
    endif

    return
  end subroutine openfile

end program
