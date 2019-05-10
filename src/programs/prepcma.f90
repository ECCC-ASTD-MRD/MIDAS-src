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

program midas_prepcma
  !
  ! :Purpose: Read observation files that are in bgckalt format (i.e. output from
  !           the background check) and transform them to the ObsSpaceData
  !           format. All observations will collected in a single ObsSpaceData
  !           structure. The structure is output in binary format.
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
  NAMELIST /NAMPREPCMA/ cmahdr,cmabdy,cmadim,obsout,brpform
  character(len=256) :: cmahdr,cmabdy,cmadim,obsout,brpform

  integer :: fnom, fclos, get_max_rss, nulnam, ierr, datestamp
  type(struct_obs),       target  :: obsSpaceData

  real(kind=8) :: hx_dummy(1,1)
  integer, parameter :: zero_ens=0
  integer :: ncmahdr,ncmahx,ncmabdy,ncmadim,nobsout,nbrpform

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
  end if

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
  write(*,*) '> midas-prepcma: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  call obs_class_initialize('ENKF')
  call obs_initialize( obsSpaceData, mpi_local=obsf_filesSplit() )

  call tmg_start(11,'READ_OBS')
  call obsf_readFiles( obsSpaceData )
  call tmg_stop(11)

  write(*,*) 'midas_prepcma: obs_numheader(obsSpaceData)', obs_numheader(obsSpaceData)
  write(*,*) 'midas_prepcma: obs_numbody(obsSpaceData)  ', obs_numbody  (obsSpaceData)

  if (mpi_nprocs /= 1) then
    call obs_expandToMpiGlobal( obsSpaceData )
  end if

  !- Write the results
  write(*,*)
  write(*,*) '> midas-prepcma: writing to files'

  nobsout = 20
  call openfile(nobsout,obsout,'NEW','FORMATTED')
  call obs_print(obsSpaceData,nobsout)
  close(nobsout)

  ncmahdr = 21
  call openfile(ncmahdr,cmahdr,'NEW','UNFORMATTED')

  ncmabdy = 22
  call openfile(ncmabdy,cmabdy,'NEW','UNFORMATTED')

  ncmadim = 23
  call openfile(ncmadim,cmadim,'NEW','FORMATTED')

  !- Write the results in CMA format
  ncmahx  = -1
  call obs_write(obsSpaceData,hx_dummy,zero_ens,ncmahdr,ncmabdy,ncmahx,ncmadim)

  close(ncmahdr)
  close(ncmabdy)
  close(ncmadim)

  nbrpform = 0
  call openfile(nbrpform,brpform,'NEW','FORMATTED')
  !! This used to contain a .true. or .false. value indicating if observations passed the QCVar
  !! Since, this is not the case, we can write .false.
  write(nbrpform,*) .false.
  close(nbrpform)

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
  end if

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

    write(*,*) 'program prepcma (openfile): filename ', trim(filename)
    open(unit=unitnum,file=trim(filename),access='sequential', &
         form=trim(fileform),status=trim(filestat),IOSTAT=ierror)
    if (ierror /= 0) then
      call utl_abort('program prepcma (openfile): could not open file '// trim(filename))
    end if
  end subroutine openfile

end program
