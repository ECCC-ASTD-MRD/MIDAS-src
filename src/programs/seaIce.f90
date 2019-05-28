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

program midas_seaIce
  !
  ! :Purpose: Main program for sea ice assimilation
  ! :Origin:  ensmanip/midas-ensmanip.ftn90
  !

  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use controlVector_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use varNameList_mod
  use ramDisk_mod
  use bmatrix_mod

  implicit none

  type(struct_gsv) :: statevector

  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_anl => null()

  real(8), pointer :: field(:,:,:)

  real(8), allocatable :: controlVector(:)

  integer :: fclos, fnom, fstopc, newdate, nstamp, ierr
  integer :: memberIndex, lonIndex, latIndex, levIndex, stepIndex, numStep
  integer :: jlon, jlat, jj, jvar, ip3
  integer :: idate, itime, nulnam
  integer :: nlons, nlats
  integer :: get_max_rss
  integer :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
  integer :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
  character(len=128) :: filename

  ! namelist variables
  integer :: date
  integer :: oneobs_lons(100), oneobs_lats(100)

  NAMELIST /NAMTESTSEAICE/ date, oneobs_lons, oneobs_lats

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-SEAICE             --",/,' //   &
        '14x,"-- Sea ice data assimilation --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !
  !- 0. MPI, tmg initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_SEAICE' )

  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !
  !- 1. Set/Read values for the namelist NAMTESTSEAICE
  !
  
  !- 1.1 Setting default values for namelist variables

  date = 2011020100
  oneobs_lons(:)=-1
  oneobs_lats(:)=-1

  !- 1.2 Read the namelist
  nulnam=0
  ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namtestseaice, iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-seaIce: Error reading namelist')
  if( mpi_myid == 0 ) write(*,nml=namtestseaice)
  ierr=fclos(nulnam)

  nlons=0
  do jj = 1, size(oneobs_lons)
     if(oneobs_lons(jj) >= 1) nlons=nlons+1  
  enddo
  nlats=0
  do jj = 1, size(oneobs_lats)
     if(oneobs_lats(jj) >= 1) nlats=nlats+1  
  enddo

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  !- 2.1 Decompose date(yyyymmddhh) into idate(YYYYMMDD) itime(HHMMSShh)
  !      calculate date-time stamp for postproc.ftn 
  idate   = date/100
  itime   = (date-idate*100)*1000000
  ierr    = newdate(nstamp, idate, itime, 3)
  if( mpi_myid == 0 ) write(*,*)' idate= ', idate, ' time= ', itime
  if( mpi_myid == 0 ) write(*,*)' date= ', date, ' stamp= ', nstamp

  !- 2.2 Initialize variables of the model states
  call gsv_setup

  !
  !- Initialize the Temporal grid
  !
  call tim_setup
  call tim_setDatestamp(nstamp)
  numstep = tim_nstepobsinc

  !- 2.3 Initialize the analysis grid
  if (mpi_myid == 0) write(*,*)
  if (mpi_myid == 0) write(*,*)' Set hco parameters for analysis grid'

  call hco_SetupFromFile( hco_anl, './analysisgrid', ' ', 'ANALYSISGRID')
  call vco_setupFromFile( vco_anl, './analysisgrid' )

  call mpivar_setup_latbands(hco_anl % nj,           & ! IN
                             latPerPE, latPerPEmax,  & ! OUT
                             myLatBeg, myLatEnd )      ! OUT
  call mpivar_setup_lonbands(hco_anl % ni,           & ! IN
                             lonPerPE, lonPerPEmax,  & ! OUT
                             myLonBeg, myLonEnd )      ! OUT


  ! DO THE TEST FOR DIFFUSION B MATRIX

  ! Allocate the statevector
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                    datestamp_opt=tim_getDatestamp(),mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)

  ! Setup the B matrix (which also setup the control vector module and cvm_nvadim)
  call bmat_setup(hco_anl,vco_anl)

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !==============================================
  !- Compute columns of B matrix
  !==============================================
  !
  if ( nlons >= 1 .and. nlats >= 1 ) then

     allocate(controlVector(cvm_nvadim))

     write(*,*) '***************************'
     write(*,*) 'Compute columns of B matrix'
     write(*,*) '***************************'

     write(*,*) 'number of longitudes =',nlons
     write(*,*) 'number of latitudes  =',nlats

     do jvar = 1, vnl_numvarmax

        if (.not. gsv_varExist(statevector_opt=statevector, varName=vnl_varNameList(jvar)) ) cycle

        filename = 'columnB_' // trim(vnl_varNameList(jvar)) // '.fst'

        ip3 = 0
        do jlon = 1, nlons
           do jlat = 1, nlats

              call gsv_zero(statevector)
              field => gsv_getField3d_r8(statevector, varName_opt=vnl_varNameList(jvar))

              if(oneobs_lats(jlat) >= statevector%myLatBeg .and. oneobs_lats(jlat) <= statevector%myLatEnd .and.  &
                   oneobs_lons(jlon) >= statevector%myLonBeg .and. oneobs_lons(jlon) <= statevector%myLonEnd) then
                 if(vnl_varLevelFromVarname(vnl_varNameList(jvar)) == 'SF') then
                    field(oneobs_lons(jlon),oneobs_lats(jlat),1) = 1.0D0
                 else
                    call utl_abort('midas-seaIce: Expecting only surface field')
                 endif
              endif

              controlVector(:) = 0.0d0
              call bmat_sqrtBT(controlVector,cvm_nvadim,statevector)
              call bmat_sqrtB(controlVector,cvm_nvadim,statevector)

              write(*,*)'seaice: writing out the column of B, jlon,jlat=',jlon,jlat
              call flush(6)

              ip3 = ip3 + 1
              call gsv_writeToFile(statevector,filename,'ONEOBS_'//trim(vnl_varNameList(jvar)),  &
                                   ip3_opt=ip3,unitConversion_opt=.true.)

           enddo
        enddo

     enddo

     deallocate(controlVector)
     
  endif ! if any oneobs selected

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- 6.  MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_SEAICE' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if( mpi_myid == 0 ) write(*,*) ' MIDAS-SEAICE ENDS'
  if( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_seaIce
