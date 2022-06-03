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

program midas_randomPert
  !
  ! :Purpose: Main program for generating an ensemble of random perturbations
  !           based on the B matrix (can be homogeneous/isotropic or
  !           ensemble-based).
  !
  use version_mod
  use midasMpi_mod
  use ramDisk_mod
  use controlVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bmatrix_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use randomNumber_mod
  use utilities_mod
  use gridVariableTransforms_mod
  implicit none

  type(struct_gsv) :: statevector

  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()

  real(8), pointer :: field(:,:,:)
  
  integer :: fclos, fnom, fstopc, newdate, nstamp, ierr, status
  integer :: memberIndex, lonIndex, latIndex, cvIndex, levIndex, nkgdim
  integer :: idate, itime, ndate, nulnam
  integer :: get_max_rss, n_grid_point, n_grid_point_glb

  integer :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
  integer :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd

  real(8), allocatable :: controlVector(:), controlVector_mpiglobal(:)
  real(8), allocatable :: gdmean(:,:,:)
  real(4), allocatable :: ensemble_r4(:,:,:,:)
  real(8), allocatable :: avg_pturb_var(:)
  real(8), allocatable :: avg_pturb_var_glb(:)
  real(8), allocatable :: pturb_var(:,:,:)

  character(len=10) :: cldate
  character(len=3)  :: clmember
  character(len=25) :: clfiname
  character(len=12) :: etiket
  character(len=12) :: out_etiket
  
  logical  :: remove_mean, smoothVariances, mpiTopoIndependent
  integer  :: nens, seed, date, numBits
  NAMELIST /NAMENKF/nens, seed, date, out_etiket, remove_mean,  &
                    smoothVariances, mpiTopoIndependent, numBits

  call ver_printNameAndVersion('randomPert','Generation of random perturbations')

  !
  !- 0. MPI, tmg initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ierr = fstopc('MSGLVL','ERRORS',0)

  !
  !- 1. Set/Read values for the namelist NAMENKF
  !
  
  !- 1.1 Setting default values
  nens  = 10
  seed  = 1
  date  = 1900120100
  remove_mean = .true.
  out_etiket='RANDOM_PERT' 
  smoothVariances = .false.
  mpiTopoIndependent = .false.
  numBits = 32
  
  !- 1.2 Read the namelist
  nulnam=0
  ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namenkf, iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-randomPert: Error reading namelist')
  if( mmpi_myid == 0 ) write(*,nml=namenkf)
  ierr=fclos(nulnam)

  ndate   = date
  write(cldate, '(I10)') ndate

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !- 2.1 Decompose ndate(yyyymmddhh) into date(YYYYMMDD) time(HHMMSShh)
  !      calculate date-time stamp for postproc.ftn 
  if (ndate <= 0) then
    nstamp = ndate
  else
    idate   = ndate/100
    itime   = (ndate-idate*100)*1000000
    ierr    = newdate(nstamp, idate, itime, 3)
  end if
  if( mmpi_myid == 0 ) write(*,*) ' idate= ', idate, ' time= ', itime
  if( mmpi_myid == 0 ) write(*,*) ' date= ', ndate, ' stamp= ', nstamp

  !- 2.2 Initialize variables of the model states
  call gsv_setup

  !
  !- Initialize the Temporal grid
  !
  call tim_setup
  call tim_setDatestamp(nstamp)

  !- 2.3 Initialize the Analysis grid
  if (mmpi_myid == 0) write(*,*) ''
  if (mmpi_myid == 0) write(*,*) ' Set hco parameters for analysis grid'
  call hco_setupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_setupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  call mmpi_setup_latbands(hco_anl % nj,                & ! IN
                           latPerPE, latPerPEmax, myLatBeg, myLatEnd ) ! OUT
  call mmpi_setup_lonbands(hco_anl % ni,                & ! IN
                           lonPerPE, lonPerPEmax, myLonBeg, myLonEnd ) ! OUT

  !- 2.4 Initialize the vertical coordinate from the analysisgrid file
  call vco_setupFromFile(vco_anl, './analysisgrid', ' ')
 
  !- 2.5 Initialize the B_hi matrix
  call bmat_setup(hco_anl, hco_core, vco_anl)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.6 Initialize the gridded variable transform module
  call gvt_setup(hco_anl,hco_core,vco_anl)
  if ( gsv_varExist(varName='HU') ) call gvt_setupRefFromTrialFiles('HU')

  !
  !- 3. Memory allocations
  !

  !- 3.1 Allocate the statevector
  call gsv_allocate(statevector, 1, hco_anl, vco_anl, &
                    dateStamp_opt=nstamp, mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)
  nkgdim = statevector%nk
  allocate(ensemble_r4(myLonBeg:myLonEnd, myLatBeg:myLatEnd, nkgdim, nEns))

  !- 3.2 Allocate auxillary variables
  allocate(gdmean(myLonBeg:myLonEnd, myLatBeg:myLatEnd, nkgdim), STAT=status)
  if ( status /= 0 ) then
    call utl_abort('midas-randomPert: PROBLEM WITH ALLOCATING OF GDMEAN')
  end if

  allocate(controlVector(cvm_nvadim))
  allocate(controlVector_mpiglobal(cvm_nvadim_mpiglobal))

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 4. Compute an ensemble of random perturbations
  !

  if( mmpi_myid == 0 ) write(*,*) '******************'
  if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: COMPUTE the mean of the random ', &
                                  'perturbations of all the members'

  if( mpiTopoIndependent ) then
    call rng_setup(abs(seed))
  else
    call rng_setup(abs(seed+mmpi_myid))
  end if

  gdmean(:,:,:) = 0.0D0
  call gsv_getField(statevector,field)

  !- 4.1 Generate a (potentially) biased ensemble
  do memberIndex = 1, NENS

    if( mmpi_myid == 0 ) write(*,*) ' computing member number= ', memberIndex

    !- 4.1.1 Create a random control vector in spectral space

    if( mpiTopoIndependent ) then
      !- Global vector (for testing different mpi topology, less efficient)
      do cvIndex = 1, cvm_nvadim_mpiglobal
        controlVector_mpiglobal(cvIndex) = rng_gaussian()
      end do
      call bmat_reduceToMPILocal( controlVector, controlVector_mpiglobal )
    else
      !- Local vector (different seed for each processor, more efficient)
      do cvIndex = 1, cvm_nvadim
        controlVector(cvIndex) = rng_gaussian()
      end do
    end if

    !- 4.1.2 Transform to control variables in physical space
    call bmat_sqrtB(controlVector, cvm_nvadim, & ! IN
                    statevector               )  ! OUT

    !- 4.1.3 Running ensemble sum
    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) = real(field(lonIndex, latIndex, levIndex), 4)
          gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) + field(lonIndex, latIndex, levIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end do

  call gsv_deallocate(statevector)
  
  !- 4.2 Remove the ensemble mean
  if ( REMOVE_MEAN ) then
    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          gdmean(lonIndex, latIndex, levIndex) = gdmean(lonIndex, latIndex, levIndex) / real(NENS, 8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (lonIndex, memberIndex, latIndex, levIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) =  &
                ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) -  &
                real(gdmean(lonIndex, latIndex, levIndex), 4)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if
  
  !- 4.3 Smooth variances to horizontally constant values
  if ( smoothVariances ) then
  
    allocate(pturb_var(myLonBeg:myLonEnd, myLatBeg:myLatEnd, statevector%nk))
    allocate(avg_pturb_var(statevector%nk), avg_pturb_var_glb(statevector%nk))
    pturb_var(:,:,:) = 0.0D0
    avg_pturb_var(:) = 0.0D0
    avg_pturb_var_glb(:) = 0.0D0
  
    do memberIndex = 1, NENS  
      !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)  
      do levIndex = 1, nkgdim
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
               pturb_var(lonIndex, latIndex, levIndex) = pturb_var(lonIndex, latIndex, levIndex) +  &
                   real(ensemble_r4(lonIndex, latIndex, levIndex, memberIndex), 8)**2
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do
  
    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)  
    do levIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          pturb_var(lonIndex, latIndex, levIndex) = pturb_var(lonIndex, latIndex, levIndex) / real(NENS, 8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  
    do latIndex = myLatBeg, myLatEnd
      do lonIndex = myLonBeg, myLonEnd
        !$OMP PARALLEL DO PRIVATE (levIndex)  
        do levIndex = 1, nkgdim
          avg_pturb_var(levIndex) = avg_pturb_var(levIndex) + pturb_var(lonIndex, latIndex, levIndex)
        end do
        !$OMP END PARALLEL DO
      end do
    end do
  
    n_grid_point=(lonPerPE)*(latPerPE)
    call rpn_comm_allreduce(n_grid_point, n_grid_point_glb, 1,  &
                            "mpi_double_precision", "mpi_sum", "GRID", ierr)
    call rpn_comm_allreduce(avg_pturb_var, avg_pturb_var_glb, nkgdim,  &
                            "mpi_double_precision", "mpi_sum", "GRID", ierr)
  
    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)  
    do levIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          if( pturb_var(lonIndex, latIndex, levIndex) > 0.0d0 ) then
            pturb_var(lonIndex, latIndex, levIndex) = sqrt( pturb_var(lonIndex, latIndex, levIndex))
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  
    !$OMP PARALLEL DO PRIVATE (levIndex)  
    do levIndex = 1, nkgdim
      if( avg_pturb_var_glb(levIndex) > 0.0d0 ) then
        avg_pturb_var_glb(levIndex) = sqrt( avg_pturb_var_glb(levIndex) / real(n_grid_point_glb, 8) )
      end if
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex, memberIndex)    
    do memberIndex = 1, NENS
      do levIndex = 1, nkgdim
        if( avg_pturb_var_glb(levIndex) > 0.0d0 ) then
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              if( pturb_var(lonIndex, latIndex, levIndex) > 0.0d0 ) then
                ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) =   &
                    ensemble_r4(lonIndex, latIndex, levIndex, memberIndex) *  &
                    real(avg_pturb_var_glb(levIndex)/pturb_var(lonIndex, latIndex, levIndex), 4)
              end if
            end do
          end do
        end if
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  !- 4.4 Write the perturbations
  do memberIndex = 1, NENS
    if( mmpi_myid == 0 ) write(*,*)
    if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: pre-processing for writing member number= ', memberIndex

    call gsv_allocate(statevector, 1, hco_anl, vco_anl, &
                      dateStamp_opt=nstamp, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_getField(statevector,field)

    !$OMP PARALLEL DO PRIVATE (lonIndex, latIndex, levIndex)    
    do levIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          field(lonIndex, latIndex, levIndex) = ensemble_r4(lonIndex, latIndex, levIndex, memberIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    write(clmember, '(I3.3)') memberIndex
    if (ndate > 0) then
      clfiname = './pert_'//trim(cldate)//'_'//trim(clmember)
    else
      clfiname = './pert_'//trim(clmember)
    end if
    if( mmpi_myid == 0 ) write(*,*) 'midas-randomPert: processing clfiname= ', clfiname

    call gio_writeToFile(statevector, clfiname, out_etiket,      & ! IN
                         numBits_opt=numBits, unitConversion_opt=.true.)  ! IN

    call gsv_deallocate(statevector)

  end do

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 5.  Memory deallocations
  !
  deallocate(gdmean, STAT=status)
  if ( status /= 0 ) then
    call utl_abort('midas-randomPert: PROBLEM WITH DEALLOCATE OF GDMEAN')
  end if

  deallocate(ensemble_r4)
  deallocate(controlVector)  
  deallocate(controlVector_mpiglobal)  

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  !
  !- 6.  MPI, tmg finalize
  !  
  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if( mmpi_myid == 0 ) write(*,*) ' --------------------------------'
  if( mmpi_myid == 0 ) write(*,*) ' MIDAS-RANDOMPERT ENDS'
  if( mmpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_randomPert
