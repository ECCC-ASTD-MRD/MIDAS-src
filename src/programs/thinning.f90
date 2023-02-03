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

program midas_thinning
  !
  !:Purpose: The thinning program reduces the density of observations for the purpose of
  !          assimilation.
  !
  !
  !:Algorithm: After setting up the ``obsSpace_data`` object, the thinning program calls
  !            the ``filt_suprep`` routine which rejects certain observations based on
  !            blacklists and other checks. 
  !
  !            -- 
  !
  !            Specific routines are then called for thinning each observation types. 
  !            These routines are found in ``thinning_mod`` and are controlled by the 
  !            following namelists:
  !
  !            |
  !
  !            ======================= ====== 
  !             Namelist                    
  !            ======================= ====== 
  !             ``thin_surface``          
  !             ``thin_raobs``            
  !             ``thin_aircraft``         
  !             ``thin_satwind``          
  !             ``thin_gpsro``           
  !             ``thin_gbgps``            
  !             ``thin_aladin``           
  !             ``thin_csr``              
  !             ``thin_scat``             
  !             ``thin_tovs``             
  !             ``thin_hyper``            
  !            ======================= ====== 
  !
  !            |
  !
  !            Observations that are rejected by the thinning routines have their 11th bit *flag* set. 
  !            In a subsequent step, these observations are removed from observation files. 
  !
  !            |
  !
  !
  !:File I/O: 
  !
  !           ============================================== ================================================================
  !            Input and Output Files (NWP applicaton)        Description of file
  !           ============================================== ================================================================
  !            ``flnml``                                      In - Main namelist file with parameters user may modify
  !            ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state -> necessary for observations in pressure
  !                                                           coordinates. 
  !            Observation Files                              In/Out - Sqlite/Burp observation files for different families. 
  !           ============================================== ================================================================
  !
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use tovs_nl_mod
  use obsFilter_mod
  use thinning_mod
  use fSQLite

  implicit none

  type(struct_obs) :: obsSpaceData
  integer :: istamp,exdb,exfin
  integer :: ierr, dateStampFromObs
  integer :: get_max_rss

  istamp = exdb('THINNING', 'DEBUT', 'NON')

  call ver_printNameAndVersion('thinning', 'Observation thinning')

  ! MPI initilization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0, 'Main')

  ! 1. Top level setup

  call ram_setup

  ! 2. configuration the job

  call utl_tmg_start(10,'--Observations')

  !- Initialize observation file names
  call obsf_setup(dateStampFromObs, 'thinning')

  !- Allocate obsSpaceData
  call obs_class_initialize('ALL')
  call obs_initialize(obsSpaceData, mpi_local = .true.)

  !- Setup obsFilter_mod
  call filt_setup('bgck')

  !- Read observations
  call utl_tmg_start(11,'----ReadObsFiles')
  call obsf_readFiles(obsSpaceData)
  call utl_tmg_stop(11)

  write(*,*) 'Memory Used: ', get_max_rss() / 1024,'Mb'

  !- Initialize TOVS processing
  if (obs_famExist(obsSpaceData, 'TO')) call tvs_setup

  !- Select the elements to "assimilate" and apply rejection flags
  call filt_suprep(obsSpaceData)

  call utl_tmg_stop(10)

  !- Setup timeCoord module and set dateStamp from env variable
  call tim_setup()
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      ! use dateStamp from obs if not set by env variable
      call tim_setDateStamp(dateStampFromObs)
    else
      call utl_abort('midas-thinning: DateStamp was not set')
    end if
  end if

  ! 3. Do the Thinning

  !- Set bit 11 of flag, process one observation type at a time
  call thn_thinHyper(obsSpaceData)
  call thn_thinTovs(obsSpaceData)
  call thn_thinCSR(obsSpaceData)
  call thn_thinScat(obsSpaceData)
  call thn_thinSatWinds(obsSpaceData)
  call thn_thinAircraft(obsSpaceData)
  call thn_thinSurface(obsSpaceData, 'SF') ! surface data thinning    
  if (obs_famExist(obsSpaceData, 'TM')) then
    call thn_thinSurface(obsSpaceData, 'TM') ! insitu SST thinning
    call thn_thinSatSST(obsSpaceData)        ! satellite SST thinning
  end if
  call thn_thinGbGps(obsSpaceData)
  call thn_thinGpsRo(obsSpaceData)
  call thn_thinAladin(obsSpaceData)
  write(*,*) 'Memory Used: ',get_max_rss() / 1024,'Mb'

  if (obs_famExist( obsSpaceData, 'UA' )) then
    write(*,*) 'midas-thinning: WARNING: radiosonde observations found in obsSpaceData!'
    write(*,*) '                These observations cannot be thinned using the stand-alone'
    write(*,*) '                MIDAS thinning program. Instead they should be thinned'
    write(*,*) '                in combination with doing the background check in the'
    write(*,*) '                obsSelection program.'
  end if

  !- Update obs files and remove rejected obs (bit 11) from file (obsFileClean)
  call obsf_writeFiles(obsSpaceData)
  call obsf_cleanObsFiles()
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! 4. Job termination

  istamp = exfin('THINNING', 'FIN', 'NON')

  !- deallocate obsSpaceData
  call obs_finalize(obsSpaceData)

  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

end program midas_thinning
