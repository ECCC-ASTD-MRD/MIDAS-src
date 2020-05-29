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

module thinning_mod
  ! MODULE thinning (prefix='thn' category='1. High-level functionality')
  !
  ! :Purpose: Using observation-type-specific algorithms, set bit 11 of 'flag'
  !           on data that are not to be assimilated.
  !
  ! :Note:    This module is intended to group all of the thinning methods in a
  !           single fortran module.
  !
  !           So far, only aladin wind data are treated.
  !
  use mpi_mod
  use bufr_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use obsSpaceData_mod
  use timeCoord_mod
  use codtyp_mod
  use physicsFunctions_mod
  use utilities_mod
  use kdtree2_mod
  implicit none
  private

  public :: thn_thinAladin, thn_thinHyper, thn_thinTovs

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  ! thn_thinAladin
  !--------------------------------------------------------------------------
  subroutine thn_thinAladin(obsdat)
    implicit none

    ! ARGUMENTS
    type(struct_obs), intent(inout) :: obsdat

    ! NAMELIST VARIABLES
    integer :: keepNthVertical ! keep every nth vertical datum

    namelist /thin_aladin/keepNthVertical

    ! LOCAL VARIABLES
    integer :: nulnam
    integer :: fnom, fclos, ierr

    keepNthVertical=-1

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_aladin','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error opening file flnml')
      read(nulnam,nml=thin_aladin,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_aladin)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAladin: Namelist block thin_aladin is missing in the namelist.'
      write(*,*) '                The default value will be taken.'
    end if

    if (keepNthVertical > 0) then
      call thn_keepNthObs(obsdat, 'AL', keepNthVertical)
    end if

  end subroutine thn_thinAladin

  !--------------------------------------------------------------------------
  ! thn_thinTovs
  !--------------------------------------------------------------------------
  subroutine thn_thinTovs(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Namelist variables
    integer :: delta    ! 

    namelist /thin_tovs/delta

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Default namelist values
    delta = 100

    ! Read the namelist for TOVS observations (if it exists)
    if (utl_isNamelistPresent('thin_tovs','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error opening file flnml')
      read(nulnam,nml=thin_tovs,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error reading namelist')
      write(*,nml=thin_tovs)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinTovs: Namelist block thin_tovs is missing in the namelist.'
      write(*,*) '              The default value will be taken.'
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, 'TO', codtyp_get_codtyp('amsua'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, 'TO', codtyp_get_codtyp('amsub'), &
                      codtyp2_opt=codtyp_get_codtyp('mhs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, 'TO', codtyp_get_codtyp('atms'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, 'TO', codtyp_get_codtyp('mwhs2'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinTovs

  !--------------------------------------------------------------------------
  ! thn_thinHyper
  !--------------------------------------------------------------------------
  subroutine thn_thinHyper(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Namelist variables
    logical :: removeUnCorrected ! indicate if obs without bias correction should be removed
    integer :: deltmax           ! time window by bin (from bin center to bin edge) (in minutes)
    integer :: deltax            ! thinning (dimension of box sides) (in km)
    integer :: deltrad           ! radius around box center for chosen obs (in km)
    namelist /thin_hyper/removeUnCorrected, deltmax, deltax, deltrad

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Default namelist values
    removeUnCorrected = .true.
    deltmax = 22
    deltax  = 150
    deltrad = 45

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_hyper','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error opening file flnml')
      read(nulnam,nml=thin_hyper,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error reading namelist')
      write(*,nml=thin_hyper)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinHyper: Namelist block thin_hyper is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('airs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('iasi'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('cris'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('crisfsr'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinHyper

!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


  !--------------------------------------------------------------------------
  ! thn_keepNthObs
  !--------------------------------------------------------------------------
  subroutine thn_keepNthObs(obsdat, familyType, keepNthVertical)
    !
    ! :Purpose: Of the observations in a column that have not already been
    !           rejected, keep every nth observation and throw out the rest.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    ! Locals:
    integer, parameter :: PROFILE_NOT_FOUND=-1
    integer :: headerIndex, bodyIndex, bodyIndex2, bodyIndexStart, bodyIndexEnd
    integer :: flag
    integer :: countKeepN ! count to keep every Nth observation in the column
    integer :: newProfileId

    countKeepN=0

    ! Loop over all body indices (columns) of the family of interest and
    ! thin each column independently of the others
    call obs_set_current_body_list(obsdat, familyType)
    BODY: do 
      bodyIndex = obs_getBodyIndex(obsdat)
      if (bodyIndex < 0) exit BODY

      ! If the datum is not being assimilated, ignore it
      if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle BODY

      ! If datum already rejected, ignore it
      flag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
      if (btest(flag,9)) cycle BODY

      headerIndex  = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
      newProfileId = obs_headElem_i(obsdat, OBS_PRFL, headerIndex)

      countKeepN=countKeepN + 1
      if ( countKeepN == keepNthVertical .or. &
           new_column() ) then
        ! Reset the counter and keep this observation
        countKeepN=0
      else
        ! Reject this observation
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(flag,11))
      end if

    end do BODY

  contains
    function new_column()
      ! Determine whether the current observation begins a new vertical column
      ! (Assume that observations are in chronological order)
      !
      ! Note:  This method has been written with aladin observations in mind.
      !        It might be necessary to generalize the method.
      implicit none
      logical :: new_column

      integer, save :: previousProfileId=huge(previousProfileId)

      if (newProfileId == PROFILE_NOT_FOUND) then
        ! The profile ID for this element is missing.
        ! Assume that it is the same as the previous element
        newProfileId = previousProfileId
      end if

      if (newProfileId /= previousProfileId) then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column

  end subroutine thn_keepNthObs

  !--------------------------------------------------------------------------
  ! thn_tovsFilt
  !--------------------------------------------------------------------------
  subroutine thn_tovsFilt(obsdat, delta, familyType, codtyp, codtyp2_opt)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: delta
    character(len=*), intent(in)    :: familyType
    integer, intent(in)             :: codtyp
    integer, optional, intent(in)   :: codtyp2_opt

    ! Locals:
    integer :: nblat, nblon, headerIndex, headerIndexKeep, latIndex, lonIndex, latIndex2
    integer :: gridIndex, ngrdtot, obsTime, obsDate, numHeader, ierr
    integer :: bodyIndex, headerIndex1, headerIndex2, stepIndex, obsCount, obsIndex
    integer :: loscan, hiscan, icntqc, obsFlag, obsFov, icntptg, nbstn2
    integer :: icntkept, icntothr, procIndex, procIndexKeep, minLonBurpFile
    integer :: allMinLonBurpFile(mpi_nprocs)
    logical :: global1, global2
    real(4) :: zlength, zlatrad, zdlon, zlatchk, zlonchk, zlontmp
    real(4) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(4) :: obsLatInDegrees, obsLonInDegrees
    real(4) :: zlat1, zlon1, zlat2, zlon2, distance, minDistance
    real(4) :: allMinDistance(mpi_nprocs)
    real(4) :: zdatflgs, zlatmid, zlonmid, zlatobs, zlonobs
    real(4) :: pcttot, pctqc, pctoth, pctkep
    real(8) :: dlhours
    real(8), allocatable :: stepObsIndex(:)
    real(4), allocatable :: zlatdeg(:), zlatg(:), zlong(:), zdobs(:)
    logical, allocatable :: valid(:), indexs(:)
    integer, allocatable :: ngrd(:), nbstn(:), igrds(:), zcentregen(:)
    integer, allocatable :: obsLonBurpFile(:), obsLatBurpFile(:), nassim(:)
    integer, allocatable :: obsDateStamp(:), bufref(:), link(:)
    integer, allocatable :: headerIndexList(:), headerIndexList2(:)
    integer, external :: newdate
    ! Locals related to kdtree2
    type(kdtree2), pointer            :: tree
    integer, parameter                :: maxNumSearch = 100
    integer                           :: numFoundSearch, resultIndex
    type(kdtree2_result)              :: searchResults(maxNumSearch)
    real(kdkind)                      :: maxRadius = 100.d0
    real(kdkind)                      :: refPosition(3)
    real(kdkind), allocatable         :: obsPosition3d(:,:)

    ! Local parameters:
    integer, parameter :: xlat=10000, xlon=40000, xdelt=50
    integer, parameter :: centre_gen_global(3)=(/53, 74, 160/)
    integer, parameter :: mxscanamsua=30
    integer, parameter :: mxscanamsub=90
    integer, parameter :: mxscanatms =96
    integer, parameter :: mxscanmwhs2=98

    write(*,*) 'thn_tovsFilt: Starting, ', trim(codtyp_get_name(codtyp))

    numHeader = obs_numHeader(obsdat)

    ! Check if we have any observations to process
    allocate(valid(numHeader))
    valid(:) = .false.
    do headerIndex = 1, numHeader
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp) then
        valid(headerIndex) = .true.
      else if (present(codtyp2_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp2_opt) then
          valid(headerIndex) = .true.
        end if
      end if
    end do
    if (count(valid(:)) == 0) then
      write(*,*) 'thn_tovsFilt: no observations for this instrument'
      return
    end if

    write(*,*) 'obsCount initial = ', count(valid(:))

    nblat = 2*xlat/delta
    nblon = xlon/delta

    ! Allocations
    allocate(zlatdeg(nblat))
    allocate(zlatg(nblat*nblon))
    allocate(zlong(nblat*nblon))
    allocate(nbstn(nblat*nblon))
    allocate(bufref(nblat*nblon))
    allocate(link(numHeader))
    allocate(headerIndexList(numHeader))
    allocate(headerIndexList2(numHeader))
    allocate(ngrd(nblat))
    allocate(igrds(numHeader))
    allocate(indexs(numHeader))
    allocate(zcentregen(numHeader))
    allocate(obsDateStamp(numHeader))
    allocate(obsLonBurpFile(numHeader))
    allocate(obsLatBurpFile(numHeader))
    allocate(nassim(numHeader))
    allocate(zdobs(numHeader))
    allocate(stepObsIndex(numHeader))
    allocate(obsPosition3d(3,numHeader))

    write(*,*) 'thn_tovsFilt: after alloction'

    ! Initialize some arrays
    nbstn(:) = 0
    zcentregen(:) = 0
    stepObsIndex(:) = 0.0d0
    obsPosition3d(:,:) = 0.0

    ! Set up the grid used for thinning
    ngrdtot = 0
    do latIndex = 1, nblat
      zlatdeg(latIndex) = (latIndex*180./nblat) - 90.
      zlatrad           = zlatdeg(latIndex) * mpc_pi_r4 / 180.
      zlength           = xlon * cos(zlatrad)
      ngrd(latIndex)    = nint(zlength/delta)
      ngrd(latIndex)    = max(ngrd(latIndex),1)
      ngrdtot           = ngrdtot + ngrd(latIndex)
    end do

    gridIndex = 0
    do latIndex = 1, nblat
      zdlon = 360./ngrd(latIndex)
      do lonIndex = 1, ngrd(latIndex)
        zlatg(gridIndex+lonIndex) = zlatdeg(latIndex)
        zlong(gridIndex+lonIndex) = (lonIndex-1)*zdlon
      end do
      gridIndex = gridIndex + ngrd(latIndex)
    end do

    write(*,*) 'thn_tovsFilt: after grid setup'

    ! Loop over all observation locations
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      ! Originating centre of data
      zcentregen(headerIndex) = obs_headElem_i(obsdat, OBS_ORI, headerIndex)

      ! Date stamp for each observation
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      ierr = newdate(obsDateStamp(headerIndex), obsDate, obsTime*10000+2900, 3)

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile(headerIndex) = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile(headerIndex) < 0) then
        obsLonBurpFile(headerIndex) = obsLonBurpFile(headerIndex) + 36000
      end if
      obsLatBurpFile(headerIndex) = 9000+nint(100.0*obsLatInDegrees)

      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(stepObsIndex(headerIndex), tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      ! location array for kdtree
      obsPosition3d(1,headerIndex) = RA * sin(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(2,headerIndex) = RA * cos(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(3,headerIndex) = RA *                    sin(obsLatInRad)

      ! Associate each observation to a grid point
      zlatchk = (obsLatBurpFile(headerIndex) - 9000.) / 100.
      zlonchk = obsLonBurpFile(headerIndex) / 100.
      zlontmp = zlonchk
      if (zlonchk > 180.) zlonchk = zlonchk - 360.
      do latIndex = 1,nblat-1
        if (zlatchk <  zlatdeg(1))     zlatchk = zlatdeg(1)
        if (zlatchk >= zlatdeg(nblat)) zlatchk = zlatdeg(nblat) - 0.5
        if (zlatchk >= zlatdeg(latIndex) .and. zlatchk < zlatdeg(latIndex+1)) then
          gridIndex = 1
          do latIndex2 = 1, latIndex-1
            gridIndex = gridIndex + ngrd(latIndex2)
          end do
          zdlon = 360./ngrd(latIndex)
          gridIndex = gridIndex + ifix(zlontmp/zdlon)
          exit
        end if
      end do
      igrds(headerIndex) = gridIndex
      nbstn(gridIndex) = nbstn(gridIndex) + 1

    end do ! headerIndex

    write(*,*) 'thn_tovsFilt: after first pass'

    ! start of "passe2"
    nullify(tree)
    tree => kdtree2_create(obsPosition3d, sort=.true., rearrange=.true.)
    HEADER1: do headerIndex1 = 1, (numHeader-1)
        
      if ( .not. valid(headerIndex1) ) cycle HEADER1

      ! Find all obs within 10km
      refPosition(:) = obsPosition3d(:,headerIndex1)
      call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numFoundSearch,&
                             nalloc=maxNumSearch, results=searchResults)
      if (numFoundSearch >= maxNumSearch) then
        call utl_abort('thn_tovsFilt: the parameter maxNumSearch must be increased')
      end if
      if (numFoundSearch == 0) then
        call utl_abort('thn_tovsFilt: no match found. This should not happen!!!')
      end if

      ! Loop over all of these nearby locations
      HEADER2: do resultIndex = 1, numFoundSearch
        headerIndex2 = searchResults(resultIndex)%idx
          
        if ( .not. valid(headerIndex2) ) cycle HEADER2

        ! Certaines stations locales nous envoient 
        ! le mauvais numero d'orbite. On ne peut donc
        ! pas s'y fier.
        ! Il faut comparer le temps de la reception des
        ! donnees
        if ( zcentregen(headerIndex1) /= zcentregen(headerIndex2) ) then
          if ( (obs_headElem_i(obsdat,OBS_FOV,headerIndex1) ==  &
                obs_headElem_i(obsdat,OBS_FOV,headerIndex2)) ) then
            if ( (obs_elem_c(obsdat,'STID',headerIndex1) ==  &
                  obs_elem_c(obsdat,'STID',headerIndex2)) ) then
            
              ! Difference (in hours) between obs time
              call difdatr(obsDateStamp(headerIndex1),obsDateStamp(headerIndex2),dlhours)

              ! Si la difference est moins de 6 minutes,
              ! on peut avoir affaire a un rars

              if ( abs(dlhours) <= 0.1 ) then

                ! ON RAMENE LES LATS ENTRE -90 ET 90
                ! ET LES LONS ENTRE 0 ET 360
                zlat1 = ( obsLatBurpFile(headerIndex1) - 9000.) / 100.
                zlon1 =   obsLonBurpFile(headerIndex1) / 100.
                zlat2 = ( obsLatBurpFile(headerIndex2) - 9000.) / 100.
                zlon2 =   obsLonBurpFile(headerIndex2) / 100.

                ! SI LA DISTANCE ENTRE DEUX STNS < MINIMUM
                ! ETABLI ICI A 10KM (PAR EXPERIENCE)
                distance = separa(zlon1,zlat1,zlon2,zlat2) * float(xlat) / 90.
                if (distance <= 10.) then

                  ! si l'element_i est global, on doit le garder et rejeter l'element_j
                  global1 = any(centre_gen_global(:) == zcentregen(headerIndex1))
                  if (global1) then 
                    valid(headerIndex2) = .false.
                  else
                    ! toutefois, ca ne signifie pas que l'element_j est un rars
                    ! VERIFIER SI LA STATION 2 EST RARS
                    global2 = any(centre_gen_global(:) == zcentregen(headerIndex2))

                    ! Si l'element_j est global, rejeter l'element_i
                    ! Si les 2 elements sont rars, garder le 1er
                    if (global2) then 
                      valid(headerIndex1) = .false.
                      cycle HEADER1
                    else
                      valid(headerIndex2) = .false.
                    end if

                  end if

                end if ! distance <= 10

              end if ! abs(dlhours) <= 0.1

            end if
          end if
        end if

      end do HEADER2
    end do HEADER1
    call kdtree2_destroy(tree)

    write(*,*) 'thn_tovsFilt: obsCount after second pass = ', count(valid(:))
    ! end of passe2

    ! start of passe1
    if      ( codtyp == codtyp_get_codtyp('amsua') ) then
      loscan   = 1
      hiscan   = mxscanamsua
    else if ( codtyp == codtyp_get_codtyp('amsub') ) then
      loscan   = 1
      hiscan   = mxscanamsub
    else if ( codtyp == codtyp_get_codtyp('atms') ) then
      loscan   = 2
      hiscan   = mxscanatms - 1
    else if ( codtyp == codtyp_get_codtyp('mwhs2') ) then
      loscan   = 1
      hiscan   = mxscanmwhs2
    end if

    icntqc = 0
    nassim(:) = 0
    do headerIndex = 1, numHeader

      if ( .not. valid(headerIndex) ) cycle

      ! Look at the obs flags
      zdatflgs = 0.

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY
        
        ! If not a blacklisted channel (note that bit 11 is set in 
        ! satqc_amsu*.f for blacklisted channels)
        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        if ( .not. btest(obsFlag,11) ) then
          nassim(headerIndex) = nassim(headerIndex) + 1
          if ( btest(obsFlag,9) ) then
            zdatflgs = zdatflgs + 1.0
          end if
        end if
      end do BODY

      ! fixer le % de rejets a 100% si aucun canal n'est assimilable         
      if ( zdatflgs == 0. .and. nassim(headerIndex) == 0 ) then
        zdatflgs = 1.
      else
        zdatflgs = zdatflgs / max(nassim(headerIndex),1)  
      end if

      obsFov = obs_headElem_i(obsdat, OBS_FOV, headerIndex)
      if ( zdatflgs >= 0.80 ) then
        icntqc = icntqc + 1
        valid(headerIndex) = .false.
      else if (obsFov < loscan .or.  &
               obsFov > hiscan ) then
        icntqc = icntqc + 1
        valid(headerIndex) = .false.
      end if

    end do
    ! end of passe1

    write(*,*) 'thn_tovsFilt: obsCount after third pass = ', count(valid(:))

    ! start of passe3
    
    ! Calculate distance of obs. from center of its grid box center
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = igrds(headerIndex)
      if (nbstn(gridIndex) /= 0) then
        latIndex = (zlatg(gridIndex)+90.)/(180./nblat)
        zdlon = 360./ngrd(latIndex)
        zlatmid = zlatg(gridIndex) + 0.5*(180./nblat)
        zlonmid = zlong(gridIndex) + 0.5*zdlon
        zlatobs = (obsLatBurpFile(headerIndex) - 9000.) / 100.
        zlonobs = obsLonBurpFile(headerIndex) / 100.
        zdobs(headerIndex) = separa(zlonobs,zlatobs,zlonmid,zlatmid) * &
             float(xlat) / 90.
      end if
    end do

    ! Create a linked list of observations (link to grid point)
    bufref(:) = 0
    link(:) = 0
      
    obsCount = 0
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = igrds(headerIndex)
      if (nbstn(gridIndex) /= 0) then
        obsCount = obsCount + 1
        headerIndexList(obsCount) = headerIndex
        link(obsCount) = bufref(gridIndex)
        bufref(gridIndex) = obsCount
      end if
    end do

    ! BOUCLE SUR LES STEPOBS
    do stepIndex = 1, tim_nstepobs

      ! BOUCLE SUR TOUTES LES BOITES
      icntptg = 0
      do gridIndex = 1, ngrdtot
        if (nbstn(gridIndex) /= 0) then
          icntptg = icntptg + 1
        end if
        nbstn2 = 0
        obsIndex = bufref(gridIndex)
        do
          if (obsIndex == 0) exit
          headerIndex = headerIndexList(obsIndex)
          if ( igrds(headerIndex) == gridIndex  .and. &
               valid(headerIndex)               .and. &
               nint(stepObsIndex(headerIndex)) == stepIndex ) then
            nbstn2 = nbstn2 + 1
            headerIndexList2(nbstn2) = headerIndex
          end if

          obsIndex = link(obsIndex)
        end do

        minDistance = 1000000.             

        ! CHOISIR LA PLUS PRES DU CENTRE DE LA BOITE, EN DECA DE 75KM.
        do obsIndex = 1, nbstn2
          if (zdobs(headerIndexList2(obsIndex)) < minDistance) then
            minDistance = zdobs(headerIndexList2(obsIndex))
            minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
            headerIndexKeep = headerIndexList2(obsIndex)
          end if
        end do

        ! Check for multiple obs with same distance to grid point
        if (nbstn2 > 0) then
          if ( count(zdobs(headerIndexList2(1:nbstn2)) == minDistance) > 1 ) then
            write(*,*) 'multiple obs with same minDistance!', &
                 zdobs(headerIndexList2(1:nbstn2))
            ! resolve ambiguity by choosing obs with min value of lon
            minLonBurpFile = 10000000
            do obsIndex = 1, nbstn2
              if (zdobs(headerIndexList2(obsIndex)) == minDistance) then
                if (obsLonBurpFile(headerIndexList2(obsIndex)) < minLonBurpFile) then
                  minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
                  headerIndexKeep = headerIndexList2(obsIndex)
                end if
              end if
            end do
          end if
        end if

        do obsIndex = 1, nbstn2
          valid(headerIndexList2(obsIndex)) = .false.
        end do
        if (nbstn2 > 0 .and. minDistance <= 75. ) then
          valid(headerIndexKeep) = .true.
        end if

        ! Communicate the distance of chosen observation among all mpi tasks
        call rpn_comm_allgather(minDistance,    1, 'mpi_real4',  &
                                allMinDistance, 1, 'mpi_real4', 'grid', ierr)

        ! Choose the closest to the center of the box among all mpi tasks
        minDistance = 1000000.
        do procIndex = 1, mpi_nprocs
          if (allMinDistance(procIndex) < minDistance) then
            minDistance = allMinDistance(procIndex)
            procIndexKeep = procIndex
          end if
        end do

        ! Adjust flags to only keep 1 observation among all mpi tasks
        if (minDistance < 1000000.) then
          if ( count(allMinDistance(:) == minDistance) > 1 ) then
            write(*,*) 'multiple obs with same minDistance! (MPI) ', allMinDistance(:)
            ! resolve ambiguity by choosing obs with min value of lon
            call rpn_comm_allgather(minLonBurpFile,    1, 'mpi_integer',  &
                                    allMinLonBurpFile, 1, 'mpi_integer', 'grid', ierr)
            minLonBurpFile = 10000000
            do procIndex = 1, mpi_nprocs
              if (allMinDistance(procIndex) == minDistance) then
                if (allMinLonBurpFile(procIndex) < minLonBurpFile) then
                  minLonBurpFile = allMinLonBurpFile(procIndex)
                  procIndexKeep = procIndex
                end if
              end if
            end do            
          end if
          if (mpi_myid /= (procIndexKeep-1)) then
            if (nbstn2 > 0) then
              valid(headerIndexKeep) = .false.
            end if
          end if
        end if

      end do ! gridIndex
    end do ! stepIndex

    ! end of passe3

    write(*,*) 'thn_tovsFilt: obsCount after thinning = ', count(valid(:))

    ! start of modflg
    obsCount = 0
    do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) then
        if (present(codtyp2_opt)) then
          if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp2_opt) then
            cycle
          end if
        else
          cycle
        end if
      end if
     
      obsCount = obsCount + 1

      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY2: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY2
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY2
      end if
    end do

    ! end of modflg

    ! start of sommair
    icntkept = count(valid)
    icntothr = obsCount - count(valid) - icntqc

    pcttot = 100.0
    pctqc  = (float(icntqc)   / float(obsCount)) * 100.0
    pctoth = (float(icntothr) / float(obsCount)) * 100.0
    pctkep = (float(icntkept) / float(obsCount)) * 100.0
         
    write(*,100)
100 format(/,' SOMMAIRE DES RESULTATS',/)
    write(*,200) obsCount,pcttot,icntqc,pctqc,icntothr,pctoth, &
                 icntkept,pctkep,delta,ngrdtot,icntptg
200 format(' NB.STNS TOTAL AU DEBUT EN ENTREE        =    ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES AU CONTROLE QUALITATIF =  - ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES POUR LES AUTRES RAISONS=  - ',I7, &
           '  =  ',F6.1,'%',/, &
           '                                              ', &
           '-------- = -------',/, &
           ' NB.STNS GARDEES A LA FIN                =    ',I7, &
           '  =  ',F6.1,'%',/////, &
           ' NB.PTS GRILLE AVEC LA RESOLUTION ',I5,' KM   =    ', &
           I7,/, &
           ' NB.PTS GRILLE TRAITES                       =    ',I7)

    ! end of sommair

    deallocate(valid)
    deallocate(zlatdeg)
    deallocate(zlatg)
    deallocate(zlong)
    deallocate(nbstn)
    deallocate(bufref)
    deallocate(link)
    deallocate(headerIndexList)
    deallocate(headerIndexList2)
    deallocate(ngrd)
    deallocate(igrds)
    deallocate(indexs)
    deallocate(zcentregen)
    deallocate(obsDateStamp)
    deallocate(obsLonBurpFile)
    deallocate(obsLatBurpFile)
    deallocate(nassim)
    deallocate(zdobs)
    deallocate(stepObsIndex)
    deallocate(obsPosition3d)

    write(*,*) 'thn_tovsFilt: finished'

  end subroutine thn_tovsFilt


  FUNCTION SEPARA(XLON1,XLAT1,XLON2,XLAT2)
    IMPLICIT NONE

    REAL XLAT1,XLAT2,XLON1,XLON2,SEPARA,COSVAL,DEGRAD,RADDEG

    RADDEG=180.0/3.14159265358979
    DEGRAD=1.0/RADDEG
    COSVAL=SIN(XLAT1*DEGRAD)*SIN(XLAT2*DEGRAD)+ &
         COS(XLAT1*DEGRAD)*COS(XLAT2*DEGRAD)* &
         COS((XLON1-XLON2)*DEGRAD)

    IF (COSVAL.LT.-1.0D0) THEN
      COSVAL=-1.0D0
    ELSE IF (COSVAL.GT.1.0D0) THEN
      COSVAL=1.0D0
    ENDIF
    SEPARA=ACOS(COSVAL)*RADDEG

  END FUNCTION SEPARA


  !--------------------------------------------------------------------------
  ! thn_thinByLatLonBoxes
  !--------------------------------------------------------------------------
  subroutine thn_thinByLatLonBoxes(obsdat, removeUnCorrected, &
                                   deltmax, deltax, deltrad,  &
                                   familyType, codtyp)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical, intent(in)             :: removeUnCorrected
    integer, intent(in)             :: deltmax
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltrad
    character(len=*), intent(in)    :: familyType
    integer, intent(in)             :: codtyp

    ! Locals:
    integer :: headerIndex, bodyIndex, obsDate, obsTime, obsFlag
    integer :: nblat, nblon, latIndex, numChannels, delMinutes
    integer :: lonBinIndex, latBinIndex, timeBinIndex
    integer :: ierr, nsize, procIndex, countHeader
    real(4) :: latr, length, distance
    real(8) :: lonBoxCenterInDegrees, latBoxCenterInDegrees
    real(8) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, stepObsIndex
    integer, parameter :: LAT_LENGTH = 10000
    integer, parameter :: LON_LENGTH = 40000
    real(4), allocatable :: latdeg(:)
    integer, allocatable :: ngrd(:)
    integer, allocatable :: headerIndexKeep(:,:,:), numChannelsKeep(:,:,:)
    integer, allocatable :: allHeaderIndex(:,:,:,:), allNumChannels(:,:,:,:)
    integer, allocatable :: delMinutesKeep(:,:,:)
    integer, allocatable :: allDelMinutes(:,:,:,:)
    integer, allocatable :: procIndexKeep(:,:,:)
    real(4), allocatable :: distanceKeep(:,:,:)
    real(4), allocatable :: allDistance(:,:,:,:)
    logical, allocatable :: rejectThisHeader(:)
    logical :: keepThisObs
    integer :: obsLonBurpFile, obsLatBurpFile
    character(len=9) :: stnid

    write(*,*) 'thn_thinByLatLonBoxes: Starting, ', trim(codtyp_get_name(codtyp))

    ! Initial setup
    nblat = nint( 2. * real(LAT_LENGTH) / real(deltax) )
    nblon = nint(      real(LON_LENGTH) / real(deltax) )
    allocate(headerIndexKeep(nblat,nblon,tim_nstepobs))
    allocate(numChannelsKeep(nblat,nblon,tim_nstepobs))
    allocate(distanceKeep(nblat,nblon,tim_nstepobs))
    allocate(delMinutesKeep(nblat,nblon,tim_nstepobs))
    allocate(latdeg(nblat))
    allocate(ngrd(nblat))
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax
    latdeg(:)              = 0.0
    ngrd(:)                = 0

    allocate(allHeaderIndex(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allNumChannels(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allDistance(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allDelMinutes(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(procIndexKeep(nblat,nblon,tim_nstepobs))
    procIndexKeep(:,:,:) = -1

    ! set spatial boxes properties
    ! latdeg(:) : latitude (deg) of northern side of the box
    ! ngrd(:)   : number of longitudinal boxes at this latitude
    do latIndex = 1, nblat
      latdeg(latIndex) = (latIndex*180./nblat) - 90.
      if ( latdeg(latIndex) <= 0.0 ) then
        latr = latdeg(latIndex) * MPC_PI_R8 / 180.
      else
        latr = latdeg(latIndex-1) * MPC_PI_R8 / 180.
      end if
      length = LON_LENGTH * cos(latr)
      ngrd(latIndex)   = nint(length/deltax)
    end do

    countHeader = 0

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER

      countHeader = countHeader + 1

      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile < 0) obsLonBurpFile = obsLonBurpFile + 36000
      obsLatBurpFile = 9000+nint(100.0*obsLatInDegrees)

      numChannels = 0

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        ! mark for rejection if not bias corrected (bit 6 not set)
        if (removeUnCorrected) then
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          if (.not. btest(obsFlag,6)) then
            call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
            call obs_bodySet_i(obsdat, OBS_ASS, bodyIndex, obs_notAssimilated)
          end if
        end if

        ! count the number of accepted channels
        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) then
          numChannels = numChannels + 1
        end if

      end do BODY

      ! Determine the lat and lon bin indexes
      do latIndex = 1, nblat
        if ( obsLatInDegrees <= (latdeg(latIndex)+0.000001) ) then
          latBinIndex = latIndex
          exit
        end if
      end do
      lonBinIndex = int( obsLonBurpFile/(36000.0/ngrd(latBinIndex)) ) + 1
      if ( lonBinIndex > ngrd(latBinIndex) ) lonBinIndex = ngrd(latBinIndex)

      ! Determine the time bin index
      call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      timeBinIndex = nint(stepObsIndex)
      delMinutes = nint(60.0 * tim_dstepobs * abs(real(timeBinIndex) - stepObsIndex))

      ! Determine distance from box center
      latBoxCenterInDegrees = latdeg(latBinIndex) - 0.5 * (180./nblat)
      lonBoxCenterInDegrees = (360. / ngrd(latBinIndex)) * (lonBinIndex - 0.5)
      obsLat = (obsLatBurpFile - 9000.) / 100.
      obsLon = obsLonBurpFile / 100.
      distance = 1.0d-3 * phf_calcDistance(MPC_RADIANS_PER_DEGREE_R8 * latBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * lonBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLat, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLon )

      stnid = obs_elem_c(obsdat,'STID',headerIndex)

      ! Apply thinning criteria
      keepThisObs = .false.

      ! keep if distance to box center smaller than limit and 
      ! time under limit fixed at input and maximise number of channels 
      if ( (distance < deltrad)                                              .and. &
           (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
           (delMinutes <= deltmax) ) keepThisObs = .true.
     
      ! keep the closest to bin central time
      if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
        if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
          keepThisObs = .false.
        else
          if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            keepThisObs = .false.
          end if
         
          if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              keepThisObs = .false.
            end if
          end if

        end if
      end if

      ! save the observation if thinning criteria is satisfied
      if ( keepThisObs ) then
        headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
        distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
        delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
        numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
      end if

    end do HEADER

    ! return if no observations for this instrument
    if (countHeader == 0) then
      call deallocLocals()
      write(*,*) 'thn_thinByLatLonBoxes: no observations for this instrument'
      return
    end if

    ! communicate results to all other mpi tasks
    nsize = nblat * nblon * tim_nstepobs
    call rpn_comm_allgather(distanceKeep, nsize, 'mpi_real4',  &
                            allDistance,  nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(delMinutesKeep, nsize, 'mpi_integer',  &
                            allDelMinutes,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(numChannelsKeep, nsize, 'mpi_integer',  &
                            allNumChannels,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(headerIndexKeep, nsize, 'mpi_integer',  &
                            allheaderIndex,  nsize, 'mpi_integer', 'grid', ierr)

    ! reset arrays that store info about kept obs
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax

    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, nblon
        do latBinIndex = 1, nblat

          ! Apply thinning criteria to results from all mpi tasks
          do procIndex = 1, mpi_nprocs

            headerIndex = allHeaderIndex(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            distance    = allDistance(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            delMinutes  = allDelMinutes(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            numChannels = allNumChannels(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            
            keepThisObs = .false.

            ! keep if distance to box center smaller than limit and 
            ! time under limit fixed at input and maximise number of channels 
            if ( (distance < deltrad)                                              .and. &
                 (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
                 (delMinutes <= deltmax) ) keepThisObs = .true.
     
            ! keep the closest to bin central time
            if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
                keepThisObs = .false.
              else
                if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  keepThisObs = .false.
                end if
         
                if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                    keepThisObs = .false.
                  end if
                end if

              end if
            end if

            ! save the observation if thinning criteria is satisfied
            if ( keepThisObs ) then
              procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)   = procIndex
              headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
              distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
              delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
              numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
            end if

          end do ! procIndex

        end do
      end do
    end do

    ! determine which headerIndex values are rejected
    allocate(rejectThisHeader(obs_numheader(obsdat)))
    rejectThisHeader(:) = .true.
    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, nblon
        do latBinIndex = 1, nblat
          if (procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == mpi_myid+1) then
            headerIndex = headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)
            rejectThisHeader(headerIndex) = .false.
          end if
        end do
      end do
    end do

    ! modify flags, by setting bit 11 for those rejected
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER2

      if (.not. rejectThisHeader(headerIndex)) cycle HEADER2

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
      end do BODY2

    end do HEADER2
    deallocate(rejectThisHeader)

    call deallocLocals()

    write(*,*) 'thn_thinByLatLonBoxes: Finished'

  contains

    subroutine deallocLocals()
      implicit none

      deallocate(headerIndexKeep)
      deallocate(numChannelsKeep)
      deallocate(distanceKeep)
      deallocate(delMinutesKeep)
      deallocate(latdeg)
      deallocate(ngrd)
      deallocate(allHeaderIndex)
      deallocate(allNumChannels)
      deallocate(allDistance)
      deallocate(allDelMinutes)
      deallocate(procIndexKeep)
      
    end subroutine deallocLocals

  end subroutine thn_thinByLatLonBoxes

end module thinning_mod
