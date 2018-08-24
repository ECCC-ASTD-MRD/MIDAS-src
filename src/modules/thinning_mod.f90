
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

!--------------------------------------------------------------------------
!! MODULE thinning (prefix="thn")
!!
!! *Purpose*: To group all the thinning methods in a single fortran module
!!
!--------------------------------------------------------------------------
module thinning_mod
!!$  use timeCoord_mod
!!$  use obsFiles_mod
  use bufr_mod
  use obsSpaceData_mod
    use utilities_mod
    use fSQLite
  implicit none
  private
  public thn_thinAladin

contains

!!$!tim_getStepObsIndex() ! replaces stepobs(), 
!!$                       ! but last arg is now numStep rather than dddt
!!$
!!$get_analysis_time() ! in BGCKALT record, used by iasi, scat, satwinds
!!$
!!$element_list = (/ values according to family type /)
!!$
!!$IdentifyObsBoxForAllObservations()
!!$  {! Loop over observations once and do this for each obs:
!!$      identify_time_bin_for_the_obs() --> call tim_getStepObsIndex()
!!$      identify_layer_for_the_obs()
!!$      identify_latlon_bin_for_the_obs()
!!$      separa() de Peter H. - for calculating deltrad (iasi)
!!$      [also calculate the time, layer, distance differences?]
!!$   return
!!$  }
!!$
!!$set_spatial_grid_boxes()  ! used by iasi, scat

  subroutine thn_thinAladin(obsdat)
    implicit none

    ! ARGUMENTS
    type(struct_obs), intent(inout) :: obsdat


    ! deltrad:  iasi also uses this in addition to thinningSteps,
    !           thinningDistance.  The radial distance from box centre.
    ! scat:  uses deltax1 and deltax2, the latter being for the resolution of
    !        circumpolar orbits (???)


    ! NAMELIST VARIABLES
    integer :: nulnam
    real :: stepHours ! assimilation step size
    integer :: keepNthVertical ! keep every nth vertical datum

    ! Distance that limits to one observation per box.
    ! Depending on the strategy, this distance is either between two observations
    ! or between an observation and the thinning-box centre.
    real :: thinningDistance

    namelist /thin_aladin/stepHours, thinningDistance, keepNthVertical


    ! LOCAL VARIABLES
    integer :: fnom, fclos, ierr, flag
!!$    integer :: numStep ! number of assimilation steps in the assimilation window
!!$    integer :: obsDateStamp



    ! Read the namelist for Aladin observations
    nulnam = 0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    if(ierr.ne.0) call utl_abort('midas-obsSelect: Error opening file flnml')
    read(nulnam,nml=thin_aladin,iostat=ierr)
    if(ierr.ne.0) call utl_abort('midas-obsSelect: Error reading namelist')
    write(*,nml=thin_aladin)
    ierr=fclos(nulnam)

    ! For now, Aladin needs no thinning.  Do nothing.
    ! For now, Aladin needs no thinning.  Do nothing.
    ! For now, Aladin needs no thinning.  Do nothing.

    call keepNthObs(obsdat, 'AL', keepNthVertical)


!!$    ! Calculate the number of time bins in the assimilation window
!!$    numStep = 2*nint((3.d0 - stepHours/2.d0)/stepHours) + 1
!!$    ! Calculate the number of horizontal boxes (for iasi, at least)
!!$    ! Calculate the number of vertical boxes (for iasi, at least)
!!$
!!$    obsf_cfilnam(1)='obsal_0001_0001'
!!$    call obsf_setup(obsDateStamp, 'thinning')
!!$    call obsf_readFiles(obsSpaceData)
!!$
!!$    Create outFileName
!!$
!!$    Scan le fichier d'entree pour chercher:
!!$      - la date et heure de l'enregistrement resume (date, temps)
!!$      - le nombre de rapports (nb_rpts)
!!$      - le nombre maximum d''observations par rapport (nb_obs_max)
!!$
!!$    ! ALGORITHM DISCREPANCY, depending on observation type:
!!$    ! =====================================================
!!$    Allocate arrays based on nb_rpts, nb_obs_max ! i.e. obs reports
!!$    (otherwise de- and re-allocate on each pass through the loop)
!!$    ! NOTE:  scat allocates based on horizontal, vertical, temporal boxes;
!!$    !        i.e. observation cubes.  But it does not attempt to save all
!!$    !        observation attributes for each box.
!!$    ! whereas satwinds, radiosondes, iasi allocate based on observation reports
!!$    ! (but iasi re-allocates for every record, not for all obs in the file)
!!$    !         --> during selection, SW loops over the observations and rejects
!!$    !             the obs if it is within a radius of distance_thin of another
!!$    !             selected obs
!!$    !         --> during selection, radiosondes loops over the observations
!!$    !             applying selection criteria, keeping only the best in a box
!!$    !         --> during selection, iasi loops over observations and keeps all
!!$    !             obs that satisfy its criteria
!!$    !         --> for each BURP report, scat loops over the observations and
!!$    !             calls scat_select().  Scat_select identifies the observation
!!$    !             cube that fits the obs and checks to see whether the
!!$    !             observation is best so far for that cube
!!$    !
!!$    ! Because there are different approaches, I cannot do the allocation
!!$    ! automatically in a service of this module.
!!$
!!$    ierr = newdate(middleDateStamp, date, temps*10000, 3)
!!$    call tim_getstamplist(dateStampList, numStep, middleDateStamp)
!!$
!!$    do records in inFileName
!!$      Ignore some records, based on idtyp
!!$      read pertinent values (different for each obs type)
!!$      adjust longitude value
!!$      call identify_time_bin_for_the_obs()
!!$      call identify_layer_for_the_obs()
!!$      record all of these values in arrays
!!$    end do
!!$
!!$    call thinIt_aladin()
!!$    ! satwinds:  calls it 'BRUTE'
!!$    ! radiosondes:  calls it simply 'thinning_model', 'thinning_es', 
!!$    !                               'backlisting_ecmwf'
!!$    ! iasi:  called simply 'thinning' in one monolithic routine
!!$    !        (thinning rules given in comments at beginning of file)
!!$    ! scat:  called 'scat_select'
!!$
!!$    Calculate the number of observations per report for the output file
!!$
!!$    do selected reports
!!$      Transfer data (from inFileName [or from memory?] to outFileName
!!$                    (different data for each obs type)
!!$    end do
!!$
!!$    deallocate
  end subroutine thn_thinAladin


!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


!--------------------------------------------------------------------------
!! SUBROUTINE keepNthObs of module thinning_mod
!!
!! *Purpose*: Of the observations in a column that have not already been
!!            rejected, keep every nth observation and throw out the rest.
!!
!!            Set bit 9 of OBS_FLG on observations that are to be rejected.
!!
!--------------------------------------------------------------------------
  subroutine keepNthObs(obsdat, familyType, keepNthVertical)
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    integer, parameter :: BIT9=int(Z'200')
    integer :: headerIndex, bodyIndex, bodyIndex2, bodyIndexStart, bodyIndexEnd
    integer :: flag
    integer :: countKeepN ! count to keep every Nth observation in the column
    integer :: newProfileId
    real :: level

    countKeepN=0

    ! Loop over all header indices of the family of interest
    call obs_set_current_header_list(obsdat, familyType)
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      ! Loop over all body indices (still in the same family)
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        ! if datum already rejected, ignore it
        flag           = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex  )
        if(btest(flag,BIT9))cycle BODY

        headerIndex    = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
        bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN , headerIndex)
        bodyIndexEnd =   obs_headElem_i(obsdat, OBS_NLV , headerIndex) &
                       + bodyIndexStart - 1
        level          = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex  )

        ! Obtain supplementary parameters that are stored as observations
        BODY_SUPP: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
          if (       obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex2 ) == BUFR_NEPR &
              .and.  obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex2 ) == level &
             ) then
            newProfileId = obs_bodyElem_i(obsdat, OBS_VAR, bodyIndex2 )
          end if
        end do BODY_SUPP

        countKeepN=countKeepN + 1

        if(     countKeepN == keepNthVertical &
           .OR. new_column() &
          )then
          ! Reset the counter and keep this observation
          countKeepN=0
        else
          ! Reject this observation
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(flag, BIT9))
        end if

      end do BODY
    end do HEADER


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

      if(newProfileId /= previousProfileId)then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column
  end subroutine keepNthObs

END module thinning_mod
