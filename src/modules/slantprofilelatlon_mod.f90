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
!! MODULE slantprofilelatlon (prefix="slp" category='7. Low-level data objects and utilities')
!!
!! *Purpose*: calculation of latitudes/longitudes on slant-path based on ColumnData.
!!
!! @Author M. Bani Shahabadi, June 2018 
!
!--------------------------------------------------------------------------
module slantprofilelatlon_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod 
  use mpivar_mod

  implicit none
  save
  private

  ! public procedures
  !public :: slp_calcLatLonTrial, slp_calcLatLonAnal
  public :: slp_calcLatLon

contains 

!  subroutine slp_calcLatLonTrial
!  end subroutine slp_calcLatLonTrial
!
!  subroutine slp_calcLatLonAnal
!  end subroutine slp_calcLatLonAnal

  subroutine slp_calcLatLon(column,obsSpaceData)
    !
    !**s/r slp_calcLatLon - Computation of lat/lon on the slant path
    !                 for radiance observations
    !
    !*    Purpose:  -To replace the vertical fields in column
    !                with line-of-sight fields.
    !
    implicit none
    type(struct_columnData) :: column
    type(struct_obs) :: obsSpaceData

    integer :: indexHeader, varLevelIndex, levelIndex
    integer :: numLevels
    integer :: ioout
    real(8) :: obsLat, obsLatRad, obsLon, obsLonRad, altGeometric
    real(8) :: satZen, satZenRad, satAzim, satAzimRad, elevAngRad
    real(8) :: latSlantPathRad, lonSlantPathRad, distAlongPath
    real(8) :: obsCordGlb(3), slantPathCordGlb(3), unitx(3), unity(3), unitz(3), unitSatLoc(3), unitSatGlb(3)
    real(8), allocatable :: latSlantPath(:), lonSlantPath(:), gzColInMetres(:)
    real(8), pointer :: GZ_column(:)
    character(len=2) :: varLevel
    character(len=*) :: varName

    print *, 'start slp_calcLatLon:'

    ioout = 100 + mpi_myid

    ! Loop over all header indices of the 'TO' family:
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      indexHeader = obs_getHeaderIndex(obsSpaceData)
      if (indexHeader < 0) exit HEADER

      ! read obsLat/obsLon/angles from obsSpaceData header
      obsLatRad = obs_headElem_r(obsSpaceData,OBS_LAT,indexHeader)
      obsLonRad = obs_headElem_r(obsSpaceData,OBS_LON,indexHeader)
      satAzim = obs_headElem_r(obsSpaceData,OBS_AZA,indexHeader)
      satZen = obs_headElem_r(obsSpaceData,OBS_SZA,indexHeader)

      ! convert angles to radian unit
      satAzimRad = satAzim * MPC_RADIANS_PER_DEGREE_R8
      satZenRad = satZen * MPC_RADIANS_PER_DEGREE_R8
      elevAngRad = 0.5d0 * MPC_PI_R8 - satZenRad

      obsCordGlb  = RA * (/ cos(obsLatRad)*cos(obsLonRad) , cos(obsLatRad)*sin(obsLonRad) , sin(obsLatRad) /)
      unitz = (/  cos(obsLatRad)*cos(obsLonRad) , cos(obsLatRad)*sin(obsLonRad)  , sin(obsLatRad) /)
      unitx = (/ -sin(obsLonRad) , cos(obsLonRad) , 0.d0 /)
      unity = (/ -sin(obsLatRad)*cos(obsLonRad) , -sin(obsLatRad)*sin(obsLonRad) , cos(obsLatRad) /)

      ! unit vector towards satellite in local coordinate
      unitSatLoc = (/ cos(elevAngRad)*sin(satAzimRad) , cos(elevAngRad)*cos(satAzimRad) , sin(elevAngRad) /)
      ! unit vector towards satellite in global coordinate
      unitSatGlb = unitSatLoc(1) * unitx + unitSatLoc(2) * unity + unitSatLoc(3) * unitz
      ! loop through thermo/momentum levels
      do varLevelIndex = 1, 2
        if (varLevelIndex == 1) then
          varLevel = 'TH'
          varName = 'GZ_T'
        else
          varLevel = 'MM'
          varName = 'GZ_M'
        end if

        numLevels = col_getNumLev(column,varLevel)

        allocate(latSlantPath(numLevels))
        allocate(lonSlantPath(numLevels))
        allocate(gzColInMetres(numLevels))
        latSlantPath(:) = 0.0
        lonSlantPath(:) = 0.0
        gzColInMetres(:) = 0.0

        ! unit: m
        GZ_column  => col_getColumn(column,indexHeader,varName,varLevel)
        gzColInMetres(:) = GZ_column(:)

        do levelIndex = 1, numLevels
          ! Geometric altitude (m)
          altGeometric = RA * gzColInMetres(levelIndex) / (RA - gzColInMetres(levelIndex))
          ! distance along line of sight (m)
          distAlongPath = altGeometric / cos(satZenRad) 

          ! unit: m
          slantPathCordGlb(:) = obsCordGlb(:) + distAlongPath * unitSatGlb(:) 

          latSlantPathRad = atan(slantPathCordGlb(3)/sqrt(slantPathCordGlb(1)**2+slantPathCordGlb(2)**2))
          lonSlantPathRad = atan2(slantPathCordGlb(2),slantPathCordGlb(1))

          ! convert to deg
          latSlantPath(levelIndex) = latSlantPathRad * MPC_DEGREES_PER_RADIAN_R8
          lonSlantPath(levelIndex) = lonSlantPathRad * MPC_DEGREES_PER_RADIAN_R8
        end do

        ! print output angles/lat/lon/GZ
        if (varLevelIndex == 1) then
          write(ioout,'(4(f9.4,1x))') obsLatRad, obsLonRad, satAzimRad, satZenRad
        end if

        write(ioout,'(a2)') varLevel

        do levelIndex = 1, numLevels
          write(ioout,'(f9.4,1x)',advance='no') latSlantPath(levelIndex)
        end do
        write(ioout,*)
        do levelIndex = 1, numLevels
          write(ioout,'(f9.4,1x)',advance='no') lonSlantPath(levelIndex)
        end do
        write(ioout,*)

        do levelIndex = 1, numLevels
          write(ioout,'(f12.4,1x)',advance='no') gzColInMetres(levelIndex)
        end do
        write(ioout,*)

        deallocate(latSlantPath)
        deallocate(lonSlantPath)
        deallocate(gzColInMetres)
      end do

      write(ioout,*)

    end do HEADER

  end subroutine slp_calcLatLon

end module slantprofilelatlon_mod
