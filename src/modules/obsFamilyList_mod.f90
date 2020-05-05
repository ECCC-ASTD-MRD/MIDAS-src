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

module obsFamilyList_mod
  ! MODULE varNameList (prefix='ofl' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: Contains a list of all recognizable observation families along with 
  !           additional information and procedures regarding these families.
  !
  use utilities_mod
  
  implicit none
  save
  private

  ! public variables (parameters)
  public :: ofl_numFamily
  public :: ofl_familyList

  ! public procedures
  public :: ofl_isFamilyTypeInList

  integer,          parameter :: ofl_numFamily = 15
  character(len=2), parameter :: ofl_familyList(ofl_numFamily)= (/ &
                                 'UA','AI','SF','SC','SW','PR','RO','GP','RA', &
                                 'TO','CH','TM','AL','GL','HY'/)

  ! Description of obs family types
  ! -------------------------------
  ! UA - radiosondes (RAOBS and SURFACE)
  ! AI - aircraft measurements (AIREPS)
  ! SF - surface obs (SURFACE)
  ! SC - scatterometer measurements (SURFACE)
  ! SW - winds from satellite measurements (SATWINDS) 
  ! PR - ground-based profiler data 
  ! RO - radio occultation data (GPS-RO; TT/HU/P0)
  ! GP - ground-based GPS data (GB-GPS)
  ! RA - radar precipitation data
  ! TO - brightness temperatures (TT/HU/P0/constituents)
  ! CH - retrieved chemical constituent data
  ! TM - Sea-surface temperature data (SST)
  ! AL - Aladin lidar horizontal wind data
  ! GL - Sea-ice concentration data 
  ! HY - Hydrological data
  !
  ! NWP obs families having data as function of pressure level: UA,AI,SW 
  ! NWP obs families contributing to the surface NWP obs group: UA,SF,SC,GP,RA
  ! NWP obs families having data as a function of altitude: PR,AL
                                                   
  contains

    !--------------------------------------------------------------------------
    ! ofl_isFamilyTypeInList
    !--------------------------------------------------------------------------
    function ofl_isFamilyTypeInList(familyName) result(familyFound)
      !
      ! :Purpose: To identify if input obs family is part of the available list 
      !

      implicit none

      !Arguments:
      character(len=2) :: familyName
      logical          :: familyFound

      if ( any(ofl_familyList(:) == familyName) ) then
         familyFound = .true.
      else
         familyFound = .false.
      end if
      
    end function ofl_isFamilyTypeInList

end module obsFamilyList_mod