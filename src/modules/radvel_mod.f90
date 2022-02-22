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

module radvel_mod
  ! MODULE radvel_mod (prefix='rdv' category='4. Observation operators')
  !
  ! :Purpose: Containing commonly used functions for the assimilation of Doppler velocity
  !
  !
  use mpi_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: rdv_getlatlonHRfromRange, rdv_getRangefromH, rdv_getHfromRange 


contains 


  subroutine rdv_getlatlonHRfromRange(antennaLat, antennaLon, beamElevation, beamAzimuth, &
                     radarAltitude, beamRange, latSlant,lonSlant, beamHeight, beamDistance)
    !
    ! :Purpose: Computation of  lat-lon , height  of the trajectory
    !           along the radar beam from range of the radar beam
    !
    implicit none

    ! Argument 
    real(8), intent(in)  :: antennaLat
    real(8), intent(in)  :: antennaLon
    real(8), intent(in)  :: beamElevation
    real(8), intent(in)  :: beamAzimuth
    real(8), intent(in)  :: radarAltitude
    real(8), intent(in)  :: beamRange
    real(8), intent(out) :: LatSlant
    real(8), intent(out) :: lonSlant
    real(8), intent(out) :: beamHeight
    real(8), intent(out) :: beamDistance
    ! Local
    real(8)              :: Re 
    
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2 * (4./3.)

    !compute height of radar observation
    call rdv_getHfromRange(beamRange, radarAltitude, beamElevation, beamHeight)

    ! distance following surface of the earth from Doviak and Zrnic (2.28c)
    beamDistance = atan(beamRange*cos(beamElevation)/(beamRange*sin(beamElevation)+Re+radarAltitude))*Re

    ! lat lon of the path along the radar beam  
    latSlant = asin( sin(antennaLat)*cos(beamDistance/ec_wgs_R2) + cos(antennaLat)*sin(beamDistance/ec_wgs_R2)*cos(beamAzimuth))
    lonSlant = antennaLon + atan2(sin(beamAzimuth)*sin(beamDistance/ec_wgs_R2)*cos(antennaLat), cos(beamDistance/ec_wgs_R2)-sin(antennaLat)*sin(latSlant))

  end subroutine rdv_getlatlonHRfromRange

  subroutine rdv_getHfromRange(beamRange, radarAltitude, beamElevation, beamHeight)
    !
    ! :Purpose: Computation of height of the radar beam
    !           from range of the radar beam
    !
    implicit none

    ! Argument 
    real(8) , intent(in)  :: beamRange
    real(8) , intent(in)  :: radarAltitude
    real(8) , intent(in)  :: beamElevation
    real(8) , intent(out) :: beamHeight
    ! Local
    real(8)               :: Re 
   
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2*(4./3.)
    ! height of radar beam  from range at beamElevation and radarAltitude 
    beamHeight = sqrt(beamRange**2.+(Re+radarAltitude)**2.+2.*beamRange*(Re+radarAltitude)*sin(beamElevation))-(Re)

  end subroutine rdv_getHfromRange

  subroutine rdv_getRangefromH(beamHeight, radarAltitude, beamElevation, beamRange)
    !
    ! :Purpose: Computation of range of the radar beam from height of the radar beam
    !
    implicit none
    
    ! Argument
    real(8) , intent(in)  :: beamHeight
    real(8) , intent(in)  :: radarAltitude
    real(8) , intent(in)  :: beamElevation
    real(8) , intent(out) :: beamRange 
    ! Local
    real(8)               :: a, b, c, Re
  
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2*(4./3.)

    a = 1.
    b = 2.*(Re + radarAltitude)*sin(beamElevation)
    c = -(Re + beamHeight)**2. + (Re + radarAltitude)**2.
    ! range of radar beam from height and elevation of the radar beam 
    beamRange  = (-b + sqrt( b**2. - 4.*a*c )) / (2.*a)
     
  end subroutine rdv_getRangefromH

end module radvel_mod
