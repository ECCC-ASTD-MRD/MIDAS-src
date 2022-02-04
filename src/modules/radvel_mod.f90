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
  ! :Note: prefix not used for all public variables
  !
  use mpi_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: radvel_getlatlonHRfromRange, radvel_getRangefromH, radvel_getHfromRange 


contains 


  subroutine radvel_getlatlonHRfromRange(beamLat, beamLon, beamElevation, beamAzimuth, radarAltitude, & 
                                            beamRange, latSlant,lonSlant, beamHeight, beamDistance)
    !
    ! :Purpose: Computation of  lat-lon , height  of the trajectory
    !            along the radar beam from range of the radar beam
    !
    
    implicit none
    ! Argument 
    real(8), intent(in)  :: beamLat, beamLon, beamElevation, beamAzimuth, radarAltitude
    real(8), intent(in)  :: beamRange
    real(8), intent(out) :: LatSlant, lonSlant, beamHeight, beamDistance
    ! Local
    real(8)              :: Re 
    
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2 * (4./3.)

    !compute height of radar observation
    call radvel_getHfromRange(beamRange, radarAltitude, beamElevation, beamHeight)

    ! distance following surface of the earth from Doviak and Zrnic (2.28c)
    beamDistance = atan(beamRange*cos(beamElevation)/(beamRange*sin(beamElevation)+Re+radarAltitude))*Re

    ! lat lon of the path along the radar beam  
    latSlant = asin( sin(beamLat)*cos(beamDistance/ec_wgs_R2) + cos(beamLat)*sin(beamDistance/ec_wgs_R2)*cos(beamAzimuth))
    lonSlant = beamLon + atan2(sin(beamAzimuth)*sin(beamDistance/ec_wgs_R2)*cos(beamLat), cos(beamDistance/ec_wgs_R2)-sin(beamLat)*sin(latSlant))

  end subroutine radvel_getlatlonHRfromRange

  subroutine radvel_getHfromRange(beamRange, radarAltitude, beamElevation, beamHeight)
    !
    ! :Purpose: Computation of height of the radar beam
    !            from range of the radar beam
    !
    implicit none
    ! Argument 
    real(8) , intent(in)  :: beamRange, radarAltitude, beamElevation
    real(8) , intent(out) :: beamHeight
    ! Local
    real(8)               :: Re 
   
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2*(4./3.)
    ! height of radar beam  from range at beamElevation and radarAltitude 
    beamHeight = sqrt(beamRange**2.+(Re+radarAltitude)**2.+2.*beamRange*(Re+radarAltitude)*sin(beamElevation))-(Re)

  end subroutine radvel_getHfromRange

  subroutine radvel_getRangefromH(beamHeight, radarAltitude, beamElevation, beamRange)
    !
    ! :Purpose: Computation of range of the radar beam from height of the radar beam
    !
    implicit none
    ! Argument
    real(8) , intent(in)  :: beamHeight,  radarAltitude, beamElevation
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
     
  end subroutine radvel_getRangefromH

end module radvel_mod
