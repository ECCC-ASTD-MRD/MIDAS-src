
module radialVelocity_mod
  ! MODULE radialVelocity_mod (prefix='rdv' category='5. Observation operators')
  !
  ! :Purpose: Containing commonly used functions for the assimilation of Doppler velocity
  !
  !
  use codePrecision_mod
  use earthConstants_mod
  implicit none
  save
  private

  ! public procedures
  public :: rdv_getlatlonHRfromRange, rdv_getRangefromH


contains 


  subroutine rdv_getlatlonHRfromRange(antennaLat, antennaLon, beamElevation, beamAzimuth, & !in
                                      radarAltitude, beamRange,                           & !in
                                      latSlant, lonSlant, beamHeight, beamDistance)         !out
    !
    ! :Purpose: Computation of  lat-lon , height  of the trajectory
    !           along the radar beam from range of the radar beam
    !
    ! Angles are expressed in radians
    ! Height and distances in meters
    !
    implicit none

    ! Argument 
    real(pre_obsReal), intent(in)  :: antennaLat
    real(pre_obsReal), intent(in)  :: antennaLon
    real(pre_obsReal), intent(in)  :: beamElevation
    real(pre_obsReal), intent(in)  :: beamAzimuth
    real(pre_obsReal), intent(in)  :: radarAltitude
    real(pre_obsReal), intent(in)  :: beamRange
    real(pre_obsReal), intent(out) :: LatSlant
    real(pre_obsReal), intent(out) :: lonSlant
    real(pre_obsReal), intent(out) :: beamHeight
    real(pre_obsReal), intent(out) :: beamDistance
    ! Local
    real(pre_obsReal)              :: Re 
    
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! ec_wgs_R2 = 6371007.1809
    ! effective radius of the earth
    Re = ec_wgs_R2 * (4./3.)

    !compute height of radar observation
    beamHeight = sqrt(beamRange**2.+(Re+radarAltitude)**2.+2.*beamRange*(Re+radarAltitude)*sin(beamElevation))-(Re)

    ! distance following surface of the earth from Doviak and Zrnic (2.28c)
    beamDistance = atan(beamRange*cos(beamElevation)/(beamRange*sin(beamElevation)+Re+radarAltitude))*Re

    ! lat lon of the path along the radar beam  
    latSlant = asin( sin(antennaLat)*cos(beamDistance/ec_wgs_R2) + cos(antennaLat)*sin(beamDistance/ec_wgs_R2)*cos(beamAzimuth))
    lonSlant = antennaLon + atan2(sin(beamAzimuth)*sin(beamDistance/ec_wgs_R2)*cos(antennaLat), cos(beamDistance/ec_wgs_R2)-sin(antennaLat)*sin(latSlant))

  end subroutine rdv_getlatlonHRfromRange

  subroutine rdv_getRangefromH(beamHeight, radarAltitude, beamElevation, beamRange)
    !
    ! :Purpose: Computation of range of the radar beam from height of the radar beam
    !
    implicit none
    
    ! Argument
    real(pre_obsReal) , intent(in)  :: beamHeight
    real(pre_obsReal) , intent(in)  :: radarAltitude
    real(pre_obsReal) , intent(in)  :: beamElevation
    real(pre_obsReal) , intent(out) :: beamRange 
    ! Local
    real(pre_obsReal)               :: a, b, c, Re

    if ( radarAltitude > beamHeight ) then 
      !beamHeight is below radar antenna which may cause the equation below to return garbage
      !this happens in a few edge cases where its okay to return zero
      beamRange = 0.0
    else

      ! Radius of sphere of equal area from earthconstants_mod.f90
      ! ec_wgs_R2 = 6371007.1809
      ! effective radius of the earth
      Re = ec_wgs_R2*(4./3.)

      a = 1.
      b = 2.*(Re + radarAltitude)*sin(beamElevation)
      c = -(Re + beamHeight)**2. + (Re + radarAltitude)**2.
      ! range of radar beam from height and elevation of the radar beam 
      beamRange  = (-b + sqrt( b**2. - 4.*a*c )) / (2.*a)

    end if
  
     
  end subroutine rdv_getRangefromH

end module radialVelocity_mod
