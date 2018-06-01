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
!! MODULE aladin_mod (prefix="ala")
!!
!! *Purpose*: Derived type and procedures related to the aladin observations of
!!            wind.
!!
!--------------------------------------------------------------------------
module aladin_mod
  ! The instrument used to make the measurement is called "Atmospheric LAser 
  ! Doppler INstrument", or ALADIN.
  use MathPhysConstants_mod
  implicit none
  private
  
  ! public procedures
  public :: ala_aladin, ala_aladin_tl, ala_aladin_ad


  contains

  real(8) function ala_aladin(uu, vv, azimuth)
    !
    ! Purpose:
    !          To calculate the wind speed observed along the horizontal line of
    !          sight (HLOS) from the model components of the wind, uu and vv.
    !
    ! Note:
    !          The azimuth is the direction of the HLOS observation, from target
    !          to satellite, and is not that of the trajectory of the satellite.
    !
    implicit none
    real(8), intent(in) :: uu, vv
    real(8), intent(in) :: azimuth ! in degrees

    ! HLOS target-to-satellite direction, CW from true north, in radians
    real(8) :: hlos

    hlos = azimuth * MPC_RADIANS_PER_DEGREE_R8
    ala_aladin = vv*cos(hlos) + uu*sin(hlos)

  end function ala_aladin


  real(8) function ala_aladin_tl(del_uu, del_vv, azimuth)
    !
    ! Purpose: TLM VERSION
    !          To calculate the delta of the wind speed observed along the
    !          horizontal line of sight (HLOS) from the deltas of the model
    !          components of the wind, d-uu and d-vv.
    !
    ! Note:
    !          The azimuth is the direction of the HLOS observation, from target
    !          to satellite, and is not that of the trajectory of the satellite.
    !
    implicit none
    real(8), intent(in) :: del_uu, del_vv, azimuth

    real(8) :: hlos ! HLOS direction CW from true north

    ! The satellite supplies the azimuth in degrees
    hlos = azimuth * MPC_RADIANS_PER_DEGREE_R8
    ala_aladin_tl = del_vv*cos(hlos) + del_uu*sin(hlos)

  end function ala_aladin_tl


  subroutine ala_aladin_ad(del_uu, del_vv, del_aladin, azimuth)
    !
    ! Purpose: ADJOINT VERSION
    !          From the delta of the wind speed observed along the horizontal
    !          line of sight (HLOS) calculate the deltas of the HLOS wind speed
    !          as well as of the model components of the wind, uu and vv.
    !
    ! Note:
    !          The azimuth is the direction of the HLOS observation, from target
    !          to satellite, and is not that of the trajectory of the satellite.
    !
    implicit none
    real(8), intent(inout) :: del_uu, del_vv, del_aladin
    real(8), intent(in)    :: azimuth

    real(8) :: hlos ! HLOS direction CW from true north

    ! The satellite supplies the azimuth in degrees
    hlos = azimuth * MPC_RADIANS_PER_DEGREE_R8

    del_uu     = del_uu + del_aladin*sin(hlos)
    del_vv     = del_vv + del_aladin*cos(hlos)
    del_aladin = 0
  end subroutine ala_aladin_ad


end module aladin_mod
