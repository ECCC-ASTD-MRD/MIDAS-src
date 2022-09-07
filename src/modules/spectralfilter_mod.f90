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

module spectralFilter_mod
  ! MODULE spectralFilter_mod (prefix='spf' category='8. Low-level utilities and constants')
  !
  ! :Purpose: For computing spectral filter functions
  !
  implicit none
  save
  private

  ! public procedures
  public :: spf_FilterResponseFunction

CONTAINS

!--------------------------------------------------------------------------
! spf_FilterResponseFunction
!--------------------------------------------------------------------------
  function spf_FilterResponseFunction(totalWaveNumber, waveBandIndex, waveBandPeaks, nWaveBand ) result(ResponseFunction) 
    implicit none

    real(8) :: totalWaveNumber
    real(8) :: ResponseFunction 
    integer :: waveBandIndex
    integer :: nWaveBand
    integer :: waveBandPeaks(:)

    real(8) :: linearResponse, lowerLimit, center, upperLimit
    real(8), parameter :: pi = 2.d0*asin(1.d0)

    if (waveBandIndex == nWaveBand ) then
       ! This wave band contains the largest scales.
       !
       ! The response function is 1 total wave number <= waveBandPeaks(nWaveBand)
       ! and decreases to 0 at waveBandPeaks(nWaveBand-1)
       !
       !                    response=1 |---
       !                               |   \
       !                               |    \
       !                    response=0 |------------
       !       waveBandPeaks(nWaveBand) <-|  |-> waveBandPeaks(nWaveBand-1)
       !
       lowerlimit = real(waveBandPeaks(waveBandIndex  ),8)
       upperlimit = real(waveBandPeaks(waveBandIndex-1),8)

       if ( totalWaveNumber < lowerlimit ) then
          ResponseFunction = 1.d0
       else if ( totalWaveNumber <= upperlimit ) then
          linearResponse = (upperlimit-totalWaveNumber) / (upperlimit-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    else if ( waveBandIndex /= 1 ) then
       ! This wave band contains intermediate scales (i.e., not the largest or the smallest).
       !
       ! The response function is 1 (only) for the total wave number = waveBandPeaks(waveBandIndex)
       ! and decreases to 0 at both waveBandPeaks(waveBandIndex+1) and waveBandPeaks(waveBandIndex-1)
       !
       !                    response=1 |      -
       !                               |     / \
       !                               |    /   \
       !                    response=0 |------------
       !  waveBandPeaks(waveBandIndex+1) <-|     |-> waveBandPeaks(waveBandIndex-1)
       !                                      |-> waveBandPeaks(waveBandIndex)
       !
       center     = real(waveBandPeaks(waveBandIndex  ),8)
       upperlimit = real(waveBandPeaks(waveBandIndex-1),8)
       lowerlimit = real(waveBandPeaks(waveBandIndex+1),8)

       if (      totalWaveNumber >  lowerlimit .and. &
                 totalWaveNumber <= center     ) then
          linearResponse = (totalWaveNumber-lowerlimit) / (center-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else if ( totalWaveNumber >  center      .and. &
                 totalWaveNumber <  upperlimit  ) then
          linearResponse = (upperlimit-totalWaveNumber) / (upperlimit-center)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    else
       !
       ! This wave band contains the smallest scales.
       !
       ! The response function is 1 total wave number >= waveBandPeaks(nWaveBand)
       ! and decreases to 0 at waveBandPeaks(nWaveBand-1)
       !
       !                    response=1 |    ---
       !                               |   /
       !                               |  / 
       !                    response=0 |------------
       !waveBandPeaks(waveBandIndex+1) <-|  |-> waveBandPeaks(1)
       !
       upperlimit = real(waveBandPeaks(waveBandIndex  ),8)
       lowerlimit = real(waveBandPeaks(waveBandIndex+1),8)

       if      ( totalWaveNumber > upperlimit ) then
          ResponseFunction = 1.d0
       else if ( totalWaveNumber > lowerlimit ) then
          linearResponse = (totalWaveNumber-lowerlimit) / (upperlimit-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    end if

  end function spf_FilterResponseFunction

END MODULE SpectralFilter_mod
