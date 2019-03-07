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
!! MODULE randomNumber_mod (prefix="rng" category='7. Low-level data objects and utilities')
!!
!! *Purpose*: A Gaussian random number generator (RNG) module
!!
!--------------------------------------------------------------------------
module randomnumber_mod
  implicit none
  save
  private
 
  ! public procedures
  public :: rng_setup, rng_gaussian

  real(8) :: gset

  integer :: idum_local
  integer :: iset
  
  logical :: initialized = .false.

contains

  !--------------------------------------------------------------------------
  ! rng_Setup
  !--------------------------------------------------------------------------
  subroutine rng_setup(seed)
    implicit none

    integer, intent(in) :: seed

    REAL(8) :: RANDOMNUMBER,V1,V2,FAC,RR
    INTEGER :: J

    if (initialized) then
       write(*,*) 'rng_setup: WARNING: you are re-initializing the module!!!'
    end if

    idum_local = -seed
    iset = 0

1   V1=2.D0*RANDOM()-1.D0
    V2=2.D0*RANDOM()-1.D0
    RR=V1**2+V2**2
    IF(RR>=1.D0) GO TO 1
    FAC=SQRT(-2.D0*LOG(RR)/RR)
    GSET=V1*FAC
    randomNumber=V2*FAC
    ISET=1

    initialized = .true.

    write(*,*) 'rng_setup: done using seed = ', seed
 
  end subroutine rng_setup
  
  !--------------------------------------------------------------------------
  ! rng_Gaussian
  !--------------------------------------------------------------------------
  FUNCTION rng_gaussian() result(randomNumberGaussian)
    implicit none

    ! ADAPTED from the book:
    !
    ! Book-title  NUMERICAL RECIPES in FORTRAN
    !             The Art of Scientific Computing
    !             (First Edition)
    !
    ! Authors     Press, Flannery, Teukolsky, Vetterling
    !
    ! OBJECT      Returns a normally distributed deviate
    !             with zero mean and unit variance, using
    !             RANDOM as the source of uniform
    !             deviates
    
    REAL(8) :: randomNumberGaussian,V1,V2,FAC,RR
    INTEGER :: J

    IF (ISET.EQ.0.OR.IDUM_LOCAL.eq.999) THEN
1      V1=2.D0*RANDOM()-1.D0
       V2=2.D0*RANDOM()-1.D0
       RR=V1**2+V2**2
       IF(RR.GE.1.D0)GO TO 1
       FAC=SQRT(-2.D0*LOG(RR)/RR)
       GSET=V1*FAC
       randomNumberGaussian=V2*FAC
       ISET=1
    ELSE
       randomNumberGaussian=GSET
       ISET=0
    ENDIF
  
  END FUNCTION RNG_GAUSSIAN
  
  !--------------------------------------------------------------------------
  ! random
  !--------------------------------------------------------------------------
  FUNCTION RANDOM() result(randomnumber)
    implicit none
    
    ! FUNCTION RANDOM
    !
    ! ADAPTED from the book:
    !
    ! Book-title  NUMERICAL RECIPES in FORTRAN
    !             The Art of Scientific Computing
    !             (First Edition)
    !
    ! Authors     PRESS, FLANNERY, TEUKOLSKY, VETTERLING
    !
    ! OBJECT      Returns a random deviate between 0.0 and 1.0.

    REAL(8), save :: RRAND(97)
    INTEGER, save :: IX1,IX2,IX3
    INTEGER, save :: IFF=1
    
    REAL(8) :: RM1,RM2,RANDOMNUMBER
    INTEGER :: J,M1,IA1,IC1,M2,IA2,IC2,M3,IA3,IC3
    PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247D-6)
    PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773D-6)
    PARAMETER (M3=243000,IA3=4561,IC3=51349)
    
    IF (IDUM_LOCAL.LT.0.OR.IFF.EQ.0) THEN
       IFF=1
       IX1=MOD(IC1-IDUM_LOCAL,M1)
       IX1=MOD(IA1*IX1+IC1,M1)
       IX2=MOD(IX1,M2)
       IX1=MOD(IA1*IX1+IC1,M1)
       IX3=MOD(IX1,M3)
       DO J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          RRAND(J)=(real(IX1,8)+real(IX2,8)*RM2)*RM1
       END DO
       IDUM_LOCAL=1
    ENDIF
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IA2*IX2+IC2,M2)
    IX3=MOD(IA3*IX3+IC3,M3)
    J=1+(97*IX3)/M3
    IF(J.GT.97.OR.J.LT.1) then
       write(6,*) 'Input error in RANDOM for  J = ',J
       stop
    endif
    RANDOMNUMBER=RRAND(J)
    RRAND(J)=(real(IX1,8)+real(IX2,8)*RM2)*RM1
    
  end FUNCTION RANDOM
  
end module randomnumber_mod
