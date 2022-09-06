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

module tovs_extrap_mod
  ! MODULE tovs_extrap (prefix='' category='5. Observation operators')
  !
  ! :Purpose: Old code for extrapolation of temperature, when model top was
  !           still at 10hPa, to extend up to 0.1hPa for calls to RTTOV.
  !
      implicit none
      save
      private

      ! public procedures
      public :: extrap, lextrap, aextrap

!     Contents previously in comtovxt
!     -------------------------------
! === Extrapolation of temperature up to 0.1 mbs ===
!     Purpose: to define the vertical extrapolation parameters
!              required by the TOVS radiative transfer model. These
!              are used to extrapolate a temperature profile from
!              approximately to 20 mbs to 0.1 mb, using 18 levels
!              from approximately 400mbs to 20 mbs.
!
!     Revision 001 J. Halle  Fev 2000
!              extrapolation to 7 upper levels only.
!
!
!     JPXTLVIN                   : number of temperature levels used
!                                  for extrapolation
!     JPXTLVOU                   : number of temperature levels to be
!                                  extrapolated
!     JPXTOUMX                   : maximum number of temperature levels to be
!                                  extrapolated
!     MLVXTIN (JPXTLVIN)         : levels used for extrapolation
!     COEFF   (JPXTLVIN,JPXTOUMX): vertical extrapolation coefficients
!     TREFIN  (JPXTLVIN)         : reference temperature (input  levels)
!     TREFOU  (JPXTOUMX)         : reference temperature (output levels)

      INTEGER, PARAMETER :: JPXTLVIN= 18
      INTEGER, PARAMETER :: JPXTLVOU=  7
      INTEGER, PARAMETER :: JPXTOUMX=  9

      INTEGER            :: MLVXTIN(JPXTLVIN)

      REAL*8             :: COEFF (JPXTLVIN,JPXTOUMX)
      REAL*8             :: TREFIN(JPXTLVIN)
      REAL*8             :: TREFOU(JPXTOUMX)

      DATA MLVXTIN /                                                 &
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,  &
        25, 26, 27 /

      DATA COEFF /                                                       &
         -0.5439D0,  0.0986D0,  0.2883D0,  0.2490D0,  0.0137D0, -0.3474D0, -0.2183D0,  &
         -0.0916D0,  0.0580D0,  0.2026D0,  0.2813D0,  0.3344D0, -0.0961D0, -0.4738D0,  &
          0.2504D0,  0.3186D0,  0.0006D0, -0.3001D0,                             &
          0.0880D0,  0.0314D0, -0.0469D0, -0.1341D0, -0.0834D0,  0.0611D0, -0.0304D0,  &
         -0.0843D0,  0.0998D0,  0.2731D0,  0.2598D0,  0.1852D0, -0.0118D0, -0.1806D0,  &
          0.0151D0,  0.1037D0,  0.1282D0,  0.1532D0,                             &
          0.6220D0, -0.0697D0, -0.3213D0, -0.3509D0, -0.1239D0,  0.2714D0,  0.0909D0,  &
         -0.0374D0,  0.1594D0,  0.3408D0,  0.2351D0,  0.0390D0,  0.0536D0,  0.0784D0,  &
         -0.1803D0, -0.0724D0,  0.2489D0,  0.5551D0,                             &
          1.2860D0, -0.3116D0, -0.6370D0, -0.3330D0, -0.0400D0,  0.2315D0,  0.1073D0,  &
          0.0161D0,  0.1668D0,  0.3014D0,  0.1870D0, -0.0052D0,  0.1122D0,  0.2133D0,  &
         -0.3026D0, -0.1904D0,  0.2886D0,  0.7436D0,                             &
          1.7806D0, -0.4868D0, -0.8326D0, -0.2319D0,  0.0540D0,  0.1063D0,  0.0949D0,  &
          0.0808D0,  0.1200D0,  0.1491D0,  0.1030D0,  0.0336D0,  0.1411D0,  0.2033D0,  &
         -0.3388D0, -0.2399D0,  0.2338D0,  0.6844D0,                             &
          1.8949D0, -0.4038D0, -0.7716D0, -0.1868D0,  0.0578D0,  0.0515D0,  0.0826D0,  &
          0.0956D0,  0.0677D0,  0.0351D0,  0.0108D0, -0.0102D0,  0.1004D0,  0.1710D0,  &
         -0.2467D0, -0.1952D0,  0.1315D0,  0.4446D0,                             &
          1.7767D0, -0.2038D0, -0.6043D0, -0.2209D0,  0.0272D0,  0.1665D0,  0.0840D0,  &
          0.0053D0, -0.0204D0, -0.0484D0, -0.0228D0,  0.0193D0,  0.0729D0,  0.0912D0,  &
         -0.1542D0, -0.1243D0,  0.0682D0,  0.2558D0,                             &
          1.5776D0, -0.1120D0, -0.4432D0, -0.1009D0,  0.0476D0,  0.0507D0,  0.0471D0,  &
          0.0312D0, -0.0225D0, -0.0755D0, -0.0393D0,  0.0248D0,  0.0461D0,  0.0397D0,  &
         -0.0892D0, -0.0694D0,  0.0387D0,  0.1460D0,                             &
          1.2832D0, -0.0175D0, -0.2829D0, -0.0348D0,  0.0495D0,  0.0123D0,  0.0183D0,  &
          0.0155D0, -0.0098D0, -0.0354D0, -0.0177D0,  0.0135D0,  0.0220D0,  0.0174D0,  &
         -0.0423D0, -0.0334D0,  0.0168D0,  0.0685D0 /

      DATA TREFIN  /                                                     &
          223.73D0,  221.25D0,  219.14D0,  217.28D0,  215.27D0,  213.20D0,  212.16D0,  &
          211.52D0,  213.15D0,  214.65D0,  216.64D0,  218.70D0,  221.39D0,  224.25D0,  &
          228.07D0,  232.72D0,  237.76D0,  242.52D0 /

      DATA TREFOU  /                                                     &
          242.25D0,  256.61D0,  264.15D0,  262.14D0,  254.40D0,  244.88D0,  237.22D0,  &
          231.38D0,  227.33D0 /

  CONTAINS


    SUBROUTINE EXTRAP ( PROFIN, PROFOUT, JPMOLEV, JPLEV, KNPF )
      !
      ! :Purpose: Extrapolate temperature profile above 20mb
      !           on RTTOV levels (up to 0.1 mbs).
      !
      ! :Arguments:
      !     i   PROFIN (JPMOLEV,KNPF) : Temperature profile (to be extrapolated)
      !     o   PROFOUT(JPLEV,KNPF)   : Temperature profile (      extrapolated)
      !     i   JPMOLEV               : number of levels (RT model) from NWP
      !     i   JPLEV                 : number of pressure levels
      !     i   KNPF                  : Number of profiles
      !
      IMPLICIT NONE
      INTEGER JI, JJ, JK, KNPF, ILEV, JPMOLEV, JPLEV, NLVLS_XTRAP, JPMOTOP
      REAL*8 PROFIN(JPMOLEV,KNPF), PROFOUT(JPLEV,KNPF)


!*    1.  Initialize output temperature profile
!     .   -------------------------------------

      JPMOTOP = JPLEV - JPMOLEV + 1
      NLVLS_XTRAP = JPLEV - JPMOLEV
      DO JI= NLVLS_XTRAP+1, JPLEV
         DO JK = 1, KNPF
            PROFOUT(JI,JK) = PROFIN(JI-JPMOTOP+1,JK)
         ENDDO
      ENDDO

!*    2.  Extrapolation of temperatures
!     .   -----------------------------

      DO JJ= 1,NLVLS_XTRAP
         DO JK = 1, KNPF
            PROFOUT(JJ,JK) = TREFOU(JJ)
         ENDDO
      ENDDO

      DO JJ=1,NLVLS_XTRAP
         DO JI=1,JPXTLVIN
            ILEV = MLVXTIN(JI)
            DO JK = 1, KNPF
               PROFOUT(JJ,JK) = PROFOUT(JJ,JK) +  &
               COEFF(JI,JJ)*(PROFIN(ILEV-JPMOTOP+1,JK)-TREFIN(JI))
            ENDDO
         ENDDO
      ENDDO

    END SUBROUTINE EXTRAP


    SUBROUTINE LEXTRAP ( PROFIN, PROFOUT, JPMOLEV, JPLEV, KNPF )
      !
      ! :Purpose: Tangent linear of extrapolation of temperature profile
      !           above 20mb on RTTOV levels (up to 0.1 mbs).
      !
      ! :Arguments:
      !     i   PROFIN (JPMOLEV,KNPF) : Temperature profile (to be extrapolated)
      !     o   PROFOUT(JPLEV,KNPF)   : Temperature profile (      extrapolated)
      !     i   JPMOLEV               : number of levels (RT model) from NWP
      !     i   JPLEV                 : number of pressure levels
      !     i   KNPF                  : Number of profiles
      !
      IMPLICIT NONE
      INTEGER JI, JJ, JK, KNPF, ILEV, JPMOLEV, JPLEV, NLVLS_XTRAP, JPMOTOP
      REAL*8 PROFIN(JPMOLEV,KNPF), PROFOUT(JPLEV,KNPF)

!*    1.  Initialize output temperature profile
!     .   -------------------------------------

      JPMOTOP = JPLEV - JPMOLEV + 1
      NLVLS_XTRAP = JPLEV - JPMOLEV
      DO JI= NLVLS_XTRAP+1, JPLEV
         DO JK = 1, KNPF
            PROFOUT(JI,JK) = PROFIN(JI-JPMOTOP+1,JK)
         ENDDO
      ENDDO

!*    2.  Extrapolation of temperatures
!     .   -----------------------------

      DO JJ= 1,NLVLS_XTRAP
         DO JK = 1, KNPF
            PROFOUT(JJ,JK) = 0.0D0
         ENDDO
      ENDDO

      DO JJ=1,NLVLS_XTRAP
         DO JI=1,JPXTLVIN
            ILEV = MLVXTIN(JI)
            DO JK = 1, KNPF
               PROFOUT(JJ,JK) = PROFOUT(JJ,JK) +  &
               COEFF(JI,JJ)*PROFIN(ILEV-JPMOTOP+1,JK)
            ENDDO
         ENDDO
      ENDDO

    END SUBROUTINE LEXTRAP


    SUBROUTINE AEXTRAP( PROFIN, PROFOUT, JPMOLEV, JPLEV, KNPF )
      !
      ! :Purpose: Adjoint of extrapolation of temperature profile above 20mb
      !           on RTTOV levels (up to 0.1 mbs).
      !
      ! :Arguments:
      !     o   PROFIN (JPMOLEV,KNPF) : output adjoint of temperature profile
      !     i   PROFOUT(JPLEV,KNPF)   : input adjoint of temperature profile
      !     i   JPMOLEV               : number of levels (RT model) from NWP
      !     i   JPLEV                 : number of pressure levels
      !     i   KNPF                  : Number of profiles
      !
      IMPLICIT NONE
      INTEGER JI, JJ, JK, KNPF, ILEV, JPMOLEV, JPLEV, NLVLS_XTRAP, JPMOTOP
      REAL*8 PROFIN(JPMOLEV,KNPF), PROFOUT(JPLEV,KNPF)

!*    1.  Initialize output adjoint of temperature profile
!     .   ------------------------------------------------

      JPMOTOP = JPLEV - JPMOLEV + 1
      NLVLS_XTRAP = JPLEV - JPMOLEV
      DO JI= 1, JPMOLEV
         DO JK = 1, KNPF
            PROFIN(JI,JK) = 0.0D0
         ENDDO
      ENDDO

!*    2.  Adjoint of extrapolation of temperatures
!     .   ----------------------------------------

      DO JJ = NLVLS_XTRAP, 1, -1
        DO JI = JPXTLVIN, 1, -1
          ILEV = MLVXTIN(JI)
          DO JK = KNPF, 1, -1
            PROFIN(ILEV-JPMOTOP+1,JK) = PROFIN(ILEV-JPMOTOP+1,JK) +  &
                                        PROFOUT(JJ,JK)*COEFF(JI,JJ)
          ENDDO
        ENDDO
      ENDDO

      DO JJ = NLVLS_XTRAP, 1, -1
        DO JK = KNPF, 1, -1
          PROFOUT(JJ,JK) = 0.0D0
        ENDDO
      ENDDO

      DO JI = JPLEV, NLVLS_XTRAP+1, -1
        DO JK = KNPF, 1, -1
          PROFIN(JI-JPMOTOP+1,JK) = PROFIN(JI-JPMOTOP+1,JK) +  &
                                    PROFOUT(JI,JK)
          PROFOUT(JI,JK) = 0.0D0
        ENDDO
      ENDDO

    END SUBROUTINE AEXTRAP

end module tovs_extrap_mod
