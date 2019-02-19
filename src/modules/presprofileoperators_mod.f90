!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!! MODULE presProfileOperators (prefix="ppo" category='7. Low-level data objects and utilities')
!!
!! *Purpose*: Vertical interpolation subroutines, including the special routines
!!            designed to interpolate to the (widely spaced) RTTOV pressure levels.
!!
!--------------------------------------------------------------------------
module presProfileOperators_mod
  implicit none
  save
  private

  ! public procedures
  public :: ppo_intAvg, ppo_intAvgTl, ppo_intAvgTl_v2, ppo_intAvgAd, ppo_lintv

  contains

    SUBROUTINE ppo_INTAVG (PVLEV,PVI,KNI,KNPROF,KNO,PPO,PVO)
!--------------------------------------------------------------------
!
!**s/r INTAVG   - Forward interpolator based on piecewise 
!                 weighted averaging.
!
!Author  : Y.J. Rochon *ARQX/EC Nov 2005
!          Starting points: LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!
!Revisions: 
!          Saroja Polavarapu *ARMA/EC Nov 2005
!          - Completed split into three routines for consistency with
!            LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!Arguments
!     i   PVLEV(KNI,KNPROF)    : Vertical levels, pressure (source)
!    i/o  PVI(KNI,KNPROF)      : Input vector to be interpolated (source)
!     i   KNI                  : Number of input levels (source)
!     i   KNPROF               : Number of profiles
!     i   KNO                  : Number of output levels (destination)
!     i   PPO(KNO)             : Vertical levels, pressure (destination)
!    o/i  PVO(KNO,KNPROF)      : Output interpolated profiles (destination) 
!
!    -------------------
!
! Purpose: Forward interpolator based on piecewise
!          weighted averaging in log of pressure   
!          of one-dimensional vectors.
!
! Comments:
!
!     1) vlev in eta coordinate associated to pvlev in pressure.
!                  Y = sum(H*PVI)
!                    = PVO
!     
!               with ZPZ=H on output of LAYERAVG
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,intent(in)  :: KNI, KNO,KNPROF
      REAL(8) ,intent(in)  :: PVLEV(KNI,KNPROF)
      REAL(8) ,intent(in)  :: PVI(KNI,KNPROF)
      REAL(8) ,intent(in)  :: PPO(KNO)
      REAL(8) ,intent(out) ::  PVO(KNO,KNPROF)
!*********************
      INTEGER  JO, JN
      REAL(8)  ZLNPI (KNI),ZPS(KNI)
      REAL(8)  ZLNPO (KNO)

! --- Apply weighted averaging 

      ZLNPO(:)=LOG(PPO(:))

      DO JN = 1, KNPROF
        ZLNPI(:) = LOG( PVLEV(:,JN) )
        PVO(:,JN) = 0.0D0

        DO JO=1,KNO
          CALL LAYERAVG(ZLNPO,ZLNPI,PVI(:,JN),KNO,KNI,  &
               JO,PVO(JO,JN) )
        END DO
      END DO

    END SUBROUTINE PPO_INTAVG
    

    SUBROUTINE ppo_INTAVGTL(PVLEV,dP_dPsfc,PVI,PVIG,PPS,KNI,KNPROF,KNO,PPO,PVO)
!--------------------------------------------------------------------
!
!!!s/r INTAVGTL - Application of tangent Linear of piecewise weighted averaging.
!
!
!Author  : Y.J. Rochon *ARQX/EC Nov 2005
!          Starting points: LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!
!Revisions:
!          Saroja Polavarapu *ARMA/EC Nov 2005
!          - Completed split into three routines for consistency with
!!           LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!          S. Pellerin, ARMA, August 2008
!          - Optimisation, call to LAYERAVG2
!
!Arguments
!     i   PVLEV(KNI,KNPROF)    : Vertical levels, pressure (source)
!     i   dP_dPsfc(KNI,KNPROF) : derivative of pressure wrt sfc pressure for vertical coordinate (source)
!     i   PVI(KNI,KNPROF)      : Increments of vector on input levels
!     i   PVIG(KNI,KNPROF)     : Vector on input levels (source)
!     i   PPS(KNPROF)             : Input surface pressure increments
!     i   KNI                  : Dimension of input levels (source)
!     i   KNI                     : Number of input levels (source)
!     i   KNPROF                  : Number of profiles
!     i   KNO                     : Number of output levels (destination)
!     i   PPO(KNO)                : Vertical levels, pressure (destination)
!     o   PVO(KNO,KNPROF)         : Increments of vector on output levels
!
!!   -------------------
!
! Purpose: Application of tangent linear
!          of piecewise weighted averaging
!          in log of pressure of one-dimensional vectors.
!
! Comments:
!
!     1) vlev in eta coordinate associated to pvlev in pressure.
!
!
!                  dY = sum(H*PVI) + PPS*sum(PVIG*dH/dPs))
!                     =     ZPVO    +       ZPVOPS
!                     = PVO
!
!               where 
!
!                     dPZ(i)/dPs = sum_k (dH(i)/dPVLEV(k) * dPVLEV(k)/dPs)
!                                = sum_k (dH(i)/dZLNPI(k) * dP_dPsfc(k)/PVLEV(k))
!
!               with ZPZ(k)=zpresb(k)/PVLEV(k) on input to LAYERAVG and
!                    ZPZ(i)=H(i) on output of LAYERAVG
! 
!!------------------------------------------------------------------- 
!
      IMPLICIT NONE
      INTEGER,intent(in)  :: KNI, KNO,KNPROF
      REAL(8),intent(in)   :: PVLEV(KNI,KNPROF)
      REAL(8),intent(in)   :: dP_dPsfc(KNI,KNPROF)
      REAL(8),intent(in)   :: PPO(KNO)
      REAL(8),intent(in)   :: PVI(KNI,KNPROF)
      REAL(8),intent(in)   :: PVIG(KNI,KNPROF)
      REAL(8),intent(in)   :: PPS(KNPROF)
      REAL(8),intent(out)  :: PVO(KNO,KNPROF)
!********************************************
      REAL(8)  ZLNPI (KNI,knprof),ZPZ(KNI), ZPS(KNI,knprof)
      REAL(8)  ZLNPO (KNO),ZPZPS,ZPVOPS(kno,knprof),PSS,ZPVO(kno,knprof)

! --- Apply weighted averaging 

      ZLNPO(:) = LOG(PPO(:))
      PSS=0.0D0

  
      ZPS(:,:) = dP_dPsfc(:,:) / PVLEV(:,:)
      ZLNPI(:,:) = LOG( PVLEV(:,:) )
      

      CALL LAYERAVG_TL(ZLNPO,ZLNPI,pvig,PVI,KNO,  &
                     KNI,PPS,ZPVO,ZPZ,ZPS,ZPVOPS,knprof,pvo)

    END SUBROUTINE PPO_INTAVGTL


    SUBROUTINE ppo_INTAVGTL_v2(PVLEV,dP_dPsfc,PVI,PVIG,PP,KNI,KNPROF,KNO,PPO,PVO)
      IMPLICIT NONE
      INTEGER,intent(in)  :: KNI, KNO,KNPROF
      REAL(8),intent(in)   :: PVLEV(KNI,KNPROF)
      REAL(8),intent(in)   :: dP_dPsfc(KNI,KNPROF)
      REAL(8),intent(in)   :: PPO(KNO)
      REAL(8),intent(in)   :: PVI(KNI,KNPROF)
      REAL(8),intent(in)   :: PVIG(KNI,KNPROF)
      REAL(8),intent(in)   :: PP(KNI,KNPROF)
      REAL(8),intent(out)  :: PVO(KNO,KNPROF)
      REAL(8)  ZLNPI (KNI,knprof),ZPZ(KNI), ZPS(KNI,knprof)
      REAL(8)  ZLNPO (KNO),ZPZPS,ZPVOPS(kno,knprof),PSS,ZPVO(kno,knprof)


      ZLNPO(:) = LOG(PPO(:))
      PSS=0.0D0

  
      ZPS(:,:) = dP_dPsfc(:,:) / PVLEV(:,:)
      ZLNPI(:,:) = LOG( PVLEV(:,:) )
      

      CALL LAYERAVG_TL_v2(ZLNPO,ZLNPI,pvig,PVI,KNO,  &
                     KNI,PP,ZPVO,ZPZ,ZPS,ZPVOPS,knprof,pvo)

    END SUBROUTINE PPO_INTAVGTL_v2


    SUBROUTINE ppo_INTAVGAD(PVLEV,dP_dPsfc,PVI,PVIG,PPS,KNI,KNPROF,KNO,PPO,PVO)
!--------------------------------------------------------------------
!
!!!s/r INTAVGAD - Adjoint of piecewise weighted averaging.
!
!
!Author  : Y.J. Rochon *ARQX/EC Nov 2005
!          Starting points: LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!
!Revisions:
!          Saroja Polavarapu *ARMA/EC Nov 2005
!          - Completed split into three routines for consistency with
!!           LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!          S. Pellerin, ARMA, August 2008
!          - Optimisation, call to LAYERAVG2
!
!Arguments
!     i   PVLEV(KNI,KNPROF)    : Vertical levels, pressure (source)
!     o   PVI(KNI,KNPROF)      : Adjoint of increments on input levels (output)
!     i   PVIG(KNI,KNPROF)     : Vector to be interpolated (source)
!     o   PPS(KNPROF)             : Output adjoint of sfc P incr
!     i   KNI                  : Dimension of input levels (source)
!     i   KNI                     : Number of input levels (source)
!     i   KNPROF                  : Number of profiles
!     i   KNO                     : Number of output levels
!     i   PPO(KNO)                : Vertical levels, pressure
!     i   PVO(KNO,KNPROF)         : Input adjoint of interpolated profiles
!
!!   -------------------
!
!!Purpose: Performs calculation for adjoint of
!          piecewise weighted averaging
!          in log of pressure of one-dimensional vectors.
!
! Comments:
!
!     1) vlev in eta coordinate associated to pvlev in pressure.
!
!                  PVI(i,jn) = PVI(i,jn) + sum_jo (PVO(jo,jn)*H(i))
!                            = PVI(i,jn) + sum_jo (ZPZ(jo,i))
!
!                  PPS(jn) = PPS(jn) + sum_jo (PVO(jo,jn)*sum_i (PVIG(i,jn)*dH(i)/dPs))
!                          = PPS(jn) + sum_jo (ZPZPS(jo))
!
!!-------------------------------------------------------------------  
      IMPLICIT NONE
     
      INTEGER,intent(in) :: KNI, KNO, KNPROF
      REAL(8),intent(in) :: PVLEV(KNI,KNPROF)
      REAL(8),intent(in) :: dP_dPsfc(KNI,KNPROF)
      REAL(8),intent(in) :: PPO(KNO), PVO(KNO,KNPROF)
      REAL(8),intent(in) :: PVIG(KNI,KNPROF)
      REAL(8),intent(out):: PVI(KNI,KNPROF)
      REAL(8),intent(out):: PPS(KNPROF)
!****************************************************
      REAL(8) ZLNPI (KNI,knprof),ZPZ(KNI)
      REAL(8) ZPS(KNI,knprof)
      REAL(8) ZLNPO (KNO)


! --- Apply weighted averaging 

      ZLNPO(:)=LOG(PPO(:))

      ZPS(:,:) = dP_dPsfc(:,:) / PVLEV(:,:)
      ZLNPI(:,:) = LOG( PVLEV(:,:) )

      CALL LAYERAVG_AD(ZLNPO,ZLNPI,pvig,pvi,  &
                     KNO,KNI,PPS,PVO,ZPZ,ZPS,knprof)

    END SUBROUTINE PPO_INTAVGAD


    SUBROUTINE LAYERAVG(PX1,PX2,PY2,KN1,KN2,KI,PMEAN)

      IMPLICIT NONE

      INTEGER,intent(in) :: KN1,KN2,KI
      REAL(8),intent(in)  :: PX1(KN1),PX2(KN2),PY2(KN2)
      REAL(8),intent(out) :: PMEAN
!---------------------------------------------------------
!
!!!s/r LAYERAVG - Perform integration between points PX1(KI-1) and
!                 PX1(KI+1) with  piecewise linear weighting having 
!                 weights of zero at  ki-1 and ki+1 and max weight at ki.
!
!                  Output: Weighted mean value, its contribution to
!                          the TLM*increment or to the adjoint.
!
!Author  : Y.J. Rochon *ARQX/EC Nov. 2005
!
!!Revisions: 
!
!    -------------------
!
!Purpose: Perform integration between points PX1(KI-1) and
!         PX1(KI+1) with  piecewise linear weighting having 
!         weights of zero at  ki-1 and ki+1 and max weight at ki.
!  
!         Output: Weighted mean value, its contribution to
!                 the TLM*increment or to the adjoint.
!
!Arguments:
!
!   INPUT
!
!
!     PX1:      Reference levels (e.g. lnP; in increasing values)
!     PX2:      Available levels (e.g. lnP; in increasing values)
!     PY2:      Values at levels (e.g. temperature) 
!     KN1:      Dimension of PX1.
!     KN2:      Dimension of other arrays.
!     KI:       Identifies region of interest: PX1(KI-1) to PX1(KI+1)
!
!   OUTPUT:
!
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) using
!               background values
!-----------------------------------------------------------
      REAL(8) :: PZ(KN2)
      REAL(8) Z1,Z2,Z3,ZW1,ZW2
      REAL(8) zsum
      INTEGER J
      logical :: SKIP, test

! --- Identify boundary points

      z2=px1(ki)

      if (ki == 1) then 
         z1=2.0D0*z2-px1(ki+1)
      else
         z1=px1(ki-1)
      end if   

      if (ki == kn1) then
         z3=2.0D0*z2-z1
      else   
         z3=px1(ki+1)
      end if
      if (z3 > px2(kn2)) z3=px2(kn2)

      skip= .false.
      if (z2 >= px2(kn2)) then
         z3=px2(kn2)
         z2=px2(kn2)
         skip = .true.
      end if

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

      pz(1:kn2)=0.0D0
      test = .false.
      loop1:do j=1,kn2-1
         if (px2(j) >= z3) exit loop1

         if (px2(j) <= z2.and.px2(j+1) > z1) then 

            call sublayer(z1,z2,z3,px2(j),px2(j+1),  &
                          py2(j),py2(j+1),zw1,zw2 )
            pz(j)=pz(j)+zw1
            pz(j+1)=pz(j+1)+zw2
            test = .true.
         end if

         if (px2(j) < z3.and.px2(j+1) >= z2.and. .not. skip) then

            call sublayer(z3,z2,z1,px2(j),px2(j+1),  &
                          py2(j),py2(j+1),zw1,zw2 )
            pz(j)=pz(j)+zw1
            pz(j+1)=pz(j+1)+zw2
            test = .true.
         end if
       end do loop1

       if ( .not. test) pz(j)=1.0D0

! --- Apply forward interpolator
 
       zsum=0.0D0
       pmean=0.0D0
       do j=1,kn2
         pmean=pmean+pz(j)*py2(j)
         zsum=zsum+pz(j)
       end do
       pmean=pmean/zsum
     

     END SUBROUTINE LAYERAVG


    SUBROUTINE LAYERAVG_TL(PX1,PX2,PY2,PY2INCR,KNO,KNI,PPS,PMEAN,PZ,PZS,PZPS,knprof,pvo)

      IMPLICIT NONE
      INTEGER,intent(in) :: KNO,KNI,knprof
      REAL(8),intent(out) :: pvo(kno,knprof)
      REAL(8),intent(out) :: PMEAN(kno,knprof)
      REAL(8),intent(in)  :: PX1(KNO),PX2(KNI,knprof),PY2(KNI,knprof),PY2INCR(KNI,knprof)
      REAL(8),intent(in)  :: PPS(knprof)
      REAL(8),intent(in)  :: PZS(KNI,knprof)
      REAL(8),intent(out) :: PZ(KNI),PZPS(kno,knprof)      
!---------------------------------------------------------
!
! s/r LAYERAVG - Perform integration between points PX1(KI-1) and
!                 PX1(KI+1) with  piecewise linear weighting having 
!                 weights of zero at  ki-1 and ki+1 and max weight at ki
!
!                  Output: Weighted mean value, its contribution to
!                          the TLM*increment or to the adjoint.
!
!Author  : Y.J. Rochon *ARQX/EC Nov. 2005
!
!Revisions: 
!          S. Pellerin, ARMA, August 2008
!              - Introduction of version 2
!              - Introduction of profile and kno loops to reduce the
!                number of calls (optimisation)
!
!    -------------------
!
!Purpose: Perform integration between points PX1(KI-1) and
!         PX1(KI+1) with  piecewise linear weighting having 
!         weights of zero at  ki-1 and ki+1 and max weight at ki.
!  
!         Output: Weighted mean value, its contribution to
!                 the TLM*increment or to the adjoint.
!
!Arguments:
!
!   INPUT
!     PX1:      Reference levels (e.g. lnP; in increasing values)
!     PX2:      Available levels (e.g. lnP; in increasing values)
!     PY2:      Values at levels (e.g. temperature) 
!     PY2INCR:  Increments
!     PZ:       Extra array related to gradients w.r.t. Ps 
!     KNO:      Dimension of PX1.
!     KNI:      Dimension of other arrays.
!     PZS:      dlnP/dPs
!     KI:       Identifies region of interest: PX1(KI-1) to PX1(KI+1)
!     PPS:      Surface pressure increment
!
!   OUTPUT:
!
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) increments
!     PZ:       Resultant accumulated contribution factors TLM 
!     PZPS:     Surface pressure related TLM*increment term
!
!-----------------------------------------------------------
      REAL(8) Z1,Z2,Z3,ZW1,ZW2
      REAL(8) zsum
      INTEGER J,jo,levIndex
      logical :: skip, test
! --- Identify boundary points

      do levIndex = 1, knprof
        do jo = 1, kno
          z2=px1(jo)

          if (jo == 1) then 
            z1=2.0D0*z2-px1(jo+1)
          else
            z1=px1(jo-1)
          end if   

          if (jo == kno) then
            z3=2.0D0*z2-z1
          else   
            z3=px1(jo+1)
          end if
          if (z3 > px2(kni,levIndex)) z3=px2(kni,levIndex)

          skip = .false.
          if (z2 >= px2(kni,levIndex)) then
            z3=px2(kni,levIndex)
            z2=px2(kni,levIndex)
            skip = .true.
          end if

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

          pzps(jo,levIndex)=0.0D0
          pz(1:kni)=0.0D0
          test = .false.

          loop2:do j=1,kni-1
            if (px2(j,levIndex) >= z3) exit loop2

            if (px2(j,levIndex) <= z2.and.px2(j+1,levIndex) > z1) then 

              call sublayer(z1,z2,z3,px2(j,levIndex),px2(j+1,levIndex), &
                            py2(j,levIndex),py2(j+1,levIndex),zw1,zw2,  &
                            pzs(j,levIndex),pzs(j+1,levIndex),pzps(jo,levIndex))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if

            if (px2(j,levIndex) < z3.and.px2(j+1,levIndex) >= z2.and. .not. skip) then

              call sublayer(z3,z2,z1,px2(j,levIndex),px2(j+1,levIndex), &
                            py2(j,levIndex),py2(j+1,levIndex),zw1,zw2,  &
                            pzs(j,levIndex),pzs(j+1,levIndex),pzps(jo,levIndex))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if
          end do loop2

          if (.not. test) pz(j)=1.0D0
! --- Apply forward interpolator (kflag=0), determine TLM*increment (kflag=1)
!     or use TLM for adjoint calc (kflag=2)

          zsum=0.0D0
          pmean(jo,levIndex)=0.0D0
          do j=1,kni      
            pmean(jo,levIndex)=pmean(jo,levIndex)+pz(j)*py2incr(j,levIndex)
            zsum=zsum+pz(j)
          end do
          pmean(jo,levIndex)=pmean(jo,levIndex)/zsum
          if (test) pzps(jo,levIndex)=pzps(jo,levIndex)*pps(levIndex)/zsum
          PVO(JO,LEVINDEX)=pmean(jo,levIndex)+PZPS(jo,levIndex)
         
        end do
      end do

    END SUBROUTINE LAYERAVG_TL


    SUBROUTINE LAYERAVG_TL_v2(PX1,PX2,PY2,PY2INCR,KNO,KNI,PP,PMEAN,PZ,PZS,PZPS,knprof,pvo)
      IMPLICIT NONE
      INTEGER,intent(in) :: KNO,KNI,knprof
      REAL(8),intent(out) :: pvo(kno,knprof)
      REAL(8),intent(out) :: PMEAN(kno,knprof)
      REAL(8),intent(in)  :: PX1(KNO),PX2(KNI,knprof),PY2(KNI,knprof),PY2INCR(KNI,knprof)
      REAL(8),intent(in)  :: PP(kni,knprof)
      REAL(8),intent(in)  :: PZS(KNI,knprof)
      REAL(8),intent(out) :: PZ(KNI),PZPS(kno,knprof)      
      REAL(8) ::  PZP(KNI)
      REAL(8) Z1,Z2,Z3,ZW1,ZW2,zp1,zp2
      REAL(8) zsum, zsum2
      INTEGER J,jo,levIndex
      logical :: skip, test
! --- Identify boundary points

      do levIndex = 1, knprof
        do jo = 1, kno
          z2=px1(jo)

          if (jo == 1) then 
            z1=2.0D0*z2-px1(jo+1)
          else
            z1=px1(jo-1)
          end if   

          if (jo == kno) then
            z3=2.0D0*z2-z1
          else   
            z3=px1(jo+1)
          end if
          if (z3 > px2(kni,levIndex)) z3=px2(kni,levIndex)

          skip = .false.
          if (z2 >= px2(kni,levIndex)) then
            z3=px2(kni,levIndex)
            z2=px2(kni,levIndex)
            skip = .true.
          end if

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

          pzps(jo,levIndex)=0.0D0
          pz(1:kni)=0.0D0
          pzp(1:kni)=0.0D0
          test = .false.

          loop2:do j=1,kni-1
            if (px2(j,levIndex) >= z3) exit loop2

            if (px2(j,levIndex) <= z2.and.px2(j+1,levIndex) > z1) then 

              call sublayer_v2(z1,z2,z3,px2(j,levIndex),px2(j+1,levIndex), &
                            py2(j,levIndex),py2(j+1,levIndex),zw1,zw2,  &
                            pzs(j,levIndex),pzs(j+1,levIndex),zp1,zp2)
              pzp(j) = pzp(j) + zp1
              pzp(j+1) = pzp(j+1) + zp2
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if

            if (px2(j,levIndex) < z3.and.px2(j+1,levIndex) >= z2.and. .not. skip) then

              call sublayer_v2(z3,z2,z1,px2(j,levIndex),px2(j+1,levIndex), &
                            py2(j,levIndex),py2(j+1,levIndex),zw1,zw2,  &
                            pzs(j,levIndex),pzs(j+1,levIndex),zp1,zp2)
              pzp(j) = pzp(j) + zp1
              pzp(j+1) = pzp(j+1) + zp2
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if
          end do loop2

          if (.not. test) then
            pz(j)=1.0D0
            !pzp(j)=1.0D0
          end if
! --- Apply forward interpolator (kflag=0), determine TLM*increment (kflag=1)
!     or use TLM for adjoint calc (kflag=2)

          zsum=0.0D0
          zsum2=0.0D0
          pmean(jo,levIndex)=0.0D0
          do j=1,kni      
            pmean(jo,levIndex)=pmean(jo,levIndex)+pz(j)*py2incr(j,levIndex)
            pzps(jo,levIndex)=pzps(jo,levIndex)+pzp(j)*pp(j,levIndex)
            zsum=zsum+pz(j)
            zsum2=zsum2+pzp(j)
          end do
          pmean(jo,levIndex)=pmean(jo,levIndex)/zsum

          !if (test) pzps(jo,levIndex)=pzps(jo,levIndex)/zsum2
          if (test) pzps(jo,levIndex)=pzps(jo,levIndex)/zsum
          PVO(JO,LEVINDEX)=pmean(jo,levIndex)+pzps(jo,levIndex)
         
        end do
      end do

    END SUBROUTINE LAYERAVG_TL_v2


    SUBROUTINE LAYERAVG_AD(PX1,PX2,PY2,PY2INCR,KNO,KNI,PPS,PMEAN,PZ,PZS,knprof)

      IMPLICIT NONE
      INTEGER,intent(in) ::  KNO,KNI,knprof
      REAL(8),intent(in) :: PMEAN(kno,knprof)
      REAL(8),intent(in) ::  PX1(KNO),PX2(KNI,knprof),PY2(KNI,knprof)
      REAL(8),intent(inout) :: PY2INCR(KNI,knprof)
      REAL(8),intent(in) ::  PZS(KNI,knprof)
      REAL(8),intent(inout) ::  PPS(knprof)
      REAL(8),intent(out) ::  PZ(KNI)
!---------------------------------------------------------
!
! s/r LAYERAVG - Perform integration between points PX1(KI-1) and
!                 PX1(KI+1) with  piecewise linear weighting having 
!                 weights of zero at  ki-1 and ki+1 and max weight at ki
!
!                  Output: Weighted mean value, its contribution to
!                          the TLM*increment or to the adjoint.
!
!Author  : Y.J. Rochon *ARQX/EC Nov. 2005
!
!Revisions: 
!          S. Pellerin, ARMA, August 2008
!!             - Introduction of version 2
!!             - Introduction of profile and kno loops to reduce the
!!               number of calls (optimisation)
!
!    -------------------
!
!Purpose: Perform integration between points PX1(KI-1) and
!         PX1(KI+1) with  piecewise linear weighting having 
!         weights of zero at  ki-1 and ki+1 and max weight at ki.
!  
!         Output: Weighted mean value, its contribution to
!                 the TLM*increment or to the adjoint.
!
!Arguments:
!
!   INPUT
!
!     PX1:      Reference levels (e.g. lnP; in increasing values)
!     PX2:      Available levels (e.g. lnP; in increasing values)
!     PY2:      Values at levels (e.g. temperature) 
!     PY2INCR:  Increments
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1)
!     KNO:      Dimension of PX1.
!     KNI:      Dimension of other arrays.
!     PZS:      dlnP/dPs
!     KI:       Identifies region of interest: PX1(KI-1) to PX1(KI+1)
!     PPS:      Surface pressure adjoint (in/out)
!
!   OUTPUT:
!
!     PZ:       Resultant accumulated contribution factors Adjoint when KFLAG=2 
!
!-----------------------------------------------------------
      REAL(8) Z1,Z2,Z3,ZW1,ZW2,PZPS(kno,knprof)
      REAL(8) zsum
      INTEGER J,jo,jn
      logical test,skip
! --- Identify boundary points

      do jn = 1, knprof
        do jo = 1, kno
          z2=px1(jo)

          if (jo == 1) then 
            z1=2.0D0*z2-px1(jo+1)
          else
            z1=px1(jo-1)
          end if   

          if (jo == kno) then
            z3=2.0D0*z2-z1
          else   
            z3=px1(jo+1)
          end if
          if (z3 > px2(kni,jn)) z3=px2(kni,jn)

          skip = .false.
          if (z2 >= px2(kni,jn)) then
            z3=px2(kni,jn)
            z2=px2(kni,jn)
            skip = .true.
          end if

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

          pzps(jo,jn)=0.0D0
          pz(1:kni)=0.0D0
          test = .false.

          do j=1,kni-1
            if (px2(j,jn) >= z3) exit

            if (px2(j,jn) <= z2.and.px2(j+1,jn) > z1) then 

              call sublayer(z1,z2,z3,px2(j,jn),px2(j+1,jn), &
                            py2(j,jn),py2(j+1,jn),zw1,zw2,  &
                            pzs(j,jn),pzs(j+1,jn),pzps(jo,jn))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if

            if (px2(j,jn) < z3.and.px2(j+1,jn) >= z2.and. .not. skip) then

              call sublayer(z3,z2,z1,px2(j,jn),px2(j+1,jn), &
                            py2(j,jn),py2(j+1,jn),zw1,zw2,  &
                            pzs(j,jn),pzs(j+1,jn),pzps(jo,jn))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              test = .true.
            end if
          end do

          if (.not. test) pz(j)=1.0D0
! --- Apply forward interpolator (kflag=0), determine TLM*increment (kflag=1)
!     or use TLM for adjoint calc (kflag=2)

          zsum=0.0D0
          do j=1,kni
            zsum=zsum+pz(j)         
          end do
          pz(1:kni)=pz(1:kni)*pmean(jo,jn)/zsum
          py2incr(1:KNI,JN)=py2incr(1:KNI,JN)+PZ(1:KNI)
          if (test) pzps(jo,jn)=pzps(jo,jn)*pmean(jo,jn)/zsum
          PPS(JN)=PPS(JN)+PZPS(jo,jn)
          
        end do
      end do

    END SUBROUTINE LAYERAVG_AD


    SUBROUTINE SUBLAYER(z1,z2,z3,x1,x2,t1,t2,w1,w2,pzs1,pzs2,pzps)
      IMPLICIT NONE
      REAL(8),intent(in) :: z1,z2,z3,x1,x2,t1,t2
      REAL(8),intent(out) :: w1,w2
      REAL(8),intent(inout) ,optional :: pzps
      REAL(8),intent(in)    ,optional :: pzs1,pzs2
!-----------------------------------------------------------------
!
!  Written by Yves J. Rochon *ARQX/EC, Nov. 2005
!
!  Revisions:
!            Saroja Polavarapu *ARMA/EC, Nov 2005
!            - Setting of y1 and y2 moved from calling routine.           
!
!  Revisions:
!
!  Purpose: Determine weight coefficients to assign to NWP variables
!           at x1 and x2. Weights are determined from integration over
!           the range (y1,y2), which is within the ranges of the
!           NWP layer (x1,x2) and of the RTM layer (z1,z2). Intergrals
!           are approximated via the trapezoidal rule:
!
!              integral of f(x) from y1 to y2 = (f(y1)+f(y2))/2*abs(y1-y2) 
!
!           This is synonomous to having an integrand linear in x.
!
!           Normalization done in calling routine.
!
!  Input:
!
!       z1.........Outer boundary of RTM level (could be above are below z2)
!       z2.........Inner boundary of RTM level 
!                      (position of RTM reference level)
!       z3.........Second outer boundary
!       x1.........Upper boundary of NWP layer (x1<x2)
!       x2.........Lower boundary of NWP layer
!       t1.........Variable value at upper NWP level.
!       t2.........Variable value at lower NWP level.
!       pzs1.......dlnP/dPs = dx1/dPs
!       pzs2.......dlnP/dPs = dx2/dPs
!       pzps.......Current gradient contribution for weights*variable w.r.t Ps
!
!  Output:
!
!       w1.........Weight assigned to variable at upper NWP level
!       w2.........Weight assigned to variable at lower NWP level
!       pzps.......Updated gradient contribution for weights*variable w.r.t Ps
!
!  Other:
!
!       tot........Evaluated integral
!       g1.........Gradient of weights*variables w.r.t. x1
!       g2.........Gradient of weights*variables w.r.t. x2
!
!-----------------------------------------------------------------
      REAL(8) y1,y2,tot,d,w10,w20,dz,dx,dy,dzd,dxd,g1,g2
      REAL(8) a1,a2,aa1,aa2
      logical bot,top
!
! --- Identify and set upper and lower boundaries of
!     integration/averaging layers.
!
!     y1.........Upper boundary of integral range (y1<y2)
!     y2.........Lower boundary of integral range

      top = .false.
      bot = .false.
      if (z1 < z3) then
         y1=z1
         if (x1 > z1) then
            y1=x1
            top = .true.
         end if
         y2=z2
         if (x2 < z2) then
            y2=x2
            bot = .true.
         end if
      else
         y1=z2
         if (x1 > z2) then
            y1=x1
            top = .true.
         end if
         y2=z1
         if (x2 < z1) then
            y2=x2
            bot = .true.
         end if
      end if

! --- Set weights for forward interpolator and linear model contribution to TLM

      dy=y2-y1
      dz=z1-z2

      if (abs(dz) < 1D-14) then
         write(*,*) 'SUBLAYER: ERROR: dz is zero. dz = ',dz
         write(*,*) 'z1,z2,z3 = ',z1,z2,z3
         write(*,*) 'x1,x2    = ',x1,x2
         write(*,*) 't1,t2    = ',t1,t2
!bue         stop
         w1=0.0D0
         w2=0.0D0
         return
      else
         dzd=1.0D0/dz
      end if
      w1=(z1-y1)*dzd*dy
      w2=(z1-y2)*dzd*dy
      w10=w1
      w20=w2
      dx=(x2-x1)
      if (abs(dx) < 1D-14) then
         write(*,*) 'SUBLAYER: ERROR: dx is zero. dx = ',dx
         write(*,*) 'z1,z2,z3 = ',z1,z2,z3
         write(*,*) 'x1,x2    = ',x1,x2
         write(*,*) 't1,t2    = ',t1,t2
!bue         stop
         w1=0.0D0
         w2=0.0D0
         return
      else
         dxd=1.0D0/dx
      end if

      if (z1 < z3 .and. .not. bot) then
         d=(x2-z2)*dxd
         w1=w1+w2*d
         w2=w2*(1.0D0-d)
      else if (z1 > z3 .and. .not. top) then
         d=(x2-z2)*dxd
         w2=w2+w1*(1.0D0-d)
         w1=w1*d
      end if
      tot=t1*w1+t2*w2
     
! --- Provide NLM contribution to TLM (gradients w.r.t. Ps)

!        Determine gradient of 'tot' w.r.t. x1

      if (present(pzs1) .and. present(pzs2) .and. present(pzps)) then
        aa1=0.0D0
        aa2=0.0D0
        a1=0.0D0
        a2=0.0D0
        if (top) then
          a1=-(dy+(z1-y1))*dzd 
          a2=-(z1-y2)*dzd
        else if (z1 > z3) then
          a1=(x2-z2)*dxd*dxd*w10
          a2=-a1
        end if
        if (z1 < z3.and. .not. bot) then
          aa2=(x2-z2)*dxd*dxd*w20
          aa1=-aa2
          if (top) then
            a1=a1+a2*d+aa1
            a2=a2*(1.0D0-d)+aa2
          end if
        end if
        g1=a1*t1+a2*t2
                      
!        Determine gradient of 'tot' w.r.t. x2

        aa1=0.0D0
        aa2=0.0D0
        a1=0.0D0
        a2=0.0D0
        if (bot) then
          a1=dzd*(z1-y1)
          a2=((z1-y2)-dy)*dzd
        else if (z1 < z3) then
          a1=((x2-z2)*dxd+1.0D0)*dxd*w20 
          a2=-a1
        end if
        if (z1 > z3 .and. .not. top) then
          aa1=(1.0D0-(x2-z2)*dxd)*dxd*w10
          aa2=-aa1
          if (bot) then
            a2=a2+a1*(1.0D0-d)+aa2
            a1=a1*d+aa1
          end if
        end if
        g2=a1*t1+a2*t2

!        Accumulate for gradient w.r.t. Ps

        pzps=pzps+g1*pzs1+g2*pzs2

      end if

  END SUBROUTINE SUBLAYER

  SUBROUTINE SUBLAYER_v2(z1,z2,z3,x1,x2,t1,t2,w1,w2,pzs1,pzs2,zp1,zp2)
    IMPLICIT NONE
    REAL(8),intent(in) :: z1,z2,z3,x1,x2,t1,t2
    REAL(8),intent(out) :: w1,w2
    REAL(8),intent(out) :: zp1, zp2
    REAL(8),intent(in)  :: pzs1,pzs2

    REAL(8) y1,y2,tot,d,w10,w20,dz,dx,dy,dzd,dxd,g1,g2
    REAL(8) a1,a2,aa1,aa2
    logical bot,top

    top = .false.
    bot = .false.
    if (z1 < z3) then
       y1=z1
       if (x1 > z1) then
          y1=x1
          top = .true.
       end if
       y2=z2
       if (x2 < z2) then
          y2=x2
          bot = .true.
       end if
    else
       y1=z2
       if (x1 > z2) then
          y1=x1
          top = .true.
       end if
       y2=z1
       if (x2 < z1) then
          y2=x2
          bot = .true.
       end if
    end if

! --- Set weights for forward interpolator and linear model contribution to TLM

    dy=y2-y1
    dz=z1-z2

    if (abs(dz) < 1D-14) then
       write(*,*) 'SUBLAYER_v2: ERROR: dz is zero. dz = ',dz
       write(*,*) 'z1,z2,z3 = ',z1,z2,z3
       write(*,*) 'x1,x2    = ',x1,x2
       write(*,*) 't1,t2    = ',t1,t2
!bue         stop
       w1=0.0D0
       w2=0.0D0
       return
    else
       dzd=1.0D0/dz
    end if
    w1=(z1-y1)*dzd*dy
    w2=(z1-y2)*dzd*dy
    w10=w1
    w20=w2
    dx=(x2-x1)
    if (abs(dx) < 1D-14) then
       write(*,*) 'SUBLAYER_v2: ERROR: dx is zero. dx = ',dx
       write(*,*) 'z1,z2,z3 = ',z1,z2,z3
       write(*,*) 'x1,x2    = ',x1,x2
       write(*,*) 't1,t2    = ',t1,t2
!bue         stop
       w1=0.0D0
       w2=0.0D0
       return
    else
       dxd=1.0D0/dx
    end if

    if (z1 < z3 .and. .not. bot) then
       d=(x2-z2)*dxd
       w1=w1+w2*d
       w2=w2*(1.0D0-d)
    else if (z1 > z3 .and. .not. top) then
       d=(x2-z2)*dxd
       w2=w2+w1*(1.0D0-d)
       w1=w1*d
    end if
    tot=t1*w1+t2*w2
     
! --- Provide NLM contribution to TLM (gradients w.r.t. Ps)

!        Determine gradient of 'tot' w.r.t. x1

    aa1=0.0D0
    aa2=0.0D0
    a1=0.0D0
    a2=0.0D0
    if (top) then
      a1=-(dy+(z1-y1))*dzd 
      a2=-(z1-y2)*dzd
    else if (z1 > z3) then
      a1=(x2-z2)*dxd*dxd*w10
      a2=-a1
    end if
    if (z1 < z3.and. .not. bot) then
      aa2=(x2-z2)*dxd*dxd*w20
      aa1=-aa2
      if (top) then
        a1=a1+a2*d+aa1
        a2=a2*(1.0D0-d)+aa2
      end if
    end if
    g1=a1*t1+a2*t2
                  
!        Determine gradient of 'tot' w.r.t. x2

    aa1=0.0D0
    aa2=0.0D0
    a1=0.0D0
    a2=0.0D0
    if (bot) then
      a1=dzd*(z1-y1)
      a2=((z1-y2)-dy)*dzd
    else if (z1 < z3) then
      a1=((x2-z2)*dxd+1.0D0)*dxd*w20 
      a2=-a1
    end if
    if (z1 > z3 .and. .not. top) then
      aa1=(1.0D0-(x2-z2)*dxd)*dxd*w10
      aa2=-aa1
      if (bot) then
        a2=a2+a1*(1.0D0-d)+aa2
        a1=a1*d+aa1
      end if
    end if
    g2=a1*t1+a2*t2

!        Accumulate for gradient w.r.t. Ps

    zp1=g1*pzs1
    zp2=g2*pzs2

  END SUBROUTINE SUBLAYER_v2

  SUBROUTINE PPO_LINTV (PVLEV,PVI,KNI, KNPROF,KNO,PPO,PVO)

      !*
      !***s/r PPO_LINTV  - Linear interpolation and constant value extrapolation.
      !*                Input pressure levels can vary for each profile.
      !*
      !*Arguments
      !*     i   PVLEV(KNI,KNPROF)    : Vertical levels, pressure (source)
      !*     i   PVI(KNI,KNPROF)      : Vector to be interpolated (source)
      !*     i   KNI                     : Number of input levels (source)
      !*     i   KNPROF                  : Number of profiles
      !*     i   KNO                     : Number of output levels (destination)
      !*     i   PPO(KNO)                : Vertical levels, pressure (destination)
      !*     o   PVO(KNO,KNPROF)         : Interpolated profiles (destination)
      !*
      !!*    -------------------
      !**    Purpose: Performs the vertical interpolation in log of pressure
      !*              and constant value extrapolation of one-dimensional vectors.
      
      IMPLICIT NONE
      INTEGER,intent(in)  :: KNI, KNO,KNPROF
      REAL(8) ,intent(in)  :: PVLEV(KNI,KNPROF)
      REAL(8) ,intent(in)  :: PVI(KNI,KNPROF)
      REAL(8) ,intent(in)  :: PPO(KNO)
      REAL(8) ,intent(out) ::  PVO(KNO,KNPROF)
!********************************************
      INTEGER  JI, JK, JO, profileIndex, IK, IORDER
      
      REAL(8)     ZPI (0:KNI+1,KNPROF)
      REAL(8)     ZPVI(0:KNI+1,KNPROF)
      INTEGER  IL  (KNO    ,KNPROF)
      
      REAL(8) ZW1, ZW2
      REAL(8) ZP, XI, ZRT, ZP1, ZP2

      !
      !**   1. Initialization for vertical extrapolation (extra dummy levels)
      !     .  --------------------------------------------------------------
      !
      
      ZPI(0,:)=2000.D0
      ZPI(KNI+1,:)=2000.D0
      
      !
      !**      1.1 Determine if input pressure levels are in ascending or
      !     .      descending order.
      !     .     -------------------------------------------------------
      !
      IF ( PVLEV(1,1)  <  PVLEV(KNI,1) ) THEN
         IORDER = 1
      ELSE
         IORDER = -1
      END IF
      !
      !**   2. Compute pressure levels pressure
      !     .  ------------------------------------------------
      !

      !
      !**   2.1 Source levels
      !     .   -------------
      !
          
      ZPI(1:KNI,:) = PVLEV(1:KNI,:)

      !
      !
      !*    3.  Interpolate in log of pressure or extrapolate with constant value
      !*    .   for each destination pressure level
      !     .   -----------------------------------------------------------------
      !

      !
      !
      !*    .  3.1  Find the adjacent level below
      !     .       -----------------------------
      !
      !
      
      IL(:,:)=0
      !
      DO JI=1,KNI
         DO profileIndex = 1, KNPROF
            DO JO=1,KNO
              ZRT = PPO(JO)
               ZP = ZPI(JI,profileIndex)
               XI = SIGN(1.0D0,IORDER*(ZRT-ZP))
               IL(JO,profileIndex) = IL(JO,profileIndex) + MAX(0.0D0,XI)
            END DO
         END DO
      END DO
      !
      !
      !*    .  3.2  Fill extra levels, for constant value extrapolation
      !     .       ---------------------------------------------------
      !
      
      DO profileIndex = 1, KNPROF
         DO JK = 1, KNI
            ZPVI(JK,profileIndex) = PVI(JK,profileIndex)
         END DO
      END DO
      DO profileIndex = 1, KNPROF
         ZPVI(0    ,profileIndex) = PVI(1  ,profileIndex)
         ZPVI(KNI+1,profileIndex) = PVI(KNI,profileIndex)
      END DO
      !
      !
      !*    .  3.3  Interpolation/extrapolation
      !     .       ---------------------------
      !
      
      DO profileIndex = 1, KNPROF
        DO JO=1,KNO
          ZP = PPO(JO)
          IK = IL(JO,profileIndex)
          ZP1 = ZPI(IK  ,profileIndex)
          ZP2 = ZPI(IK+1,profileIndex)
          ZW1 = LOG(ZP/ZP2)/LOG(ZP1/ZP2)
          ZW2 = 1.D0 - ZW1
          PVO(JO,profileIndex) = ZW1*ZPVI(IK,profileIndex) +  ZW2*ZPVI(IK+1,profileIndex)
        END DO
      END DO
     
      
    END SUBROUTINE PPO_LINTV
    
  end module presProfileOperators_mod
