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
!! MODULE presProfileOperators (prefix="ppo" category='7')
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
  public :: ppo_intAvg, ppo_intAvgTl, ppo_intAvgAd, ppo_lintv

  contains

    SUBROUTINE ppo_INTAVG (PVLEV,PVI,KNIDIM,KNI,KNPROF,KNO,PPO,PVO)
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
!     i   KFLAG                   : Indicates purpose of calc
!                                   0 - Application of TLM as forward model
!                                   1 - Application as TLM with increments
!                                   2 - Setting of adjoint elements
!     i   PVLEV(KNIDIM,KNPROF)    : Vertical levels, pressure (source)
!    i/o  PVI(KNIDIM,KNPROF)      : kflag=0, Input vector to be interpolated (source)
!                                   kflag=1, Input increments of above 
!                                   kflag=2, Output adjoint of above
!     i   PVIG(KNIDIM,KNPROF)     : kflag>0, Vector to be interpolated (source)
!                                   Not used otherwise 
!    i/o  PPS(KNPROF)             : kflag=0, Not needed (source)
!                                   kflag=1, Input surface pressure increments
!                                   kflag=2, Output adjoint of above
!     i   KNIDIM                  : Dimension of input levels (source)
!     i   KNI                     : Number of input levels (source)
!     i   KNPROF                  : Number of profiles
!     i   KNO                     : Number of output levels (destination)
!     i   PPO(KNO)                : Vertical levels, pressure (destination)
!    o/i  PVO(KNO,KNPROF)         : kflag=0, Output interpolated profiles (destination)
!                                   kflag=1, Output increments of above         
!                                   kflag=2, Input adjoint of above            
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
!
!     2) Cases:
!
!               a) KFLAG=0 for application as forward interpolator
!
!                  Y = sum(H*PVI)
!                    = PVO
!     
!               with ZPZ=H on output of LAYERAVG
!
!               b) KFLAG=1 for TLM case:
!
!                  dY = sum(H*PVI) + PPS*sum(PVIG*dH/dPs))
!                     =     ZPVO    +       ZPVOPS
!                     = PVO
!
!               where
!
!                     dPZ(i)/dPs = sum_k (dH(i)/dPVLEV(k) * dPVLEV(k)/dPs)
!                                = sum_k (dH(i)/dZLNPI(k) * zpresb(k)/PVLEV(k))
!
!               with ZPZ(k)=zpresb(k)/PVLEV(k) on input to LAYERAVG and
!                    ZPZ(i)=H(i) on output of LAYERAVG
!            
!               c) KFLAG=2 for adjoint case:
!               
!                  PVI(i,jn) = PVI(i,jn) + sum_jo (PVO(jo,jn)*H(i))
!                            = PVI(i,jn) + sum_jo (ZPZ(jo,i))
!                    
!                  PPS(jn) = PPS(jn) + sum_jo (PVO(jo,jn)*sum_i (PVIG(i,jn)*dH(i)/dPs))
!                          = PPS(jn) + sum_jo (ZPZPS(jo))
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  JO, JN
      INTEGER  KNIDIM, KNI, KNO, KNPROF, KFLAG

      REAL*8 PVLEV(KNIDIM,KNPROF)
      REAL*8 PPO(KNO), PVO(KNO,KNPROF)
      REAL*8 PVI(KNIDIM,KNPROF)

      REAL*8  ZLNPI (KNI),ZPZ(KNI),ZPG(KNI),ZPS(KNI)
      REAL*8  ZLNPO (KNO),ZPZPS,PSS

! --- KFLAG defines type of calculation: 0=Forward, 1=TLM, 2=ADJ

      KFLAG=0

! --- Apply weighted averaging 

      ZLNPO(1:KNO)=LOG(PPO(1:KNO))
      PSS=0.0D0
      ZPG(1:KNI)=0.0D0
      ZPZ(1:KNI)=0.0D0
      ZPS(1:KNI)=0.0D0
      ZPZPS=0.0D0

      DO JN = 1, KNPROF
         ZLNPI(1:KNI)=LOG(PVLEV(1:KNI,JN))
         PVO(1:KNO,JN)=0.0D0

         DO JO=1,KNO
            CALL LAYERAVG(KFLAG,ZLNPO,ZLNPI,PVI(1:KNI,JN),ZPG,KNO,KNI,  &
                          JO,PSS,PVO(JO,JN),ZPZ,ZPS,ZPZPS)
         ENDDO
      END DO

    END SUBROUTINE PPO_INTAVG
    

    SUBROUTINE ppo_INTAVGTL(PVLEV,dP_dPsfc,PVI,PVIG,PPS,KNIDIM,KNI,KNPROF,KNO,PPO,PVO)
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
!            LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!          S. Pellerin, ARMA, August 2008
!          - Optimisation, call to LAYERAVG2
!
!Arguments
!     i   PVLEV(KNIDIM,KNPROF)    : Vertical levels, pressure (source)
!     i   dP_dPsfc(KNIDIM,KNPROF) : derivative of pressure wrt sfc pressure for vertical coordinate (source)
!     i   PVI(KNIDIM,KNPROF)      : Increments of vector on input levels
!     i   PVIG(KNIDIM,KNPROF)     : Vector on input levels (source)
!     i   PPS(KNPROF)             : Input surface pressure increments
!     i   KNIDIM                  : Dimension of input levels (source)
!     i   KNI                     : Number of input levels (source)
!     i   KNPROF                  : Number of profiles
!     i   KNO                     : Number of output levels (destination)
!     i   PPO(KNO)                : Vertical levels, pressure (destination)
!     o   PVO(KNO,KNPROF)         : Increments of vector on output levels
!
!    -------------------
!
! Purpose: Application of tangent linear
!          of piecewise weighted averaging
!          in log of pressure of one-dimensional vectors.
!
! Comments:
!
!     1) vlev in eta coordinate associated to pvlev in pressure.
!
!     2) Cases:
!
!               a) KFLAG=0 for application as forward interpolator
!
!                  Y = sum(H*PVI)
!                    = PVO
!
!               with ZPZ=H on output of LAYERAVG
!
!               b) KFLAG=1 for TLM case:
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
!               c) KFLAG=2 for adjoint case:
!
!                  PVI(i,jn) = PVI(i,jn) + sum_jo (PVO(jo,jn)*H(i))
!                            = PVI(i,jn) + sum_jo (ZPZ(jo,i))
!
!                  PPS(jn) = PPS(jn) + sum_jo (PVO(jo,jn)*sum_i (PVIG(i,jn)*dH(i)/dPs))
!                          = PPS(jn) + sum_jo (ZPZPS(jo))
! 
!-------------------------------------------------------------------- 
!
      IMPLICIT NONE
      INTEGER  JI, JN
      INTEGER  KNIDIM, KNI, KNO, KNPROF, KFLAG

      REAL*8 PVLEV(KNIDIM,KNPROF)
      REAL*8 dP_dPsfc(KNIDIM,KNPROF)
      REAL*8 PPO(KNO), PVO(KNO,KNPROF)
      REAL*8 PVI(KNIDIM,KNPROF)
      REAL*8 PVIG(KNIDIM,KNPROF)
      REAL*8 PPS(KNPROF)

      REAL*8  ZLNPI (KNI,knprof),ZPZ(KNI), ZPS(KNI,knprof)
      REAL*8  ZLNPO (KNO),ZPZPS,ZPVOPS(kno,knprof),PSS,ZPVO(kno,knprof)

! --- KFLAG defines type of calculation: 0=Forward, 1=TLM, 2=ADJ

      KFLAG=1

! --- Apply weighted averaging 

      ZLNPO(1:KNO)=LOG(PPO(1:KNO))
      PSS=0.0D0

      DO JN = 1, KNPROF
        DO JI = 1, KNI
          ZPS(JI,jn)=dP_dPsfc(ji,jn)/PVLEV(JI,JN)
          ZLNPI(ji,JN)=LOG(PVLEV(ji,JN))
        END DO
      enddo

      CALL LAYERAVG2(KFLAG,ZLNPO,ZLNPI,pvig(1:kni,:),PVI(1:kni,:),KNO,  &
                     KNI,PPS,ZPVO,ZPZ,ZPS,ZPVOPS,knprof,pvo)

    END SUBROUTINE PPO_INTAVGTL


    SUBROUTINE ppo_INTAVGAD(PVLEV,dP_dPsfc,PVI,PVIG,PPS,KNIDIM,KNI,KNPROF,KNO,PPO,PVO)
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
!            LINTV2, LLINTV2, and ALINTV2 by J. Halle et al.
!          S. Pellerin, ARMA, August 2008
!          - Optimisation, call to LAYERAVG2
!
!Arguments
!     i   PVLEV(KNIDIM,KNPROF)    : Vertical levels, pressure (source)
!     o   PVI(KNIDIM,KNPROF)      : Adjoint of increments on input levels (output)
!     i   PVIG(KNIDIM,KNPROF)     : Vector to be interpolated (source)
!     o   PPS(KNPROF)             : Output adjoint of sfc P incr
!     i   KNIDIM                  : Dimension of input levels (source)
!     i   KNI                     : Number of input levels (source)
!     i   KNPROF                  : Number of profiles
!     i   KNO                     : Number of output levels
!     i   PPO(KNO)                : Vertical levels, pressure
!     i   PVO(KNO,KNPROF)         : Input adjoint of interpolated profiles
!
!    -------------------
!
!!Purpose: Performs calculation for adjoint of
!          piecewise weighted averaging
!          in log of pressure of one-dimensional vectors.
!
! Comments:
!
!     1) vlev in eta coordinate associated to pvlev in pressure.
!
!     2) Cases:
!
!               a) KFLAG=0 for application as forward interpolator  
!
!                  Y = sum(H*PVI)
!                    = PVO
!
!               with ZPZ=H on output of LAYERAVG
!
!               b) KFLAG=1 for TLM case:
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
!               with ZPZ(k)=dP_dPsfc(k)/PVLEV(k) on input to LAYERAVG and
!                    ZPZ(i)=H(i) on output of LAYERAVG
!
!               c) KFLAG=2 for adjoint case:
!
!                  PVI(i,jn) = PVI(i,jn) + sum_jo (PVO(jo,jn)*H(i))
!                            = PVI(i,jn) + sum_jo (ZPZ(jo,i))
!
!                  PPS(jn) = PPS(jn) + sum_jo (PVO(jo,jn)*sum_i (PVIG(i,jn)*dH(i)/dPs))
!                          = PPS(jn) + sum_jo (ZPZPS(jo))
!
!--------------------------------------------------------------------  
      IMPLICIT NONE
      INTEGER  JI, JO, JN
      INTEGER  KNIDIM, KNI, KNO, KNPROF, KFLAG

      REAL*8 PVLEV(KNIDIM,KNPROF)
      REAL*8 dP_dPsfc(KNIDIM,KNPROF)
      REAL*8 PPO(KNO), PVO(KNO,KNPROF)
      REAL*8 PVI(KNIDIM,KNPROF)
      REAL*8 PVIG(KNIDIM,KNPROF)
      REAL*8 PPS(KNPROF)

      REAL*8 ZLNPI (KNI,knprof),ZPZ(KNI),ZPG(KNI,knprof)
      REAL*8 ZPS(KNI,knprof),zpvo(kno,knprof)
      REAL*8 ZLNPO (KNO),ZPZPS(kno,knprof),ZPVOPS(kno,knprof)

! --- KFLAG defines type of calculation: 0=Forward, 1=TLM, 2=ADJ

      KFLAG=2

! --- Apply weighted averaging 

      ZLNPO(1:KNO)=LOG(PPO(1:KNO))
      ZPG(1:KNI,:)=0.0D0
      PVI(:,:)=0.0D0

      DO JN = 1, KNPROF
        DO JI = 1, KNI
          ZPS(JI,jn)=dP_dPsfc(ji,jn)/PVLEV(JI,JN)
          ZLNPI(ji,jn)=LOG(PVLEV(ji,JN))
        END DO
      enddo

      CALL LAYERAVG2(KFLAG,ZLNPO,ZLNPI,pvig(1:kni,:),pvi(1:kni,:),  &
                     KNO,KNI,PPS,PVO,ZPZ,ZPS,ZPZPS,knprof,zpvo)

    END SUBROUTINE PPO_INTAVGAD


    SUBROUTINE LAYERAVG(KFLAG,PX1,PX2,PY2,PY2INCR,KN1,KN2,KI,PPS,PMEAN,PZ,PZS,PZPS)

      IMPLICIT NONE

      INTEGER KN1,KN2,KI,KFLAG
      REAL*8 PX1(KN1),PX2(KN2),PY2(KN2),PY2INCR(KN2)
      REAL*8 PPS,PMEAN,PZ(KN2),PZPS,PZS(KN2)
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
!     KFLAG:    Indicates purpose of calc
!               0 - Application as forward interpolator
!               1 - Application for TLM* increments
!               2 - Setting of adjoint elements
!
!     PX1:      Reference levels (e.g. lnP; in increasing values)
!     PX2:      Available levels (e.g. lnP; in increasing values)
!     PY2:      Values at levels (e.g. temperature) 
!     PY2INCR:  Not relevant when kflag=0 or 2
!               Increments when kflag=1
!     PZ:       Extra array related to gradients w.r.t. Ps for KFLAG.gt.0
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) when KFLAG=2.
!               Output otherwise.
!     KN1:      Dimension of PX1.
!     KN2:      Dimension of other arrays.
!     PZS:      dlnP/dPs
!     KI:       Identifies region of interest: PX1(KI-1) to PX1(KI+1)
!     PPS:      Surface pressure increment when KFLAG=1, not needed otherwise.
!
!   OUTPUT:
!
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) using
!               background values when kflag=0 and increments when kflag=1.
!     PZ:       Resultant accumulated contribution factors
!               Adjoint when KFLAG=2 
!               TLM when KFLAG=1 or 0.
!     PZPS:     Surface pressure related adjoint term when KFLAG=2
!               Surface pressure related TLM*increment term when KFLAG=1
!
!-----------------------------------------------------------
      REAL*8 Z1,Z2,Z3,ZW1,ZW2
      REAL*8 zsum
      INTEGER J,IC,ISKIP

! --- Identify boundary points

      z2=px1(ki)

      if (ki.eq.1) then 
         z1=2.0D0*z2-px1(ki+1)
      else
         z1=px1(ki-1)
      endif   

      if (ki.eq.kn1) then
         z3=2.0D0*z2-z1
      else   
         z3=px1(ki+1)
      endif
      if (z3.gt.px2(kn2)) z3=px2(kn2)

      iskip=0
      if (z2.ge.px2(kn2)) then
         z3=px2(kn2)
         z2=px2(kn2)
         iskip=1
      endif

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

      pzps=0.0D0
      pz(1:kn2)=0.0D0
      ic=0
      do j=1,kn2-1
         if (px2(j).ge.z3) go to 1000

         if (px2(j).le.z2.and.px2(j+1).gt.z1) then 

            call sublayer(z1,z2,z3,px2(j),px2(j+1),kflag,  &
                          py2(j),py2(j+1),zw1,zw2,  &
                          pzs(j),pzs(j+1),pzps)
            pz(j)=pz(j)+zw1
            pz(j+1)=pz(j+1)+zw2
            ic=1
         endif

         if (px2(j).lt.z3.and.px2(j+1).ge.z2.and.iskip.eq.0) then

            call sublayer(z3,z2,z1,px2(j),px2(j+1),kflag,  &
                          py2(j),py2(j+1),zw1,zw2,  &
                          pzs(j),pzs(j+1),pzps)
            pz(j)=pz(j)+zw1
            pz(j+1)=pz(j+1)+zw2
            ic=1
         endif
      enddo
 1000 continue
      if (ic.eq.0) pz(j)=1.0D0

! --- Apply forward interpolator (kflag=0), determine TLM*increment (kflag=1)
!     or use TLM for adjoint calc (kflag=2)
 
      zsum=0.0D0
      if (kflag.eq.0) then
         pmean=0.0D0
         do j=1,kn2
            pmean=pmean+pz(j)*py2(j)
            zsum=zsum+pz(j)
         end do
         pmean=pmean/zsum
      else if (kflag.eq.1) then
         pmean=0.0D0
         do j=1,kn2      
            pmean=pmean+pz(j)*py2incr(j)
            zsum=zsum+pz(j)
         end do
         pmean=pmean/zsum
         if (ic.ne.0) pzps=pzps*pps/zsum
      else if (kflag.eq.2) then
         do j=1,kn2
            zsum=zsum+pz(j)         
         end do
         pz(1:kn2)=pz(1:kn2)*pmean/zsum
         if (ic.ne.0) pzps=pzps*pmean/zsum
      endif

    END SUBROUTINE LAYERAVG


    SUBROUTINE LAYERAVG2(KFLAG,PX1,PX2,PY2,PY2INCR,KNO,KNI,PPS,PMEAN,PZ,PZS,PZPS,knprof,pvo)

      IMPLICIT NONE
      INTEGER KNO,KNI,KFLAG,knprof
      REAL*8 PX1(KNO),PX2(KNI,knprof),PY2(KNI,knprof),PY2INCR(KNI,knprof)
      REAL*8 PPS(knprof),PMEAN(kno,knprof),PZ(KNI),PZPS(kno,knprof)
      REAL*8 PZS(KNI,knprof),pvo(kno,knprof)
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
!
!     KFLAG:    Indicates purpose of calc
!               0 - Application as forward interpolator
!               1 - Application for TLM* increments
!               2 - Setting of adjoint elements
!
!     PX1:      Reference levels (e.g. lnP; in increasing values)
!     PX2:      Available levels (e.g. lnP; in increasing values)
!     PY2:      Values at levels (e.g. temperature) 
!     PY2INCR:  Not relevant when kflag=0 or 2
!               Increments when kflag=1
!     PZ:       Extra array related to gradients w.r.t. Ps for KFLAG.gt.0
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) when KFLAG=2.
!               Output otherwise.
!     KNO:      Dimension of PX1.
!     KNI:      Dimension of other arrays.
!     PZS:      dlnP/dPs
!     KI:       Identifies region of interest: PX1(KI-1) to PX1(KI+1)
!     PPS:      Surface pressure increment when KFLAG=1, not needed otherwise.
!
!   OUTPUT:
!
!     PMEAN:    Mean weigthed value for PX1(KI-1) to PX1(KI+1) using
!               background values when kflag=0 and increments when kflag=1.
!     PZ:       Resultant accumulated contribution factors
!               Adjoint when KFLAG=2 
!               TLM when KFLAG=1 or 0.
!     PZPS:     Surface pressure related adjoint term when KFLAG=2
!               Surface pressure related TLM*increment term when KFLAG=1
!
!-----------------------------------------------------------
      REAL*8 Z1,Z2,Z3,ZW1,ZW2
      REAL*8 zsum
      INTEGER J,IC,ISKIP,jo,jn

! --- Identify boundary points

      do jn = 1, knprof
        do jo = 1, kno
          z2=px1(jo)

          if (jo.eq.1) then 
            z1=2.0D0*z2-px1(jo+1)
          else
            z1=px1(jo-1)
          endif   

          if (jo.eq.kno) then
            z3=2.0D0*z2-z1
          else   
            z3=px1(jo+1)
          endif
          if (z3.gt.px2(kni,jn)) z3=px2(kni,jn)

          iskip=0
          if (z2.ge.px2(kni,jn)) then
            z3=px2(kni,jn)
            z2=px2(kni,jn)
            iskip=1
          endif

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

          pzps(jo,jn)=0.0D0
          pz(1:kni)=0.0D0
          ic=0

!          call sublayer2(z1,z2,z3,px2(:,jn),kflag,  &
!                         py2(:,jn),pzs(:,jn),pzps(jo,jn),pz,kni,ic,iskip)
          do j=1,kni-1
            if (px2(j,jn).ge.z3) go to 1000

            if (px2(j,jn).le.z2.and.px2(j+1,jn).gt.z1) then 

              call sublayer(z1,z2,z3,px2(j,jn),px2(j+1,jn),kflag,  &
                            py2(j,jn),py2(j+1,jn),zw1,zw2,  &
                            pzs(j,jn),pzs(j+1,jn),pzps(jo,jn))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              ic=1
            endif

            if (px2(j,jn).lt.z3.and.px2(j+1,jn).ge.z2.and.iskip.eq.0) then

              call sublayer(z3,z2,z1,px2(j,jn),px2(j+1,jn),kflag,  &
                            py2(j,jn),py2(j+1,jn),zw1,zw2,  &
                            pzs(j,jn),pzs(j+1,jn),pzps(jo,jn))
              pz(j)=pz(j)+zw1
              pz(j+1)=pz(j+1)+zw2
              ic=1
            endif
          enddo
 1000     continue
          if (ic.eq.0) pz(j)=1.0D0
! --- Apply forward interpolator (kflag=0), determine TLM*increment (kflag=1)
!     or use TLM for adjoint calc (kflag=2)

          zsum=0.0D0
          if (kflag.eq.0) then
            pmean(jo,jn)=0.0D0
            do j=1,kni
              pmean(jo,jn)=pmean(jo,jn)+pz(j)*py2(j,jn)
              zsum=zsum+pz(j)
            end do
            pmean(jo,jn)=pmean(jo,jn)/zsum
          else if (kflag.eq.1) then
            pmean(jo,jn)=0.0D0
            do j=1,kni      
              pmean(jo,jn)=pmean(jo,jn)+pz(j)*py2incr(j,jn)
              zsum=zsum+pz(j)
            end do
            pmean(jo,jn)=pmean(jo,jn)/zsum
            if (ic.ne.0) pzps(jo,jn)=pzps(jo,jn)*pps(jn)/zsum
            PVO(JO,JN)=pmean(jo,jn)+PZPS(jo,jn)
          else if (kflag.eq.2) then
            do j=1,kni
              zsum=zsum+pz(j)         
            end do
            pz(1:kni)=pz(1:kni)*pmean(jo,jn)/zsum
            py2incr(1:KNI,JN)=py2incr(1:KNI,JN)+PZ(1:KNI)
            if (ic.ne.0) pzps(jo,jn)=pzps(jo,jn)*pmean(jo,jn)/zsum
            PPS(JN)=PPS(JN)+PZPS(jo,jn)
          endif
        enddo
      enddo

    END SUBROUTINE LAYERAVG2


    SUBROUTINE SUBLAYER(z1,z2,z3,x1,x2,imode,t1,t2,w1,w2,pzs1,pzs2,pzps)
      IMPLICIT NONE
      REAL*8 z1,z2,z3,x1,x2,t1,t2,w1,w2,pzs1,pzs2,pzps
      INTEGER imode
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
!       imode......0:  Forward interpolator and linear model contribution
!                      to TLM only
!                  >0: Also provide NLM contribution to TLM (gradients 
!                      w.r.t. Ps)
!       t1.........Variable value at upper NWP level.
!       t2.........Variable value at lower NWP level.
!       dzs1.......dlnP/dPs = dx1/dPs
!       dzs2.......dlnP/dPs = dx2/dPs
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
      REAL*8 y1,y2,tot,d,w10,w20,dz,dx,dy,dzd,dxd,g1,g2
      REAL*8 a1,a2,aa1,aa2
      integer ibot,itop
!
! --- Identify and set upper and lower boundaries of
!     integration/averaging layers.
!
!     y1.........Upper boundary of integral range (y1<y2)
!     y2.........Lower boundary of integral range

      itop=0
      ibot=0
      if (z1.lt.z3) then
         y1=z1
         if (x1.gt.z1) then
            y1=x1
            itop=1
         endif
         y2=z2
         if (x2.lt.z2) then
            y2=x2
            ibot=1
         endif
      else
         y1=z2
         if (x1.gt.z2) then
            y1=x1
            itop=1
         endif
         y2=z1
         if (x2.lt.z1) then
            y2=x2
            ibot=1
         endif
      endif

! --- Set weights for forward interpolator and linear model contribution to TLM

      dy=y2-y1
      dz=z1-z2
!     dzd=1.0/dz
      if (abs(dz).lt.1D-14) then
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
      endif
      w1=(z1-y1)*dzd*dy
      w2=(z1-y2)*dzd*dy
      w10=w1
      w20=w2
      dx=(x2-x1)
!     dxd=1.0/dx
      if (abs(dx).lt.1D-14) then
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
      endif

      if (z1.lt.z3.and.ibot.eq.0) then
         d=(x2-z2)*dxd
         w1=w1+w2*d
         w2=w2*(1.0D0-d)
      else if (z1.gt.z3.and.itop.eq.0) then
         d=(x2-z2)*dxd
         w2=w2+w1*(1.0D0-d)
         w1=w1*d
      end if
      tot=t1*w1+t2*w2
     
! --- Provide NLM contribution to TLM (gradients w.r.t. Ps)

      IF (imode.gt.0) THEN

!        Determine gradient of 'tot' w.r.t. x1

         aa1=0.0D0
         aa2=0.0D0
         a1=0.0D0
         a2=0.0D0
         if (itop.eq.1) then
            a1=-(dy+(z1-y1))*dzd 
            a2=-(z1-y2)*dzd
         else if (z1.gt.z3) then
            a1=(x2-z2)*dxd*dxd*w10
            a2=-a1
         end if
         if (z1.lt.z3.and.ibot.eq.0) then
            aa2=(x2-z2)*dxd*dxd*w20
            aa1=-aa2
            if (itop.eq.1) then
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
         if (ibot.eq.1) then
            a1=dzd*(z1-y1)
            a2=((z1-y2)-dy)*dzd
         else if (z1.lt.z3) then
            a1=((x2-z2)*dxd+1.0D0)*dxd*w20 
            a2=-a1
         end if
         if (z1.gt.z3.and.itop.eq.0) then
            aa1=(1.0D0-(x2-z2)*dxd)*dxd*w10
            aa2=-aa1
            if (ibot.eq.1) then
               a2=a2+a1*(1.0D0-d)+aa2
               a1=a1*d+aa1
            end if
         end if
         g2=a1*t1+a2*t2

!        Accumulate for gradient w.r.t. Ps

         pzps=pzps+g1*pzs1+g2*pzs2

      ENDIF

    END SUBROUTINE SUBLAYER 

    SUBROUTINE PPO_LINTV (PVLEV,PVI,KNIDIM,KNI, KNPROF,KNO,PPO,PVO)

      !*
      !***s/r PPO_LINTV  - Linear interpolation and constant value extrapolation.
      !*                Input pressure levels can vary for each profile.
      !*
      !*Arguments
      !*     i   PVLEV(KNIDIM,KNPROF)    : Vertical levels, pressure (source)
      !*     i   PVI(KNIDIM,KNPROF)      : Vector to be interpolated (source)
      !*     i   KNIDIM                  : Dimension of input levels (source)
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
      
      INTEGER  JI, JK, JO, JN, IK, IORDER
      INTEGER  KNIDIM, KNI, KNO, KNPROF, ILEN, IERR
      
      REAL*8 PVLEV(KNIDIM,KNPROF)
      REAL*8 PPO(KNO), PVO(KNO,KNPROF)
      REAL*8 PVI(KNIDIM,KNPROF)
      
      REAL*8     ZPI (0:KNI+1,KNPROF)
      REAL*8     ZPO (KNO    ,KNPROF)
      REAL*8     ZPVI(0:KNI+1,KNPROF)
      INTEGER  IL  (KNO    ,KNPROF)
      
      REAL*8 ZW1, ZW2
      REAL*8 ZP, XI, ZRT, ZP1, ZP2
      !
      !**   0. Dynamic memory allocation for temporary vectors
      !     .  -----------------------------------------------
      !
050   CONTINUE
      !
      !     ... removed and replaced by automatic arrays, jh feb 2003 ......
      !
      !**   1. Initialization for vertical extrapolation (extra dummy levels)
      !     .  --------------------------------------------------------------
      !
100   CONTINUE
      
      ZPI(0,:)=2000.D0
      ZPI(KNI+1,:)=2000.D0
      
      !
      !**      1.1 Determine if input pressure levels are in ascending or
      !     .      descending order.
      !     .     -------------------------------------------------------
      !
      IF ( PVLEV(1,1) .LT. PVLEV(KNI,1) ) THEN
         IORDER = 1
      ELSE
         IORDER = -1
      ENDIF
      !
      !**   2. Compute pressure levels pressure
      !     .  ------------------------------------------------
      !
200   CONTINUE
      !
      !**   2.1 Source levels
      !     .   -------------
      !
      
      DO JN = 1, KNPROF
         DO JK = 1, KNI
            ZPI(JK,JN) = PVLEV(JK,JN)
         ENDDO
      ENDDO
      !
      !**   2.2 Destination levels
      !     .   ------------------
      !
      
      DO JN = 1, KNPROF
         DO JK = 1, KNO
            ZPO(JK,JN) = PPO(JK)
         ENDDO
      ENDDO
      !
      !*    3.  Interpolate in log of pressure or extrapolate with constant value
      !*    .   for each destination pressure level
      !     .   -----------------------------------------------------------------
      !
300   CONTINUE
      !
      !
      !*    .  3.1  Find the adjacent level below
      !     .       -----------------------------
      !
310   CONTINUE
      !
      
      IL(:,:)=0
      !
      DO JI=1,KNI
         DO JN = 1, KNPROF
            DO JO=1,KNO
               ZRT = ZPO(JO,JN)
               ZP = ZPI(JI,JN)
               XI = SIGN(1.0D0,IORDER*(ZRT-ZP))
               IL(JO,JN) = IL(JO,JN) + MAX(0.0D0,XI)
            ENDDO
         ENDDO
      ENDDO
      !
      !
      !*    .  3.2  Fill extra levels, for constant value extrapolation
      !     .       ---------------------------------------------------
      !
320   CONTINUE
      
      DO JN = 1, KNPROF
         DO JK = 1, KNI
            ZPVI(JK,JN) = PVI(JK,JN)
         ENDDO
      ENDDO
      DO JN = 1, KNPROF
         ZPVI(0    ,JN) = PVI(1  ,JN)
         ZPVI(KNI+1,JN) = PVI(KNI,JN)
      ENDDO
      !
      !
      !*    .  3.3  Interpolation/extrapolation
      !     .       ---------------------------
      !
330   CONTINUE
      
      DO JN = 1, KNPROF
         DO JO=1,KNO
            IK = IL(JO,JN)
            ZP = ZPO(JO,JN)
            ZP1 = ZPI(IK  ,JN)
            ZP2 = ZPI(IK+1,JN)
            ZW1 = LOG(ZP/ZP2)/LOG(ZP1/ZP2)
            ZW2 = 1.D0 - ZW1
            PVO(JO,JN) = ZW1*ZPVI(IK,JN) +  ZW2*ZPVI(IK+1,JN)
         ENDDO
      ENDDO
      !
      !*    4.  Deallocate memory
      !     .   -----------------
      !
400   CONTINUE
      !
      !     ... removed and replaced by automatic arrays, jh feb 2003 ......
      !
      
    END SUBROUTINE PPO_LINTV
    
  end module presProfileOperators_mod
