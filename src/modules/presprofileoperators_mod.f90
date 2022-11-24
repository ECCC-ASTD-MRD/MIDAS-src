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

module presProfileOperators_mod
  ! MODULE presProfileOperators_mod (prefix='ppo' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Vertical interpolation, integration, and layer averaging subroutines.
  !           Includes the special routines designed to interpolate to the 
  !           (widely spaced) RTTOV pressure levels.
  !
  use utilities_mod

  implicit none
  save
  private

  ! public procedures
  
  ! Linear interpolation routine
  public :: ppo_lintv
  
  ! Stand-alone interpolator calling routine providing interpolation weights 
  ! For weights W(:,:) and initial vector X(:), the interpolated vector would be W*X
  ! Will call one of the following prodedures:
  !    ppo_layeravgInterpWgts     - weights for linear piecewise weighted averaging interpolation
  !    ppo_piecewiseLinearWgts    - weights for piecewise linear interpolation
  ! The first of routines is similar in content to those used with RTTOV pressure levels.
  public :: ppo_vertInterpWgts    ! Main call routine

  ! Routines common to ppo_vertIntegWgts and ppo_vertAvgWgts indicated below.
  !     ppo_getLevelIndex - get the vertical input level index for level
  !                         within target layer and nearest specified layer boundary
  public :: ppo_vertLayersSetup
  
  ! Work arrays for ppo_vertInteg* and ppo_vertAvg*
  real(8), allocatable :: boundaries(:), weights(:,:)
  
  ! Stand-alone integration routines providing integration weights 
  ! For weights W(:,:) and initial vector X(:), the integrated values would be W*X  
  public :: ppo_vertIntegWgts

  ! Stand-alone layer averaging routine (in lp(P)) providing weights 
  ! For weights W(:,:) and initial vector X(:), the average values would be W*X  
  public :: ppo_vertAvgWgts
  
  contains

  !--------------------------------------------------------------------------
  ! PPO_LINTV
  !--------------------------------------------------------------------------
  SUBROUTINE PPO_LINTV (PVLEV,PVI,KNI, KNPROF,KNO,PPO,PVO)
      !
      !:Purpose: To perform the vertical interpolation in log of pressure and
      !          constant-value extrapolation of one-dimensional vectors. Input
      !          pressure levels can vary for each profile.
      !
      implicit none

      ! Arguments:
      real(8) ,intent(in)  :: pvlev(kni,knprof) ! Vertical levels, pressure (source)
      real(8) ,intent(in)  :: pvi(kni,knprof)   ! Vector to be interpolated (source)
      integer ,intent(in)  :: kni               ! Number of input levels (source)
      integer ,intent(in)  :: knprof            ! Number of profiles
      integer ,intent(in)  :: kno               ! Number of output levels (destination)
      real(8) ,intent(in)  :: ppo(kno,knprof)   ! Vertical levels, pressure (destination)
      real(8) ,intent(out) :: pvo(kno,knprof)   ! Interpolated profiles (destination)

      ! Locals:
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
              ZRT = PPO(JO,profileIndex)
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
          ZP = PPO(JO,profileIndex)
          IK = IL(JO,profileIndex)
          ZP1 = ZPI(IK  ,profileIndex)
          ZP2 = ZPI(IK+1,profileIndex)
          ZW1 = LOG(ZP/ZP2)/LOG(ZP1/ZP2)
          ZW2 = 1.D0 - ZW1
          PVO(JO,profileIndex) = ZW1*ZPVI(IK,profileIndex) +  ZW2*ZPVI(IK+1,profileIndex)
        END DO
      END DO
           
  END SUBROUTINE PPO_LINTV

  !==========================================================================
  !---- Stand-alone interpolation routine providing interpolation weights ---
  ! For weights W(:,:) and initial vector X(:), the integrated values would be W*X.  

  !--------------------------------------------------------------------------
  ! ppo_vertInterpWgts
  !--------------------------------------------------------------------------
  subroutine ppo_vertInterpWgts(pressInput,pressTarget,numInputLevs,numTargetLevs, &
                                 wgts,kstart,kend,method_opt,skipType_opt,outbound_opt,success_opt)
    !
    !:Purpose: Determination of interpolation weights for interpolation to points 
    !          in a profile. Applies interpolation in log(Pressure).
    !
    !:Input:
    !       :pressInput:       pressure on reference column levels assumed to be in ascending order
    !       :pressTarget:      target pressure levels assumed to be in ascending order
    !       :numInputLevs:     number of input/reference levels
    !       :numTargetLevs:    number of target levels
    !       :method_opt:       Specified interpolation method
    !       :skipType_opt:     Skipping processing of specific target levels depending on case:
    !                          'default'  - extrapolation allowed and skipping application via input success_opt only
    !                          'noExtrap' - no extrapolation as well as additional skipping via input success_opt
    !                          'doAll&noExtrap'  - no extrapolation only (all other levels processed)
    !:InOut:
    !       :outbound_opt:     Flag set when beyond range of reference levels
    !       :success_opt:      LOgical indicating a valid target level
    !
    !:Output:
    !
    !       :wgts:             Interpolation coefficients/weights
    !       :kstart:           Index of first relevant referenence/input level for each target level
    !       :kend:             Index of last relevant referenence/input level for each target levbel
    ! 
    implicit none

    ! Arguments:    
    integer, intent(in) :: numInputLevs  ! # of original/input vertical levels
    integer, intent(in) :: numTargetLevs ! # of target vertical levels
    real(8), intent(in) :: pressInput(numInputLevs)  ! pressure on reference column levels assumed to be in ascending order
    real(8), intent(in) :: pressTarget(numTargetLevs) ! target pressure levels assumed to be in ascending order
    character(len=*), intent(in), optional :: method_opt   ! Specified interpolation method
    character(len=*), intent(in), optional :: skipType_opt ! Skipping processing of specific target levels depending on case

    integer, intent(inout), optional :: outbound_opt(numTargetLevs) ! Flag indicating if obs outside model vertical range (0 for no)
    logical, intent(inout), optional :: success_opt(numTargetLevs)  ! success of interpolation
    real(8), intent(out) :: wgts(numTargetLevs,numInputLevs) ! Averaging weights
    integer, intent(out) :: kstart(numTargetLevs) ! Index of first relevant original/input level for each target level
    integer, intent(out) :: kend(numTargetLevs)   ! Index of last relevant original/input level for each target level

    ! Locals:
    integer :: TargetIndex 
    real(8) :: logPressInput(numInputLevs),logPressTarget(numTargetLevs)
    real(8) :: pieceLinearWgts(numTargetLevs,2)
    logical :: success(numTargetLevs)
    character(len=20) :: method,skipType

    if (present(method_opt)) then
      method = method_opt
    else
      method = 'default'
    end if
        
    if (present(skipType_opt)) then
      skipType = skipType_opt
    else
      skipType = 'default'
    end if 
 
    if (present(success_opt)) then
      if ( trim(skipType) == 'doAll&noExtrap' ) then
        success(:)=.true.
      else
        success(:) = success_opt(:)
      end if 
    else
      success(:) = .true.
    end if
    
    ! Flag to prevent extrapolation beyond the range of the reference column levels
        
    if ( trim(skipType) == 'noExtrap' .or. trim(skipType) == 'doAll&noExtrap' ) then

      do TargetIndex=1,numTargetLevs

        ! Check if target level is outside reference column range
        if ( PressTarget(TargetIndex) < PressInput(1) .or. &
             PressTarget(TargetIndex) > PressInput(numInputLevs) ) then
          success(TargetIndex)=.false.
          success_opt(TargetIndex)=.false.
          if (.not.present(outbound_opt)) cycle
          if (PressTarget(TargetIndex) < PressInput(1)) then
            outbound_opt(TargetIndex)=1
          else
            outbound_opt(TargetIndex)=2
          end if
        end if
      end do
    end if
        
    if ( trim(method) == 'default' ) then

      ! Piecewise log-linear interpolation
      
      logPressInput(:) = log(PressInput(:))
      logPressTarget(:) = log(PressTarget(:))

      call ppo_piecewiseLinearWgts(logPressTarget,logPressInput,numTargetLevs,numInputLevs, &
                                 pieceLinearWgts,kstart,success)
                                              
      kend(:) = kstart(:) + 1      
      do TargetIndex=1,numTargetLevs
        if (success(TargetIndex)) then
          kend(TargetIndex) = kstart(TargetIndex) + 1      
          wgts(TargetIndex,kstart(TargetIndex))=pieceLinearWgts(TargetIndex,1)
          wgts(TargetIndex,kend(TargetIndex))=pieceLinearWgts(TargetIndex,2)
        else
          kstart(TargetIndex) = 1
          kend(TargetIndex) = 1
          wgts(TargetIndex,1)=0.0d0
        end if
      end do

    else if ( trim(method) == 'wgtAvg' .and. numTargetLevs > 1) then
    
      ! Piecewise weighted averaging according to distance 
      ! Involves all model levels within the profile range
      
      logPressInput(:) = log(PressInput(:))
      logPressTarget(:) = log(PressTarget(:))
      
      call ppo_layeravgInterpWgts(logPressTarget,logPressInput,numTargetLevs,numInputLevs, &
                                  wgts,kstart,kend,validLevel_opt=success)

    else
      call utl_abort('ppo_vertInterpWgts: This interpolation observation operator is not recognized')
    end if
    
  end subroutine ppo_vertInterpWgts
  
  !------------ routines for interface ppo_linearInterpWgts ----------------
  
  !--------------------------------------------------------------------------
  ! ppo_piecewiseLinearWgts
  !--------------------------------------------------------------------------
  subroutine ppo_piecewiseLinearWgts(pvo,pvi,kno,kni,wgts,kstart,validLevel_opt)
    !
    !:Purpose: To obtain peacewise linear interpolation weigths.
    !          Assumes pv*(i) < pv*(i+1). Constant values extrapolation is applied.
    !
    implicit none

    ! Arguments:
    real(8) ,intent(in)  :: pvi(kni)     ! Vertical levels, pressure (source)
    integer ,intent(in)  :: kni          ! Number of input levels (source)
    integer ,intent(in)  :: kno          ! Number of output levels (destination)
    real(8) ,intent(in)  :: pvo(kno)     ! Vertical levels, pressure (destination)
    real(8) ,intent(out) :: wgts(kno,2)  ! Interpolation weights (destination)
    integer, intent(out) :: kstart(kno)  ! Index i of pvlev level associated to  pvo(j,1)
                                         ! pvo(j,2) is for pvlev level i+1
    logical, intent(out), optional :: validLevel_opt(kno)
    
    ! Locals:
    integer :: pviIndex, pvoIndex, ik
    integer :: lowerlevel(0:KNO)     
    real(8) :: zw1, zp2
    logical :: validLevel(kno)
      
    if (present(validLevel_opt)) then
      validLevel(:) = validLevel_opt(:)
    else
      validLevel(:) = .true.
    end if
    
    ! Find the adjacent level below
    
    lowerlevel(:)=1
    lowerlevel(1:kno) = 0
    do pvoIndex=1,kno
      if (.not.validLevel(pvoIndex)) cycle 
      do pviIndex=max(1,lowerlevel(pvoIndex-1)),kni
        if ( pvo(pvoIndex) < pvi(pviIndex) ) then
          lowerlevel(pvoIndex) = pviIndex
          exit 
        end if
      end do
      if (pviIndex > kni .or. lowerlevel(pvoIndex) == 0 ) then
        lowerlevel(pvoIndex:kno) = kni+1
        exit
      end if
    end do
    
    ! Determine interpolation/extrapolation weights
  
    !$OMP PARALLEL DO PRIVATE(pvoIndex,ik,zp2,zw1)    
    do pvoIndex=1,kno
      if (.not.validLevel(pvoIndex)) then    
        kstart(pvoIndex)=1
        wgts(pvoIndex,1) = 0.0d0
        wgts(pvoIndex,2) = 0.0d0
        cycle
      end if
      
      ik = lowerlevel(pvoIndex)
      if (ik <=1 ) then
        kstart(pvoIndex)=1
        wgts(pvoIndex,1) = 1.0d0
        wgts(pvoIndex,2) = 0.0d0
      else if ( ik > kni ) then
        kstart(pvoIndex)=kni-1
        wgts(pvoIndex,1)=0.0d0
        wgts(pvoIndex,2)=1.0d0
      else
        zp2 = pvi(ik-1)
        zw1 = (pvo(pvoIndex)-zp2)/(pvi(ik)-zp2)
        wgts(pvoIndex,2) =zw1
        wgts(pvoIndex,1) = 1.0d0 - zw1
        kstart(pvoIndex) = ik - 1
      end if
    end do
    !$OMP END PARALLEL DO
           
  end subroutine ppo_piecewiseLinearWgts

  !--------------------------------------------------------------------------
  ! ppo_layeravgInterpWgts
  !--------------------------------------------------------------------------
  subroutine ppo_layeravgInterpWgts(PX1,PX2,KN1,KN2,PZ,kstart,kend,PZS1_opt, &
                                     PZS2_opt,PZDPS_opt,rttov_opt,validLevel_opt)
    !
    !:Purpose: Determine profile interpolation weights by considering integrations 
    !          over of series of  segments [PX1(KI-1),PX1(KI+1)] using  piecewise 
    !          linear weighting with weights of zero at  KI-1 and KI+1 and 
    !          max weight at KI. KI ranges from 1 to KN1.
    !
    !          Can also provide gradient contributions from both 
    !          linear and non-linear components of interpolator. The non-linear 
    !          components (case PZS* is present) stem from vertical coordinate 
    !          independent variables (e.g. dependency on Ps). The gradients 
    !          contributions from the linear components are the interpolation weights.
    !
    !          For the interpolation model f(x) where
    !
    !                f(v,x) = F(v)*x
    !                       ~ F(vo)*xo + F(vo)*(x-xo) + (dF/dv)*xo*(v-vo)
    !                       = F(vo)*x + (dF/dv)*xo*(v-vo)                (eqn 1)
    ! 
    !          
    !                 F(vo):  array of interpolation weights
    !                         = array of gradients from the linear component
    !                 (dF/dv)*xo:  
    !                         array of gradients from linearized component.
    !                         (dF/dv) or (dF/dv)*(v-vo) provided when pzs* is present
    !
    !          Method:
    !             - Piecewise weighted interpolation in ln(P).
    !             - Journal reference:
    !               Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
    !               Jacobian mapping between coordinate systems in data assimilation,
    !               Q. J. of the Royal Met. Soc., vol 133, 1547-1558, 2007.
    !               (www.interscience.wiley.com) DOI:10.1002/qj.117
    !
    !               URL:
    !               http://www3.interscience.wiley.com/cgi-bin/fulltext/116320452/PDFSTART
    !  
    !:Comments:
    !                    
    !     Assumption: PX1(i)<PX1(i+1) & PX2(i)<PX2(i+1)
    !
    !     1) The input profile is now extrapolated at constant value.
    !
    !     The impact is of most practical significance for instruments where 
    !     the weighting function peaks at or near the surface such as SSMI.
    !
    !     This approach increases the weights of contribution from
    !     the lowest and highest input domain levels for output
    !     layers intersecting these input domain boundaries. It consists
    !     of applying contant value extrapolation by introducing  
    !     'fake' or 'virtual' layers. For the lowest level of the input domain, 
    !     as example, this implies creating a virtual surface layer which 
    !     extends to the lower boundary of the output domain layer which 
    !     contains the input domain surface. This increases the
    !     contributing weight of the surface which would be otherwise 
    !     understimated in the original code due to the interpolator 
    !     actually doing piecewise weighted averaging.
    !
    !     2) COmment out use of 'zb' for consistency with RTTOV-9 when 
    !     rttov_opt = .true.
    !     See the four lines ending with !C1 and version 7 comment above.
    !
    !     3) A major reduction in computational time results from only
    !     assigning values to the non-zero ranges of the 2D output
    !     arrays. These ranges of the 2D areas are identified by 'kstart'
    !     and 'kend'. Initialization to zero for values within these ranges 
    !     is done using 1D work arrays 'zpz' and 'zpzd', with the resulting
    !     values then being assigned to the related elements of 'PZ' and 
    !     'PZDPS'.
    !
    !     Therefore, elements of 'PZ' and 'PZDPS' outside these ranges 
    !     could be undefined (i.e. NaN if not 0.0) and should not be used.
    !
    !     Other notable reductions in computational time stem (a) from inlining
    !     of 'sublayer' code (applied to a reduced degree in 'layeravg' as
    !     compared to 'rttov_layeravg_*) and (b) from updating the start 
    !     position 'istart' of the loops over J. The latter was faciliated 
    !     by moving the loop over KI inside 'layeravg'.
    !
    !     These improvements were originally devised and implemented by 
    !     Deborah Salmond and Mats Hamrud (ECMWF) in the RTTOV-9 routines 
    !     'rttov_layeravg*'.
    !
    !     4) Contributors to improvements and changes to the original version 
    !     of 'layeravg' and to 'rttov_layeravg*': members of the RTTOV9 
    !     development team, namely Niels Bormann, Alan Geer, Deborah Salmond,
    !     and Mats Hamrud of ECMWF, Peter Rayer and Roger Saunders of the 
    !     Met Office and Pascal Brunel of Meteo-France, and Y.J. Rochon (EC).
    !                                
    implicit none
     
    !Arguments:
    integer, intent(in)   :: KN1          ! Dimension of PX1
    integer, intent(in)   :: KN2          ! Dimension of other arrays
    real(8), intent(in)   :: PX1(KN1)     ! Levels of output domain (e.g. lnP; in increasing values)
    real(8), intent(in)   :: PX2(KN2)     ! Levels of input domain (e.g. lnP; in increasing values)
    real(8), intent(out)  :: PZ(KN1,KN2)  ! F(vo): Resultant accumulated weighting factors 
                                          ! for interpolation from input to output domain
    integer, intent(out)  :: KSTART(KN1)  ! Start index for relevant PZ row
    integer, intent(out)  :: KEND(KN1)    ! End index of relevant PZ row

    real(8), intent(in), optional   :: PZS1_opt(KN1) ! dPX1/dv or perturbation dPX1/dv * delta(v) where PX1(v), e.g. v=Ps    
    real(8), intent(in), optional   :: PZS2_opt(KN2) ! dPX2/dv or perturbation dPX2/dv * delta(v) where PX2(v), e.g. v=Ps
    real(8), intent(out), optional  :: PZDPS_opt(KN1,KN2) ! dF/dv or perturbations (dF/dv)*(v-vo): 
    !                                                     ! Resultant accumulated factors for the
    !                                                     ! gradients w.r.t. v (or perturbations)
    !                                                     ! associated to coordinates
    logical, intent(in), optional   :: RTTOV_OPT  ! Commented out use of 'zb' when .true. for consistency with RTTOV-9
    logical, intent(in), optional   :: validLevel_opt(kn1) ! Logical indicating validity of each output level
                  
    !Locals:                
    logical  :: LGRADP,lgradp1,lgradp2,validLevel(kn1)
    integer  :: J,IC,ISTART,KI
    real(8)  :: Z1(0:KN1+1),ZW1,ZW2,ZSUM,ZB,ZBPS,ZWPS1,ZWPS2,ZDXD
    real(8)  :: PZ1(0:KN1+1),PZS2(KN2)
    real(8)  :: DZ(KN1+1),DZD(KN1+1),DXD(KN2+1)
    real(8)  :: ZPZ(KN2),ZPZD(KN2)
    real(8)  :: WX1,WX2,Y1,Y2,DY,DAD

    real(8), parameter :: WGT1=0.0d0, WGT2=1.0d0, WGT3=0.0d0

    !- Set integration/averaging range boundaries for output domain.
    !  Range of integration for each layer ki is z1(ki-1) to z1(ki+1).
    !  Weighting function is linear with weights of 0.0 at z1(ki-1) and 
    !  z1(ki+1) and a weight of 1.0 at z1(ki).

    z1(1:kn1)=px1(1:kn1)
    z1(0)=2.0*px1(1)-px1(2)
    z1(kn1+1)=2.0*px1(kn1)-px1(kn1-1)

    if (present(pzs1_opt)) then 
      lgradp1=.true.
      lgradp=.true.     
      pz1(1:kn1)=pzs1_opt(1:kn1)
      pz1(0)=2.0*pzs1_opt(1)-pzs1_opt(2)
      pz1(kn1+1)=2.0*pzs1_opt(kn1)-pzs1_opt(kn1-1)
    else
      lgradp1=.false.
      lgradp=.false.
      pz1(0:kn1+1)=0.0d0
    end if
    if (present(pzs2_opt)) then 
      lgradp2=.true. 
      lgradp=.true.     
      pzs2(:)=pzs2_opt(:)
    else
      lgradp2=.false.
      pzs2(:)=0.0d0
    end if       
    if (present(validLevel_opt)) then
      validLevel(:) = validLevel_opt(:)
    else
      validLevel(:) = .true.
    end if
     
    !- Pre-calculate values (dzd and dxd) used by subroutine sublayer

    dz(1:kn1+1)=z1(1:kn1+1)-z1(0:kn1)
    dzd(1:kn1+1)=1.0d0/dz(1:kn1+1)

    dxd(2:kn2)=1.0d0/(px2(2:kn2)-px2(1:kn2-1))
    dxd(1)=dxd(2)
    dxd(kn2+1)=dxd(kn2)

    !- Determine forward interpolator or TLM coefficients
       
    !- Loop over output domain levels for determining 
    !  contributing weights of input domain levels over 
    !  segments [PX1(KI-1),PX1(KI+1)]

    istart=1
    do ki=1,kn1

      if (.not.validLevel(ki)) then
         kstart(ki) = 1
         kend(ki) = 1
         pz(ki,:) = 0.0d0
         if (lgradp) pzdps_opt(ki,:)=0.0d0
         cycle
      end if         
      
      ! -- Consider constant value extrapolations cases for output domain 
      !    layers entirely outside the input domain.

      if (z1(ki+1).le.px2(1)) then

        ! pz(ki,1) is set to 1.0 to force constant value extrapolation
        ! of field (e.g. temperature) above highest input level.

        pz(ki,1)=1.0d0
        if (lgradp) pzdps_opt(ki,1)=0.0
        kstart(ki)=1
        kend(ki)=1
        CYCLE

      else if (z1(ki-1).ge.px2(kn2)) then

        ! pz(ki:kni,kn2) is set to 1.0 to force constant value extrapolation 
        ! of field (e.g. temperature) below lowest input level.
 
        pz(ki:kn1,kn2)=1.0d0
        if (lgradp) pzdps_opt(ki:kn1,kn2)=0.0
        kstart(ki:kn1)=kn2
        kend(ki:kn1)=kn2
        return
      end if   

      ! -- Consider piecewise averaging interpolation for output domain
      !    layers and layer segments within the input domain.
      !
      ! -- Loop over input layers within the (z1(ki-1),z1(ki)) and 
      !    (z1(ki),z1(ki+1)) integration ranges.
      !
      ! Accumulate contributions to integration components over the 
      !     different segments.

      ic=0
      zpz(istart)=0.0d0
      if (lgradp) then
        zpzd(istart)=0.0d0
        do j=istart,kn2-1
          zpz(j+1)=0.0d0
          zpzd(j+1)=0.0d0

          if (px2(j).ge.z1(ki+1)) exit

          if (px2(j).lt.z1(ki).and.px2(j+1).gt.z1(ki-1)) then

            ! Integration over the segment of the range (z1(ki-1),z1(ki))
            ! intersecting with the range (px2(j),px2(j+1))=(x1,x2)

            call ppo_sublayerInterpWgts(z1(ki-1),z1(ki),dzd(ki),wgt1,wgt2, &
                 px2(j),px2(j+1),dxd(j+1),zpz(j),zpz(j+1), &
                 pz1(ki-1),pz1(ki), &
                 pzs2(j),pzs2(j+1),zpzd(j),zpzd(j+1),lgradp1,lgradp2)
            if (ic.eq.0) ic=j
          endif

          if (px2(j+1).gt.z1(ki)) then

            ! Integration over the segment of the range (z1(ki),z1(ki+1))
            ! intersecting with the range (px2(j),px2(j+1))=(x1,x2).

            call ppo_sublayerInterpWgts(z1(ki),z1(ki+1),dzd(ki+1),wgt2,wgt3, &
                 px2(j),px2(j+1),dxd(j+1),zpz(j),zpz(j+1), &
                 pz1(ki),pz1(ki+1), &
                 pzs2(j),pzs2(j+1),zpzd(j),zpzd(j+1),lgradp1,lgradp2)
            if (ic.eq.0) ic=j
          endif
        enddo

        if (ic.eq.0) then
          j=kn2    
          pz(ki,j)=1.0d0
          pzdps_opt(ki,j)=0.0
          kstart(ki)=j
          kend(ki)=j
          CYCLE
	else
          istart=ic
          kstart(ki)=istart
          kend(ki)=j	
        end if

      else

        ! Same as above but with a compressed subset of 'sublayer' inlined 
        ! for improved speed at least in setting interplation weights.
        ! This follows the corresponding change in 'rttov_layeravg' 
        ! by Deborah Salmond and Mats Hamrud (ECMWF) and takes advantage
        ! of the known wgt* values in each case.
        !
        ! See 'ppo_sublayerInterpWgts' for information on applied equations.

        do j=istart,kn2-1
          zpz(j+1)=0.0d0

          if (px2(j).ge.z1(ki+1)) exit
	  
          if (px2(j).lt.z1(ki).and.px2(j+1).gt.z1(ki-1)) then

            ! Integration over the segment of the range (z1(ki-1),z1(ki))
            ! intersecting with the range (px2(j),px2(j+1))=(x1,x2)

            if (z1(ki-1).lt.px2(j)) then
              y1=px2(j)
              wx1=(y1-z1(ki-1))
            else
              y1=z1(ki-1)
              wx1=0.0d0
            end if
            if (z1(ki).gt.px2(j+1)) then
              wx2=(px2(j+1)-z1(ki-1))
              dy=(px2(j+1)-y1)*dzd(ki)
              zpz(j)=zpz(j)+dy*wx1
              zpz(j+1)=zpz(j+1)+dy*wx2
            else
              dad=dxd(j+1)*dz(ki)*(z1(ki)-px2(j+1))
              dy=(z1(ki)-y1)*dzd(ki)
              zpz(j)=zpz(j)+dy*(wx1-dad)
              zpz(j+1)=zpz(j+1)+dy*(dz(ki)+dad)
            end if

            if (ic.eq.0) ic=j
          endif

          if (px2(j+1).gt.z1(ki)) then

            ! Integration over the segment of the range (z1(ki),z1(ki+1))
            ! intersecting with the range (px2(j),px2(j+1))=(x1,x2).

            if (z1(ki+1).gt.px2(j+1)) then
              y2=px2(j+1)
              wx2=z1(ki+1)-y2
            else
              y2=z1(ki+1)
              wx2=0.0d0
            end if
            if (z1(ki).lt.px2(j)) then
              wx1=z1(ki+1)-px2(j)
              dy=(y2-px2(j))*dzd(ki+1)
              zpz(j)=zpz(j)+dy*wx1
              zpz(j+1)=zpz(j+1)+dy*wx2
            else  
              dad=dxd(j+1)*dz(ki+1)*(z1(ki)-px2(j))
              dy=(y2-z1(ki))*dzd(ki+1)
              zpz(j)=zpz(j)+dy*(dz(ki+1)-dad)
              zpz(j+1)=zpz(j+1)+dy*(wx2+dad)
            end if

            if (ic.eq.0) ic=j
          endif
        enddo

        if (ic.eq.0) then
          j=kn2    
          pz(ki,j)=1.0d0
          kstart(ki)=j
          kend(ki)=j
          CYCLE
	else
          istart=ic
          kstart(ki)=istart
          kend(ki)=j	
        end if

      end if

      ! -- Consider constant value extrapolation contribution for output domain
      !    layers crossing lower or upper input domain boundaries.
      !
      ! -- Extend to below (and/or above) lowest (and/or highest) input levels
      !    if output layer intersects the corresponding input domain boundary.
      !
      !    A virtual layer covering the region between the input domain boundary 
      !    and the boundary of the output layer is added. This increases the
      !    contributing in weight of the input domain boundary roughly by the 
      !    relative thickness of the virtual layer to the entire output layer 
      !    thickness.
      !
      !    This follows updates by Alan Geer (ECMWF) in rttov_layeravg 
      !    of RTTOV-9. It is done for two reasons:
      !
      !       1) Piecewise weighted averaging using only the region
      !          within the input domain will tend to underestimate the impact
      !          of the input level boundary.
      !
      !       2) The regression coefficients of the output domain
      !          model are based on complete output domain layers, this being
      !          consistent with RTTOV convention. 
    
      if (z1(ki+1).gt.px2(kn2).and.z1(ki-1).lt.px2(kn2)) then

        ! Increase contribution from lowest input level.

        zw1=0.0d0
        zw2=0.0d0
        zwps1=0.0d0
        zwps2=0.0d0
        zb=z1(ki+1)
        zbps=pz1(ki+1)
        if (present(rttov_opt)) then
          if (rttov_opt) then
            if (z1(ki+1).gt.2*px2(kn2)-px2(kn2-1).and.ki.eq.kn1) then     !C1
              zb=2*px2(kn2)-px2(kn2-1)                                    !C1
              zbps=2*pzs2(kn2)-pzs2(kn2-1)                                !C1
            end if                                                        !C1
          end if                                                       
        end if
        zdxd=1.0d0/(zb-px2(kn2))
        if (z1(ki).lt.zb) then
          call ppo_sublayerInterpWgts(z1(ki),z1(ki+1),dzd(ki+1),wgt2,wgt3, &
                          px2(kn2),zb,zdxd,zw1,zw2,  &
                          pz1(ki),pz1(ki+1), &
                          pzs2(kn2),zbps,zwps1,zwps2,lgradp1,lgradp2)
        end if
        if (z1(ki).gt.px2(kn2)) then
          call ppo_sublayerInterpWgts(z1(ki-1),z1(ki),dzd(ki),wgt1,wgt2, &
                          px2(kn2),zb,zdxd,zw1,zw2, &
                          pz1(ki-1),pz1(ki), &
                          pzs2(kn2),zbps,zwps1,zwps2,lgradp1,lgradp2)
        end if
        zpz(kn2)=zpz(kn2)+zw1+zw2
        if (lgradp) zpzd(kn2)=zpzd(kn2)+zwps1+zwps2
      end if       
      if (z1(ki-1).lt.px2(1).and.z1(ki+1).gt.px2(1)) then

        ! Increase contribution from highest input level.

        zw1=0.0d0
        zw2=0.0d0
        zwps1=0.0d0
        zwps2=0.0d0
        zb=z1(ki-1)
        zbps=pz1(ki-1)
        if (present(rttov_opt)) then
          if (rttov_opt) then
            if (z1(ki-1).lt.2*px2(1)-px2(2).and.ki.eq.1) then          !C1
              zb=2*px2(1)-px2(2)                                       !C1
              zbps=2*pzs2(1)-pzs2(2)                                   !C1
            end if                                                     !C1
          end if
        end if
        zdxd=1.0d0/(px2(1)-zb)
        if (z1(ki).gt.zb) then
          call ppo_sublayerInterpWgts(z1(ki-1),z1(ki),dzd(ki),wgt1,wgt2, &
                          zb,px2(1),zdxd,zw1,zw2, &
                          pz1(ki-1),pz1(ki), &
                          zbps,pzs2(1),zwps1,zwps2,lgradp1,lgradp2)
        end if
        if (z1(ki).lt.px2(1)) then
          call ppo_sublayerInterpWgts(z1(ki),z1(ki+1),dzd(ki+1),wgt2,wgt3, &
                          zb,px2(1),zdxd,zw1,zw2, &
                          pz1(ki),pz1(ki+1), &
                          zbps,pzs2(1),zwps1,zwps2,lgradp1,lgradp2)
        end if
        zpz(1)=zpz(1)+zw1+zw2
        if (lgradp) zpzd(1)=zpzd(1)+zwps1+zwps2
      end if

      ! -- Normalize sum to unity (instead of calculating and dividing by
      !    weighting denominator)

      zsum=1.0d0/sum(zpz(kstart(ki):kend(ki)))
      pz(ki,kstart(ki):kend(ki))=zpz(kstart(ki):kend(ki))*zsum

      if (lgradp) then

        ! Normalize and account for denominator gradients, i.e.

        ! d[sum1(w*t)/sum2(w)]/dv = sum1[t*(dw/dv)]/sum2(w) 
        !                          -sum1[t*w]*sum2[(dw/dv)]/sum2(w)^2       

        pzdps_opt(ki,kstart(ki):kend(ki))=(zpzd(kstart(ki):kend(ki)) & 
                             -pz(ki,kstart(ki):kend(ki)) &
                             *sum(zpzd(kstart(ki):kend(ki))))*zsum
      end if

    end do

  end subroutine ppo_layeravgInterpWgts

  !--------------------------------------------------------------------------
  ! ppo_sublayerInterpWgts
  !--------------------------------------------------------------------------
  subroutine ppo_sublayerInterpWgts(z1,z2,dzd,wgt1,wgt2,x1,x2,dxd,w1,w2, &
                                     pzs1,pzs2,pxs1,pxs2,wps1,wps2,lgradpx,lgradpz)
    !
    !:Purpose: Determine weight coefficient contributions to w1 and w2 to assign
    !          to input domain (e.g. NWP model) variables at x1 and x2. Weights 
    !          are determined from integration over the intersecting segment (y1,y2) 
    !          of the ranges (x1,x2) for the input domain and (z1,z2) for the 
    !          output domain. Integrals are approximated via the trapezoidal rule:
    !
    !          integral of f(x)=w(x)*t(x) from y1 to y2
    !
    !                                      = (f(y1)+f(y2))/2*(y2-y1)
    !                                      = w(y1)*t(y1)+w(y2)*t(y2)
    !                                      = w1*t(x1)+w2*t(x2)
    !                                      = w1*t1+w2*t2     
    !
    !          This is synonomous to having an integrand linear in x.
    !
    !          In the above (and below) equation(s), w1 and w2 are contributions to 
    !          the input values.
    !
    !          The weights for linearized contributions of non-linear interpolator
    !          components, i.e. gradient w.r.t. the vertical coordinate 
    !          independent variable (e.g. v*=Ps), are calculated 
    !          when LGRADP* is .true.:
    !
    !                pzps = pzps + (df/dx1)*(dx1/dvx1)+(df/dx2)*(dx2/dvx2)
    !                            + (df/dz1)*(dz1/dvz1)+(df/dz2)*(dz2/dvz2)
    !
    !                     = pzpz + (dw1/dx1*t1+dw2/dx1*t2)*pxs1
    !                            + (dw1/dx2*t1+dw2/dx2*t2)*pxs2
    !                            + (dw1/dz1*t1+dw2/dz1*t2)*pzs1
    !                            + (dw1/dz2*t1+dw2/dz2*t2)*pzs2
    !
    !          This routine provides terms on the right-hand-side.
    !
    !          Note: pxs* and pzs* can be provided either as gradients or 
    !          perturbations.
    ! 
    !          Method:
    !          - Piecewise weighted interpolation in ln(P).
    !          - Journal reference:
    !            Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
    !            Jacobian mapping between coordinate systems in data assimilation,
    !            Q. J. of the Royal Met. Soc., vol 133, 1547-1558, 2007. 
    !            (www.interscience.wiley.com) DOI:10.1002/qj.117
    !
    !           URL: 
    !           http://www3.interscience.wiley.com/cgi-bin/fulltext/116320452/PDFSTART 
    !
    !:Comments:
    !
    !     Assumptions:
    !     
    !        1) x1<x2
    !
    !        2) z1<z2
    !
    !        3) The ranges (z1,z2) and (x1,x2) overlap. The overlap region
    !        will be identified as (y1,y2) with y1<y2.
    !
    !     1) w(y1) and w(y2) are obtained by linear interpolation of the linear
    !     weighting function w with w(z1)=wgt1 and w(z2)=wgt2.
    !                                      
    !     2) The w1 and w2 above are determined by expanding t(y1) and t(y2)
    !     as linear functions of t(x1)=t1 and t(x2)=t2.
    !
    !     3) The factor of 1/2 in
    !
    !                          (f(y1)+f(y2))/2*(y2-y1)
    !                           = w(y1)*t(y1)+w(y2)*t(y2)
    !
    !     is omitted as normalization is performed in the calling routine
    !     LAYERAVG.
    !
    !     4) The code version of the interpolator part of 'int_sublayerInterpWgts'
    !     provided for RTTOV-9 assumed the following conditions:
    !
    !        (wgt1,wgt2)=(0,1),  d1=(y1-x1)=0 from y1=x1 or
    !                            wy1=0 from wgt1=0 and y1=z1,
    !     or
    !
    !        (wgt1,wgt2)=(1,0),  d2=(y2-x2)=0 when y2=x2 or
    !                            wy2=0 from wgt2=0 and y2=z2
    !
    !     and took account of the implications on d* and wy*.
    !
    !     The version presented here has each step accompanied by
    !     related equations. It does not assume the above restrictions on  
    !     wgt1 and wgt2. This version is provided in the 
    !     comments section of the RTTOV-9 module 'rttov_sublayer'.
    !
    !     5)  When LGRADP*=.true., this routine provides terms needed for 
    !     the gradients w.r.t.the vertical coordinate independent variable, 
    !     e.g. Ps. 
    !                           
    implicit none

    !Arguments:                           
    logical, intent(in) :: LGRADPz ! Outpout domain logical indicating if gradient w.r.t. vertical coordinate 
                                   ! independent variable if required
				   ! i.e. d/dv where P(v) (e.g. v=Ps).
                                   ! True for yes.

    logical, intent(in) :: LGRADPx ! Input domain logical indicating if gradient w.r.t. vertical coordinate 
                                   ! independent variable if required
				   ! i.e. d/dv where P(v) (e.g. v=Ps).
                                   ! True for yes.
    real(8), intent(in)  :: z1   ! Upper level of output domain half-layer (z1<z2)
    real(8), intent(in)  :: z2   ! Lower level of output domain half-layer
    real(8), intent(in)  :: x1   ! Upper boundary of input layer (x1<x2)
    real(8), intent(in)  :: x2   ! Lower boundary of input layer
    real(8), intent(in)  :: wgt1 ! Weight at z1 (0.0 or 1.0 or ...)
    real(8), intent(in)  :: wgt2 ! Weight at z2 (1.0 or 0.0 or ...)
    real(8), intent(in)  :: pzs1 ! dz1/dvz or perturbation dz1/dvz * delta(vz)
    real(8), intent(in)  :: pzs2 ! dz2/dvz or perturbation dz2/dvz * delta(vz)
    real(8), intent(in)  :: pxs1 ! dx1/dvx or perturbation dx1/dvx * delta(vx)
                                 ! (required when LGRADPx=.true.)
    real(8), intent(in)  :: pxs2 ! dx2/dvx or perturbation dx2/dvx * delta(vx)
    real(8), intent(in)  :: dxd  ! 1.0/(x2-x1)=1.0/dx
    real(8), intent(in)  :: dzd  ! 1.0/(z2-z1)=1.0/dz
    real(8), intent(inout) :: w1 ! Starting (in) and updated (out) weight 
    !                            ! assigned to variable at upper level x1
    real(8), intent(inout) :: w2 ! Starting (in) and updated (out) weight 
    !                            ! assigned to variable at upper level x2
    real(8), intent(inout) :: wps1 ! Starting (in) and updated (out) value of
    !                              ! (pxs1*dw1/dx1 + pxs2*dw1/dx2
    !                              ! +pzs1*dw1/dz1 + pzs2*dw1/dz2)
    real(8), intent(inout) :: wps2 ! Starting (in) and updated (out) value of
    !                              ! pxs1*dw2/dx1 + pxs2*dw2/dx2
    !                              ! +pzs1*dw2/dz1 + pzs2*dw2/dz2)
    !                              ! (required when LGRADP*=.true.)
     
    !Locals:           
    real(8)  :: y1  ! Upper boundary of integral range (y1<y2)
    real(8)  :: y2  ! Lower boundary of integral range
    real (8) :: d1,d2,wx1,wx2,wy1,wy2,dy,dzdd
    real(8)  :: a1,a2,dxy1,dxy2,dxyd1,dxyd2,c1,c2
    real(8)  :: dydzdd,dy1z1,dy2z2,d1d1,d1d2
    integer  :: ibot,itop
    real(8) :: zthreshold

    !  1. Initialization

    !- Set upper and lower boundaries of integration layer (y1,y2):
    !        (y1,y2) = ( max(x1,z1), min(x2,z2) )
    itop=0
    ibot=0
    y1=z1
    y2=z2
    if (y1.lt.x1) then
      y1=x1
      itop=1
    end if
    if (y2.gt.x2) then
      y2=x2
      ibot=1
    end if
    dy=y2-y1

    !- Verify for negative and zero (and near-zero) y2-y1 values.
     
    zthreshold=epsilon(1.0d0)*100.0
    c1=abs(y2)*zthreshold

    if (dy.lt.-c1) then
      write(*,*) 'y1,y2 = ',y1,y2 
      write(*,*) 'z1,z2 = ',z1,z2 
      write(*,*) 'x1,x2 = ',x1,x2
      call utl_abort("ppo_sublayerInterpWgts: dy is negative")
    else if (dy.lt.c1) then
      ! dy~0; weights can be taken as having values of zero.
      ! w1 and w2 unchanged.
      return
    end if
    
    !  2. Interpolation weights (equivalent to gradient contributions from
    !     linear component of interpolator)

    !- Determine w(y1) and w(y2) of the integral f = w(y1)*t(y1)+w(y2)*t(y2)
    !                                              = wy1*t(y1)+wy2*t(y2)

    !  Weights are obtained by linear interpolation of the linear
    !  weighting function w with w(z1)=wgt1 and w(z2)=wgt2.

    dzdd=dzd*(wgt2-wgt1)
    dy1z1=y1-z1
    wx1=wgt1+dy1z1*dzdd
     
    dydzdd=dy*dzdd
    wx2=wx1+dydzdd    !  wx2=wgt1+(y2-z1)*dzdd=wgt2+(y2-z2)*dzdd
     
    wy1=dy*wx1
    wy2=dy*wx2

   !- Determine contribution of w1 and w2 for f=w1*t(x1)+w2*t(x2). 

   !  The w1 and w2 contributions above are determined by expanding t(y1)
   !  and t(y2) in f = w(y1)*t(y1)+w(y2)*t(y2) as linear functions of 
   !  t(x1)=t1 and t(x2)=t2.
   !
   !        t(y1) = t1+(y1-x1)*(t2-t1)/(x2-x1) = t1+(y1-x1)*dxd*(t2-t1)
   !        t(y2) = t2+(y2-x2)*(t2-t1)/(x2-x1) = t2+(y2-x2)*dxd*(t2-t1)
   !
   !  Therefore,
   !
   !        f = wy1*[t1+(y1-x1)*dxd*(t2-t1)]
   !           +wy2*[t2+(y2-x2)*dxd*(t2-t1)]
   !
   !          = t1*[wy1-wy1*(y1-x1)*dxd-wy2*(y2-x2)*dxd]
   !           +t2*[wy1*(y1-x1)*dxd+wy2+wy2*(y2-x2)*dxd]
   !
   !  Aside: Contribution to w1+w2 = wy1+wy2 = 2 * w[(y2+y1)/2]
 
    dxy1=y1-x1
    d1=dxy1*dxd

    dxy2=y2-x2
    d2=dxy2*dxd
    
    d1d1=1.0d0-d1
    d1d2=1.0d0+d2  

    w1=w1+wy1*d1d1-wy2*d2
    w2=w2+wy1*d1+wy2*d1d2  

    !  Aside used below: wy1*d1d1-wy2*d2 = wy1 + wy2 - (wy1*d1d1-wy2*d2)

    !  3. Provide gradient contributions from linearized terms of non-linear
    !     component of interpolator (i.e. vertical coordinate dependency)

    if (LGRADPx) then
     
      !  -- 3.1 Provide linearized contribution of non-linear interpolator
      !     component to TLM - this part in relation to the input vertical  
      !     coordinate:
      !
      !                    (df/dx1)*(dx1/dv)+(df/dx2)*(dx2/dv)
      !                        =  (dw1/dx1*t1+dw2/dx1*t2)*pxs1
      !                         + (dw1/dx2*t1+dw2/dx2*t2)*pxs2
      !
      !     a1 and a2 are used below to denote dw*/dx1 and dw*/dx2.
      !
      !     Gradients of w1 and w2 with respect to x1 and x2 are done in two parts:
      !
      !          - Part I:  The gradients calc excludes the cases where one 
      !                     might have y1=x1 and or y2=x2.
      !          - Part II: The gradient terms related to y1=x1 and or y2=x2 are
      !                     added.
      !
      !  -- Gradient of f w.r.t x1
      !    
      !  -- Gradients of w1 and w2 w.r.t. x1:
      !
      !     PART I:
      !
      !         For this case:
      !
      !            d(d1)/dx1 = dxd*(-1+dxy1*dxd) 
      !            d(d2)/dx1 = dxd*dxy2*dxd
      !
      !            d(wy1)/dx1=0, d(wy2)/dx1=0
      !
      !         Note: d(dj)/dx* = 0.0 when dj=0.0 due to y1=x1 or y2=x2
      
      dxyd1=d1*dxd
      dxyd2=d2*dxd                                
      a1=wy1*(dxd-dxyd1)-wy2*dxyd2
      a2=-a1        ! a2=wy1*(-dxd+dxyd1)+wy2*dxyd2

      !     PART II:
      !
      !        For the extra terms: (gradients w.r.t. y1=x1)
      !
      !            d(d1)/dy1 = dxd
      !            d(d2)/dy1 = 0
      !
      !            d(wy1)/dy1 = -wx1+dzdd*dy
      !            d(wy2)/dy1 = -wx2
      
      if (itop.eq.1) then

        ! Case y1=x1 (d1=0)

        c1=wy1*dxd
        c2=-wx1+dydzdd
        a1=a1-c1+c2+wx2*d2
        a2=a2+c1-wx2*d1d2
      end if

      !  -- Add to accumulated gradient terms

      wps1=wps1+a1*pxs1
      wps2=wps2+a2*pxs1

      !  -- Gradient of f w.r.t. x2
      !
      !  -- Gradients of w1 and w2 w.r.t. x2:
      !
      !     PART I:
      !
      !         For this case:
      !
      !            d(d1)/dx2 = -dxd*dxy1*dxd
      !            d(d2)/dx2 = -dxd*(1+dxy2*dxd)
      !
      !            d(wy1)/dx2=0,  d(wy2)/dx2=0
      
      a1=wy1*dxyd1+wy2*(dxd+dxyd2)
      a2=-a1         ! a2=-wy1*dxyd1-wy2*(dxd+dxyd2)
                    
      !     PART II:
      !     
      !        For the extra terms: (gradients w.r.t. y2=x2)
      !
      !            d(d1)/dy2 = 0
      !            d(d2)/dy2 = dxd
      !            
      !            d(wy1)/dy2 = wx1
      !            d(wy2)/dy2 = wx2+dzdd*dy
                
      if (ibot.eq.1) then

        ! Case y2=x2 (d2=0)

        c1=wy2*dxd
        c2=wx2+dydzdd
        a1=a1-c1+wx1*d1d1
        a2=a2+c1+c2+wx1*d1
      end if

      !  -- Add to accumulated gradient terms

      wps1=wps1+a1*pxs2
      wps2=wps2+a2*pxs2

    end if

    if (LGRADPz) then
     
      !  -- 3.2 Provide linearized contribution of non-linear interpolator
      !     component to the TLM - this part in relation to the output vertical 
      !     coordinate:
      !
      !                    (df/dz1)*(dz1/dv)+(df/dz2)*(dz2/dv)
      !
      !                        =  (dw1/dz1*t1+dw2/dz1*t2)*pzs1
      !                         + (dw1/dz2*t1+dw2/dz2*t2)*pzs2
      !
      !     a1 and a2 are used below to denote dw*/dz1 and dw*/dz2.
      !
      !     Gradients of w1 and w2 with respect to z1 and z2 are done in two parts:
      !
      !          - Part I:  The gradients calc excludes the cases where one 
      !                     might have y1=z1 and or y2=z2.
      !          - Part II: The gradient terms related to y1=z1 and or y2=z2 are
      !                     added.
      !
      !  -- Gradient of f w.r.t z1
      !
      !  -- Gradients of w1 and w2 w.r.t. z1:
      !
      !     PART I: Excludes y1=z1 consideration
      !
      !        d(wy1)/dz1 = dy*dzdd*((y1-z1)*dzd-1)=dy*dzdd*(y1-z2)*dzd
      !        d(wy2)/dz1 = dy*dzdd*(y2-z2)*dzd
      
      dy2z2=y2-z2
      dxyd2=dydzdd*dzd
      dxyd1=dxyd2*dy1z1
      dxyd2=dxyd2*dy2z2
      a1=(dxyd1-dydzdd)*d1d1-dxyd2*d2
      a2=dxyd1+dxyd2-(dydzdd+a1) 

      if (itop.eq.0) then 

        ! PART II: y1=z1
        !
        !          d(d1)/dy1 = dxd
        !          d(d2)/dy1 = 0
        !
        !          d(wy1)/dy1 = -wx1+dzdd*dy 
        !          d(wy2)/dy1 = -wx2

        c1=wy1*dxd
        c2=-wx1+dydzdd
        a1=a1-c1+c2*d1d1+wx2*d2
        a2=a2+c1+c2*d1-wx2*d1d2
      end if

      !  -- Add to accumulated gradient terms

      wps1=wps1+a1*pzs1
      wps2=wps2+a2*pzs1

      !  -- Gradient of f w.r.t. z2
      !
      !  -- Gradients of w1 and w2 w.r.t. z2:
      !
      !     PART I: Excludes y2=z2 consideration
      !
      !         d(wy1)/dz2 = -dy*dzdd*(y1-z1)*dzd
      !         d(wy2)/dz2 = -dy*dzdd*(1+(y2-z2)*dzd)=-dy*dzdd*(y2-z1)*dzd
      
      a1=-dxyd1*d1d1+(dxyd2+dydzdd)*d2
      a2=-(dxyd1+dxyd2+dydzdd+a1)

      if (ibot.eq.0) then

        ! PART II: y2=z2
        !     
        !          d(d1)/dy2 = 0
        !          d(d2)/dy2 = dxd
        !            
        !          d(wy1)/dy2 = wx1
        !          d(wy2)/dy2 = wx2+dzdd*dy
        
        c1=wy2*dxd
        c2=wx2+dydzdd
        a1=a1-c1-c2*d2+wx1*d1d1
        a2=a2+c1+c2*d1d2+wx1*d1 
      end if

      !  -- Add to accumulated gradient terms

      wps1=wps1+a1*pzs2
      wps2=wps2+a2*pzs2

    end if

  end subroutine ppo_sublayerInterpWgts

  !==========================================================================
  !--- Stand-alone integration (and averaging) routines providing weights ---
  ! For weights W(:,:) and initial vector X(:), the integrated (or average) 
  ! values would be W*X.  

  !--------------------------------------------------------------------------
  ! ppo_vertLayersSetup
  !--------------------------------------------------------------------------
  subroutine ppo_vertLayersSetup(operatorType,pressInput,numInputLevs)
    !
    !:Purpose: Preliminary calculations for producing components required for
    !          vertical integration (or averaging) w.r.t. pressure for the full  
    !          vertical rangeor a set of target layers. To be called before 
    !          routine ppo_vertIntegWgts or ppo_vertAvgWgts.
    !
    !          Integration calculations are performed appling quadratic interpolation 
    !          in P between level. 
    !
    !:Input:
    !         :operatorType:       'integ' for intergration; 'avg' for averaging
    !         :pressInput:         Reference input levels
    !         :numInputLevs:       # of model vertical levels
    !
    !:Output:
    !         :boundaries(numInputLevs+1):  Input layer boundaries assuming
    !                                       provided input levels can be taken as
    !                                       mid-layer values.
    !         :weights:                     Second order Lagrangian
    !                                       interp integration weights or 
    !                                       unity for averaging weights
    !
    !:Comments:
    !         - This subroutine does the following:
    !
    !           - Setting of layer boundaries
    !           - If integration, determining integration weights associated 
    !             to second order Lagrangian interpolation. Otherwise, initialize
    !             weights to unity.
    !
    !         - Layer boundaries are taken as mid-point between provided levels in
    !           lnP coordinate. Layer values are set to be the values
    !           interpolated to the mid-point in P within the various layers.
    !           Interpolation in P is done quadratically. 
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: operatorType       ! 'integ' for integration; 'avg' for averaging
    integer, intent(in)    :: numInputLevs             ! # of model vertical levels
    real(8), intent(in)    :: pressInput(numInputLevs) ! Reference input levels

    ! Locals:

    integer   :: levelIndex
    real(8)   :: zp, zp1, zp2, zp3, zr1, zr2, zr3

    ! Determine P boundaries of layers and save weights for
    ! use in setting integration (or averaging) weights.
    ! N.B.: Boundaries of layers set to mid-point between input levels
      
    ! Calculate layer boundaries

    if (allocated(boundaries)) deallocate(boundaries)
    allocate(boundaries(numInputLevs+1))
    
    boundaries(1)=pressInput(1)
    boundaries(numInputLevs+1)= pressInput(numInputLevs)

    !$OMP PARALLEL DO PRIVATE(levelIndex)    
    DO levelIndex = 2, numInputLevs
       boundaries(levelIndex)=sqrt(pressInput(levelIndex-1)*pressInput(levelIndex))
    END DO
    !$OMP END PARALLEL DO 

    if (allocated(weights)) deallocate(weights)
    
    if ( trim(operatorType) == 'avg' ) then
    
      ! Initialize weights as scaling factors of unity.

      allocate(weights(numInputLevs,1))      
      weights(:,:) = 1.0d0
    
    else

      ! Set second degree Lagrangian interpolator weights

      allocate(weights(numInputLevs,numInputLevs))            
      weights(:,:) = 0.0d0
    
      ! Interpolation to mid-layer level in P using
      ! second degree Lagrangian interpolator.
      ! N.B.: Integration is w.r.t. P
    
      ! Calculating for levelIndex=1
    
      zp1= pressInput(1)
      zp2= pressInput(2)
      zp3= pressInput(3)
      zp = (boundaries(2)+boundaries(1))/2.0
      zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
      zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
      zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
      weights(1,1)=zr1
      weights(2,1)=zr2
      weights(3,1)=zr3

      !$OMP PARALLEL DO PRIVATE(levelIndex,zp1,zp2,zp3,zp,zr1,zr2,zr3)    
      DO levelIndex=2,numInputLevs-1
        zp1=pressInput(levelIndex-1)
        zp2=pressInput(levelIndex)
        zp3=pressInput(levelIndex+1)
        zp=(boundaries(levelIndex+1)+boundaries(levelIndex))/2.0
        zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
        zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
        zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
        weights(levelIndex-1,levelIndex)=zr1
        weights(levelIndex,levelIndex)=zr2
        weights(levelIndex+1,levelIndex)=zr3
      ENDDO
      !$OMP END PARALLEL DO 
    
      ! Calculating  for levelIndex=numInputLevs
    
      zp1= pressInput(numInputLevs-2)
      zp2= pressInput(numInputLevs-1)
      zp3= pressInput(numInputLevs)
      zp = (boundaries(numInputLevs+1)+boundaries(numInputLevs))/2.0
      zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
      zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
      zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
      weights(numInputLevs-2,numInputLevs)=zr1
      weights(numInputLevs-1,numInputLevs)=zr2
      weights(numInputLevs,numInputLevs)=zr3

    end if

  end subroutine ppo_vertLayersSetup

  !--------------------------------------------------------------------------
  ! ppo_vertIntegWgts
  !--------------------------------------------------------------------------
  subroutine ppo_vertIntegWgts(targetLayersTop,targetLayersBot,numInputLevs, &
                                numTargetLevs,kstart,kend,wgts,wgts_opt,skipType_opt, &
                                outbound_opt,success_opt,dealloc_opt)   
    !
    !:Purpose: To calculate integration weights "wgts" required for vertical integration w.r.t.
    !          pressure for the full vertical range or a set of target layers. 
    !          Given the calculated weights and a user intergrand array vector X, the integral 
    !          for a given layer i would be given by sum(wgts(i,:)*X(:))
    !
    !          Integration calculations are performed applying quadratic interpolation 
    !          in P between level.
    !
    !          Routine ppo_vertLayersSetup to be called beforehand to generate Lagrangian weights
    !          and related layer boundaries (arrays 'weights' and 'boundaries')
    ! 
    !:Input:
    !         :targetLayersTop:    top of target layers
    !         :targetLayersBot:    bottom of target layers
    !         :numInputLevs:       # of original/input vertical levels
    !         :numTargetLevs:      # of target vertical levels
    !         :kstart:             Index of first relevant original/input level for each target level
    !         :kend:               Index of last relevant original/input level for each target level
    !                              If kstart and kend are non-zero on input, 
    !                              the input are initial estimates of the values.
    !         :weights:            See routine ppo_vertLayersSetup
    !         :boundaries:         Boundaries of input layers
    !         :skipType_opt:       Skipping processing of specific target layers depending on case:
    !                              'default' - skipping application via input success_opt only
    !                              'doAll&noExtrap' - application of both success_opt and outbound_opt
    !         :outbound_opt:       Flag set when beyond range of reference/input levels
    !         :success_opt:        Logical indicating a valid target layer
    !         :dealloc_opt:        Logical indicating if deallocation is desired when done. (default: .false.)
    !
    !:Output:
    !         :wgts(numTargetLevs,numInputLevs):     Integration weights
    !         :wgts_opt(numTargetLevs,numInputLevs): Part of integrtation weights not related to resolution
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: numInputLevs  ! # of original/input vertical levels
    integer, intent(in) :: numTargetLevs ! # of target vertical levels
    real(8), intent(in) :: targetLayersTop(numTargetLevs) ! top of target layers
    real(8), intent(in) :: targetLayersBot(numTargetLevs) ! bottom of target layers

    real(8), intent(out) :: wgts(numTargetLevs,numInputLevs) ! Averaging weights
    integer, intent(inout) :: kstart(numTargetLevs) ! Index of first relevant original/input level for each target level
    integer, intent(inout) :: kend(numTargetLevs)   ! Index of last relevant original/input level for each target level

    character(len=*), intent(in), optional :: skipType_opt ! Skipping processing of specific target layers depending on case
    logical, intent(in), optional  :: dealloc_opt ! Logical indicating if deallocation is desired when done
    real(8), intent(out), optional :: wgts_opt(numTargetLevs,numInputLevs) ! Part of averaging weights not related to resolution
    integer, intent(inout), optional :: outbound_opt(numTargetLevs) ! Flag indicating if obs outside input vertical range (0 for no)
    logical, intent(inout), optional :: success_opt(numTargetLevs)  ! success of interpolation

    ! Locals:
    integer :: TargetIndex   
    logical :: success(numTargetLevs)
    character(len=20) :: skipType
    integer, parameter :: ivweights=2  ! Order of Lagrangian interpolation.

    integer :: levelIndex,JK,ILMAX2,ILMIN2
    integer :: ILMIN, ILMAX
    real(8) :: zp, zp1, zp2, zp3, zr1, zr2, zr3, ptop, pbtm
    
    if (present(skipType_opt)) then
      skipType = skipType_opt
    else
      skipType = 'default'
    end if 
 
    if (present(success_opt)) then
      if ( trim(skipType) == 'doAll&noExtrap') then
        if (present(outbound_opt)) then
          where (outbound_opt(:) == 0)
            success(:) = .true.
          elsewhere
            success(:) = .false.
          end where
        else
          success(:) = .true.
        end if
      else
        success(:) = success_opt(:)
      end if 
    else
      success(:) = .true.
    end if
        
    do TargetIndex=1,numTargetLevs

      if ( .not.success(TargetIndex) ) then
        wgts(TargetIndex,:) = 0.0D0
        wgts_opt(TargetIndex,:) = 0.0D0          
        kstart(TargetIndex)=1
        kend(TargetIndex)=1
        cycle
      end if

      ptop = targetLayersTop(TargetIndex)
      pbtm = targetLayersBot(TargetIndex)
         
      ! Find the range of vertical levels over which to perform the integration
      ! and set the integration weights over this range.
          
      ilmin=1
      ilmax=numInputLevs
      if (ptop <= boundaries(1)*1.01 .and. &
          pbtm >= boundaries(numInputLevs+1)*0.99) then

        ! Total column integration part

        !$OMP PARALLEL DO PRIVATE(jk,levelIndex)                  
        do jk = 1,numInputLevs
          do levelIndex=max(1,jk-ivweights),min(numInputLevs,jk+ivweights)
            wgts(TargetIndex,jk)=wgts(TargetIndex,jk) &
	          +(boundaries(levelIndex+1) &
                  -boundaries(levelIndex))*weights(jk,levelIndex)
            if (present(wgts_opt)) &
              wgts_opt(TargetIndex,jk)=wgts_opt(TargetIndex,jk)+ &
	                               weights(jk,levelIndex)
          end do
        end do
        !$OMP END PARALLEL DO                 
             
      else

        ! Partial column integration part (special treatment at boundaries)
     
        ! Identify input layer boundaries just within the target layer.
             
        ilmin = ppo_getLevelIndex(ptop, boundaries, 'top', numInputLevs+1)
        ilmax = ppo_getLevelIndex(pbtm, boundaries, 'btm', numInputLevs+1)
               
        if (ilmin == ilmax+1) then

          ! Entire target layer within one input layer
                
          levelIndex=ilmax
          if (levelIndex < 3) levelIndex=3
          if (levelIndex > numInputLevs) levelIndex=numInputLevs
          zp1=boundaries(levelIndex-2)
          zp2=boundaries(levelIndex-1)
          zp3=boundaries(levelIndex)
          zp=(ptop+pbtm)/2.0
          zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
          zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
          zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                
          wgts(TargetIndex,levelIndex-2)=(pbtm-ptop)*zr1
          wgts(TargetIndex,levelIndex-1)=(pbtm-ptop)*zr2
          wgts(TargetIndex,levelIndex)=(pbtm-ptop)*zr3
          if (present(wgts_opt)) then
            wgts_opt(TargetIndex,levelIndex-2)=zr1
            wgts_opt(TargetIndex,levelIndex-1)=zr2
            wgts_opt(TargetIndex,levelIndex)=zr3
          end if
          ilmin=levelIndex-2
          ilmax=levelIndex
                  
        else
                
          ! Determine terms from the inner layers (excluding the lower and upper
          ! boundary layers when these layers not covering entire input layers)
                
          if (pbtm >= boundaries(numInputLevs)*0.99) then
            ilmax2=numInputLevs
          else
            ilmax2=ilmax-1
          end if
          if (ptop <= boundaries(1)*1.01) then
            ilmin=1
            ilmin2=ilmin
          else
            ilmin2=ilmin
          end if
          if (ilmin2 <= ilmax2) then
            !$OMP PARALLEL DO PRIVATE(jk,levelIndex)                  
            do jk = ilmin2,ilmax2
              do levelIndex=max(1,jk-ivweights),min(numInputLevs,jk+ivweights)
                wgts(TargetIndex,jk)=wgts(TargetIndex,jk)+(boundaries(levelIndex+1) &
                         -boundaries(levelIndex))*weights(jk,levelIndex)
                if (present(wgts_opt)) &
                  wgts_opt(TargetIndex,jk)=wgts_opt(TargetIndex,jk)+weights(jk,levelIndex)
              end do
            end do
            !$OMP END PARALLEL DO               
          end if
                
          ! Determine terms from the lower and upper boundary layers
          ! when these layers do not cover entire input layers.
                
          if (pbtm < boundaries(numInputLevs)*0.99) then
                     
            levelIndex=ilmax+1
            if (levelIndex > numInputLevs) levelIndex=numInputLevs
            if (levelIndex < 3) levelIndex=3
            zp1=boundaries(levelIndex-2)
            zp2=boundaries(levelIndex-1)
            zp3=boundaries(levelIndex)
            zp=(boundaries(ilmax)+pbtm)/2.0
            zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
            zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
            zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
            wgts(TargetIndex,levelIndex-2)=wgts(TargetIndex,levelIndex-2)+(pbtm - boundaries(ilmax))*zr1
            wgts(TargetIndex,levelIndex-1)=wgts(TargetIndex,levelIndex-1)+(pbtm - boundaries(ilmax))*zr2
            wgts(TargetIndex,levelIndex)=wgts(TargetIndex,levelIndex)+(pbtm - boundaries(ilmax))*zr3
           
            if (present(wgts_opt)) then
              wgts_opt(TargetIndex,levelIndex-2)=wgts_opt(TargetIndex,levelIndex-2)+zr1
              wgts_opt(TargetIndex,levelIndex-1)=wgts_opt(TargetIndex,levelIndex-1)+zr2
              wgts_opt(TargetIndex,levelIndex)=wgts_opt(TargetIndex,levelIndex)+zr3
            end if
            ilmax=levelIndex
                  
          end if
                  
          if (ptop > boundaries(1)*1.01) then
                     
            levelIndex=ilmin-1
            if (levelIndex < 1) levelIndex=1
            if (levelIndex > numInputLevs-2) levelIndex=numInputLevs-2
            zp1= boundaries(levelIndex)
            zp2= boundaries(levelIndex+1)
            zp3= boundaries(levelIndex+2)
            zp = (boundaries(ilmin)+ptop)/2.0
            zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
            zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
            zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
            wgts(TargetIndex,levelIndex)=wgts(TargetIndex,levelIndex)+(boundaries(ilmin)-ptop)*zr1
            wgts(TargetIndex,levelIndex+1)=wgts(TargetIndex,levelIndex+1)+(boundaries(ilmin)-ptop)*zr2
            wgts(TargetIndex,levelIndex+2)=wgts(TargetIndex,levelIndex+2)+(boundaries(ilmin)-ptop)*zr3
                   
            if (present(wgts_opt)) then
              wgts_opt(TargetIndex,levelIndex)=wgts_opt(TargetIndex,levelIndex)+zr1
              wgts_opt(TargetIndex,levelIndex+1)=wgts_opt(TargetIndex,levelIndex+1)+zr2
              wgts_opt(TargetIndex,levelIndex+2)=wgts_opt(TargetIndex,levelIndex+2)+zr3
 	    end if
            ilmin=levelIndex
            if (ilmax < levelIndex+2) ilmax=levelIndex+2
                   
          end if
          if (ilmin > ilmax-2) ilmin=ilmax-2
        end if
      end if

      if (kstart(TargetIndex) > 0 .and. kend(TargetIndex) > 0) then
        if (abs(kstart(TargetIndex)-ilmin) > 1 .or. &
	    abs(kend(TargetIndex)-ilmax) > 1) then
	    
          write(*,*) 'ppo_vertIntegWgts: Suspected error in layer', &
	        ' identification: ',TargetIndex,kstart(TargetIndex),ilmin, &
	        kend(TargetIndex),ilmax
        end if
      end if

      kstart(TargetIndex)=ilmin
      kend(TargetIndex)=ilmax
       
    end do

    if (present(dealloc_opt)) then
      if (dealloc_opt) deallocate(weights,boundaries)
    end if
    
  end subroutine ppo_vertIntegWgts

  !--------------------------------------------------------------------------
  ! ppo_getLevelIndex
  !--------------------------------------------------------------------------
  integer function ppo_getLevelIndex(level, layerBoundaryLevels, topbtm, numBoundaries)
    !
    !:Purpose: To get the vertical input level index for level
    !          within target layer and nearest specified layer boundary.  
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: numBoundaries  ! Number of layer boundaries       
    real(8), intent(in) :: level          ! Target layer index      
    real(8), intent(in) :: layerBoundaryLevels(numBoundaries) ! Layer boundaries
    character(len=*), intent(in) :: topbtm      ! indicating whether we are looking for top or bottom level

    ! Locals
    integer     :: ilev1, ilev2
    integer     :: levelIndex

    ! Find the model levels adjacent to pressure level 

    ! Default values
    
    if (level < 0.0d0) then
      if ((topbtm == 'btm') .or. (topbtm == 'BTM')) then
        ppo_getLevelIndex = numBoundaries
      endif
      if ((topbtm == 'top') .or. (topbtm == 'TOP')) then
        ppo_getLevelIndex = 1
      endif                                                  
    endif
      
    ilev1=0
    ilev2=1
    do levelIndex=1,numBoundaries
      if (level > layerBoundaryLevels(levelIndex)) then
        ilev1=levelIndex
        ilev2=levelIndex+1
      else
        exit
      endif
    enddo

    ! Find the input level index

    ! If we are looking for top level, the index is the level immediately 
    ! below. if looking for bottom level, the index is the one immediately 
    ! above.
    
    if ((topbtm == 'btm') .or. (topbtm == 'BTM')) then
      ppo_getLevelIndex=ilev1
    else if ((topbtm == 'top') .or. (topbtm == 'TOP')) then
      ppo_getLevelIndex=ilev2
    endif

    if (ppo_getLevelIndex < 1) ppo_getLevelIndex=1
    if (ppo_getLevelIndex > numBoundaries) ppo_getLevelIndex=numBoundaries
  
  end function ppo_getLevelIndex

  !--------------------------------------------------------------------------
  ! ppo_vertAvgWgts
  !--------------------------------------------------------------------------
  subroutine ppo_vertAvgWgts(targetLayersTop,targetLayersBot,numInputLevs, &
                             numTargetLevs,kstart,kend,wgts,wgts_opt,skipType_opt, &
                             outbound_opt,success_opt,dealloc_opt)   
    !
    !:Purpose: To calculate averaging weights "wgts" required for vertical averaging 
    !          w.r.t. ln(pressure) for the full vertical range or a set of target layers. 
    !          Given the calculated weights and a user input array vector X, the average 
    !          for a given layer i would be given by sum(wgts(i,:)*X(:))
    !
    !          Routine ppo_vertLayersSetup to be called beforehand to initial weigths
    !          and related layer boundaries (arrays 'weights' and 'boundaries')
    ! 
    !:Input:
    !         :targetLayersTop:    top of target layers
    !         :targetLayersBot:    bottom of target layers
    !         :numInputLevs:       # of original/input vertical levels
    !         :numTargetLevs:      # of target vertical levels
    !         :kstart:             Index of first relevant original/input level for each target level
    !         :kend:               Index of last relevant original/input level for each target level
    !                              If kstart and kend are non-zero on input, 
    !                              the input are initial estimates of the values.
    ! 
    !         :weights:            See routine ppo_vertLayersSetup
    !         :boundaries:         Boundaries of input layers
    !         :skipType_opt:       Skipping processing of specific target layers depending on case:
    !                              'default' - skipping application via input success_opt only
    !                              'doAll&noExtrap' - application of both success_opt and outbound_opt
    !         :outbound_opt:       Flag set when beyond range of reference/input levels
    !         :success_opt:        Logical indicating a valid target layer
    !         :dealloc_opt:        Logical indicating if deallocation is desired when done. (default: .false.)
    !
    !:Output:
    !         :wgts(numTargetLevs,numInputLevs):     Averaging weights
    !         :wgts_opt(numTargetLevs,numInputLevs): Part of averaging weights not related to resolution
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: numInputLevs  ! # of original/input vertical levels
    integer, intent(in) :: numTargetLevs ! # of target vertical levels
    real(8), intent(in) :: targetLayersTop(numTargetLevs) ! top of target layers
    real(8), intent(in) :: targetLayersBot(numTargetLevs) ! bottom of target layers

    real(8), intent(out)   :: wgts(numTargetLevs,numInputLevs) ! Averaging weights
    integer, intent(inout) :: kstart(numTargetLevs) ! Index of first relevant original/input level for each target level
    integer, intent(inout) :: kend(numTargetLevs)   ! Index of last relevant original/input level for each target level

    character(len=*), intent(in), optional :: skipType_opt ! Skipping processing of specific target layers depending on case
    logical, intent(in), optional  :: dealloc_opt ! Logical indicating if deallocation is desired when done
    real(8), intent(out), optional :: wgts_opt(numTargetLevs,numInputLevs) ! Part of averaging weights not related to resolution
    integer, intent(inout), optional :: outbound_opt(numTargetLevs) ! Flag indicating if obs outside input vertical range (0 for no)
    logical, intent(inout), optional :: success_opt(numTargetLevs)  ! success of interpolation

    ! Locals:
    integer :: TargetIndex   
    logical :: success(numTargetLevs)
    character(len=20) :: skipType

    integer :: levelIndex,ILMAX2,ILMIN2
    integer :: ILMIN, ILMAX
    real(8) :: SumWeights, TargetLayerThickWgt, ptop, pbtm
    
    if (present(skipType_opt)) then
      skipType = skipType_opt
    else
      skipType = 'default'
    end if 
 
    if (present(success_opt)) then
      if ( trim(skipType) == 'doAll&noExtrap') then
        if (present(outbound_opt)) then
          where (outbound_opt(:) == 0)
             success(:) = .true.
          elsewhere
             success(:) = .false.
          end where
        else
           success(:) = .true.
        end if
      else
        success(:) = success_opt(:)
      end if 
    else
      success(:) = .true.
    end if
        
    do TargetIndex=1,numTargetLevs

      if ( .not.success(TargetIndex) ) then
        wgts(TargetIndex,:) = 0.0D0
        wgts_opt(TargetIndex,:) = 0.0D0          
        kstart(TargetIndex)=1
        kend(TargetIndex)=1
        cycle
      end if

      ptop = targetLayersTop(TargetIndex)
      pbtm = targetLayersBot(TargetIndex)
      TargetLayerThickWgt=1.0D0/(min(pbtm,boundaries(numInputLevs+1))-max(ptop,boundaries(1)))
         
      ! Find the range of vertical levels over which to perform the averaging
      ! and set the averaging weights over this range.
          
      ilmin=1
      ilmax=numInputLevs
      if (ptop <= boundaries(1)*1.01 .and. pbtm >= boundaries(numInputLevs+1)*0.99) then

        ! Total column averaging part

        SumWeights=1.0D0/sum(weights(1:numInputLevs,1))
        !$OMP PARALLEL DO PRIVATE(levelIndex)                  
        do levelIndex = 1,numInputLevs
          wgts(TargetIndex,levelIndex)=(boundaries(levelIndex+1) &
                -boundaries(levelIndex))*TargetLayerThickWgt
          if (present(wgts_opt)) &
            wgts_opt(TargetIndex,levelIndex)=SumWeights
        end do
        !$OMP END PARALLEL DO                 
             
      else

        ! Partial column averaging part (special treatment at boundaries)
     
        ! Identify the vertical input level indices for levels
        ! within target layer and nearest specified layer boundary. 

        ilmin = ppo_getLevelIndex(ptop, boundaries, 'top', numInputLevs+1)
        ilmax = ppo_getLevelIndex(pbtm, boundaries, 'btm', numInputLevs+1)
               
        if (ilmin == ilmax+1) then

          ! Entire target layer within one input layer
                
          levelIndex=ilmin
          if (levelIndex < 1) levelIndex=1
          if (levelIndex > numInputLevs) levelIndex=numInputLevs

          wgts(TargetIndex,levelIndex)=1.0D0
          if (present(wgts_opt)) wgts_opt(TargetIndex,levelIndex)=1.0D0
	  
          ilmin=levelIndex
          ilmax=levelIndex+1
                  
        else
                
          ! Determine terms from the inner layers (excluding the lower and upper
          ! boundary layers when these layers not covering entire input layers)
                
          if (pbtm >= boundaries(numInputLevs)*0.99) then
            ilmax2=numInputLevs
          else
            ilmax2=ilmax-1
          end if
          if (ptop <= boundaries(1)*1.01) then
            ilmin=1
            ilmin2=ilmin
          else
            ilmin2=ilmin
          end if	
	    
          SumWeights=1.0D0/sum(weights(ilmin:ilmax,1))
	  
          if (ilmin2 <= ilmax2) then
            !$OMP PARALLEL DO PRIVATE(levelIndex)                  
            do levelIndex = ilmin2,ilmax2
              wgts(TargetIndex,levelIndex)= &
	        (boundaries(levelIndex+1)-boundaries(levelIndex))*TargetLayerThickWgt
              if (present(wgts_opt)) wgts_opt(TargetIndex,levelIndex)=SumWeights
            end do
            !$OMP END PARALLEL DO               
          end if
                
          ! Determine terms from the lower and upper boundary layers
          ! when these layers do not cover entire input layers.
                
          if (pbtm < boundaries(numInputLevs)*0.99) then
                     
            levelIndex=ilmax+1
            if (levelIndex > numInputLevs) levelIndex=numInputLevs
            if (levelIndex < 1) levelIndex=1
                   
            wgts(TargetIndex,levelIndex)= &
	      (pbtm - boundaries(ilmax))*TargetLayerThickWgt
            
            if (present(wgts_opt)) wgts_opt(TargetIndex,levelIndex)=SumWeights

            ilmax=levelIndex
                  
          end if
                  
          if (ptop > boundaries(1)*1.01) then
                     
            levelIndex=ilmin-1
            if (levelIndex < 1) levelIndex=1
            if (levelIndex > numInputLevs) levelIndex=numInputLevs
                   
            wgts(TargetIndex,levelIndex)= &
	      (boundaries(ilmin)-ptop)*TargetLayerThickWgt

            if (present(wgts_opt)) wgts_opt(TargetIndex,levelIndex)=SumWeights
	      
            ilmin=levelIndex
            if (ilmax < levelIndex+1) ilmax=levelIndex+1
                   
          end if
          if (ilmin > ilmax-1) ilmin=ilmax-1
        end if
      end if

      kstart(TargetIndex)=ilmin
      kend(TargetIndex)=ilmax
       
    end do

    if (present(dealloc_opt)) then
      if (dealloc_opt) deallocate(weights,boundaries)
    end if
    
  end subroutine ppo_vertAvgWgts
   
end module presProfileOperators_mod
