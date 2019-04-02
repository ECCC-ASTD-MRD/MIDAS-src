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
!! MODULE physicsFunctions (prefix='phf' category='7. Low-level data objects and utilities')
!!
!! *Purpose*: A collection of basic functions for various purposes 
!!            (e.g. computing saturation vapour pressure)
!!
!--------------------------------------------------------------------------
module physicsFunctions_mod
  use MathPhysConstants_mod
  use earthconstants_mod
  use utilities_mod
  
  implicit none
  private

  ! public procedures
  public :: FOEW8, FODLE8, FOQST8, FODQS8, FOEFQ8, FOQFE8, FOTVT8, FOTTV8
  public :: FOHR8, FOEWA8, FODLA8, FOQSA8, FODQA8, FOHRA8, FOTW8, FOTI8
  public :: FODTW8, FODTI8, FOTWI8, FODTWI8, FOEW8_CMAM, FOEI8_CMAM, FOERAT8_CMAM
  public :: FOEWI8_CMAM, FODLE8_CMAM, FOQST8_CMAM, FOTW8_CMAM, FOTI8_CMAM, FODTW8_CMAM
  public :: FODTI8_CMAM, FOTWI8_CMAM, FODTWI8_CMAM, FQBRANCH, FOEFQL, fotvvl, FOEFQA
  public :: FOEFQPSA, fottva, folnqva
  public :: phf_convert_z_to_pressure,phf_convert_z_to_gz
  public :: phf_get_tropopause, phf_get_pbl, phf_calcDistance, phf_calcDistanceFast
  public :: phf_alt2geopotential, phf_gravityalt, phf_gravitysrf

  LOGICAL :: initialized = .false.
  LOGICAL :: NEW_TETENS_COEFS

!**s/r physicsFunctions  - REAL*8 thermodyanmic statement functions.
!
!     Author: JM Belanger CMDA/SMC  Aug. 2000
!
!     Revisions: 
!                Yves J. Rochon *ARQX/SMC, Sept 2004
!                   - Added functions of sat. vapour pressure
!                     over water (FOEW8_CMAM) and over ice (FOEI8_CMAM),
!                     resultant sat. specific humidity (FOQST8_CMAM), and
!                     others.
!                Yves Rochon, ARQI/AQRD, Feb 2017
!                   - Added phf_convert* and phf_get_*
!
!                Stephane Laroche, Sept 2017
!                   - Added AERK thermodynamic fonctions for FOTW8 and FODTW8
!
!     REAL*8 version of thermodynamic functions based on
!     fintern.cdk in the physics library.
!

  contains

!
! ==============================================================================
!

     subroutine tetens_coefs_switch

      INTEGER*4      :: NULNAM,IER,FNOM,FCLOS
      CHARACTER *256 :: NAMFILE

      NEW_TETENS_COEFS = .false.
      NAMELIST /NAMPHY/NEW_TETENS_COEFS
      NAMFILE=trim("flnml")
      nulnam=0
      IER=FNOM(NULNAM,NAMFILE,'R/O',0)

      READ(NULNAM,NML=NAMPHY,IOSTAT=IER)
      if(IER.ne.0) then
        write(*,*) 'No valid namelist NAMPHY found'
      endif

      iER=FCLOS(NULNAM)

      write(*,*) 'new_tetens_coefs = ',new_tetens_coefs
      initialized = .true.

     end subroutine tetens_coefs_switch

!
! ==============================================================================
!
!     FONCTION DE TENSION DE VAPEUR SATURANTE (TETENS) - EW OU EI SELON TT
      real*8 function FOEW8(TTT)
        implicit none
        real*8 TTT
        FOEW8 = 610.78D0*EXP( MIN(SIGN(17.269D0,TTT-MPC_TRIPLE_POINT_R8),SIGN &
          (21.875D0,TTT-MPC_TRIPLE_POINT_R8))*ABS(TTT-MPC_TRIPLE_POINT_R8)/ &
          (TTT-35.86D0+MAX(0.D0,SIGN(28.2D0,MPC_TRIPLE_POINT_R8-TTT))))
      end function foew8
!
!     FONCTION CALCULANT LA DERIVEE SELON T DE  LN EW (OU LN EI)
      real*8 function FODLE8(TTT)
        implicit none
        real*8 TTT
        FODLE8=(4097.93D0+MAX(0.D0,SIGN(1709.88D0,MPC_TRIPLE_POINT_R8-TTT))) &
         /((TTT-35.86D0+MAX(0.D0,SIGN(28.2D0,MPC_TRIPLE_POINT_R8-TTT)))*  &
         (TTT-35.86D0+MAX(0.D0,SIGN(28.2D0,MPC_TRIPLE_POINT_R8-TTT))))
      end function FODLE8
!
!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT)
      real*8 function FOQST8(TTT,PRS) 
        implicit none
        real*8 TTT,PRS
        FOQST8=MPC_EPS1_R8/(MAX(1.D0,PRS/FOEW8(TTT))-MPC_EPS2_R8)
      end function FOQST8
!
!     FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
      real*8 function FODQS8(QST,TTT) 
        implicit none
        real*8 TTT,QST
        FODQS8=QST*(1.D0+MPC_DELTA_R8*QST)*FODLE8(TTT)
      end function FODQS8
!     QST EST LA SORTIE DE FOQST
!
!     FONCTION CALCULANT TENSION VAP (EEE) FN DE HUM SP (QQQ) ET PRS
      real*8 function FOEFQ8(QQQ,PRS)  
        implicit none
        real*8 QQQ,PRS
        FOEFQ8= MIN(PRS,(QQQ*PRS) / (MPC_EPS1_R8 + MPC_EPS2_R8*QQQ))
      end function FOEFQ8
!
!     FONCTION CALCULANT HUM SP (QQQ) DE TENS. VAP (EEE) ET PRES (PRS)
      real*8 function FOQFE8(EEE,PRS)  
        implicit none
        real*8 EEE,PRS
        FOQFE8= MIN(1.D0,MPC_EPS1_R8*EEE / (PRS-MPC_EPS2_R8*EEE))
      end function FOQFE8
!
!     FONCTION CALCULANT TEMP VIRT. (TVI) DE TEMP (TTT) ET HUM SP (QQQ)
      real*8 function FOTVT8(TTT,QQQ)  
        implicit none
        real*8 TTT,QQQ
        FOTVT8= TTT * (1.0D0 + MPC_DELTA_R8*QQQ)
      end function FOTVT8
!
!     FONCTION CALCULANT TTT DE TEMP VIRT. (TVI) ET HUM SP (QQQ)
      real*8 function FOTTV8(TVI,QQQ)  
        implicit none
        real*8 TVI,QQQ
        FOTTV8= TVI / (1.0D0 + MPC_DELTA_R8*QQQ)
      end function FOTTV8
!
!     FONCTION CALCULANT HUM REL DE HUM SP (QQQ), TEMP (TTT) ET PRES (PRS)
!     HR = E/ESAT
      real*8 function FOHR8(QQQ,TTT,PRS) 
        implicit none
        real*8 QQQ,TTT,PRS
        FOHR8 = MIN(PRS,FOEFQ8(QQQ,PRS)) / FOEW8(TTT)
      end function FOHR8
!
!     LES 5 FONCTIONS SUIVANTES SONT VALIDES DANS LE CONTEXTE OU ON
!     NE DESIRE PAS TENIR COMPTE DE LA PHASE GLACE DANS LES CALCULS
!     DE SATURATION.
!
!     FONCTION DE VAPEUR SATURANTE (TETENS)
      real*8 function FOEWA8(TTT) 
        implicit none
        real*8 TTT
        FOEWA8=610.78D0*EXP(17.269D0*(TTT-MPC_TRIPLE_POINT_R8)/(TTT-35.86D0))
      end function FOEWA8
!
!     FONCTION CALCULANT LA DERIVEE SELON T DE LN EW
      real*8 function FODLA8(TTT) 
        implicit none
        real*8 TTT
        FODLA8=17.269D0*(MPC_TRIPLE_POINT_R8-35.86D0)/(TTT-35.86D0)**2
      end function FODLA8
!
!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE
      real*8 function FOQSA8(TTT,PRS) 
        implicit none
        real*8 TTT,PRS
        FOQSA8=MPC_EPS1_R8/(MAX(1.D0,PRS/FOEWA8(TTT))-MPC_EPS2_R8)
      end function FOQSA8
!
!     FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
      real*8 function FODQA8(QST,TTT) 
        implicit none
        real*8 QST,TTT
        FODQA8=QST*(1.D0+MPC_DELTA_R8*QST)*FODLA8(TTT)
      end function FODQA8
!
!     FONCTION CALCULANT L'HUMIDITE RELATIVE
      real*8 function FOHRA8(QQQ,TTT,PRS) 
        implicit none
        real*8 QQQ,TTT,PRS
        FOHRA8=MIN(PRS,FOEFQ8(QQQ,PRS))/FOEWA8(TTT)
      end function FOHRA8
!
!     LES 6 FONCTIONS SUIVANTES SONT REQUISES POUR LA TEMPERATURE
!     EN FONCTION DE LA TENSION DE VAPEUR SATURANTE.
!     (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)
!
!     FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
!     SATURANTE PAR RAPPORT A EW.
      real*8 function FOTW8(EEE) 
        implicit none
        real*8 EEE

       if(.not.initialized) call tetens_coefs_switch
       if(new_tetens_coefs) then
         FOTW8=(30.11D0*LOG(EEE/610.94D0)-17.625D0*MPC_TRIPLE_POINT_R8)/ &
              (LOG(EEE/610.94D0)-17.625D0)
        else
         FOTW8=(35.86D0*LOG(EEE/610.78D0)-17.269D0*MPC_TRIPLE_POINT_R8)/ &
              (LOG(EEE/610.78D0)-17.269D0)
        endif

      end function FOTW8
!
!     FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
!     SATURANTE PAR RAPPORT A EI.
      real*8 function FOTI8(EEE) 
        implicit none
        real*8 EEE
        FOTI8=(7.66D0*LOG(EEE/610.78D0)-21.875D0*MPC_TRIPLE_POINT_R8)/ &
             (LOG(EEE/610.78D0)-21.875D0)
      end function FOTI8
!
!     FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
!     VAPEUR SATURANTE (EW).
      real*8 function FODTW8(TTT,EEE) 
        implicit none
        real*8 TTT,EEE

        if(.not.initialized) call tetens_coefs_switch
        if(new_tetens_coefs) then
         FODTW8=(30.11D0-TTT)/EEE/(LOG(EEE/610.94D0)-17.625D0)
        else
         FODTW8=(35.86D0-TTT)/EEE/(LOG(EEE/610.78D0)-17.269D0)
        endif

      end function FODTW8
!
!     FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
!     VAPEUR SATURANTE (EI).
      real*8 function FODTI8(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FODTI8=(7.66D0-TTT)/EEE  &
                      /(LOG(EEE/610.78D0)-21.875D0)
      end function FODTI8
!
!     FONCTION DE L'AJUSTEMENT DE LA TEMPERATURE.
      real*8 function FOTWI8(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FOTWI8=MAX(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*FOTW8(EEE)  &
                     -MIN(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*FOTI8(EEE)
      end function FOTWI8
!
!     FONCTION DE L'AJUSTEMENT DE LA DERIVEE DE LA TEMPERATURE.
      real*8 function FODTWI8(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FODTWI8=MAX(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*FODTW8(TTT,EEE)  &
                      -MIN(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*FODTI8(TTT,EEE)
      end function FODTWI8
!
!     LES 7 FONCTIONS SUIVANTES POUR EW-EI SONT REQUISES POUR LES MODELES
!     CMAM ET CGCM. (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)
!
!     FONCTION DE TENSION DE VAPEUR SATURANTE - EW
      real*8 function FOEW8_CMAM(TTT)  
        implicit none
        real*8 TTT
        FOEW8_CMAM= 100.D0*EXP(21.656D0-5418.D0/TTT)
      end function FOEW8_CMAM

!     FONCTION DE TENSION DE VAPEUR SATURANTE - EI
      real*8 function FOEI8_CMAM(TTT)  
        implicit none
        real*8 TTT
        FOEI8_CMAM= 100.D0*EXP(24.292D0-6141.D0/TTT)
      end function FOEI8_CMAM
!
!     FONCTION DE LA PROPORTION DE LA CONTRIBUTION DE EW VS EI
      real*8 function FOERAT8_CMAM(TTT) 
        implicit none
        real*8 TTT
        FOERAT8_CMAM=MIN(1.0D0,0.0059D0+0.9941D0*EXP(-0.003102D0*  &
            MIN(0.0D0,TTT-MPC_TRIPLE_POINT_R8)**2))
      end function FOERAT8_CMAM
!
!     FONCTION DE TENSION DE VAPEUR SATURANTE RESULTANTE - EW et EI
      real*8 function FOEWI8_CMAM(TTT)  
        implicit none
        real*8 TTT
        FOEWI8_CMAM= FOEW8_CMAM(TTT)*FOERAT8_CMAM(TTT)  &
         +(1.0D0-FOERAT8_CMAM(TTT))*FOEI8_CMAM(TTT)
      end function FOEWI8_CMAM
!
!     FONCTION DE LA DERIVE DE LN(E) PAR RAPPORT A LA TEMPERATURE
      real*8 function FODLE8_CMAM(TTT)  
        implicit none
        real*8 TTT
        FODLE8_CMAM= FOERAT8_CMAM(TTT)*5418.D0/TTT/TTT  &
         +(1.0D0-FOERAT8_CMAM(TTT))*6141.D0/TTT/TTT
      end function FODLE8_CMAM
!
!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT).
!     PREND EN COMPTE LES PHASES GLACE ET EAU.
      real*8 function FOQST8_CMAM(TTT,PRS) 
        implicit none
        real*8 TTT,PRS
        FOQST8_CMAM=MPC_EPS1_R8/(MAX(1.0D0,PRS/  &
                           FOEWI8_CMAM(TTT))-MPC_EPS2_R8)
      end function FOQST8_CMAM
!
!     LES 6 FONCTIONS SUIVANTES SONT REQUISES POUR LA TEMPERATURE
!     EN FONCTION DE LA TENSION DE VAPEUR SATURANTE POUR CMAM/CGCM
!     (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)
!
!     FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
!     SATURANTE PAR RAPPORT A EW.
      real*8 function FOTW8_CMAM(EEE) 
        implicit none
        real*8 EEE
        FOTW8_CMAM=5418.D0/(21.656D0-LOG(EEE/100.0D0))
      end function FOTW8_CMAM
!
!     FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
!     SATURANTE PAR RAPPORT A EI.
      real*8 function FOTI8_CMAM(EEE) 
        implicit none
        real*8 EEE
        FOTI8_CMAM=6141.D0/(24.292D0-LOG(EEE/100.0D0))
      end function FOTI8_CMAM
!
!     FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
!     VAPEUR SATURANTE (EW).
      real*8 function FODTW8_CMAM(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FODTW8_CMAM=TTT/EEE/(21.656D0-LOG(EEE/100.0D0))
      end function FODTW8_CMAM
!
!     FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
!     VAPEUR SATURANTE (EI).
      real*8 function FODTI8_CMAM(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FODTI8_CMAM=TTT/EEE/(24.292D0-LOG(EEE/100.0D0))
      end function FODTI8_CMAM
!
!     FONCTION DE L'AJUSTEMENT DE LA TEMPERATURE.
      real*8 function FOTWI8_CMAM(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FOTWI8_CMAM= FOERAT8_CMAM(TTT)*FOTW8_CMAM(EEE)+ &
                   (1.0D0-FOERAT8_CMAM(TTT))*FOTI8_CMAM(EEE)
      end function FOTWI8_CMAM
!
!     FONCTION DE L'AJUSTEMENT DE LA DERIVEE DE LA TEMPERATURE.
      real*8 function FODTWI8_CMAM(TTT,EEE) 
        implicit none
        real*8 TTT,EEE
        FODTWI8_CMAM= FOERAT8_CMAM(TTT)*FODTW8_CMAM(TTT,EEE)+  &
                    (1.0D0-FOERAT8_CMAM(TTT))*FODTI8_CMAM(TTT,EEE)
      end function FODTWI8_CMAM
!

!
!     function returning 0/1 depending on the minimum q branch condition
!     as discussed by Brunet (1996) to prevent getting a vapour pressure that exceeds
!     the total pressure p when q exceeds 1.
!
      real*8 function FQBRANCH(QQQ)  
        implicit none
        real*8 QQQ
        FQBRANCH= 0.5D0+SIGN(0.5D0,1.D0-(QQQ))
      end function FQBRANCH

!=============================================================================
!
!     TLM  of  THERMODYNAMIC FUNCTIONS USED IN 3DVAR
!     CONSTANTS FROM COMMON /CTESDYN/
!     NOTE: ALL FUNCTIONS WORK IN  S.I. UNITS
!           I.E PRS IN PA, QQQ IN KG/KG
!
!          ***C. Chouinard August 1998 ***
! Revision:
!          S. Pellerin *ARMA/AES - Sept. 1998
!                      -Tangent-linear operator of Tv
!
!
!      TLM  OF FUNCTION CALCULATING VAPOUR PRESSURE
!             - INPUT:  QQL ,  PERTURBATION OF LN SPECIFIC HUM
!                       PRSL ,   PERTURBATION OF SURFACE PRESSURE
!                       QQQ   ,  SPECIFIC HUMIDITY
!                       PRS   , PRESSURE
!                       PNETA   , VALUE OF ETA LEVEL
!             - OUTPUT: FOEFQL,  PERTURBATION  OF VAPOUR PRESSURE
!
      real*8 function FOEFQL(QQL,PRSL,QQQ,PRS,PNETA)  
        implicit none
        real*8 QQL,PRSL,QQQ,PRS,PNETA
        FOEFQL= FQBRANCH(QQQ)  &
           * ((QQL*MPC_EPS1_R8*PRS*QQQ/(MPC_EPS1_R8+MPC_EPS2_R8*QQQ)  &
           +  PRSL*PNETA*QQQ)/(MPC_EPS1_R8+MPC_EPS2_R8*QQQ))  &
           + (1.0D0 - FQBRANCH(QQQ))*PRSL*PNETA
      end function FOEFQL
!
!
!---- Tangent-linear operator of virtual temperature -----
!
!     qqq: backgroud specific humidity
!     ttt: backgroud temperature
!     ttl: temperature increment
!     plnql: increment of logarithm specific humidity  (del(ln q))
!
      real*8 function fotvvl(qqq,ttt,ttl,plnql) 
        implicit none
        real*8 qqq,ttt,ttl,plnql
        fotvvl=(1 + MPC_DELTA_R8*qqq)*ttl + MPC_DELTA_R8*qqq*ttt*plnql
      end function fotvvl

!=============================================================================
!
!   DEFINITION OF ADJOINTS OF THERMODYNAMIC  FUNCTIONS
!   CONSTANTS AS IN COMMON /CTESDYN/
!     NOTE: ALL  UNITS S.I.
!           I.E. PADES IN DEG K, PRS EN PA, QQQ EN KG/KG
!
!      ADJOINT OF LN SPECIFIC  HUM (QQQ) DUE TO DEWPOINT DEPRESSION CORRECTIONS
!             - INPUT : PADES ,  ADJOINT OF DEWPOINT DEPRESSION
!                       PGAMMA,  ADOINT OF VAPOUR PRESSURE RELATIONSHIP
!                       QQQ   , SPECIFIC HUMIDITY
!                       PRS   , PRESSURE
!             - OUTPUT: FOEFQA, ADJOINT OF LN SPECIFIC HUMIDITY
!
      real*8 function FOEFQA(PADES,PGAMMA,QQQ,PRS)  
        implicit none
        real*8 PADES,PGAMMA,QQQ,PRS
        FOEFQA= PADES*PGAMMA*MPC_EPS1_R8*PRS*QQQ/  &
                ((MPC_EPS1_R8+MPC_EPS2_R8*QQQ)*(MPC_EPS1_R8+MPC_EPS2_R8*QQQ))
      end function FOEFQA

!
!      ADJOINT OF SURFACE PRESSURE  DUE TO DEWPOINT DEPRESSION CORRECTIONS
!             - INPUT:  PADES ,  ADJOINT OF DEWPOINT DEPRESSION
!                       PGAMMA,  ADOINT OF VAPOUR PRESSURE RELATIONSHIP
!                       QQQ   , SPECIFIC HUMIDITY
!                       PNETA   , VALUE OF NETA
!             - OUTPUT: FOEFQPSA, ADJOINT OF SURFACE PRESSURE
!
      real*8 function FOEFQPSA(PADES,PGAMMA,QQQ,PNETA) 
        implicit none
        real*8 PADES,PGAMMA,QQQ,PNETA
        FOEFQPSA = PADES*PGAMMA*QQQ*PNETA/  &
                   (MPC_EPS1_R8+MPC_EPS2_R8*QQQ)
      end function FOEFQPSA

!
!--------------------- Adjoint of virtual temperature operator -------------------------
!
!     fottva: Adjoint of temperature due to virtual temperature correction
!     qqq:   background specific humidity
!     tva:   adjoint variable of virtual temperature
!
      real*8 function fottva(qqq,tva)  
        implicit none
        real*8 qqq,tva
        fottva= (1D0 + MPC_DELTA_R8*qqq)*tva
      end function fottva

!
!     folnqva: Adjoint of logarithm of specific humidity due to virtual temperature correction
!     qqq:   background specific humidity
!     ttt:   background temperature
!     tva:   adjoint variable of virtual temperature
!
      real*8 function folnqva(qqq,ttt,tva) 
        implicit none
        real*8 qqq,ttt,tva
        folnqva = MPC_DELTA_R8*qqq*ttt*tva
      end function folnqva

  !-------------------------------------------------------------------------------------------

  function  phf_convert_z_to_pressure(altitude,rgz_mod,press_mod,nlev,nlev_mod,lat,success) result(press)
  ! 
  ! Author   : M. Sitwell, May 2015
  !          
  ! Purpose: Converts an array of (geometric) altitudes to pressures. Uses linear interpolation
  !          in log(p).  
  !
  ! Arguments:
  !
  !   Input
  !     altitude      altitudes to convert to pressures (m)
  !     rgz_mod       geopotential heights on model levels (m), assumed to be in decending order
  !     press_mod     pressure on model levels, assumed to be in ascending order
  !     nlev          length of altitude array
  !     nlev_mod      number of model levels
  !     lat           latitude (rad)
  !
  !   Output
  !     press         converted pressures 
  !
  !!---------------------------------------------------------------------------------------
    
    implicit none

    real(8), intent(in) :: altitude(nlev),rgz_mod(nlev_mod),press_mod(nlev_mod),lat
    integer, intent(in) :: nlev,nlev_mod
    logical, intent(inout) :: success(nlev)
    real(8) :: press(nlev),rgz(nlev)
    integer :: ilev,ilev_mod

    ! Check model geopotential heights and pressures
    if (any(rgz_mod.lt.0.0D0)) call utl_abort("phf_compute_z_to_pressure: Invalid model GZ.")
    if (any(press_mod.lt.0.0D0)) call utl_abort("phf_compute_z_to_pressure: Invalid model pressure.")
    
    ! Convert altitudes to geopotential heights
    rgz = phf_convert_z_to_gz(altitude,lat,nlev)

    do ilev=1,nlev

       ! Check if height is above or below model boundaries
       if ( rgz(ilev).gt.rgz_mod(1) .or. rgz(ilev).lt.rgz_mod(nlev_mod) ) then
          success(ilev)=.false.
       end if

       if (success(ilev)) then

          ! Find model layers directly above and below rgz(ilev).
          ! After exit of loop we will have 
          ! rgz_mod(ilev_mod) >= rgz(ilev) > rgz_mod(ilev_mod+1)
          do ilev_mod=1,nlev_mod-1
             if ( rgz(ilev).le.rgz_mod(ilev_mod) .and. &
                  rgz(ilev).gt.rgz_mod(ilev_mod+1) ) exit
          end do
          
          ! Linear interpolation in gz,log(p)
          press(ilev) = press_mod(ilev_mod+1) * (press_mod(ilev_mod)/press_mod(ilev_mod+1))**( &
               (rgz(ilev)-rgz_mod(ilev_mod+1))/(rgz_mod(ilev_mod)-rgz_mod(ilev_mod+1)) )

       else
          press(ilev) = 0.0
       end if

    end do

  end function phf_convert_z_to_pressure

  !----------------------------------------------------------------------------------------

  function phf_convert_z_to_gz(altitude,lat,nlev) result(rgz)
  !
  ! Author   : M. Sitwell, June 2015
  !          
  ! Purpose: Converts altitudes to geopotential heights. Uses the Helmert formula to
  !          parameterize the latitude dependence and uses analytical result of the
  !          integral of \int g(z)dz for the altitude dependence (see J.A. Dutton 1976,
  !          p.65). At an altitude of 50 km, the altitude and geopotential height
  !          differ by around 0.2-0.5 km, depending on the latitude.
  !
  ! Arguments: 
  !
  !   Input
  !     altitude      altitudes (m)
  !     lat           latitude (rad) 
  !
  !   Output
  !     rgz           geopotential heights (m)
  !
  !!---------------------------------------------------------------------------------------
  
    implicit none

    real(8), intent(in) :: altitude(nlev),lat
    integer, intent(in) :: nlev
    real(8) :: rgz(nlev)

    rgz = (RG/9.8) * (1.-2.64D-03*cos(2.*lat)+5.9D-6*cos(2.*lat)**2) * RA*altitude/(RA+altitude)

  end function phf_convert_z_to_gz
 
  !----------------------------------------------------------------------------------------

  function phf_get_tropopause(nmodlev,pressmod,tt,gz,hu_opt) result(tropo_press)
  !
  ! Author   : Y. Rochon, ARQI/AQRD Oct 2015
  !            - Following consultation with Irena Paunova for water vapour based approach
  !                                     with Sylive Gravel (and wikipedia) for temperature based approach
  !          
  ! Purpose: Determines pressure level of tropopause. 
  !          Final tropopause is taken as max pressure (lowest altitude) from the
  !          water vapour and temperature based tropopauses.
  !
  ! Revisions: 
  !          
  ! Arguments:
  !
  !   Input
  !
  !      nmodlev      Number of model levels
  !      pressmod     Model pressure array (Pa)
  !      hu           Model specific humidity 
  !      tt           Model temperature (Kelvin)
  !      gz           Model geopotential height (m)
  !
  !   Output
  !
  !     tropo_press   Tropopause level in Pa
  !             
  !!---------------------------------------------------------------------------------------
    
    implicit none

    integer, intent(in) :: nmodlev
    real(8), intent(in) :: pressmod(nmodlev),tt(nmodlev),gz(nmodlev)
    real(8), intent(in), optional :: hu_opt(nmodlev)
   
    real(8) :: tropo_press
  
    integer :: itop,i,k,ilaps
    real(8) :: hu_ppmv1,hu_ppmv2,hu_ppmv3,xlaps,tropo_press_hu
    real(8), parameter :: press_min=6000.         ! Min tropoause pressure 60 hPa.; equivalent to ~ 20km
    real(8), parameter :: gz_min=6000.0           ! Min tropopause level in meters.
    real(8), parameter :: ppmv_threshold=10.0     
    real(8), parameter :: tgrad_threshold=0.002   ! degrees/m (2 degrees/km)
    real(8), parameter :: consth=0.160754938e+07  ! conversion from mass mixing ratio to ppmv;  1.0e+06 / (18.015/28.96)

    tropo_press=-1.0
    if (all(gz.lt.0.0))  call utl_abort('phf_get_tropopause: Missing GZ for determining tropopause pressure')

    ! Initialize tropopause pressure level using temperature gradient.
    ! Thermal tropopause is defined as the lowest level (above gz_min) at which (1) the lapse rate decreases
    ! to <= 2 C/km and (2) the average lapse rate between this level and all higher levels within 2 km are <= 2 C/km. 
    ! Ref: International Meteorological Vocabulary (2nd ed.). Geneva: Secretariat of the World Meteorological
    !      Organization. 1992. p. 636. ISBN 92-63-02182-1.
    ! The second requirement, based on hu, may give levels that are to high (pressure too low) in the winter hemisphere.

    do itop=3,nmodlev
       if (pressmod(itop).ge.press_min) exit
    end do
    itop=itop-1
       
    do i=nmodlev,itop+1,-1
       if (gz(i)-gz(nmodlev).lt.gz_min) cycle
       xlaps=-(tt(i)-tt(i-1))/(gz(i)-gz(i-1))
       if (xlaps.le.tgrad_threshold) then
          ilaps=1
          do k=i-1,itop,-1
             if (gz(k)-gz(i).gt.2000.0) exit
             xlaps=xlaps-(tt(k)-tt(k-1))/(gz(k)-gz(k-1))
             ilaps=ilaps+1
          end do
          if (xlaps/ilaps.le.tgrad_threshold) exit
       end if
    enddo
    tropo_press=pressmod(i)
    
   ! Improve on tropopause pressure levels using specific humidity if available,

    if (present(hu_opt)) then
    
      !  Use water vapour
      
       hu_ppmv1=0.0
       do i=itop,nmodlev

          ! Convert specific humidity to ppmv mixing ratio.
          ! First apply r=q/(1-q) to convert to mass mixing ratio.

          if (hu_opt(i).le.0.8.and.hu_opt(i).ge.0) then
               hu_ppmv2 = consth*hu_opt(i)/(1.0-hu_opt(i))
          else if (hu_opt(i).gt.0.8) then
               hu_ppmv2 = consth*0.8/(1.0-0.8)
          else if (hu_opt(i).lt.0.0) then
               hu_ppmv2 = 0.0
          end if

          ! Check if transition point reached.
          ! Added requirement that levels below also satisfy this condition.

          if (hu_ppmv2.ge.ppmv_threshold) then
             ilaps=1
             do k=i+1,nmodlev
                if (gz(i)-gz(k).gt.5000.0) exit
                if (hu_opt(k).le.0.8.and.hu_opt(k).ge.0) then
                   hu_ppmv3 = consth*hu_opt(k)/(1.0-hu_opt(k))
                else if (hu_opt(k).gt.0.8) then
                   hu_ppmv3 = consth*0.8/(1.0-0.8)
                else
                   hu_ppmv3=0.0
                end if
                if (hu_ppmv3.lt.ppmv_threshold) ilaps=0
             end do
             if (ilaps.eq.1) exit
          end if
          hu_ppmv1=hu_ppmv2
       end do
       
       if (hu_ppmv2.ge.ppmv_threshold.and.ilaps.eq.1) then
       
            ! Interpolate between levels
      
             if (abs(hu_ppmv2-hu_ppmv1).lt.0.1) hu_ppmv1=hu_ppmv2-0.1
!             tropo_press_hu=(log(pressmod(i))*(ppmv_threshold-hu_ppmv1)+ &
!                    log(pressmod(i-1))*(hu_ppmv2-ppmv_threshold)) &
!                   /(hu_ppmv2-hu_ppmv1)
!             tropo_press_hu=exp(tropo_press_hu)
             tropo_press_hu=pressmod(i)
             
             tropo_press=min(tropo_press,tropo_press_hu)
       else
          write(*,*) 'phf_get_tropopause: Level and specific humidity: ',itop,hu_ppmv2
          call utl_abort('phf_get_tropopause: Specific humidity too small.')
       end if
                             
    end if
    
  end function phf_get_tropopause

  !----------------------------------------------------------------------------------------

  function phf_get_pbl(nmodlev,pressmod,tt,gz,hu_opt,uu_opt,vv_opt) result(pbl_press)
  !
  ! Author   : Y. Rochon, ARQI/AQRD Oct 2015
  !            - Following consultation with Amir Aliabadi, Shuzhan Ren and Saroja Polavarapu.
  !            - RiB case based on original routine 'mixing_properties' by Chris Golaz (GFDL) 
  !              and Sungsu Park(NCAR) provided by Shuzhan Ren.  
  ! 
  ! Purpose: Determines pressure level of planetary boundary layer using 
  !          a first threshold of 0.5 for the bulk Richadson number (after Mahrt, 1981; 
  !          requires availability of uu and vv). Threshold reduced to largest value
  !          between 0.25 and 0.5 if first not satisfied. 
  !
  !          If not found with this approach, applies a variant of the Heffter approach 
  !          described in Aliabadi et al (2016), with some local variation.
  !
  ! References:
  !
  !  Aliabadi A.A., R.M. Staebler, J. de Grandpre, A. Zadra, and P.A.
  !       Vaillancourt, 2016: Comparison of Estimated Atmospheric
  !       Boundary Layer Mixing Height in teh Arctic and Southern Great
  !       Plains under Statisticallt Stable Conditions: Experimental
  !       and Numerical Aspects, Submitted to Atmosphere-Ocean (2015).
  !  Mahrt, L. 1981: Modelling depth of the stable boundary-layer,
  !       Bound-Lay. Meteorol., 21, 3-19
  !  Heffter, J.L.,1980: Transport layer depth calculations, Second Joint Conference on Applications
  !       of Air Pollution Meteorology, New Orleans, LA, 24-27 March 1980. American Meteorological
  !       Society, Boston, MA.
  !       
  ! Revisions: 
  !          
  ! Arguments: 
  !
  !   Input
  !
  !      nmodlev      Number of model levels for variables other than uu and vv
  !      pressmod     Model pressure array (Pa)
  !      tt           Model temperature (Kelvin)
  !      gz           Model geopotential height (meters)
  !      hu           Specific humidity
  !      uu           Model zonal wind component (m/s)
  !      vv           Model meridional wind component (m/s)
  !
  !   Output
  !
  !     pbl_press     PBL level in Pa
  ! 
  !   Comments
  !
  !   A) Currently assumes (uu,vv) midlayer levels approximately at tt, gz, and hu levels
  !      when size(uu).ne.nmodlev.
  !
  !!---------------------------------------------------------------------------------------
        
    implicit none

    integer, intent(in) :: nmodlev
    real(8), intent(in) :: pressmod(nmodlev),tt(nmodlev),gz(nmodlev)
    real(8), optional :: uu_opt(:),vv_opt(:),hu_opt(nmodlev)
   
    real(8) :: pbl_press
  
    integer :: itop,i,id,igradmax,inv,iRiBmax
    real(8) :: RiB1,RiB2,RiBmax,zs,thetavs,thetavh(nmodlev),us,vs,uv,hus,huh,gradmax,grad
    real(8), parameter :: kappa = 287.04/1004.67  ! R/Cp
    real(8), parameter :: RiB_threshold=0.5, reduced=0.5 

    ! Imposed min presssure of PBL height of 200 hPa (extreme; PBL height should normally be under 3km)
    real(8), parameter :: press_min=20000.  

    real(8) :: huw(nmodlev)

    pbl_press=-1.0
     
    ! Set values for lowest prognostic level 

    i = nmodlev   
    
    if (all(gz.lt.0.0))  call utl_abort('phf_get_pbl: Missing GZ for determining PBL pressure')

    ! Convert hu to mass mixing ratio

    if (present(hu_opt)) then
       huw(:)=hu_opt(:)
    else
       huw(:)=0.0
    end if
    
    if (huw(i).le.0.8.and.huw(i).ge.0) then
        hus = huw(i)/(1.0-huw(i))
    else if (huw(i).gt.0.8) then
        hus = 0.8/(1.0-0.8)
    else if (huw(i).lt.0.0) then
        hus = 0.0
    end if
    zs = gz(i)*0.001
    
    ! Potential virtual temperature at lowest prognostic level

    thetavs = tt(i)*(1.D5/pressmod(i))**kappa* (1.0 + 0.61*hus )
    thetavh(nmodlev)=thetavs

    ! Set max vertical level

    do itop=2,nmodlev-1
       if (pressmod(itop).ge.press_min) exit
    end do

    RiB1=0.0
    RiB2=0.0
    RiBmax=0.0
    iRiBmax=0
    if (present(uu_opt).and.present(vv_opt)) then
       id=nmodlev-size(uu_opt)
       if (id.gt.1.or.id.lt.0) then
          call utl_abort('phf_get_pbl: Unexpected number of UV levels, nmodlev = ' // trim(utl_str(nmodlev)) // ' , size(uu) = ' // trim(utl_str(size(uu_opt))) )    
       end if
       us = uu_opt(size(uu_opt))
       vs = vv_opt(size(vv_opt))
!   us,vs set to 0.0
!      us=0.0
!      vs=0.0

       ! Calc RiB from near-surface to level attaining RiB_threshold
 
       do i=nmodlev-1,itop,-1
  
           if (huw(i).le.0.8.and.huw(i).ge.0) then
               huh = huw(i)/(1.0-huw(i))
           else if (huw(i).gt.0.8) then
               huh = 0.8/(1.0-0.8)
           else if (huw(i).lt.0.0) then
               huh = 0.0
           end if
           thetavh(i) = tt(i)*(1.D5/pressmod(i))**kappa* ( 1.0 + 0.61*huh )

           if (id.eq.0) then
               uv = max( (uu_opt(i)-us)**2 + (vv_opt(i)-vs)**2, 1.D-8 ) 
           else
              ! Take layer midpoint values
              uv = max( ((uu_opt(i)+uu_opt(i-1))/2.0-us)**2 + ((vv_opt(i)+vv_opt(i-1))/2.0-vs)**2, 1.D-8 ) 
           end if
         
           RiB2 = grav * (thetavh(i)-thetavs) * (gz(i)*0.001-zs) / (thetavs*uv)
           if (RiBmax.lt.RiB2.and.RiB2.ge.reduced*RiB_threshold) then
              RiBmax=RiB2
              iRiBmax=i
           end if
           if (RiB2.ge.RiB_threshold) exit
           RiB1=RiB2
       end do
    else
    
       ! Calc only theta
      
       do i=nmodlev-1,itop,-1
  
           if (huw(i).le.0.8.and.huw(i).ge.0) then
               huh = huw(i)/(1.0-huw(i))
           else if (huw(i).gt.0.8) then
               huh = 0.8/(1.0-0.8)
           else if (huw(i).lt.0.0) then
               huh = 0.0
           end if
           thetavh(i) = tt(i)*(1.D5/pressmod(i))**kappa* ( 1.0 + 0.61*huh )
       end do
    end if   
    
    if (RiB2.ge.RiB_threshold) then    
   
      !  Interpolate between levels

       pbl_press=(log(pressmod(i))*(RiB_threshold-RiB1)+ &
               log(pressmod(i+1))*(RiB2-RiB_threshold)) &
               /(RiB2-RiB1)
       pbl_press=exp(pbl_press)
    else if (RiBmax.ge.reduced*RiB_threshold) then
       ! Apply to level with largest RiB between reduced*RiB_threshold and RiB_threshold
       pbl_press=pressmod(iRiBmax)
    else
    
       ! Estimate PBL level using the Heffter conditions:
       ! First find lowest inversion layer where dtheta>2K.
       ! If found, assign mid of layer as PBL level. 
       ! Otherwise, assign PBL level as that with largest
       ! theta gradient.
       
       i=nmodlev-1
       do while (i.gt.itop) 
          !if (thetavh(i)-thetavh(i+1).gt.0.0) then
          if ((thetavh(i)-thetavh(i+1))/(gz(i)-gz(i+1)).ge.0.005) then
             ! Near bottom of inversion layer found
             inv=i+1
             i=i-1
             !do while (thetavh(i)-thetavh(i+1).gt.0.0.and.i.gt.itop) 
             do while ((thetavh(i)-thetavh(i+1))/(gz(i)-gz(i+1)).ge.0.005.and.i.gt.itop) 
                 i=i-1
             end do
             if ((thetavh(i+1)-thetavh(inv)).gt.2.0)  then
                ! Apply  midlayer as PBL
                pbl_press=sqrt(pressmod(i+1)*pressmod(inv))
                exit
             end if 
          else
             i=i-1
          end if
       end do
       if (pbl_press.le.0.0) then
          gradmax=-1.D30
          igradmax=nmodlev-1
          do i=nmodlev-1,itop,-1
             grad=(thetavh(i)-thetavh(i+1))/(gz(i)-gz(i+1))
             if (gradmax.lt.grad) then
                gradmax=(thetavh(i)-thetavh(i+1))/(gz(i)-gz(i+1))
                igradmax=i
                if (grad.ge.0.005) then ! Check next layer as well
                   if ((thetavh(i-1)-thetavh(i))/(gz(i-1)-gz(i)).ge.0.005) exit
                end if
            end if
          end do          
          pbl_press=pressmod(igradmax)
          ! write(*,*) 'phf_get_pbl: Warning2 - Max allowed altitude reached for. ',pbl_press,igradmax,gradmax,RiB2,iRiBmax,RiBmax
       !else
       !   write(*,*) 'phf_get_pbl: Warning1 - Max allowed altitude reached for. ',pbl_press,i,RiB2,iRiBmax,RiBmax     
       end if
    end if

  end function phf_get_pbl

  !--------------------------------------------------------------------------
  ! calcDistance
  !--------------------------------------------------------------------------
  function phf_calcDistance(lat2, lon2, lat1, lon1) result(distanceInM)
    implicit none
    ! **Purpose:**
    ! Compute the distance between two point on earth: (lat1,lon1) and (lat2,lon2).
    ! Calcul utilisant la Formule d'Haversine.
    ! Reference: R.W. Sinnott,'Virtues of Haversine',Sky and Telescope, vol.68, no.2, 1984, p.159)

    ! arguments
    real(8) :: lat1, lon1, lat2, lon2
    real(8) :: distanceInM

    ! locals
    real(8) :: dlat, dlon, a, c

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (sin(dlat/2.d0))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2.d0))**2
    c = 2.d0 * atan2(sqrt(a),sqrt(1.d0-a))
    distanceInM = RA * c

  end function phf_calcDistance

  !--------------------------------------------------------------------------
  ! calcDistanceFast
  !--------------------------------------------------------------------------
  function phf_calcDistanceFast(lat2, lon2, lat1, lon1) result(distanceInM)
    implicit none
    ! **Purpose:**
    ! Compute the distance between two point on earth: (lat1,lon1) and (lat2,lon2).
    ! Using a quick and dirty formula good for short distances not close to the pole.

    ! arguments
    real(8) :: lat1, lon1, lat2, lon2
    real(8) :: distanceInM

    ! locals
    real(8) :: dlat, dlon

    dlon = (lon2 - lon1)*cos(lat1)
    dlat = lat2 - lat1

    distanceInM = RA * sqrt(dlon*dlon + dlat*dlat)

  end function phf_calcDistanceFast


  function phf_gravitysrf(sLat)
    !  Normal gravity on ellipsoidal surface:
    !  Input:  Latitude
    !          sin(Latitude)
    !
    !  Output: Normal gravity
    !          phf_gravitysrf         : m/s2
    !
    real(8), intent(in)  :: sLat
    real(8)              :: phf_gravitysrf
    
    real(8)              :: ks2
    real(8)              :: e2s

    ks2 = WGS_TNGk * sLat*sLat
    e2s = 1.D0 - WGS_e2 * sLat*sLat
    phf_gravitysrf = WGS_GammaE * (1.D0 + ks2) / sqrt(e2s)
  end function phf_gravitysrf


  function phf_gravityalt(sLat, Altitude)
    ! Normal gravity above the ellipsoidal surface:
    ! Input:  Latitude, altitude
    !         sin(Latitude)
    !         Altitude               : m
    !
    ! Output: Normal gravity
    !         phf_gravityalt          : m/s2
    !
    real(8), intent(in)  :: sLat
    real(8), intent(in)  :: Altitude
    real(8)              :: phf_gravityalt

    real(8)              :: C1
    real(8)              :: C2

    C1 =-2.D0/WGS_a*(1.D0+WGS_f+WGS_m-2*WGS_f*sLat*sLat)
    C2 = 3.D0/WGS_a**2
    phf_gravityalt = phf_gravitysrf(sLat)*                                   &
         (1.D0 + C1 * Altitude + C2 * Altitude**2)
  end function phf_gravityalt


  subroutine phf_alt2geopotential(altitude, latitude, geopotential, printGZ)
    !**s/r phf_alt2geopotential - Geopotential energy at a given point.
    ! Result is based on the WGS84 approximate expression for the
    ! gravity acceleration as a function of latitude and altitude,
    ! integrated with the trapezoidal rule.
    !
    ! Input:  altitude(m), latitude (rad)
    ! Output: geopotential (m2/s2)
    !
    ! Author : M. Bani Shahabadi, November 2018
    !
    real(8), intent(in)   :: altitude(:)
    real(8), intent(in)   :: latitude
    real(8), intent(inout):: geopotential(:)

    integer           :: nlev, nlev500m
    real(8), allocatable :: alt500m(:), gravity500m(:)
    real(8)           :: delAlt, aveGravity, sLat, gravity, gravityM1
    integer           :: i, ilev 
    logical, optional :: printGZ
    

    nlev = size(altitude)
    sLat = sin(latitude)

    nlev500m = int(altitude(nlev) / 500.D0)
    if ( nlev500m >= 1) then
      allocate(alt500m(0:nlev500m))
      allocate(gravity500m(0:nlev500m))

      ! Calculate gravity and height of levels for 500m-layers
      do i = 0, nlev500m
        alt500m(i) = i * 500.0D0
        gravity500m(i) = phf_gravityalt(sLat, alt500m(i))
      enddo

      geopotential(nlev) = 0.0D0
      ! integrate from surface on the 500m-layers untill below the desired altitude
      do i = 1, nlev500m
        delAlt = alt500m(i) - alt500m(i-1)
        aveGravity = 0.5D0 * (gravity500m(i) + gravity500m(i-1))
        geopotential(nlev) = geopotential(nlev) + delAlt * aveGravity
      enddo

      ! Add the contribution from top of the last 500m-layer to the altitude  
      delAlt = altitude(nlev) - alt500m(nlev500m)
      aveGravity = 0.5D0 * (phf_gravityalt(sLat, altitude(nlev)) + gravity500m(nlev500m))
      geopotential(nlev) = geopotential(nlev) + delAlt * aveGravity 
      gravityM1 = phf_gravityalt(sLat, altitude(nlev))

      deallocate(alt500m)
      deallocate(gravity500m)

    else
      ! At surface, use local gravity to get GZ
      gravity = phf_gravityalt(sLat,altitude(nlev))
      geopotential(nlev) = altitude(nlev) * gravity
      gravityM1 = gravity

    endif

    ! At upper-levels, integrate on model levels to get GZ 
    do ilev = nlev-1, 1, -1
      gravity = phf_gravityalt(sLat,altitude(ilev))
      aveGravity = 0.5D0 * (gravity + gravityM1)
      delAlt = altitude(ilev) - altitude(ilev+1)
      geopotential(ilev) = geopotential(ilev+1) + delAlt * aveGravity
      gravityM1 = gravity
    enddo

    if ( present(printGZ) ) then
      if ( printGZ ) then
        write(*,*) 'phf_alt2geopotential, GZ_T:'
        write(*,*) geopotential(:)

        printGZ = .false.
      endif
    endif

  end subroutine phf_alt2geopotential

  
end module physicsFunctions_mod


