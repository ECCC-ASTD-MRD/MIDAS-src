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

module physicsFunctions_mod
  ! MODULE physicsFunctions_mod (prefix='phf' category='8. Low-level utilities and constants')
  !
  ! :Purpose: A collection of basic functions for various purposes 
  !           (e.g. computing saturation vapour pressure)
  !
  use MathPhysConstants_mod
  use earthConstants_mod
  use utilities_mod
  use midasMpi_mod
  use message_mod
  
  implicit none
  private

  ! public procedures
  public :: phf_FOEW8, phf_FODLE8, phf_FOQST8, phf_FODQS8, phf_FOEFQ8, phf_FOQFE8, phf_FOTVT8, phf_FOTTV8
  public :: phf_FOHR8, phf_FOEWA8, phf_FODLA8, phf_FOQSA8, phf_FODQA8, phf_FOHRA8, phf_FOTW8, phf_FOTI8
  public :: phf_FODTW8, phf_FODTI8, phf_FOTWI8, phf_FODTWI8, phf_FOEW8_CMAM, phf_FOEI8_CMAM, phf_FOERAT8_CMAM
  public :: phf_FOEWI8_CMAM, phf_FODLE8_CMAM, phf_FOQST8_CMAM, phf_FOTW8_CMAM, phf_FOTI8_CMAM, phf_FODTW8_CMAM
  public :: phf_FODTI8_CMAM, phf_FOTWI8_CMAM, phf_FODTWI8_CMAM, phf_FQBRANCH, phf_FOEFQL, phf_fotvvl, phf_FOEFQA
  public :: phf_FOEFQPSA, phf_fottva, phf_folnqva
  public :: phf_convert_z_to_pressure,phf_convert_z_to_gz
  public :: phf_get_tropopause, phf_get_pbl, phf_calcDistance, phf_calcDistanceFast
  public :: phf_height2geopotential, phf_gravityalt, phf_gravitysrf

  logical           :: phf_initialized = .false.

  ! namelist variables:
  character(len=20) :: saturationCurve   ! saturationCurve must be one of 'Tetens_1930', 'Tetens_2018a', 'Tetens_2018'

  contains

   !--------------------------------------------------------------------------
   ! tetens_coefs_switch
   !--------------------------------------------------------------------------
   subroutine phf_tetens_coefs_switch
     !
     ! :Purpose: Read water saturation strategy from nml. Options are:
     !
     !             - 'Tetens_1930' was active before 2018.
     !             - 'Tetens_2018a' is a partial update that missed some functions (active before IC4)
     !             - 'Tetens_2018' completes the update to the intended 2018 specification.
     !
     integer            :: NULNAM,IERR,FNOM,FCLOS
     character(len=256) :: NAMFILE
     logical            :: validOption
     NAMELIST /NAMPHY/ saturationCurve

     !$omp critical
     if (.not.phf_initialized) then

       saturationCurve = 'Tetens_2018'

       if ( .not. utl_isNamelistPresent('NAMPHY','./flnml') ) then
         call msg( 'phf_tetens_coefs_switch', &
              'NAMPHY is missing in the namelist. Default values will be taken.', mpiAll_opt=.False.)
       else
         ! Reading the namelist
         nulnam = 0
         ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
         read(nulnam, nml=namphy, iostat=ierr)
         if ( ierr /= 0) call utl_abort('tetens_coefs: Error reading namelist')
         ierr = fclos(nulnam)
       end if

       validOption = (trim(saturationCurve) == 'Tetens_1930'  .or.  &
                      trim(saturationCurve) == 'Tetens_2018a' .or.  &
                      trim(saturationCurve) == 'Tetens_2018')
       if (.not.validOption) then
         call utl_abort('phf_tetens: WV Saturation not in expected list')
       end if

       call msg( 'phf_tetens_coefs_switch ', saturationCurve )
       phf_initialized = .true.
     end if
     !$omp end critical

   end subroutine phf_tetens_coefs_switch

   !--------------------------------------------------------------------------
   ! phf_FOEW8
   !--------------------------------------------------------------------------
   real*8 function phf_FOEW8(TTT)
     !
     ! :Purpose: Water vapour saturation pressure (Tetens) - EW or EI as a fct of TT
     !
     implicit none
     real*8 TTT

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       ! Updated coefficients 2018
       if (TTT > MPC_K_C_DEGREE_OFFSET_R8) then
         phf_FOEW8 = 610.94D0*EXP(17.625D0 * (TTT-MPC_TRIPLE_POINT_R8) / (TTT-30.11D0))
       else
         phf_FOEW8 = 610.94D0*EXP(22.587D0 * (TTT-MPC_TRIPLE_POINT_R8) / (TTT+ 0.71D0))
       end if
     else
       ! Classic Tetens 1930 coefficients
       if (TTT > MPC_TRIPLE_POINT_R8) then
         phf_FOEW8 = 610.78D0*EXP(17.269D0 * (TTT-MPC_TRIPLE_POINT_R8) / (TTT-35.86D0))
       else
         phf_FOEW8 = 610.78D0*EXP(21.875D0 * (TTT-MPC_TRIPLE_POINT_R8) / (TTT- 7.66D0))
       end if
     end if
   end function phf_foew8

   !--------------------------------------------------------------------------
   ! phf_FODLE8
   !--------------------------------------------------------------------------
   real*8 function phf_FODLE8(TTT)
     !
     ! :Purpose: FONCTION CALCULANT LA DERIVEE SELON T DE  LN EW (OU LN EI)
     !
     implicit none
     real*8 TTT

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       if (TTT > MPC_TRIPLE_POINT_R8) then
         phf_FODLE8 = 4283.76D0/(TTT-30.11D0)**2
       else
         phf_FODLE8 = 6185.90D0/(TTT+ 0.71D0)**2
       end if
     else
       if (TTT > MPC_TRIPLE_POINT_R8) then
         phf_FODLE8 = 4097.93D0/(TTT-35.86D0)**2
       else
         phf_FODLE8 = 5807.81D0/(TTT- 7.66D0)**2
       end if
     end if
   end function phf_FODLE8

   !--------------------------------------------------------------------------
   ! phf_FOQST8
   !--------------------------------------------------------------------------
   real*8 function phf_FOQST8(TTT,PRS)
     !
     ! :Purpose: FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT)
     !
     implicit none
     real*8 TTT,PRS

     phf_FOQST8=MPC_EPS1_R8/(MAX(1.D0,PRS/phf_FOEW8(TTT))-MPC_EPS2_R8)
   end function phf_FOQST8

   !--------------------------------------------------------------------------
   ! phf_FODQS8
   !--------------------------------------------------------------------------
   real*8 function phf_FODQS8(QST,TTT)
     !
     ! :Purpose: FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
     !
     implicit none
     real*8 TTT,QST

     phf_FODQS8=QST*(1.D0+MPC_DELTA_R8*QST)*phf_FODLE8(TTT)
   end function phf_FODQS8
   ! QST EST LA SORTIE DE FOQST

   !--------------------------------------------------------------------------
   ! phf_FOEFQ8
   !--------------------------------------------------------------------------
   real*8 function phf_FOEFQ8(QQQ,PRS)
     !
     ! :Purpose: FONCTION CALCULANT TENSION VAP (EEE) FN DE HUM SP (QQQ) ET PRS
     !
     implicit none
     real*8 QQQ,PRS

     phf_FOEFQ8= MIN(PRS,(QQQ*PRS) / (MPC_EPS1_R8 + MPC_EPS2_R8*QQQ))
   end function phf_FOEFQ8

   !--------------------------------------------------------------------------
   ! FOQFE8
   !--------------------------------------------------------------------------
   real*8 function phf_FOQFE8(EEE,PRS)
     !
     ! :Purpose: FONCTION CALCULANT HUM SP (QQQ) DE TENS. VAP (EEE) ET PRES (PRS)
     !
     implicit none
     real*8 EEE,PRS

     phf_FOQFE8= MIN(1.D0,MPC_EPS1_R8*EEE / (PRS-MPC_EPS2_R8*EEE))
   end function phf_FOQFE8

   !--------------------------------------------------------------------------
   ! phf_FOTVT8
   !--------------------------------------------------------------------------
   real*8 function phf_FOTVT8(TTT,QQQ)
     !
     ! :Purpose: FONCTION CALCULANT TEMP VIRT. (TVI) DE TEMP (TTT) ET HUM SP (QQQ)
     !
     implicit none
     real*8 TTT,QQQ

     phf_FOTVT8= TTT * (1.0D0 + MPC_DELTA_R8*QQQ)
   end function phf_FOTVT8

   !--------------------------------------------------------------------------
   ! phf_FOTTV8
   !--------------------------------------------------------------------------
   real*8 function phf_FOTTV8(TVI,QQQ)
     !
     !:Purpose: FONCTION CALCULANT TTT DE TEMP VIRT. (TVI) ET HUM SP (QQQ)
     !
     implicit none
     real*8 TVI,QQQ

     phf_FOTTV8= TVI / (1.0D0 + MPC_DELTA_R8*QQQ)
   end function phf_FOTTV8

   !--------------------------------------------------------------------------
   ! phf_FOHR8
   !--------------------------------------------------------------------------
   real*8 function phf_FOHR8(QQQ,TTT,PRS)
     !
     !:Purpose: FONCTION CALCULANT HUM REL DE HUM SP (QQQ), TEMP (TTT) ET PRES (PRS).
     !          Where HR = E/ESAT
     !
     implicit none
     real*8 QQQ,TTT,PRS

     phf_FOHR8 = MIN(PRS,phf_FOEFQ8(QQQ,PRS)) / phf_FOEW8(TTT)
   end function phf_FOHR8

   ! LES 5 FONCTIONS SUIVANTES SONT VALIDES DANS LE CONTEXTE OU ON
   ! NE DESIRE PAS TENIR COMPTE DE LA PHASE GLACE DANS LES CALCULS
   ! DE SATURATION.
  
   !--------------------------------------------------------------------------
   ! phf_FOEWA8
   !--------------------------------------------------------------------------
   real*8 function phf_FOEWA8(TTT)
     !
     !:Purpose: FONCTION DE VAPEUR SATURANTE (TETENS)
     !
     implicit none
     real*8 TTT

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       phf_FOEWA8=610.94D0*EXP(17.625D0*(TTT-MPC_TRIPLE_POINT_R8)/(TTT-30.11D0))
     else
       phf_FOEWA8=610.78D0*EXP(17.269D0*(TTT-MPC_TRIPLE_POINT_R8)/(TTT-35.86D0))
     end if
   end function phf_FOEWA8

   !--------------------------------------------------------------------------
   ! phf_FODLA8
   !--------------------------------------------------------------------------
   real*8 function phf_FODLA8(TTT)
     !
     !:Purpose: FONCTION CALCULANT LA DERIVEE SELON T DE LN EW
     !
     implicit none
     real*8 TTT

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       phf_FODLA8=17.625D0*(MPC_TRIPLE_POINT_R8-30.11D0)/(TTT-30.11D0)**2
     else
       phf_FODLA8=17.269D0*(MPC_TRIPLE_POINT_R8-35.86D0)/(TTT-35.86D0)**2
     end if
   end function phf_FODLA8

   !--------------------------------------------------------------------------
   ! FOQSA8
   !--------------------------------------------------------------------------
   real*8 function phf_FOQSA8(TTT,PRS)
     !
     !:Purpose: FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE
     !
     implicit none
     real*8 TTT,PRS

     phf_FOQSA8=MPC_EPS1_R8/(MAX(1.D0,PRS/phf_FOEWA8(TTT))-MPC_EPS2_R8)
   end function phf_FOQSA8

   !--------------------------------------------------------------------------
   ! phf_FODQA8
   !--------------------------------------------------------------------------
   real*8 function phf_FODQA8(QST,TTT)
     !
     !:Purpose: FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
     !
     implicit none
     real*8 QST,TTT

     phf_FODQA8=QST*(1.D0+MPC_DELTA_R8*QST)*phf_FODLA8(TTT)
   end function phf_FODQA8

   !--------------------------------------------------------------------------
   ! phf_FOHRA8
   !--------------------------------------------------------------------------
   real*8 function phf_FOHRA8(QQQ,TTT,PRS)
     !
     !:Purpose: FONCTION CALCULANT L'HUMIDITE RELATIVE
     !
     implicit none
     real*8 QQQ,TTT,PRS

     phf_FOHRA8=MIN(PRS,phf_FOEFQ8(QQQ,PRS))/phf_FOEWA8(TTT)
   end function phf_FOHRA8

   ! LES 6 FONCTIONS SUIVANTES SONT REQUISES POUR LA TEMPERATURE
   ! EN FONCTION DE LA TENSION DE VAPEUR SATURANTE.
   ! (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)

   !--------------------------------------------------------------------------
   ! phf_FOTW8
   !--------------------------------------------------------------------------
   real*8 function phf_FOTW8(EEE)
     !
     !:Purpose: FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
     !          SATURANTE PAR RAPPORT A EW.
     !
     implicit none
     real*8 EEE

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018a' .or. trim(saturationCurve) == 'Tetens_2018') then
       phf_FOTW8=(30.11D0*LOG(EEE/610.94D0)-17.625D0*MPC_TRIPLE_POINT_R8)/ &
           (LOG(EEE/610.94D0)-17.625D0)
     else
       phf_FOTW8=(35.86D0*LOG(EEE/610.78D0)-17.269D0*MPC_TRIPLE_POINT_R8)/ &
           (LOG(EEE/610.78D0)-17.269D0)
     end if
   end function phf_FOTW8

   !--------------------------------------------------------------------------
   ! phf_FOTI8
   !--------------------------------------------------------------------------
   real*8 function phf_FOTI8(EEE)
     !
     ! :Purpose: FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
     !          SATURANTE PAR RAPPORT A EI.
     !
     implicit none
     real*8 EEE

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       phf_FOTI8=(-0.71D0*LOG(EEE/610.94D0)-22.587D0*MPC_TRIPLE_POINT_R8)/ &
           (LOG(EEE/610.94D0)-22.587D0)
     else
       phf_FOTI8=( 7.66D0*LOG(EEE/610.78D0)-21.875D0*MPC_TRIPLE_POINT_R8)/ &
           (LOG(EEE/610.78D0)-21.875D0)
     end if
   end function phf_FOTI8

   !--------------------------------------------------------------------------
   ! phf_FODTH8
   !--------------------------------------------------------------------------
   real*8 function phf_FODTW8(TTT,EEE)
     !
     ! :Purpose: FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
     !           VAPEUR SATURANTE (EW).
     !
     implicit none
     real*8 TTT,EEE

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018a' .or. trim(saturationCurve) == 'Tetens_2018') then
       phf_FODTW8=(30.11D0-TTT)/EEE/(LOG(EEE/610.94D0)-17.625D0)
     else
       phf_FODTW8=(35.86D0-TTT)/EEE/(LOG(EEE/610.78D0)-17.269D0)
     endif
   end function phf_FODTW8

   !--------------------------------------------------------------------------
   ! phf_FODTI8
   !--------------------------------------------------------------------------
   real*8 function phf_FODTI8(TTT,EEE)
     !
     ! :Purpose: FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
     !           VAPEUR SATURANTE (EI).
     !
     implicit none
     real*8 TTT,EEE

     if (.not.phf_initialized) call phf_tetens_coefs_switch

     if (trim(saturationCurve) == 'Tetens_2018') then
       phf_FODTI8=(-0.71D0-TTT)/EEE / (LOG(EEE/610.94D0)-22.587D0)
     else
       phf_FODTI8=( 7.66D0-TTT)/EEE / (LOG(EEE/610.78D0)-21.875D0)
     end if
   end function phf_FODTI8

   !--------------------------------------------------------------------------
   ! phf_FOTWI8
   !--------------------------------------------------------------------------
   real*8 function phf_FOTWI8(TTT,EEE)
     !
     ! :Purpose: FONCTION DE L'AJUSTEMENT DE LA TEMPERATURE.
     !
     implicit none
     real*8 TTT,EEE

     phf_FOTWI8=MAX(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*phf_FOTW8(EEE)  &
               -MIN(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*phf_FOTI8(EEE)
   end function phf_FOTWI8

   !--------------------------------------------------------------------------
   ! phf_FODWTI8
   !--------------------------------------------------------------------------
   real*8 function phf_FODTWI8(TTT,EEE)
     !
     ! :Purpose: FONCTION DE L'AJUSTEMENT DE LA DERIVEE DE LA TEMPERATURE.
     !
     implicit none
     real*8 TTT,EEE

     phf_FODTWI8=MAX(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*phf_FODTW8(TTT,EEE)  &
                -MIN(0.0D0,SIGN(1.0D0,TTT-MPC_TRIPLE_POINT_R8))*phf_FODTI8(TTT,EEE)
   end function phf_FODTWI8

  ! LES 7 FONCTIONS SUIVANTES POUR EW-EI SONT REQUISES POUR LES MODELES
  ! CMAM ET CGCM. (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)
  
  !--------------------------------------------------------------------------
  ! phf_FOEW8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOEW8_CMAM(TTT)  
        !
        !:Purpose: FONCTION DE TENSION DE VAPEUR SATURANTE - EW
        !
        implicit none
        real*8 TTT
        phf_FOEW8_CMAM= 100.D0*EXP(21.656D0-5418.D0/TTT)
  end function phf_FOEW8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOEI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOEI8_CMAM(TTT)  
        !
        !:Purpose: FONCTION DE TENSION DE VAPEUR SATURANTE - EI
        !
        implicit none
        real*8 TTT
        phf_FOEI8_CMAM= 100.D0*EXP(24.292D0-6141.D0/TTT)
  end function phf_FOEI8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOERAT8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOERAT8_CMAM(TTT) 
        !
        !:Purpose: FONCTION DE LA PROPORTION DE LA CONTRIBUTION DE EW VS EI
        !
        implicit none
        real*8 TTT
        phf_FOERAT8_CMAM=MIN(1.0D0,0.0059D0+0.9941D0*EXP(-0.003102D0*  &
            MIN(0.0D0,TTT-MPC_TRIPLE_POINT_R8)**2))
  end function phf_FOERAT8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOEWI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOEWI8_CMAM(TTT)  
        !
        !:Purpose: FONCTION DE TENSION DE VAPEUR SATURANTE RESULTANTE - EW et EI
        !
        implicit none
        real*8 TTT
        phf_FOEWI8_CMAM = phf_FOEW8_CMAM(TTT)*phf_FOERAT8_CMAM(TTT)  &
               +(1.0D0-phf_FOERAT8_CMAM(TTT))*phf_FOEI8_CMAM(TTT)
  end function phf_FOEWI8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FODLE8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FODLE8_CMAM(TTT)  
        !
        !:Purpose: FONCTION DE LA DERIVE DE LN(E) PAR RAPPORT A LA TEMPERATURE
        !
        implicit none
        real*8 TTT
        phf_FODLE8_CMAM= phf_FOERAT8_CMAM(TTT)*5418.D0/TTT/TTT  &
             +(1.0D0-phf_FOERAT8_CMAM(TTT))*6141.D0/TTT/TTT
  end function phf_FODLE8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOQST8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOQST8_CMAM(TTT,PRS) 
        !
        !:Purpose: FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT).
        !          PREND EN COMPTE LES PHASES GLACE ET EAU.
        !
        implicit none
        real*8 TTT,PRS
        phf_FOQST8_CMAM=MPC_EPS1_R8/(MAX(1.0D0,PRS/  &
                           phf_FOEWI8_CMAM(TTT))-MPC_EPS2_R8)
  end function phf_FOQST8_CMAM

  ! LES 6 FONCTIONS SUIVANTES SONT REQUISES POUR LA TEMPERATURE
  ! EN FONCTION DE LA TENSION DE VAPEUR SATURANTE POUR CMAM/CGCM
  ! (AJOUTE PAR YVES J. ROCHON, ARQX/SMC, JUIN 2004)

  !--------------------------------------------------------------------------
  ! phf_FOTW8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOTW8_CMAM(EEE) 
        !
        !:Purpose: FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
        !          SATURANTE PAR RAPPORT A EW.
        !
        implicit none
        real*8 EEE
        phf_FOTW8_CMAM=5418.D0/(21.656D0-LOG(EEE/100.0D0))
  end function phf_FOTW8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOTI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOTI8_CMAM(EEE) 
        !
        !:Purpose: FONCTION DE LA TEMPERATURE EN FONCTION DE LA TENSION DE VAPEUR
        !          SATURANTE PAR RAPPORT A EI.
        !
        implicit none
        real*8 EEE
        phf_FOTI8_CMAM=6141.D0/(24.292D0-LOG(EEE/100.0D0))
  end function phf_FOTI8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FODTW8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FODTW8_CMAM(TTT,EEE) 
        !
        !:Purpose: FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
        !          VAPEUR SATURANTE (EW).
        !
        implicit none
        real*8 TTT,EEE
        phf_FODTW8_CMAM=TTT/EEE/(21.656D0-LOG(EEE/100.0D0))
  end function phf_FODTW8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FODTI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FODTI8_CMAM(TTT,EEE) 
        !
        !:Purpose: FONCTION DE LA DERIVE DE LA TEMPERATURE EN FONCTION DE LA TENSION DE
        !          VAPEUR SATURANTE (EI).
        !
        implicit none
        real*8 TTT,EEE
        phf_FODTI8_CMAM=TTT/EEE/(24.292D0-LOG(EEE/100.0D0))
  end function phf_FODTI8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FOTWI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FOTWI8_CMAM(TTT,EEE) 
        !
        !:Purpose: FONCTION DE L'AJUSTEMENT DE LA TEMPERATURE.
        !
        implicit none
        real*8 TTT,EEE
        phf_FOTWI8_CMAM= phf_FOERAT8_CMAM(TTT)*phf_FOTW8_CMAM(EEE)+ &
                   (1.0D0-phf_FOERAT8_CMAM(TTT))*phf_FOTI8_CMAM(EEE)
  end function phf_FOTWI8_CMAM

  !--------------------------------------------------------------------------
  ! phf_FODTWI8_CMAM
  !--------------------------------------------------------------------------
  real*8 function phf_FODTWI8_CMAM(TTT,EEE) 
        !
        !:Purpose: FONCTION DE L'AJUSTEMENT DE LA DERIVEE DE LA TEMPERATURE.
        !
        implicit none
        real*8 TTT,EEE
        phf_FODTWI8_CMAM= phf_FOERAT8_CMAM(TTT)*phf_FODTW8_CMAM(TTT,EEE)+  &
                          (1.0D0-phf_FOERAT8_CMAM(TTT))*phf_FODTI8_CMAM(TTT,EEE)
  end function phf_FODTWI8_CMAM


  !--------------------------------------------------------------------------
  ! phf_FOBRANCH
  !--------------------------------------------------------------------------
  real*8 function phf_FQBRANCH(QQQ)  
        !
        !:Purpose: function returning 0/1 depending on the minimum q branch condition
        !          as discussed by Brunet (1996) to prevent getting a vapour pressure that exceeds
        !          the total pressure p when q exceeds 1.
        !
        implicit none
        real*8 QQQ
        phf_FQBRANCH= 0.5D0+SIGN(0.5D0,1.D0-(QQQ))
  end function phf_FQBRANCH

  !=============================================================================
  !
  !     TLM  of  THERMODYNAMIC FUNCTIONS USED IN 3DVAR
  !     CONSTANTS FROM COMMON /CTESDYN/
  !     NOTE: ALL FUNCTIONS WORK IN  S.I. UNITS
  !           I.E PRS IN PA, QQQ IN KG/KG 
  
  !--------------------------------------------------------------------------
  ! phf_FOEFQL
  !--------------------------------------------------------------------------
   real*8 function phf_FOEFQL(QQL,PRSL,QQQ,PRS,PNETA)  
        !
        !:Purpose: TLM  OF FUNCTION CALCULATING VAPOUR PRESSURE
        !
        !          INPUTS:
        !  
        !          * QQL ,  PERTURBATION OF LN SPECIFIC HUM
        !          * PRSL ,   PERTURBATION OF SURFACE PRESSURE
        !          * QQQ   ,  SPECIFIC HUMIDITY
        !          * PRS   , PRESSURE
        !          * PNETA   , VALUE OF ETA LEVEL
        !
        !          OUTPUT: 
        !
        !          * FOEFQL,  PERTURBATION  OF VAPOUR PRESSURE
        !
        implicit none
        real*8 QQL,PRSL,QQQ,PRS,PNETA
        phf_FOEFQL= phf_FQBRANCH(QQQ)  &
           * ((QQL*MPC_EPS1_R8*PRS*QQQ/(MPC_EPS1_R8+MPC_EPS2_R8*QQQ)  &
           +  PRSL*PNETA*QQQ)/(MPC_EPS1_R8+MPC_EPS2_R8*QQQ))  &
           + (1.0D0 - phf_FQBRANCH(QQQ))*PRSL*PNETA
  end function phf_FOEFQL

  !--------------------------------------------------------------------------
  ! phf_FOTVVL
  !--------------------------------------------------------------------------
   real*8 function phf_fotvvl(qqq,ttt,ttl,plnql) 
        !
        !:Purpose: Tangent-linear operator of virtual temperature 
        !
        implicit none
        real*8 qqq   ! backgroud specific humidity
        real*8 ttt   ! backgroud temperature
        real*8 ttl   ! temperature increment
        real*8 plnql ! increment of logarithm specific humidity  (del(ln q))
        phf_fotvvl=(1 + MPC_DELTA_R8*qqq)*ttl + MPC_DELTA_R8*qqq*ttt*plnql
  end function phf_fotvvl

  !=============================================================================
  !
  !   DEFINITION OF ADJOINTS OF THERMODYNAMIC  FUNCTIONS
  !   CONSTANTS AS IN COMMON /CTESDYN/
  !     NOTE: ALL  UNITS S.I.
  !           I.E. PADES IN DEG K, PRS EN PA, QQQ EN KG/KG
  !

  !--------------------------------------------------------------------------
  ! phf_FOEFQA
  !--------------------------------------------------------------------------
  real*8 function phf_FOEFQA(PADES,PGAMMA,QQQ,PRS)  
        !
        !:Purpose: ADJOINT OF LN SPECIFIC  HUM (QQQ) DUE TO DEWPOINT DEPRESSION CORRECTIONS
        !
        implicit none
        real*8 PADES  ! ADJOINT OF DEWPOINT DEPRESSION
        real*8 PGAMMA ! ADOINT OF VAPOUR PRESSURE RELATIONSHIP
        real*8 QQQ    ! SPECIFIC HUMIDITY
        real*8 PRS    ! PRESSURE
        phf_FOEFQA= PADES*PGAMMA*MPC_EPS1_R8*PRS*QQQ/  &
                    ((MPC_EPS1_R8+MPC_EPS2_R8*QQQ)*(MPC_EPS1_R8+MPC_EPS2_R8*QQQ))
  end function phf_FOEFQA

  !--------------------------------------------------------------------------
  ! phf_FOEFQPSA
  !--------------------------------------------------------------------------
  real*8 function phf_FOEFQPSA(PADES,PGAMMA,QQQ,PNETA) 
        !
        !:Purpose: ADJOINT OF SURFACE PRESSURE  DUE TO DEWPOINT DEPRESSION CORRECTIONS
        !
        implicit none
        real*8 PADES  ! ADJOINT OF DEWPOINT DEPRESSION
        real*8 PGAMMA ! ADOINT OF VAPOUR PRESSURE RELATIONSHIP
        real*8 QQQ    ! SPECIFIC HUMIDITY
        real*8 PNETA  ! VALUE OF NETA
        phf_FOEFQPSA = PADES*PGAMMA*QQQ*PNETA/  &
                       (MPC_EPS1_R8+MPC_EPS2_R8*QQQ)
  end function phf_FOEFQPSA

  !--------------------------------------------------------------------------
  ! phf_FOTTVa
  !--------------------------------------------------------------------------
  real*8 function phf_fottva(qqq,tva)  
        !
        !:Purpose: Adjoint of temperature due to virtual temperature correction
        !
        implicit none
        real*8 qqq ! background specific humidity
        real*8 tva ! adjoint variable of virtual temperature
        phf_fottva= (1D0 + MPC_DELTA_R8*qqq)*tva
  end function phf_fottva

  !--------------------------------------------------------------------------
  ! phf_FOLNQVA
  !--------------------------------------------------------------------------
  real*8 function phf_folnqva(qqq,ttt,tva) 
        !
        !:Purpose:  Adjoint of logarithm of specific humidity due to virtual temperature correction
        !
        implicit none
        real*8 qqq ! background specific humidity
        real*8 ttt ! background temperature
        real*8 tva ! adjoint variable of virtual temperature
        phf_folnqva = MPC_DELTA_R8*qqq*ttt*tva
  end function phf_folnqva

  !-------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! PHF_CONVERT_Z_TO_PRESSURE
  !--------------------------------------------------------------------------
  function  phf_convert_z_to_pressure(altitude,rgz_mod,press_mod,nlev,nlev_mod,lat,success) result(press)
    !          
    !:Purpose: Converts an array of (geometric) altitudes to pressures. Uses linear interpolation
    !          in log(p).  
    implicit none

    real(8), intent(in) :: altitude(nlev)      ! altitudes to convert to pressures (m)
    real(8), intent(in) :: rgz_mod(nlev_mod)   ! geopotential heights on model levels (m), assumed to be in decending order
    real(8), intent(in) :: press_mod(nlev_mod) ! pressure on model levels, assumed to be in ascending order
    real(8), intent(in) :: lat                 ! latitude (rad)
    integer, intent(in) :: nlev                ! length of altitude array
    integer, intent(in) :: nlev_mod            ! number of model levels
    logical, intent(inout) :: success(nlev)    ! flag indicating success
    real(8) :: press(nlev)                     ! converted pressures

    real(8) :: rgz(nlev)
    integer :: ilev,ilev_mod

    ! Check model pressures
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

  !--------------------------------------------------------------------------
  ! PHF_CONVERT_Z_TO_GZ
  !--------------------------------------------------------------------------
  function phf_convert_z_to_gz(altitude,lat,nlev) result(rgz)
    !          
    !:Purpose: Converts altitudes to geopotential heights. Uses the Helmert formula to
    !          parameterize the latitude dependence and uses analytical result of the
    !          integral of \int g(z)dz for the altitude dependence (see J.A. Dutton 1976,
    !          p.65). At an altitude of 50 km, the altitude and geopotential height
    !          differ by around 0.2-0.5 km, depending on the latitude.
    !
    implicit none

    real(8), intent(in) :: altitude(nlev) ! altitudes (m)
    real(8), intent(in) :: lat            ! latitude (rad)
    integer, intent(in) :: nlev           ! number of levels
    real(8) :: rgz(nlev)                  ! geopotential heights (m)

    rgz = (ec_rg/9.8) * (1.-2.64D-03*cos(2.*lat)+5.9D-6*cos(2.*lat)**2) * ec_ra*altitude/(ec_ra+altitude)

  end function phf_convert_z_to_gz

  !--------------------------------------------------------------------------
  ! PHF_GET_TROPOPAUSE
  !--------------------------------------------------------------------------
  function phf_get_tropopause(nmodlev,pressmod,tt,height,hu_opt) result(tropo_press)
    !          
    !:Purpose: Determines pressure level of tropopause. 
    !          Final tropopause is taken as max pressure (lowest altitude) from the
    !          water vapour and temperature based tropopauses.
    !
    implicit none

    integer, intent(in) :: nmodlev                   ! Number of model levels
    real(8), intent(in) :: pressmod(nmodlev)         ! Model pressure array (Pa)
    real(8), intent(in) :: tt(nmodlev)               ! Model temperature (Kelvin)
    real(8), intent(in) :: height(nmodlev)           ! Model height (m)
    real(8), intent(in), optional :: hu_opt(nmodlev) ! Model specific humidity
    real(8) :: tropo_press                           ! Tropopause level in Pa (output)
  
    integer :: itop,i,k,ilaps
    real(8) :: hu_ppmv1,hu_ppmv2,hu_ppmv3,xlaps,tropo_press_hu
    real(8), parameter :: press_min=6000.         ! Min tropoause pressure 60 hPa.; equivalent to ~ 20km
    real(8), parameter :: height_min=6000.0           ! Min tropopause level in meters.
    real(8), parameter :: ppmv_threshold=10.0     
    real(8), parameter :: tgrad_threshold=0.002   ! degrees/m (2 degrees/km)
    real(8), parameter :: consth=0.160754938e+07  ! conversion from mass mixing ratio to ppmv;  1.0e+06 / (18.015/28.96)

    tropo_press=-1.0
    if (all(height.lt.0.0))  call utl_abort('phf_get_tropopause: Missing height for determining tropopause pressure')

    ! Initialize tropopause pressure level using temperature gradient.
    ! Thermal tropopause is defined as the lowest level (above height_min) at which (1) the lapse rate decreases
    ! to <= 2 C/km and (2) the average lapse rate between this level and all higher levels within 2 km are <= 2 C/km. 
    ! Ref: International Meteorological Vocabulary (2nd ed.). Geneva: Secretariat of the World Meteorological
    !      Organization. 1992. p. 636. ISBN 92-63-02182-1.
    ! The second requirement, based on hu, may give levels that are to high (pressure too low) in the winter hemisphere.

    do itop=3,nmodlev
       if (pressmod(itop).ge.press_min) exit
    end do
    itop=itop-1
       
    do i=nmodlev,itop+1,-1
       if (height(i)-height(nmodlev).lt.height_min) cycle
       xlaps=-(tt(i)-tt(i-1))/(height(i)-height(i-1))
       if (xlaps.le.tgrad_threshold) then
          ilaps=1
          do k=i-1,itop,-1
             if (height(k)-height(i).gt.2000.0) exit
             xlaps=xlaps-(tt(k)-tt(k-1))/(height(k)-height(k-1))
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
                if (height(i)-height(k).gt.5000.0) exit
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

  !--------------------------------------------------------------------------
  ! PHF_GET_PBL
  !--------------------------------------------------------------------------
  function phf_get_pbl(nmodlev,pressmod,tt,height,hu_opt,uu_opt,vv_opt) result(pbl_press)
    ! 
    !:Purpose: Determines pressure level of planetary boundary layer using 
    !          a first threshold of 0.5 for the bulk Richadson number (after Mahrt, 1981; 
    !          requires availability of uu and vv). Threshold reduced to largest value
    !          between 0.25 and 0.5 if first not satisfied. 
    !
    !          If not found with this approach, applies a variant of the Heffter approach 
    !          described in Aliabadi et al (2016), with some local variation.
    !
    !          References:
    !
    !          * Aliabadi A.A., R.M. Staebler, J. de Grandpre, A. Zadra, and P.A.
    !            Vaillancourt, 2016: Comparison of Estimated Atmospheric
    !            Boundary Layer Mixing Height in teh Arctic and Southern Great
    !            Plains under Statisticallt Stable Conditions: Experimental
    !            and Numerical Aspects, Submitted to Atmosphere-Ocean (2015).
    !
    !          * Mahrt, L. 1981: Modelling depth of the stable boundary-layer,
    !            Bound-Lay. Meteorol., 21, 3-19
    !
    !          * Heffter, J.L.,1980: Transport layer depth calculations, Second 
    !            Joint Conference on Applications of Air Pollution Meteorology,
    !            New Orleans, LA, 24-27 March 1980. American Meteorological
    !            Society, Boston, MA.
    !       
    !          Comments: Currently assumes (uu,vv) midlayer levels approximately 
    !          at tt, height, and hu levels when size(uu).ne.nmodlev.
    !
    implicit none

    integer, intent(in) :: nmodlev           ! Number of model levels for variables other than uu and vv
    real(8), intent(in) :: pressmod(nmodlev) ! Model pressure array (Pa)
    real(8), intent(in) :: tt(nmodlev)       ! Model temperature (Kelvin)
    real(8), intent(in) :: height(nmodlev)   ! Model height (meters)
    real(8), optional :: uu_opt(:)           ! Model zonal wind component (m/s)
    real(8), optional :: vv_opt(:)           ! Model meridional wind component (m/s)
    real(8), optional :: hu_opt(nmodlev)     ! Specific humidity
    real(8) :: pbl_press                     ! PBL level in Pa (output)
  
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
    
    if (all(height.lt.0.0))  call utl_abort('phf_get_pbl: Missing height for determining PBL pressure')

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
    zs = height(i)*0.001
    
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
       
!      us,vs set to 0.0
!       us=0.0
!       vs=0.0
 
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
         
           RiB2 = ec_rg * (thetavh(i)-thetavs) * (height(i)*0.001-zs) / (thetavs*uv)
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
          if ((thetavh(i)-thetavh(i+1))/(height(i)-height(i+1)).ge.0.005) then
             ! Near bottom of inversion layer found
             inv=i+1
             i=i-1
             !do while (thetavh(i)-thetavh(i+1).gt.0.0.and.i.gt.itop) 
             do while ((thetavh(i)-thetavh(i+1))/(height(i)-height(i+1)).ge.0.005.and.i.gt.itop) 
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
             grad=(thetavh(i)-thetavh(i+1))/(height(i)-height(i+1))
             if (gradmax.lt.grad) then
                gradmax=(thetavh(i)-thetavh(i+1))/(height(i)-height(i+1))
                igradmax=i
                if (grad.ge.0.005) then ! Check next layer as well
                   if ((thetavh(i-1)-thetavh(i))/(height(i-1)-height(i)).ge.0.005) exit
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
  ! phf_calcDistance
  !--------------------------------------------------------------------------
  function phf_calcDistance(lat2, lon2, lat1, lon1) result(distanceInM)
    !
    !:Purpose: Compute the distance between two point on earth: (lat1,lon1) 
    !          and (lat2,lon2). Calcul utilisant la Formule d'Haversine.
    !
    !          Reference: R.W. Sinnott,'Virtues of Haversine',Sky and 
    !          Telescope, vol.68, no.2, 1984, p.159)
    !
    implicit none

    ! arguments
    real(8) :: lat1, lon1, lat2, lon2
    real(8) :: distanceInM

    ! locals
    real(8) :: dlat, dlon, a, c

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (sin(dlat/2.d0))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2.d0))**2
    c = 2.d0 * atan2(sqrt(a),sqrt(1.d0-a))
    distanceInM = ec_ra * c

  end function phf_calcDistance

  !--------------------------------------------------------------------------
  ! phf_calcDistanceFast
  !--------------------------------------------------------------------------
  function phf_calcDistanceFast(lat2, lon2, lat1, lon1) result(distanceInM)
    !
    !:Purpose: Compute the distance between two point on earth: (lat1,lon1)
    !          and (lat2,lon2). Using a quick and dirty formula good for 
    !          short distances not close to the pole.
    !
    implicit none

    ! arguments
    real(8) :: lat1, lon1, lat2, lon2
    real(8) :: distanceInM

    ! locals
    real(8) :: dlat, dlon

    dlon = (lon2 - lon1)*cos(lat1)
    dlat = lat2 - lat1

    distanceInM = ec_ra * sqrt(dlon*dlon + dlat*dlat)

  end function phf_calcDistanceFast

  !--------------------------------------------------------------------------
  ! PHF_GRAVITYSRF
  !--------------------------------------------------------------------------
  function phf_gravitysrf(sLat)
    !
    !:Purpose: Normal gravity on ellipsoidal surface
    !
    implicit none
    real(8), intent(in)  :: sLat           ! sin of latitude
    real(8)              :: phf_gravitysrf ! normal gravity (m/s2)
    
    real(8)              :: ks2
    real(8)              :: e2s

    ks2 = ec_wgs_TNGk * sLat*sLat
    e2s = 1.D0 - ec_wgs_e2 * sLat*sLat
    phf_gravitysrf = ec_wgs_GammaE * (1.D0 + ks2) / sqrt(e2s)
  end function phf_gravitysrf

  !--------------------------------------------------------------------------
  ! PHF_GRAVITYALT
  !--------------------------------------------------------------------------
  function phf_gravityalt(sLat, Altitude)
    !
    !:Purpose: Normal gravity above the ellipsoidal surface
    !
    implicit none
    real(8), intent(in)  :: sLat           ! sin of latitude
    real(8), intent(in)  :: Altitude       ! altitude (m)
    real(8)              :: phf_gravityalt ! Normal gravity (m/s2)

    real(8)              :: C1
    real(8)              :: C2

    C1 =-2.D0/ec_wgs_a*(1.D0+ec_wgs_f+ec_wgs_m-2*ec_wgs_f*sLat*sLat)
    C2 = 3.D0/ec_wgs_a**2
    phf_gravityalt = phf_gravitysrf(sLat)*                                   &
         (1.D0 + C1 * Altitude + C2 * Altitude**2)
  end function phf_gravityalt

  !--------------------------------------------------------------------------
  ! PHF_HEIGHT2GEOPOTENTIAL
  !--------------------------------------------------------------------------
  subroutine phf_height2geopotential(altitude, latitude, geopotential, &
                                     printHeight)
    !
    !:Purpose: Geopotential energy at a given point.
    !          Result is based on the WGS84 approximate expression for the
    !          gravity acceleration as a function of latitude and altitude,
    !          integrated with the trapezoidal rule.
    !
    implicit none
    real(8), intent(in)   :: altitude(:)     ! altitude (m)
    real(8), intent(in)   :: latitude        ! latitude (rad)
    real(8), intent(inout):: geopotential(:) ! geopotential (m2/s2)

    integer           :: nlev, nlev500m
    real(8), allocatable :: alt500m(:), gravity500m(:)
    real(8)           :: delAlt, aveGravity, sLat, gravity, gravityM1
    integer           :: i, ilev 
    logical, optional :: printHeight
    

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
      ! At surface, use local gravity to get height
      gravity = phf_gravityalt(sLat,altitude(nlev))
      geopotential(nlev) = altitude(nlev) * gravity
      gravityM1 = gravity

    endif

    ! At upper-levels, integrate on model levels to get height 
    do ilev = nlev-1, 1, -1
      gravity = phf_gravityalt(sLat,altitude(ilev))
      aveGravity = 0.5D0 * (gravity + gravityM1)
      delAlt = altitude(ilev) - altitude(ilev+1)
      geopotential(ilev) = geopotential(ilev+1) + delAlt * aveGravity
      gravityM1 = gravity
    enddo

    if ( present(printHeight) ) then
      if ( printHeight ) then
        write(*,*) 'phf_height2geopotential, Z_T:'
        write(*,*) geopotential(:)

        printHeight = .false.
      endif
    endif

  end subroutine phf_height2geopotential
  
end module physicsFunctions_mod


