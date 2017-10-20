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
!! MODULE LQtoES (prefix="LQtoES")
!!
!! *Purpose*: To compute dewpoint depression (ES) from log(humidity) (LQ) and
!!            also tangent linear and adjoint versions of this operation.
!!
!--------------------------------------------------------------------------
module LQtoES_mod
  use MathPhysConstants_mod
  use physicsFunctions_mod
  implicit none
  save
  private
  
  ! public procedures
  public :: LQtoES, LQtoES_tl, LQtoES_ad


  contains

  real*8 function LQtoES(lq,tt,pressure)
    !
    ! Purpose:
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    real*8, intent(in) :: lq,tt,pressure
    real*8 :: husat,td
    !
    ! get the saturated vapor pressure from lq (log of specific humidity)
    !
    husat = foefq8(max(exp(lq),1.d-12),pressure)
    !
    ! now the dewpoint temperature
    !
    td = fotw8(husat)
    !
    ! finally the dewpoint depression
    !
    LQtoES = min(tt-td,MPC_MAXIMUM_ES_R8)

  end function LQtoES


  real*8 function LQtoES_tl(delLQ,delTT,delP0,LQ_g,PRES_g,dPdPsfc)
    !
    ! Purpose: TLM VERSION
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    REAL(8) :: delLQ,delTT,delP0,LQ_g,PRES_g,dPdPsfc
    REAL(8) :: ZE,ZTD,dTDdE,ZQBRANCH,HU_g
    REAL(8) :: dESdLQ,dESdTT,dESdP0

    dESdTT = 1.0d0
    !
    ! Forward calculations of saturation vapour pressure and dewpoint temperature
    ! and adjoint of vapour pressure from adjoint of dewpoint temperature
    !
    HU_g = exp(LQ_g)
    ZE = FOEFQ8(HU_g, PRES_g)

    ZTD=FOTW8(ZE)
    dTDdE=FODTW8(ZTD,ZE)
    !
    ! adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    !
    ZQBRANCH = FQBRANCH(HU_g)
    dESdLQ = - ZQBRANCH*FOEFQA(1.0d0,dTDdE,HU_g,PRES_g)

    dESdP0 = - ZQBRANCH*FOEFQPSA(1.0d0,dTDdE,HU_g,dPdPsfc)-  &
               (1.D0-ZQBRANCH)*(dTDdE*dPdPsfc)

    LQtoES_tl =  dESdLQ*delLQ + dESdP0*delP0 + dESdTT*delTT

  end function LQtoES_tl


  subroutine LQtoES_ad(delLQ,delTT,delP0,delES,LQ_g,PRES_g,dPdPsfc)
    !
    ! Purpose: ADJOINT VERSION
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    REAL*8, intent(inout) :: delLQ,delTT,delP0
    REAL*8, intent(in)  :: delES,LQ_g,PRES_g,dPdPsfc
    REAL*8 ZE,ZTD,dTDdE,ZQBRANCH,HU_g
    REAL*8 dESdLQ,dESdTT,dESdP0

    dESdTT = 1.0d0
    !
    ! Forward calculations of saturation vapour pressure and dewpoint temperature
    ! and adjoint of vapour pressure from adjoint of dewpoint temperature
    !
    HU_g = exp(LQ_g)
    ZE = FOEFQ8(HU_g, PRES_g)

    ZTD=FOTW8(ZE)
    dTDdE=FODTW8(ZTD,ZE)
    !
    ! adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    !
    ZQBRANCH = FQBRANCH(HU_g)
    dESdLQ = - ZQBRANCH*FOEFQA(1.0d0,dTDdE,HU_g,PRES_g)

    dESdP0 = - ZQBRANCH*FOEFQPSA(1.0d0,dTDdE,HU_g,dPdPsfc)-  &
               (1.D0-ZQBRANCH)*(dTDdE*dPdPsfc)

    ! TLM: delES =  dESdLQ*delLQ + dESdP0*delP0 + dESdTT*delTT
    ! ADJOINT:
    delLQ = delLQ + dESdLQ*delES
    delP0 = delP0 + dESdP0*delES
    delTT = delTT + dESdTT*delES

  end subroutine LQtoES_ad


end module LQtoES_mod