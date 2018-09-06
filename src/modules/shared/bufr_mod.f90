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
!!
!! MODULE bufr_mod
!!
!! This module contains information that is needed from the
!! Binary Universal Format for the Representation of meteorological data, 
!! maintained by the World Meteorological Organization
!!
!! Public routines:
!!v
!!v       - "buf_IsAtmosConstituent" determines if observation is
!!v         constituent/chemistrry data
!!v
!!v       - "buf_IsIntegral" determines if vertical
!!v         integral measurement.
!!v
!
! Revisions:
!           Y.J. Rochon, ARQI/AQRD, Feb 2015, March 2017
!           - Removed BUFR_NEOZ
!           - Removed BUFR_NECH
!           - Added BUFR_SCALE_EXPONENT, BUFR_UNITS_*
!           - Added public bufr_IsAtmosConstituent, bufr_IsIntegral
!           M. Sitwell, ARQI/AQRD, Feb 2015
!           - Replaced BUFR_NETR_* by BUFR_NECH_*  
!       
!------------------------------------------------------------------------
module bufr_mod

 public

 ! Table B:
 !          Universal Field-Identity Numbers
 !          (At the CMC, these have also been known as burp_id.)
 !
 ! Use the program codbuf to find out the meaning of these IDs.
 !

 integer, parameter :: BUFR_NETT=12001
 integer, parameter :: BUFR_NEUU=11003
 integer, parameter :: BUFR_NEVV=11004
 integer, parameter :: BUFR_NEGZ=10194
 integer, parameter :: BUFR_NEES=12192
 integer, parameter :: BUFR_NEFF=11002
 integer, parameter :: BUFR_NEDD=11001
 integer, parameter :: BUFR_NEUS=11215
 integer, parameter :: BUFR_NEVS=11216
 integer, parameter :: BUFR_NETS=12004
 integer, parameter :: BUFR_NESS=12203
 integer, parameter :: BUFR_NEFS=11012
 integer, parameter :: BUFR_NEDS=11011
 integer, parameter :: BUFR_NEDZ=10192
 integer, parameter :: BUFR_NEPP=07004
 integer, parameter :: BUFR_NEPS=10004
 integer, parameter :: BUFR_NEPN=10051
 integer, parameter :: BUFR_NBT1=12062
 integer, parameter :: BUFR_NBT2=12063
 integer, parameter :: BUFR_NBT3=12163
 integer, parameter :: BUFR_NEBD=15037
 integer, parameter :: BUFR_NERF=15036
 integer, parameter :: BUFR_NEHU=13210
 integer, parameter :: BUFR_NEHS=13220
 integer, parameter :: BUFR_NEZD=15031
 integer, parameter :: BUFR_NEFE=15032
 integer, parameter :: BUFR_NEZW=15035
 ! Doppler wind 
 integer, parameter :: BUFR_NEDW=11200  
 integer, parameter :: bufr_sst =22042
 
 !
 ! Table C-14: (for Code Table 08046)
 !             Atmospheric chemical or physical constituent type
 !             - Includes local values as well!
 integer, parameter :: BUFR_NECH_O3    =0   
 integer, parameter :: BUFR_NECH_H2O   =1   
 integer, parameter :: BUFR_NECH_CH4   =2   
 integer, parameter :: BUFR_NECH_CO2   =3   
 integer, parameter :: BUFR_NECH_CO    =4   
 integer, parameter :: BUFR_NECH_NO2   =5   
 integer, parameter :: BUFR_NECH_N2O   =6   
 integer, parameter :: BUFR_NECH_HCHO  =7
 integer, parameter :: BUFR_NECH_SO2   =8  
 integer, parameter :: BUFR_NECH_NH3   =9  
 integer, parameter :: BUFR_NECH_NO    =11
 integer, parameter :: BUFR_NECH_PM25  =26 
 integer, parameter :: BUFR_NECH_PM10  =27  

 ! Table B elements associated to constituents/chemistry
 ! -----------------------------------------------------
 !
 ! Element denoting exponents accompanying obs bufr element
 integer, parameter :: BUFR_SCALE_EXPONENT      = 8090
 
 ! Elements for units of constituent observations
 !
 ! When needed, used in tandem with BUFR_NECH_* identifying the consituent
 ! and BUFR_SCALE_EXPONENT when the input file contains an exponent for
 ! a power ten scale factor. 
 !
 ! For some units, multiple BUFR are available, some originally defined to
 ! be specific to certain constituents and some being local elements devised
 ! prior to official elements being assigned. While a single set would be
 ! sufficient, all are included for completeness. A sufficient set might be
 ! 15008,15026,13002,15022,15230,15023,15024,15010,15029,15198,15009,15021,15028.

 integer, parameter :: BUFR_UNIT_VMR            = 15008   ! Volume mixing ratio (vmr)
 integer, parameter :: BUFR_UNIT_VMR2           = 15208   ! Volume mixing ratio  
 integer, parameter :: BUFR_UNIT_MolePerMole    = 15026   ! Pollutant concentration (mole/mole)
 integer, parameter :: BUFR_UNIT_MolePerMole2   = 15197   ! Mixing ratio (mole/mole)
 integer, parameter :: BUFR_UNIT_MMR            = 13002   ! Mass mixing ratio (kg/kg)
 integer, parameter :: BUFR_UNIT_MMR2           = 13001   ! Humidity mass mixing ratio (kg/kg) - same as above
 integer, parameter :: BUFR_UNIT_NumberDensity  = 15022   ! Number density (1/m^3) 
 integer, parameter :: BUFR_UNIT_MolarDensity   = 15230   ! Molar density (mole/m^3)  
 integer, parameter :: BUFR_UNIT_Density        = 15023   ! Density (kg/m^3)
 integer, parameter :: BUFR_UNIT_Density2       = 15027   ! Pollutant concentration in kg/m^3 - same as above
 integer, parameter :: BUFR_UNIT_Density3       = 15223   ! Concentration = Density (kg/m^3) 
 integer, parameter :: BUFR_UNIT_AirDensity     = 15194   ! Air density (kg/m^3) 
 integer, parameter :: BUFR_UNIT_PMDensity      = 15195   ! Density of PM2.5 (kg/m^3)  
 integer, parameter :: BUFR_UNIT_PartPress      = 15010   ! Partial pressure in Pa
 integer, parameter :: BUFR_UNIT_PartPress2     = 15003   ! Ozone partial pressure in Pa (same as above)
 integer, parameter :: BUFR_UNIT_PartPress3     = 15199   ! Partial pressure in Pa (same as above)
 integer, parameter :: BUFR_UNIT_MR_NVaerosol   = 15055   ! Non-volatile aerosol mixing ratio (unitless)  
 integer, parameter :: BUFR_UNIT_ExtinctCoef    = 15029   ! Extinction coefficient (1/m)
 
 ! The following elements are all vertically integrated quantities

 integer, parameter :: BUFR_UNIT_DU             = 15198   ! Integrate number density in Dobson units DU
 integer, parameter :: BUFR_UNIT_DU2            = 15001   ! Total ozone in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_DU3            = 15005   ! Partial column for ozone in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_DU4            = 15045   ! Partial column for SO2 in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_IntegND        = 15009   ! Integrated number density (1/m^2)
 integer, parameter :: BUFR_UNIT_IntegND2       = 15012   ! Electron density per m^2 (1/m^2) 
 integer, parameter :: BUFR_UNIT_IntegDens      = 15021   ! Integrated density (kg/m^2) - same as above
 integer, parameter :: BUFR_UNIT_IntegDens2     = 15020   ! Integrated density for ozone (kg/m^2; same as above) - applicable to all
 integer, parameter :: BUFR_UNIT_IntegDens3     = 15200   ! Integrated density (kg/m^2)

 integer, parameter :: BUFR_UNIT_OptDepth       = 15024   ! Optical depth (unitless)
 integer, parameter :: BUFR_UNIT_OptDepth2      = 15196   ! Optical depth (unitless)
 integer, parameter :: BUFR_UNIT_OptDepth3      = 15062   ! Aerosol Optical depth (unitless)

 ! Additional elements associated to chemistry
 
 integer, parameter :: BUFR_UNIT_PhotoDissoc    = 15028   ! Photodissociate rate (1/sec)
 
 ! ----------------------------------------------
  
contains

!---------------------------------------------------------------
! 
! FUNCTION bufr_IsAtmosConstituent(varNumber) result(var_chm)
!!
!! *Purpose*: Determine if 'varNumber' refers to constituent data
!!            from the CH family with recognized data units. 
!!
!! @author Y.J. Rochon (ARQI/AQRD) Jan 2015, March 2017
!!v           Following recommendation by M. Buehner.
!!
!! Revisions:
!!
!! In
!!
!!v       varNumber     BUFR element number
!!
! ------------------------------------------------------------
  function bufr_IsAtmosConstituent(varNumber) result(var_chm)
  
      implicit none

      integer, intent(in)           :: varNumber
      logical                       :: var_chm
      
      if (any(varNumber.eq. (/ BUFR_UNIT_VMR, BUFR_UNIT_VMR2, BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2, &
                               BUFR_UNIT_MMR, BUFR_UNIT_MMR2, BUFR_UNIT_NumberDensity, BUFR_UNIT_MolarDensity,  &
                               BUFR_UNIT_Density, BUFR_UNIT_Density2, BUFR_UNIT_Density3, &
                               BUFR_UNIT_AirDensity, BUFR_UNIT_PMDensity, &
                               BUFR_UNIT_OptDepth, BUFR_UNIT_OptDepth2, BUFR_UNIT_OptDepth3, BUFR_UNIT_MR_NVaerosol, &
                               BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2, BUFR_UNIT_PartPress3, &
                               BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                               BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2, &
                               BUFR_UNIT_IntegDens, BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3, &
                               BUFR_UNIT_ExtinctCoef, BUFR_UNIT_PhotoDissoc /) )) then          
          var_chm=.true.
      else         
          var_chm=.false.
      end if
      
  end function bufr_IsAtmosConstituent

!------------------------------------------------------------------------------------------
!
! LOGICAL FUNCTION bufr_IsIntegral(varNumber)
!!
!! *Purpose*: Identify if obs is a vertically integrated constituent measurement.
!!            Excludes optical depths and photodissociation rates.
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2015, March 2017
!!
!! Revision: 
!!
!! In
!!
!!v       varNumber     BUFR element number
!!
! ---------------------------------------------------------------------------------------
  logical function bufr_IsIntegral(varNumber)
  
  implicit none
  integer, intent(in) :: varNumber
 
  if (any(varNumber .eq. (/ BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                            BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2, BUFR_UNIT_IntegDens, &
                            BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3 /) )) then      
      bufr_IsIntegral=.true.     
  else
      bufr_IsIntegral=.false.
  end if
  
  end function bufr_IsIntegral

end module bufr_mod
