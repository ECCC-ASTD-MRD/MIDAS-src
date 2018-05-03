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
!! MODULE hirchannels (prefix="hir")
!!
!! *PURPOSE*: HYPERSPECTRAL INFRARED CHANNEL PROPERTIES
!!
!! @author A. BEAULNE (CMDA/SMC) February 2006
!!
!!       REVISION: S. Heilliette transformation into a real module
!!                  and addition of IASI and CrIS March 2013
!!
!!       METHODS:
!!              - hir_get_assim_chan : get assimilation flag from channel number
!!              - hir_set_assim_chan : set a vector of  assimilation flags (1=assimilate)
!!
!--------------------------------------------------------------------------
Module hirchannels_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: hir_get_assim_chan,hir_set_assim_chan


  INTEGER, PARAMETER :: AIRSSNCH = 281
  INTEGER,target   :: AIRS_ASSIM(AIRSSNCH)



  INTEGER, PARAMETER :: IASISNCH = 616

  INTEGER,target   :: IASI_ASSIM(IASISNCH)

  INTEGER, PARAMETER :: CRISSNCH = 1305
  INTEGER ,target  :: CRIS_ASSIM(CRISSNCH)


contains



  subroutine hir_set_assim_chan(CINST,iselec)
    implicit none
    character (len=*),intent(in) :: CINST
    integer,intent(in) :: iselec(:)
    !*************************************************
    integer :: nch,i
    integer ,pointer :: pt(:)
    !***************************

    nch=size(iselec)

    if (nch<1) then
       Write(*,*) "Empty channel selection !"
       call utl_abort("hir_set_assim_chan from  hirchannels module")
    endif

    if (nch/=hir_get_nchan_selected(CINST)) then
       Write(*,*) "Wrong number of channels",nch
       Write(*,*) "Should be ",hir_get_nchan_selected(CINST)
       call utl_abort("hir_set_assim_chan from  hirchannels module")
    endif

    select case(trim(CINST))
    case("AIRS","airs")
       pt=>airs_assim
    case("IASI","iasi")
       pt=>iasi_assim
    case("CRIS","cris","CrIS")
       pt=>cris_assim
    case default
       Write(*,*) "Unknown instrument! ",CINST
       call utl_abort("set_assim_chan from  hirchannels module")
    end select

    do i=1,nch
       if (iselec(i)==0 .or. iselec(i)==1) then
          pt(i)=iselec(i)
       else
          Write(*,*) "Invalid selection flag ",iselec(i),i,cinst
          call utl_abort("hir_set_assim_chan from  hirchannels module")
       endif
    enddo

    nullify(pt)

  end subroutine hir_set_assim_chan


  integer function hir_get_assim_chan(CINST,ichan)
    implicit none
    character (len=*),intent(in) :: CINST
    integer,intent(in) :: ichan
    !*******************************
    integer ,pointer :: pt(:)

    select case(trim(CINST))
    case("AIRS","airs")
       pt=>airs_assim
    case("IASI","iasi")
       pt=>iasi_assim
    case("CRIS","cris","CrIS")
       pt=>cris_assim
    case default
       Write(*,*) "Unknown instrument! ",CINST
       call utl_abort("hir_get_assim_chan from  hirchannels module")
    end select

    if (ichan<1 .or. ichan>hir_get_nchan_selected(CINST)) then
       Write(*,*) "Invalid channel index",ichan
       call utl_abort("hir_get_assim_chan from  hirchannels module")
    endif

    hir_get_assim_chan = pt(ichan)

    nullify(pt)

  end function hir_get_assim_chan

End module hirchannels_mod
