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
!! *Purpose*:  Set BACKGROUND CHECK FLAGS According to values set in a table.
!!             Original values in table come from ecmwf.
!!
!! @author P. Koclas *CMC/CMSV  September 1998
!!
!! Arguments:
!!     -KVNAM= VARIABLE NAME ( BURP )
!!     -KODTYP=BURP CODE TYPE
!!     -CDFAM= FAMILY  NAME ( 'UA' , 'AI'   ...etc.. )
!!     -ZLEV = LEVEL
!!     -zbgchk=NORMALIZED BACKGROUND DEPARTURE
!!     -lmodif1020=switch to activate special criteria for backound check (*ua 10-20 mb)
!!
!--------------------------------------------------------------------------
      function isetflag(cdfam,kodtyp,kvnam,zlev,zbgchk,lmodif1020)

      use bufr_mod
      IMPLICIT NONE
      integer isetflag
      integer kvnam,kodtyp
      real*8 zlev
      real*8 zbgchk
      character*2 cdfam
      logical lmodif1020

      real*8 zgzcrit(3),zttcrit(3),zuvcrit(3),zescrit(3),zdzcrit(3)
      real*8 zpscrit(3),zpncrit(3),ztscrit(3),zswcrit(3),zzdcrit(3)
      real*8 zchcrit(3)

      isetflag=0

! ASSUME CVCORD = GEMHYB
         zttcrit(1) = 9.00D0
         zttcrit(2) = 16.00D0
         zttcrit(3) = 25.00D0

         zuvcrit(1) = 10.00D0
         zuvcrit(2) = 20.00D0
         zuvcrit(3) = 30.00D0

         zescrit(1) = 10.00D0
         zescrit(2) = 20.00D0
         zescrit(3) = 30.00D0

         zpscrit(1) = 9.00D0
         zpscrit(2) = 16.00D0
         zpscrit(3) = 25.00D0

         zpncrit(1) = 10.00D0
         zpncrit(2) = 20.00D0
         zpncrit(3) = 30.00D0

         zswcrit(1) = 10.00D0
         zswcrit(2) = 20.00D0
         zswcrit(3) = 30.00D0

         ztscrit(1) = 5.00D0
         ztscrit(2) = 25.00D0
         ztscrit(3) = 30.00D0

         zdzcrit(1) = 2.25D0
         zdzcrit(2) = 5.06D0
         zdzcrit(3) = 7.56D0

         zgzcrit(1) = 12.25D0
         zgzcrit(2) = 25.00D0
         zgzcrit(3) = 36.00D0

         zzdcrit(1) = 9.00D0
         zzdcrit(2) = 16.00D0
         zzdcrit(3) = 25.00D0

         zchcrit(1) = 9.00D0
         zchcrit(2) = 16.00D0
         zchcrit(3) = 25.00D0

         if ( kodtyp .eq. 37 ) then
           zuvcrit(2)=25.D0
         else
           zuvcrit(2)=20.D0
         endif
!C
!C     SET FLAG FOR HEIGHTS
!C
      if ( kvnam .eq. BUFR_NEGZ ) then
         if (      zbgchk .gt. zgzcrit(1) .and. zbgchk .lt. zgzcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zgzcrit(2) .and. zbgchk .lt. zgzcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zgzcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR TEMPERATURE
!C
      if ( kvnam .eq. BUFR_NETT ) then
         if (      zbgchk .gt. zttcrit(1) .and. zbgchk .lt. zttcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zttcrit(2) .and. zbgchk .lt. zttcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zttcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SATEMS
!C
      if ( kvnam .eq. BUFR_NEDZ ) then
         if (      zbgchk .gt. zdzcrit(1) .and. zbgchk .lt. zdzcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zdzcrit(2) .and. zbgchk .lt. zdzcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zdzcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR WIND COMPONENTS
!C
      if ( kvnam .eq. BUFR_NEUU .or. kvnam .eq. BUFR_NEVV ) then
         if (      zbgchk .gt. zuvcrit(1) .and. zbgchk .lt. zuvcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zuvcrit(2) .and. zbgchk .lt. zuvcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zuvcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE WIND COMPONENTS
!C
      if ( kvnam .eq. BUFR_NEUS .or. kvnam .eq. BUFR_NEVS ) then
         if (      zbgchk .gt. zswcrit(1) .and. zbgchk .lt. zswcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zuvcrit(2) .and. zbgchk .lt. zswcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zswcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR DEW POINT DEPRESSION
!C
      if ( kvnam .eq. BUFR_NEES ) then
         if (      zbgchk .gt. zescrit(1) .and. zbgchk .lt. zescrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zescrit(2) .and. zbgchk .lt. zescrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zescrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE PRESSURE
!C
      if ( kvnam .eq. BUFR_NEPS ) then
         if (      zbgchk .gt. zpscrit(1) .and. zbgchk .lt. zpscrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zpscrit(2) .and. zbgchk .lt. zpscrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zpscrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR MEAN SEA LEVEL PRESSURE
!C
      if ( kvnam .eq. BUFR_NEPN ) then
         if (      zbgchk .gt. zpncrit(1) .and. zbgchk .lt. zpncrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zpncrit(2) .and. zbgchk .lt. zpncrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zpncrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE TEMPERATURE
!C
      if ( kvnam .eq. BUFR_NETS ) then
         if (      zbgchk .gt. ztscrit(1) .and. zbgchk .lt. ztscrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. ztscrit(2) .and. zbgchk .lt. ztscrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. ztscrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR GB-GPS ZENITH DELAY
!C
      if ( kvnam .eq. BUFR_NEZD ) then
         if (      zbgchk .gt. zzdcrit(1) .and. zbgchk .lt. zzdcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zzdcrit(2) .and. zbgchk .lt. zzdcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zzdcrit(3) )then
              isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR CHEMICAL CONSTITUENTS
!C
      if ( cdfam .eq. 'CH' ) then
         if (      zbgchk .gt. zchcrit(1) .and. zbgchk .lt. zchcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zchcrit(2) .and. zbgchk .lt. zchcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zchcrit(3) )then
              isetflag =3
         endif
      endif

      return
      END FUNCTION ISETFLAG

