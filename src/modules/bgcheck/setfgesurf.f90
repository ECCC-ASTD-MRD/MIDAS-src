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
!! *Purpose*: Interpolate vertically the contents of "column" to
!!            the pressure levels of the observations.
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/AES  September 2000
!!
!--------------------------------------------------------------------------
  subroutine setFGESurf( lcolumn, lcolumng, lobsSpaceData )
  
  use EarthConstants_mod
  use MathPhysConstants_mod
  use bufr_mod
  use columnData_mod
  use obsSpaceData_mod
  implicit none
  ! arguments
  type(struct_columnData) :: lcolumn, lcolumng
  type(struct_obs)        :: lobsSpaceData
  ! locals
  integer          :: ipb, ipt, idim, headerIndex, ik, bodyIndex, ityp
  real(8)          :: zwb, zwt, zlev, zpt, zpb, zhhh
  character(len=2) :: cfam, varLevel
  logical          :: llok

  ! loop over all body rows
  BODY: do bodyIndex = 1, obs_numbody( lobsSpaceData )

    cfam = obs_getFamily( lobsSpaceData, bodyIndex = bodyIndex )
    if( cfam == 'SF'.or. cfam == 'TM' .or. cfam == 'UA' .or. cfam  == 'SC' .or. cfam == 'GP' ) then

      ! Process all data within the domain of the model (excluding GB-GPS ZTD data)
      llok = .false.

      if ( obs_bodyElem_i( lobsSpaceData, OBS_VCO, bodyIndex ) == 1 ) then

        ityp = obs_bodyElem_i( lobsSpaceData, OBS_VNM, bodyIndex )
        if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. ityp == BUFR_NESS .or. &
             ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or. ityp == BUFR_NEFS .or. ityp == BUFR_NEDS .or. &
             ityp == bufr_sst ) then

          llok = ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == 1 )

        else if ( ityp == BUFR_NEZD ) then

          ! make sure total zenith delay (from ground-based GPS) not treated
          llok=.false.

        else

          llok=(obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == 1 .and. &
                obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) >= 0)
          if ( llok ) write(*,*) 'setfgesurf: WARNING!!! unknown obs seen'
          if ( llok ) write(*,*) 'setfgesurf: ityp=',ityp,', cfam=',cfam

        end if

        if ( llok ) then

          headerIndex = obs_bodyElem_i( lobsSpaceData, OBS_HIND, bodyIndex )
          ityp         = obs_bodyElem_i( lobsSpaceData, OBS_VNM , bodyIndex )
          varLevel     = vnl_varLevelFromVarnum( ityp )
          idim = 1
          if ( varLevel == 'SF') idim = 0
          ik   = obs_bodyElem_i( lobsSpaceData, OBS_LYR, bodyIndex )
          zlev = obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex )
          zhhh = zlev

          if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. &
               ityp == BUFR_NESS .or. ityp == BUFR_NEUS .or. ityp == BUFR_NEVS ) then

            ipt  = ik + col_getOffsetFromVarno(lcolumng,ityp)
            ipb  = ipt+1
            call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex, col_getElem( lcolumn, ipb, headerIndex ) )

          else

            ipt  = ik + col_getOffsetFromVarno(lcolumng,ityp)
            ipb  = ipt+1
            zpt  = col_getHeight(lcolumng,ik,headerIndex,varLevel)
            zpb  = col_getHeight(lcolumng,ik+1,headerIndex,varLevel)
            zwb  = idim*(zpt-zhhh)/(zpt-zpb)
            zwt  = 1.d0 - zwb
    
            if ( obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) == 0 ) then

              call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
              zwb * col_getElem( lcolumn, ipb, headerIndex ) + zwt * col_getElem( lcolumn, ipt, headerIndex ))

            else

              call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
                col_getElem( lcolumn, ik + col_getOffsetFromVarno( lcolumng, ityp ), headerIndex ))

            end if

            if(obs_elem_c( lobsSpaceData, 'STID', headerIndex ) == '99999999' ) then

              write(*,*) 'setfgesurf: stn, ityp, xtr, ipt, ipb, zwt, zwb',  &
                   obs_elem_c( lobsSpaceData, 'STID', headerIndex ), ityp, &
                   obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ), ipt, ipb, zwt, zwb
              write(*,*) 'setfgesurf: gobs(ipb), gobs(ipt), fge',   &
                    col_getElem( lcolumn, ipb, headerIndex ), col_getElem( lcolumn, ipt, headerIndex ), &
                    obs_bodyElem_r( lobsSpaceData, OBS_HPHT, bodyIndex )

            endif

          end if
        end if
      end if

    end if

  end do BODY

  end subroutine setFGESurf
