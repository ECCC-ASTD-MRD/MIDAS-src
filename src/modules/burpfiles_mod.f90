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
!! MODULE burpFiles (prefix="brpf")
!!
!! *Purpose*: To store the filenames of the burp observation files and call
!!            subroutines in readBurp to read and update burp files.
!!
!--------------------------------------------------------------------------
module burpFiles_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use utilities_mod
  use mpivar_mod
  use obsSpaceData_mod
  use burpread_mod
  use ramDisk_mod
  use bufr_mod
  use utilities_mod
  use obsSubSpaceData_mod
  use burp_module
  use mpi_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: brpf_getDateStamp, brpf_readfile, brpf_updatefile
  public :: brpf_chem_read_all, brpf_chem_update_all

contains

  subroutine brpf_getDateStamp(datestamp, burpFileName)
    implicit none

    ! arguments
    integer :: dateStamp
    character(len=*), intent(in) :: burpFileName

    ! locals
    integer :: ier, inblks, nulburp, fnom, fclos, numblks
    logical :: isExist_L 
    integer :: ktime, kdate, kdate_recv, ktime_recv, ihandl, ilong
    integer :: itime, iflgs, idburp, ilat, ilon, idx, idy
    integer :: ialt, idelay, idate, irs, irunn, inblk, isup, ixaux
    integer :: insup, inxaux
    integer, allocatable :: ibuf(:)
    integer :: inrecs, mrfcls, mrfopn, mrfopc, mrbhdr, mrfloc, mrfget, mrfmxl
    integer :: istampobs, inewhh, newdate, nresume, ivals
    real(8) :: delhh
    character(len=9) :: clstnid

    !
    !- Get the date from the burp files
    !

    ier = mrfopc('MSGLVL','FATAL')

    ivals = 8
    kdate = -9999
    ktime = -9999
    nresume = 0
    nulburp = 0
    inquire(file=trim(burpFileName),exist=isExist_L)
    if ( isExist_L ) then
      ier = fnom(nulburp,trim(burpFileName),'RND+OLD',0)
      write(*,*)' Open File : ',trim(burpFileName)
      if ( ier == 0 ) then
        inblks = -1
        inblks = numblks(nulburp)
        if ( inblks > 0 ) then
          inrecs= mrfopn(nulburp,'READ')
          ilong = mrfmxl(nulburp)
          allocate(ibuf(ilong + 20))
          ibuf(1) = ilong + 20
          ihandl  = mrfloc(nulburp,0,'>>*******',-1,-1,-1,-1,-1,-1,0)
          if ( ihandl < 0 ) then
            ihandl=mrfloc(nulburp,0,'*********',-1,-1,-1,-1,-1,-1,0)
          else
            nresume=nresume+1
          end if
          if ( ihandl < 0 ) then
            write(*,*) 'AUCUN ENREGISTREMENT VALIDE DANS LE FICHIER BURP'
          else
            if ((kdate < 0.and.ktime < 0).or.nresume == 1) then 
              insup=0
              inxaux=0
              ier=mrfget(ihandl,ibuf)
              ier=mrbhdr(ibuf,itime,iflgs,clstnid,idburp,ilat,   &
                   ilon,idx,idy, ialt,idelay,idate,irs,irunn,inblk, &
                   isup,insup,ixaux,inxaux)
              ktime=itime
              kdate=idate
              if (nresume == 1) nresume=2
            end if
          end if
          deallocate(ibuf)
          ier=mrfcls(nulburp)
        end if
        end if
      ier= fclos(nulburp)
    end if

    !
    !- Set reference datestamp
    !
    ! Make sure all mpi tasks have a valid date (important for split burp files)
    call rpn_comm_allreduce(kdate,kdate_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)
    call rpn_comm_allreduce(ktime,ktime_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)
    kdate = kdate_recv
    ktime = ktime_recv
    if (nresume >= 1 ) then  
      ier = newdate(datestamp,kdate,ktime*10000,3)
    else
      ! Assumes 6-hour windows with reference times being synoptic times.
      ! Does not require kdate and ktime to be from a resume record.
      ier = newdate(istampobs,kdate,ktime*10000,3)
      delhh = 3.0d0
      call incdatr (datestamp, istampobs, delhh)
      ier = newdate(datestamp,kdate,inewhh,-3)
      ktime = ktime/100
      if (ktime >= 21 .or. ktime < 3) then
        ktime = 0
      else if(ktime >= 3 .and. ktime < 9) then
        ktime = 6
      else if(ktime >= 9 .and. ktime < 15) then
        ktime = 12
      else
        ktime = 18
      end if
      ier = newdate(datestamp,kdate,ktime*1000000,3)
      ktime = ktime*100
    end if

    write(*,*)' BURP FILES VALID DATE (YYYYMMDD) : ', kdate
    write(*,*)' BURP FILES VALID TIME     (HHMM) : ', ktime
    write(*,*)' BURP FILES DATESTAMP             : ', datestamp

  end subroutine brpf_getDateStamp


  subroutine brpf_readFile(obsdat,fileName,familyType,fileIndex)
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*) :: fileName
    character(len=*) :: familyType
    integer :: fileINdex

    ! locals
    integer :: ibeg, iend, nstn1, nstn2
    logical :: burp_chem
    integer :: jo
    real(obs_real)  :: misg

    write(*,*) ' '
    write(*,*) 'brpf_readFile: Starting'
    write(*,*) ' '
    misg = real(MPC_missingValue_R8,OBS_REAL)

    ibeg = obs_numbody(obsdat) +1
    Nstn1 = obs_numheader(obsdat)

    call brpr_readBurp(obsdat,familyType,fileName,fileIndex)
    Nstn2 = obs_numheader(obsdat)
    iend = obs_numbody(obsdat)

    burp_chem = trim(familyType) == 'CH'

    if ( trim(familyType) /= 'TO' .and. .not.burp_chem) THEN
      call FDTOUV_OBSDAT(  obsdat,Nstn1+1,Nstn2,MPC_missingValue_R4)
      call ADJUST_HUM_GZ(  obsdat,Nstn1+1,Nstn2)
      call ADJUST_SFVCOORD(obsdat,Nstn1+1,Nstn2)
    end if
    do jo = nstn1+1, nstn2
      call obs_headSet_i(obsdat,OBS_OTP,JO,fileIndex)
      call obs_setFamily(obsdat,trim(familyType),JO)
       ! For CH family, apply scaling from the element BUFR_SCALE_EXPONENT when present.
       if (burp_chem) call set_scale_chm(obsdat,JO,forward=.true.)
    end do
    !    initializations
    do jo = ibeg, iend
      if ( obs_columnActive_RB(obsdat,OBS_OMA) )  call obs_bodySet_r(obsdat,OBS_OMA ,JO,MISG)
      if ( obs_columnActive_RB(obsdat,OBS_OMP) )  call obs_bodySet_r(obsdat,OBS_OMP ,JO,MISG)
      if ( obs_columnActive_RB(obsdat,OBS_OER) )  call obs_bodySet_r(obsdat,OBS_OER ,JO,MISG)
      if ( obs_columnActive_RB(obsdat,OBS_HPHT) ) call obs_bodySet_r(obsdat,OBS_HPHT,JO,MISG)
      if ( obs_columnActive_RB(obsdat,OBS_WORK) ) call obs_bodySet_r(obsdat,OBS_WORK,JO,MISG)
    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD data (element 15031)
    if ( trim(familyType) == 'GP') then
      print * ,' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call set_err_gbgps(obsdat,Nstn1+1,Nstn2)
    end if

    write(*,*) 'brpf_readFile: obs_numheader(obsdat)', obs_numheader(obsdat)
    write(*,*) 'brpf_readFile: obs_numbody(obsdat)  ', obs_numbody  (obsdat)

  end subroutine brpf_readFile


  subroutine brpf_updateFile(obsSpaceData,fileName,familyType,fileIndex)
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*) :: fileName
    character(len=*) :: familyType
    integer :: fileIndex

    ! locals
    integer :: headerIndex

    call tmg_start(93,'POST_UPDATEBRP')

    write(*,*) 'brpf_updateFile: Starting'
      
    ! CH family: Scaling of the obs related values to be stored in the BURP files

    call obs_set_current_header_list(obsSpaceData,'CH')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      call set_scale_chm(obsSpaceData,headerIndex,forward=.false.)
    end do HEADER

    call brpr_updateBurp(obsSpaceData,familyType,fileName,fileIndex)

    write(*,*) 'brpf_updateFile: Done'

    call tmg_stop(93)

  end subroutine brpf_updateFile


  subroutine  adjust_hum_gz(obsdat,start,end)
!**s/r ADJUST_HUM_GZ  - Adjust  t-td and GZ in obsdat
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!
!*    Purpose:  - Adjust  t-td values to zesmax=30. in obsdat
!                 set Z to GZ                       in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
!
      implicit none
      integer  :: start,end

      integer  :: j,jo,rln,nlv
      integer  :: varno
      type (struct_obs), intent(inout):: obsdat

      real(obs_real)    :: zesmax,gz,obsv
      real              :: rmin

      zesmax = 30.0

      write(*,*)'   ADJUST_HUM_GZ '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,JO)
        nlv = obs_headElem_i(obsdat,OBS_NLV,JO)

        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          select case(varno)
            case(bufr_nees,bufr_ness)
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             if ( obsv > zesmax) then
                obsv = zesmax
             end if
             call obs_bodySet_r(obsdat,OBS_VAR,j, obsv )
            case(bufr_negz)
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             gz = obsv*grav
             call obs_bodySet_r(obsdat,OBS_VAR,j,gz )
          end select

        end do

      end do

      write(*,*)' DONE   ADJUST_HUM_GZ '

  end subroutine adjust_hum_gz


  subroutine  set_err_gbgps(obsdat,start,end)
!**s/r SET_ERR_GBGPS  - SET INITIAL ERROR FRO GROUND BASED GPS
!
!
!Author  : P. Koclas *CMC/CMDA  July 2013
!Revision:
!
!*    Purpose:  - PUT 15032 observation element as error of 15031 element  in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      implicit none
      real(obs_real)    :: obsv
      integer  :: start,end

      integer  :: j,jo,rln,nlv
      integer  :: varno
      type (struct_obs), intent(inout):: obsdat

      write(*,*)'   SET_ERR_GBGPS '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,jo)
        nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

        obsv = real(MPC_missingValue_R8,obs_real)
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          if ( varno == 15032 ) then
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             call obs_bodySet_i(obsdat,OBS_VNM,j,999 )
             exit
          end if

        end do
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          if ( varno == 15031 .and. obsv /= real(MPC_missingValue_R8,obs_real)) then
             call obs_bodySet_r(obsdat,OBS_OER,j,obsv)
             exit
          end if

        end do

      end do

      write(*,*)' DONE   SET_ERR_GBGPS '

  end subroutine set_err_gbgps

!--------------------------------------------------------------------------
!! *Purpose*:  Apply or unapply scaling to CH observations  by multiplying
!!             (or dividing) with 10^{exponent} where the exponent is from
!!             element BUFR_SCALE_EXPONENT if provided.           
!!
!! @author P. Du *CMC/CMDA Nov. 2014
!!
!! Revisions:
!!v     Y. Rochon, ARQI/AQRD, Feb 2015
!!v       - Generalized except for the assumption that if
!!v         exponents are present, then they must be present
!!v         for all obs and in the same sequential order.
!!v     Y. Rochon, ARQI/AQRD, March/April 2016
!!v       - Added consideration of HBHT (HPHT)
!!v       - Avoid abort when BUFR_SCALE_EXPONENT present but not the mantissa as it could 
!!v         have been filtered in brpr_readBurp.
!!v     M. Sitwell, ARQI/AQRD, March 2017
!!v       - Modified to do scaling for a single profile (i.e. one headerIndex)
!!v         instead of for all observations (loop is done outside of set_scale_chm) 
!!
!! Input:
!!v      obsdat         struct_obs instance
!!v      headerIndex    header index in obsdat to apply/unapply scaling to
!!v      forward        applies scaling if .true., unapplies scaling if .false.
!--------------------------------------------------------------------------
  subroutine set_scale_chm(obsdat,headerIndex,forward)

      implicit none
      
      type (struct_obs), intent(inout):: obsdat
      integer, intent(in) :: headerIndex
      logical, intent(in) :: forward

      integer  :: bodyIndex,rln,nlv

      real(obs_real) :: obsv
      real(obs_real) :: vomp, voma, voer, vhpht, scale
      integer        :: nexp,iobs,iexp

      real(obs_real), allocatable :: expnt(:)

      rln = obs_headElem_i(obsdat,OBS_RLN,headerIndex)
      nlv = obs_headElem_i(obsdat,OBS_NLV,headerIndex)

      allocate(expnt(nlv))

      ! Count number of power of 10 exponents
      nexp = 0
      do bodyIndex = rln, nlv + rln -1
         if (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) == bufr_scale_exponent) then
            nexp = nexp + 1
            expnt(nexp) = obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
         end if
      end do

      if (nexp == 0) then
         deallocate(expnt)
         return
      end if

      if (nexp /= nlv/2) then
         ! Skip over obs assuming mantissa was filtered out in brpr_readBurp 
         ! (not inserted in obsSpaceData) due to quality flags.
         ! Set exponent quality flag to that of a 'Suspicious element' 
         
         do bodyIndex = RLN, NLV + RLN -1
            call obs_bodySet_i(obsdat,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsdat,OBS_FLG,bodyIndex),04) )
         end do
              
         ! write(*,*) 'NLV =',nlv,' Nexp=',nexp    
         ! call utl_abort('set_scale_chm: Inconsistent number of exponents')
         deallocate(expnt)
         return
      end if

      if (forward) then
         
         ! Apply power of 10 exponents if present
         iobs = 0
         do bodyIndex = RLN, NLV + RLN -1
            if (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) /= bufr_scale_exponent) then
               iobs = iobs + 1
               obsv = obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
               call obs_bodySet_r(obsdat,OBS_VAR,bodyIndex,obsv*10**(expnt(iobs)) )
            end if
         end do
      
      else
             
         ! Unapply power of 10 exponents if present
         iobs=0
         iexp=0
         do bodyIndex = RLN, NLV + RLN -1
            if (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) == bufr_scale_exponent) then
               ! Store scaling exponents
               iexp = iexp + 1
               call obs_bodySet_r(obsdat,OBS_OMP,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_OMA,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_OER,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_HPHT,bodyIndex,expnt(iexp))
            else
               iobs=iobs+1
               obsv=obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
               vomp=obs_bodyElem_r(obsdat,OBS_OMP,bodyIndex)
               voma=obs_bodyElem_r(obsdat,OBS_OMA,bodyIndex)
               voer=obs_bodyElem_r(obsdat,OBS_OER,bodyIndex)
               vhpht=obs_bodyElem_r(obsdat,OBS_HPHT,bodyIndex)
               scale=10**(-expnt(iobs))
               call obs_bodySet_r(obsdat,OBS_VAR,bodyIndex,obsv*scale )
               call obs_bodySet_r(obsdat,OBS_OMP,bodyIndex,vomp*scale )
               call obs_bodySet_r(obsdat,OBS_OMA,bodyIndex,voma*scale )
               call obs_bodySet_r(obsdat,OBS_OER,bodyIndex,voer*scale )
               call obs_bodySet_r(obsdat,OBS_HPHT,bodyIndex,vhpht*scale )
            end if
         end do
                 
      end if

  end subroutine set_scale_chm


  subroutine  adjust_sfvcoord(obsdat,start,end)
!
!**s/r ADJUST_SFVCOORD  - Computation of HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!          S. Macpherson *ARMA  Oct 2013
!              -- add GB-GPS (GP family) element BUFR_NEZD (ele 15031)
!              -- NOTE that for GP data, ELEV = GPS Antenna Height so
!                 no adjustment is needed (SFC_VCO=0).
!
!*    Purpose:  -Compute  HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!                and INSERT INTO CMA.
!
!
!Arguments
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      implicit none
      integer  :: start,end
      integer  :: j,jo,rln,nlv
      integer  :: varno,codtyp,ity
      real     :: sfc_vco,elev
      real(obs_real) :: ppp
      type (struct_obs), intent(inout):: obsdat

      write(*,*)'   ADJUST_SFVCOORD '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,jo)
        nlv = obs_headElem_i(obsdat,OBS_NLV,jo)
        ity = obs_headElem_i(obsdat,OBS_ITY,jo)
        codtyp = ity
        elev = obs_headElem_r(obsdat,OBS_ALT,jo)
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          select case(varno)
            case(bufr_neds,bufr_nefs,bufr_neus,bufr_nevs,bufr_nets,bufr_ness,bufr_nepn,bufr_neps,bufr_nehs,bufr_nezd)

             sfc_vco= surfvcord(varno,codtyp)
             if ( varno /= bufr_nepn) then
                ppp = elev + sfc_vco
                call obs_bodySet_r(obsdat,OBS_PPP,j,ppp)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             else
                 ppp = 0.
                call obs_bodySet_r(obsdat,OBS_PPP,j,ppp)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             end if
          end select
        end do

      end do

      write(*,*)' DONE   ADJUST_SFVCOORD '

  end subroutine adjust_sfvcoord


  real function surfvcord(ilem,idtyp)
      implicit none
      integer :: ilem,idtyp,type
      real :: vcordsf2
!***********************************************************************
!
!      PURPOSE: SEt vertical coordinate for surface data.
!
!       AUTHOR:   P. KOCLAS (CMC/CMDA) December 2011
!
!       Revision : 
!
!    ARGUMENTS:
!               INPUT:
!                      -ILEMP   : BURP ELEMENT NUMBER
!                      -IDTYP   : BURP CODETYPE
!
!               OUTPUT:
!                      -SURFVCORD
!
!
!***********************************************************************
!
!
!     GENERATE TABLES TO ADJUST VERTICAL COORDINATE OF SURFACE DATA
!
!     DEFAULT VALUE 
      vcordsf2=0.

      select case(idtyp)
        case(135,136,137,138,32,34,35,37,38,159,160,161,162)
          ! UPPER AIR LAND
          type=3

        case(139,140,141,142,33,36)
          ! UPPER AIR SHIP
          type=4

        case(12,14,146)
          ! SYNOPS
          type=1

        case(13,18,145,147)
          ! SHIPS
          type=2

        case(254)
          ! SCATTEROMETER WINDS
          type=5

        case default
          type=-99

      end select

      select case(type)
         case (1)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               ! us,vs,ffs,dds
               vcordsf2=10.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees,bufr_ness)
               vcordsf2=1.5

           end select

         case (2)
           select case(ilem)
             ! us,vs,ffs,dds
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=20.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=11.5

             case (bufr_nees,bufr_ness)
               vcordsf2=11.5
           end select

         case (3)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=10.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees)
               vcordsf2=0.0

             case (bufr_ness)
               vcordsf2=1.5

           end select

         case (4)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=20.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees)
               vcordsf2=0.0

             case (bufr_ness)
               vcordsf2=1.5

           end select

         case (5)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=10.0

           end select

      end select

      surfvcord = vcordsf2

  end function  surfvcord


  subroutine fdtouv_obsdat(obsdat,start,end,ppmis)
!
!---------------------------------------------------------------
!
! Author  : P. Koclas, CMC/CMDA December  2012
!           CONVERT DD , FF  WINDS TO
!            UU (est-west),  VV (north-south) COMPONENTS
!
!    ARGUMENTS:
!                 INPUT:
!
!                       -obsdat     : CMA_table INSTANCE 
!                       -START     : FIRST OBERVATION
!                       -END       : LAST  OBERVATION
!                       -PPMIS     : MISSING VALUE  
!
!        **************************************************
!         IT IS ASSUMED THAT CMA CONTAINS ENTRIES   FOR 
!          UU AND VV  with observed values = missing value
!        **************************************************
!
!---------------------------------------------------------------
!
    implicit none
    type (struct_obs), intent(inout) :: obsdat
    real(4)        :: ppmis
    integer(4)     :: start,end

    integer(4)     :: varno,varno2,varno4
    real(4)        :: obsuv
    integer(4)     :: jo,rln,nlv,j,j2,j4,jpos,ilem
    integer(4)     :: ddflag,ffflag,newflag,uuflag,vvflag
    integer(4)     :: ilemf,ilemu,ilemv,indu_misg,indv_misg,indum,indvm
    logical        :: llmisdd,llmisff,llmis,lluv_misg,llu_misg,llv_misg
    logical        :: lluv_present,llu_present,llv_present
    real(obs_real) :: uu,vv,dd,ff
    real(obs_real) :: level_dd,level4,level,level_uu

    ffflag = 0  ! bhe 

    HEADER1: do jo = start, end
        
      rln = obs_headElem_i(obsdat,OBS_RLN,jo)
      nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

      BODY1: do j = rln, nlv + rln -1
        dd = ppmis
        ff = ppmis

        varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
        llmisdd = .true.

        if ( varno /= bufr_nedd .and. varno /= bufr_neds ) cycle BODY1

        if( varno == bufr_neds) then
          ilemf = bufr_nefs
          ilemu = bufr_neus
          ilemv = bufr_nevs
        else
          ilemf = bufr_neff
          ilemu = bufr_neuu
          ilemv = bufr_nevv
        end if

        dd       = obs_bodyElem_r(obsdat,OBS_VAR,j)
        ddflag   = obs_bodyElem_i(obsdat,OBS_FLG,j)
        level_dd = obs_bodyElem_r(obsdat,OBS_PPP,j)
          
        llu_misg = .false.
        llv_misg = .false.
        llu_present = .false.
        llv_present = .false.
        indum = -1
        indvm = -1

        ! FIND IF  U AND V ARE ALREADY IN CMA
        UVINOBSDAT: do j4 = j, nlv + rln -1

          level4 = obs_bodyElem_r(obsdat, OBS_PPP,j4)
          if (level4 == level_dd) then
            varno4 = obs_bodyElem_i(obsdat, OBS_VNM,j4)
            select case (varno4)
              case (11003,11004,11002,11001,11215,11216,11011,11012)

                obsuv = obs_bodyElem_r(obsdat, OBS_VAR,j4)
                if (  (varno4 == ilemu)     .and.  (obsuv /= ppmis) ) then
                  llu_present = .true.
                  indum = j4
                else if ( (varno4 == ilemv) .and. (obsuv /= ppmis) ) then
                  llv_present = .true.
                  indvm = j4
                end if
                
                if (  (varno4 == ilemu)     .and. (obsuv == ppmis) ) then
                  llu_misg = .true.
                  indu_misg = j4
                else if ( (varno4 == ilemv)  .and. (obsuv == ppmis) ) then
                  llv_misg = .true.
                  indv_misg = j4
                end if

            end select
          end if

        end do UVINOBSDAT

        lluv_misg = (llu_misg .and. llv_misg)
        lluv_present = (llu_present .and. llv_present)

        if ( lluv_misg) then

          CALCUV: do j2 = j, nlv + rln -1

            llmisff = .true.
            llmisdd = .true.
            llmis   = .true.
            level = obs_bodyElem_r(obsdat,OBS_PPP,j2)
            if ( level /= level_dd) cycle
            varno2 = obs_bodyElem_i(obsdat,OBS_VNM,j2)
            if ( varno2 == ilemf ) then

              ff = obs_bodyElem_r(obsdat,OBS_VAR,j2)
              ffflag = obs_bodyElem_i(obsdat,OBS_FLG,j2)
              IF ( (dd == 0.  .and. ff > 0.) .or. ( dd > 360. .or. dd  < 0.) ) then
                llmisdd = .true.
                llmisff = .true.
              else if ( dd == ppmis .or. ff == ppmis) then
                llmisdd = .true.
                llmisff = .true.
              else
                llmisdd = .false.
                llmisff = .false.
              end if

              ! IF SPEED = 0 CALM WIND IS ASSUMED.
              if (ff == 0.0) then
                dd = 0.
              end if
                   
              dd = dd + 180.
              if ( dd > 360.) dd = dd - 360.
              dd = dd * mpc_radians_per_degree_r8
                
              ! U,V COMPONENTS ARE
              uu = ff * sin(dd)
              vv = ff * cos(dd)
              if  ( ( llmisdd .eqv. .true.) .or. ( llmisff .eqv. .true. ) ) then
                llmis = .true.
                if ( indu_misg > 0 .or. indv_misg > 0 ) then
                  call obs_bodySet_i(obsdat,OBS_VNM,INDU_MISG,-1)
                  call obs_bodySet_i(obsdat,OBS_VNM,INDV_MISG,-1)
                end if
              else
                llmis = .false.
              end if

            end if
            newflag = ior(ddflag,ffflag)

            if ( indum > 0 .or. indvm > 0 ) then
              call obs_bodySet_i(obsdat,OBS_VNM,indu_misg,-1)
              call obs_bodySet_i(obsdat,OBS_VNM,indv_misg,-1)
            end if
            if (llmis) then
              if ( indum > 0 .or. indvm > 0 ) then
                call obs_bodySet_i(obsdat,OBS_FLG,induM,newflag)
                call obs_bodySet_i(obsdat,OBS_FLG,indvM,newflag)
              end if
            else if (.not.llmis) then
              call obs_bodySet_r(obsdat,OBS_VAR,indu_misg,uu)
              call obs_bodySet_i(obsdat,OBS_FLG,indu_misg,newflag)

              call obs_bodySet_r(obsdat,OBS_VAR,indv_misg,vv)
              call obs_bodySet_i(obsdat,OBS_FLG,indv_misg,newflag)
            end if

          end do CALCUV

        else                       

          if ( lluv_present ) then
            call obs_bodySet_i(obsdat,OBS_VNM,indu_misg,-1)
            call obs_bodySet_i(obsdat,OBS_VNM,indv_misg,-1)
          ELSE
            if (indum > 0) then
              call obs_bodySet_i(obsdat,OBS_VNM,indum,-1)
            end if
            if (indvm > 0) then
              call obs_bodySet_i(obsdat,OBS_VNM,indvm,-1)
            end if
          end if

        end if

      end do BODY1

    end do HEADER1


    do jo = start, end

      rln = obs_headElem_i(obsdat,OBS_RLN,jo)
      nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

      do j = rln, nlv + rln -1

        llmisdd = .true.
        varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
        level = obs_bodyElem_r(obsdat,OBS_PPP,j)

        select case (varno)

          case (bufr_neuu)
            ilem = bufr_nevv
          case (bufr_neus)
            ilem = bufr_nevs        
          case default
            cycle

        end select

        Jpos = -1

        ! TRANSFER THE FLAG BITS  FROM ONE WIND COMPONENT TO THE OTHER
        do j4 = rln, nlv + rln -1
          uu       = obs_bodyElem_r(obsdat,OBS_VAR,j4)
          level_uu = obs_bodyElem_r(obsdat,OBS_PPP,j4)
          Jpos = -1
        
          if ( level_uu == level .and. uu == ppmis ) then
            call obs_bodySet_i(obsdat,OBS_VNM,j4,-1)
          end if

          if ( level_uu == level .and. uu /= ppmis ) then
            uuflag = obs_bodyElem_i(obsdat,OBS_FLG,j4)
            varno2 = obs_bodyElem_i(obsdat,OBS_VNM,j4)

            if ( ilem == varno2 ) then
              vvflag  = obs_bodyElem_i(obsdat,OBS_FLG,j)
              newflag = ior(uuflag,vvflag)
              call obs_bodySet_i(obsdat,OBS_FLG,j, newflag)
              call obs_bodySet_i(obsdat,OBS_FLG,j4,newflag)
              jpos = j4
              exit
            end if

          end if

        end do !j4

        ! ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UU OR VV) IS MISSING
        if (jpos < 0) then
          write(*,*) ' eliminate winds for station : ', obs_elem_c(obsdat,'STID',JO),  &
            obs_bodyElem_i(obsdat,OBS_VNM,J), obs_bodyElem_r(obsdat,OBS_PPP,J)
          call obs_bodySet_i(obsdat,OBS_VNM,j,-1)
        end if

      end do !j

    end do !jo

  end subroutine fdtouv_obsdat


!--------------------------------------------------------------------------
!! *Purpose*: Retrieves information for observations from all BURP files. Currently only used for
!!            chemical constituent but can potentially be used for any family. If the BURP files
!!            are not split, then the single BURP file is read. If the BURP files are split, then
!!            each split BURP file is read and the data combined into one struct_oss_obsdata
!!            object. See the function brpf_chem_read for more details.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Revisions:
!!v       Y. Rochon, ARQI/AQRD, Nov 2016
!!v         - Added optional input argument codtyplist and option of varno <=0 (in brpf_chem_read)
!!
!! Input:
!!v           stnid         station ID of observation
!!v           varno         BUFR code 
!!v                         If <=0, to search through all codes to obtain first
!!v                         between 10000 and 16000.
!!v           nlev          number of levels in the observation
!!v           ndim          number of dimensions for the retrieved data in
!!v                         each report (e.g. ndim=1 for std, ndim=2 for
!!v                         averagine kernels) 
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           match_nlev    determines if the report matching criteria includes checking
!!v                         if the report number of levels is the same as the input
!!v                         argument nlev
!!v           codtyplist    optional CODTYP list for search (optional)
!!v           obsfam        observation family name
!!
!! Output:
!!v           burp_out      struct_oss_obsdata object
!!
!! Comment:
!!
!!v   This routine is general enough to be used by observation families
!!v   other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function brpf_chem_read_all(obsfam,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev,  &
                              codtyplist_opt, obsFilesSplit_opt) result(burp_out)

    implicit none

    character(len=9), intent(in) :: stnid
    character(len=4), intent(in) :: block_type
    integer, intent(in)          :: ndim,varno,nlev,bkstp
    logical, intent(in)          :: match_nlev
    integer, intent(in), optional :: codtyplist_opt(:)
    character(len=*), intent(in) :: obsfam
    type(struct_oss_obsdata) :: burp_out
    logical, optional        :: obsFilesSplit_opt

    character(len=500) :: filename
    logical :: found, obsFilesSplit

    if ( present(obsFilesSplit_opt) ) then
      obsFilesSplit = obsFilesSplit_opt
    else
      obsFilesSplit = .true.
    end if

    filename = brpf_get_filename(obsfam,found)

    if (found) then
       burp_out = brpf_chem_read(filename,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev, &
                                 codtyplist_opt=codtyplist_opt)
    else
       if (obsFilesSplit) then
          ! Must allocate burp_out so that it is available from ALL processors when
          ! requiring of rpn_comm_allgather via oss_obsdata_MPIallgather.
          write(*,*) "brpf_chem_read_all: Could not find/open BURP file: ",trim(filename)
          if (ndim == 1) then
             call oss_obsdata_alloc(burp_out,1,dim1=nlev)
          else
             call oss_obsdata_alloc(burp_out,1,dim1=nlev,dim2_opt=nlev)
          end if
          burp_out%nrep=0
          write(*,*) "brpf_chem_read_all: Number of reports set to ",burp_out%nrep
       else
          call utl_abort('brpf_chem_read_all: Could not find BURP file: ' // trim(filename))
       end if
    end if

    if (obsFilesSplit) call oss_obsdata_MPIallgather(burp_out)

  end function brpf_chem_read_all

!--------------------------------------------------------------------------
!! *Purpose*: Retrieve information from observation BURP file. Can retrieve
!!            either 1D or 2D data from a report. Currently only used for
!!            chemical constituent but can potentially be used for any family.
!!
!! @author M. Sitwell, ARQI/AQRD, March 2015
!!
!! Revision:
!!v    M. Sitwell, ARQI/AQRD, May 2015
!!v      - Modified to read both 1D and 2D data from a report
!!v    Y. Rochon, ARQI/AQRD, May 2016
!!v      - Updated to increment 'icount' only if varno is also found in addition
!!v        stnid, nlev, block_type and bkstp
!!v    M. Sitwell, ARQI/AQRD, June 2016
!!v      - Added match_nlev input argument
!!v    Y. Rochon, ARQI/AQRD, Oct 2016
!!v      - Added optional codtyplist argument and option of input varno<=0.
!!
!! Input:
!!v           filename      BURP file name
!!v           stnid         station ID of observation
!!v           varno         BUFR code
!!v                         If <=0, search through all codes to obtain first
!!v                         between 10000 and 16000.
!!v           nlev          number of levels in the observation
!!v           ndim          number of dimensions for the retrieved data in
!!v                         each report (e.g. ndim=1 for std, ndim=2 for
!!v                         averagine kernels) 
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           match_nlev    determines if the report matching criteria includes checking
!!v                         if the report number of levels is the same as the input
!!v                         argument nlev
!!v           codtyplist    optional CODTYP list for search
!!
!! Output: 
!!v           burp_out      struct_oss_obsdata object
!!
!! Comments:
!!
!!v   - BUFR power 10 exponent element (i.e. data with BUFR number BUFR_SCALE_EXPONENT) 
!!v     will only be applied to 1D data if present.
!!v   - As burp_out is for a specific input stnid, burp_out%code only contains the (lat/long and
!!v     time coord.) with 22 characters.
!!v   - Exponent BUFR data (i.e. data with BUFR number BUFR_SCALE_EXPONENT) will only 
!!v     be applied to 1D data.
!!v   - This routine is general enough to be used by observation families
!!v     other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function brpf_chem_read(filename,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev,  &
                          codtyplist_opt) result(burp_out)
    
    implicit none

    character(len=*), intent(in) :: filename
    character(len=9), intent(in) :: stnid
    character(len=4), intent(in) :: block_type
    integer, intent(in)          :: ndim,varno,nlev,bkstp
    logical, intent(in)          :: match_nlev
    integer, intent(in), optional :: codtyplist_opt(:)
    type(struct_oss_obsdata) :: burp_out

    character(len=9)  :: rep_stnid
    type(burp_file)   :: brp
    type(burp_rpt)    :: rep
    type(burp_block)  :: blk
    integer           :: error,ref_rpt,nrep,ref_blk,varno_ivar
    integer           :: ref_bkstp,nval,ivar,iexp,ilev,icount,icodtyp
    integer           :: date,time,ilat,ilon,iele,nele,icol
    real(8)           :: val,exponent

    ! initialize burp file, report, and block
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_READ, IOSTAT=error)
    
    if (error == 0) then
       write(*,*) "brpf_chem_read: Reading file " // trim(filename)
       write(*,*) "brpf_chem_read: Selecting STNID = ",stnid," BUFR = ",varno," block type = ",block_type
       write(*,*) "brpf_chem_read:           bkstp = ",bkstp," nlev = ",nlev," match_nlev = ",match_nlev
       if (present(codtyplist_opt)) write(*,*) "brpf_chem_read: CodeTypeList: ",codtyplist_opt(:)
    else
       call utl_abort('brpf_chem_read: Could not find/open BURP file: ' // trim(filename))
    end if

    ! get number of reports in file
    call BURP_Get_Property(brp, NRPTS=nrep)

    ! allocate memory
    if (ndim == 1) then
       call oss_obsdata_alloc(burp_out,nrep,dim1=nlev)
    else
       call oss_obsdata_alloc(burp_out,nrep,dim1=nlev,dim2_opt=nlev)
    end if
    
    icount = 0  ! counter of reports with same stnid, number of levels, and varno as input 
    ref_rpt = 0
    
    ! loop through reports    
    REPORTS: do

       ref_rpt = BURP_Find_Report(brp, REPORT=rep, SEARCH_FROM=ref_rpt, IOSTAT=error)

       if (ref_rpt<0) exit REPORTS
       
       call BURP_Get_Property(rep, STNID=rep_stnid, DATE=date, TEMPS=time, LATI=ilat, LONG=ilon, IDTYP=icodtyp) 

       if (present(codtyplist_opt)) then
          if (.not.any(codtyplist_opt(:) == icodtyp)) cycle REPORTS
       end if

       if (.not.utl_stnid_equal(stnid,rep_stnid)) cycle REPORTS

       ! loop through blocks
       ref_blk = 0
       BLOCKS: do
          
          ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
          if (ref_blk<0) exit BLOCKS
          
          call BURP_Get_Property(blk, NELE=nele, NVAL=nval, BKSTP=ref_bkstp, IOSTAT=error)

          if (.not.IS_Burp_Btyp(trim(block_type),BLOCK=blk) .or. bkstp /= ref_bkstp .or. (match_nlev.and.nval /= nlev)) cycle BLOCKS

          if (varno > 0) then
             ivar = BURP_Find_Element(blk, ELEMENT=varno, IOSTAT=error)
             if (ivar < 0) cycle BLOCKS
          else 
             ! Search for first data element within elements 10000 and 16000.
             varno_ivar=-1
             do ivar=1,nele
                varno_ivar=BURP_Get_Element(blk, INDEX=ivar, IOSTAT=error)
                if (varno_ivar >= 10000.and.varno_ivar < 16000) exit
             end do
             if (varno_ivar < 10000.or.varno_ivar >= 16000) call utl_abort('brpf_chem_read: No valid element found for STNID ' // rep_stnid )
          end if

          ! required block found if code reaches this point, retrieve data and store in burp_out
          
          if (nval > nlev) call utl_abort('brpf_chem_read: number of levels in the report (' // trim(utl_str(nval)) // &
                                         ') exceeds the specified maximum number of levels (' // trim(utl_str(nlev)) // &
                                         ') for STNID ' // rep_stnid )

          icount=icount+1
          burp_out%code(icount) = oss_obsdata_get_header_code(ilon,ilat,date,time,rep_stnid)  ! this code is a unique identifier for this report

          if (ndim == 1) then
             ! retrieve 1D data

             iexp = BURP_Find_Element(blk, ELEMENT=BUFR_SCALE_EXPONENT, IOSTAT=error)
                
             if (iexp < 0) then
                ! No exponent found in block
                do ilev=1,nval                   
                   burp_out%data1d(ilev,icount) = BURP_Get_Rval(blk, NELE_IND=ivar, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)                
                end do
             else
                ! Apply exponent
                do ilev=1,nval                   
                   val = BURP_Get_Rval(blk, NELE_IND=ivar, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)                
                   exponent = BURP_Get_Rval(blk, NELE_IND=iexp, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   burp_out%data1d(ilev,icount) = val * 10**exponent
                end do
             end if
                   
          else if (ndim == 2) then
             ! retrieve 2D data

             icol = 0
             do iele=1,nele
                ivar = BURP_Get_Element(blk, INDEX=iele, IOSTAT=error)
                if (ivar == varno) then
                   icol = icol+1
                   do ilev=1,nval                  
                      burp_out%data2d(ilev,icol,icount) = BURP_Get_Rval(blk, NELE_IND=iele, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   end do
                end if
             end do

          end if

          exit BLOCKS
          
       end do BLOCKS
      
    end do REPORTS

    ! resize first dimension of data arrays from length of nrep to icount
    call utl_resize(burp_out%code,icount)
    if (ndim == 1) then
       call utl_resize(burp_out%data1d,nlev,icount)
    else if (ndim == 2) then
       call utl_resize(burp_out%data2d,nlev,nlev,icount)
    end if

    burp_out%nrep = icount

    write(*,*) "brpf_chem_read: Reading of file complete. Number of reports found: ",burp_out%nrep
    
    ! deallocate
    Call BURP_Free(brp,iostat=error)
    Call BURP_Free(rep,iostat=error)
    Call BURP_Free(blk,iostat=error)
    
  end function brpf_chem_read

!--------------------------------------------------------------------------
!! *Purpose*: Add or modify information from all BURP files. Currently only used for
!!            chemical constituent but can potentially be used for any family.
!!            If BURP files are not split, a single files is updated. If BURP files
!!            are split, then each split file is updated. See the function brpf_chem_update
!!            for more details.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Input:
!!v           obsfam        observation family name
!!v           varno         BUFR descriptors. Number of elements must be 
!!v                         max(1,obsdata%dim2)
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           obsdata       Input struct_oss_obsdata object for varno.
!!v           multi         Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
!!v           fname_prefix  filename prefix, e.g. brpua_ (optional)
!!v
!!
!! Output: 
!!v           nrep_modified Number of modified reports
!
!! Comment:
!!
!!v   This routine is general enough to be used by observation families
!!v   other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function brpf_chem_update_all(obsfam,varno,bkstp,block_type,obsdata,multi_opt, &
                                obsFilesSplit_opt) result(nrep_modified)
    implicit none

    character(len=4), intent(in) :: block_type
    character(len=*), intent(in) :: obsfam
    type(struct_oss_obsdata), intent(inout) :: obsdata
    integer, intent(in) :: varno(:),bkstp    
    character(len=*), intent(in), optional :: multi_opt
    logical, optional :: obsFilesSplit_opt
    integer :: nrep_modified

    integer :: ierr,nrep_modified_global
    logical :: obsFilesSplit

    if ( present(obsFilesSplit_opt) ) then
      obsFilesSplit = obsFilesSplit_opt
    else
      obsFilesSplit = .true.
    end if

    if (obsFilesSplit .or. mpi_myid == 0) then
       nrep_modified = brpf_chem_update(brpf_get_filename(obsfam),varno,bkstp,block_type,obsdata,multi_opt=multi_opt)
    end if

    if (obsFilesSplit) then
       call rpn_comm_allreduce(nrep_modified,nrep_modified_global,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
       nrep_modified = nrep_modified_global
    end if

  end function brpf_chem_update_all

!--------------------------------------------------------------------------
!! *Purpose*: Add or modify information from BURP file in existing block and
!!            for specified BUFR descriptor varno(s). Provided data can be either
!!            1D or 2D data. Currently only used for chemical constituent but can
!!            potentially be used for any family.
!!
!! @author Y. Rochon, ARQI/AQRD, June 2016 (partly based on brpf_chem_read by M. Sitwell)
!!
!! Revisions:
!!v        M. Sitwell, ARQI/AQRD, Aug 2016
!!v          - Modified to preserve order of reports.
!!v        Y. Rochon, ARQI/AQRD, Jan 2017
!!v          - Added optional check for multi using 'DATA' when block_type='INFO'.
!!
!! Input:
!!v           filename      BURP file name
!!v           varno         BUFR descriptors. Number of elements must be 
!!v                         max(1,obsdata%dim2)
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           obsdata       Input struct_oss_obsdata object for varno.
!!v           multi         Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
!!
!! Output:
!!v           nrep_modified Number of modified reports
!!
!! Comments:
!!v  - Currently assumes that all elements of varno(:) are distinct from each other.
!!v  - In blocks with new data to be added/modified, if the varno already exists in the block, the
!!v    new data will overwrite the existing varno data, otherwise will append the new data
!!v    to the block.
!!v  - The settings for BURP_Write_Block should have ENCODE_BLOCK and CONVERT_BLOCK set to
!!v    .true. in all cases, including when the block has not been modified, due to problems
!!v    that can occur when writing blocks containing negative integers with datyp=4.
!!v  - This routine is general enough to be used by observation families
!!v    other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function brpf_chem_update(filename,varno,bkstp,block_type,obsdata,multi_opt) result(nrep_modified)

    implicit none

    character(len=*), intent(in) :: filename
    character(len=4), intent(in) :: block_type
    type(struct_oss_obsdata), intent(inout) :: obsdata
    integer, intent(in) :: varno(:),bkstp
    
    character(len=*), intent(in), optional :: multi_opt

    integer :: nrep_modified,ncount
    logical :: blk_found
    integer, parameter :: LNMX=100000, code_len=90

    character(len=9)  :: stnid
    character(len=code_len) :: code    ! Must be at least as large as burp_code_len
    type(burp_file)   :: brp
    type(burp_rpt)    :: rep,rep_new
    type(burp_block)  :: blk
    integer           :: error,ref_rpt,nrep,ref_blk,ndim,dim1,dim2
    integer           :: ref_bkstp,nval,ivar,ilev,istat
    integer           :: date,time,ilat,ilon,iele,nele,k
    integer, allocatable :: address(:)
    real(4), allocatable :: new_vals(:,:,:)
    logical, allocatable :: modify(:)
    
    ! Check presence of data to update
    if (obsdata%nrep <= 0) then
       write(*,*) 'brpf_chem_update: Skipped due to absence of data to update.'
       return
    end if
    
    ! Identify dimensions for the input data    
    ndim=obsdata%ndim
    dim1=obsdata%dim1
    if (ndim == 1) then
       dim2=1
    else
       dim2=obsdata%dim2
    end if
    
    if (size(varno) < dim2) call utl_abort('brpf_chem_update: Number of BUFR elements not sufficient. ' // &
                                          trim(utl_str(size(varno))) // ' vs ' // trim(utl_str(dim2)))

    if (code_len < oss_obsdata_code_len()) call utl_abort('brpf_chem_update: Length of code string' &
                                          // ' needs to be increased to ' // trim(utl_str(oss_obsdata_code_len())))
     
    ! initialize burp file, report, and block system resources
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, R2=rep_new, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file in append mode (to replace or add data in a block)
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_APPEND, IOSTAT=error)
    if (error /= 0) call utl_abort('brpf_chem_update: Could not open BURP file: ' // trim(filename))

    ! get number of reports in file
    call BURP_Get_Property(brp, NRPTS=nrep)

    allocate(address(nrep),modify(nrep),new_vals(dim1,dim2,nrep))
    address(:)=0
    modify(:)=.false.
    new_vals(:,:,:)=0.

    ! First loop through reports to identify addresses of original file as well as identify if new
    ! information should be included to that report.
    ! NOTE: The addresses of all reports have to be saved in their original order to ensure the
    !       order of the reports in the file is unchanged.
    ref_rpt=0
    ncount=0
    obsdata%irep=1
    REPORTS1: do

       ref_rpt = BURP_Find_Report(brp, REPORT=rep, SEARCH_FROM=ref_rpt, IOSTAT=error)
       if (ref_rpt<0) exit REPORTS1

       ncount=ncount+1
       address(ncount)=ref_rpt

       call BURP_Get_Property(rep, STNID=stnid, DATE=date, TEMPS=time, LATI=ilat, LONG=ilon)

       if (stnid(1:2) == '>>') cycle REPORTS1

       ! Get unique identifier for search from input data
       code = oss_obsdata_get_header_code(ilon,ilat,date,time,stnid)
       
       ! Determine if replacement/additional data likely present for this report
       if (dim1 == 1.and.dim2 == 1) then
          new_vals(1,1,ncount)=real(oss_obsdata_get_element(obsdata,trim(code),1,stat_opt=istat))
       else if (dim2 == 1) then
          new_vals(:,1,ncount)=real(oss_obsdata_get_array1d(obsdata,trim(code),stat_opt=istat))
       else 
          new_vals(:,:,ncount)=real(oss_obsdata_get_array2d(obsdata,trim(code),stat_opt=istat))
       end if

       if (istat == 0) then
          if (present(multi_opt)) then
             ! loop through blocks to find first data block
             ref_blk = 0
             BLOCKS1: do
                ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
                if (ref_blk<0) exit BLOCKS1
                if (IS_Burp_Btyp('DATA',BLOCK=blk)) then
                   if (IS_Burp_Btyp(trim(multi_opt),BLOCK=blk)) modify(ncount) = .true.
                   exit BLOCKS1
                end if
             end do BLOCKS1
          else
             modify(ncount) = .true.
          end if
       end if

    end do REPORTS1
    
    nrep_modified = count(modify)   ! number of reports with same code and, possibly, same number of obs data levels

    ! Generate new report
    Call BURP_New(rep_new, ALLOC_SPACE=10*LNMX, IOSTAT=error)

    ! second loop through reports to include the new information to the file    
    REPORTS2: do k=1,ncount
    
       call BURP_Get_Report(brp, REPORT=rep, REF=address(k), IOSTAT=error)
       
       ! Copy report header
       Call BURP_Copy_Header(TO=rep_new,FROM=rep)
       Call BURP_Init_Report_Write(brp,rep_new,IOSTAT=error)

       ! loop through blocks
       ref_blk = 0
       BLOCKS: do
          
          ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
          if (ref_blk<0) exit BLOCKS
          
          if (modify(k)) then

             call BURP_Get_Property(blk, NELE=nele, NVAL=nval, BKSTP=ref_bkstp, IOSTAT=error)

             blk_found = IS_Burp_Btyp(trim(block_type),BLOCK=blk) .and. bkstp == ref_bkstp .and. dim1 == nval
         
             if (blk_found) then
                ! Block to be modified has been found, add new data to block.
                ! If the varno is already in the block, the new data will overwrite the
                ! existing data, otherwise will append the new data to the block.

                do iele=1,dim2
                   ivar = BURP_Find_Element(blk, ELEMENT=varno(iele), IOSTAT=error)           
                   if (ivar < 0) then
                      ivar=nele+1
                      call BURP_Resize_Block(blk,ADD_NELE=1,IOSTAT=error)
                      call BURP_Set_Element(blk,NELE_IND=ivar,ELEMENT=varno(iele),IOSTAT=error)
                   end if
                
                   do ilev=1,nval 
                      call BURP_Set_Rval(blk,NELE_IND=ivar,NVAL_IND=ilev,NT_IND=1,RVAL=new_vals(ilev,iele,k),IOSTAT=error)                 
                   end do
                end do
        
             end if
          end if
             
          call BURP_Write_Block(rep_new, BLOCK=blk, ENCODE_BLOCK=.true., CONVERT_BLOCK=.true., IOSTAT=error)
         
       end do BLOCKS
       
       call BURP_Delete_Report(brp,rep,IOSTAT=error)
       call BURP_Write_Report(brp,rep_new,IOSTAT=error) 
  
    end do REPORTS2
        
    ! deallocate
    deallocate(address,modify,new_vals)
    Call BURP_Free(brp,iostat=error)
    Call BURP_Free(rep,R2=rep_new,iostat=error)
    Call BURP_Free(blk,iostat=error)
    
  end function brpf_chem_update

!--------------------------------------------------------------------------
!! *Purpose*: Returns the BURP file name assigned to the calling processor. If the
!!            input family has more than one file, the first file found will be
!!            returned. File names are assigned in the module burpFiles_mod.
!!
!! @author M. Sitwell  Sept 2016
!!
!! Input:
!!v    obsfam            observation family name
!!
!! Output:
!!v    burp_filename  file name of associated BURP file
!!v    found          logical indicating if the BURP file could be found (optional)
!--------------------------------------------------------------------------
  function brpf_get_filename(obsfam,found_opt) result(burp_filename)

    implicit none

    character(len=2), intent(in) :: obsfam
    logical, intent(out), optional :: found_opt
    character(len=128) :: burp_filename
    
    logical :: file_found
    integer :: ifile

    burp_filename = ""
    file_found = .false.
       
! THIS NEEDS TO BE MOVED, ALONG WITH THE CODE THAT CALLS IT
! PROBABLY TO obsFiles_mod.f90 (Mark Buehner)
!    do ifile=1,burp_nfiles
!       if (obsfam == burp_cfamtyp(ifile)) then
!          burp_filename = burp_cfilnam(ifile)
!          inquire(file=trim(burp_filename), exist=file_found)
!          exit
!       end if
!    end do

    if (.not.file_found) write(*,*) "brpf_get_filename: File not found for observation family " // trim(obsfam)

    if (present(found_opt)) found_opt = file_found

  end function brpf_get_filename

end module burpFiles_mod
