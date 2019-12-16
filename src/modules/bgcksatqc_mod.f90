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

module bgcksatqc_mod
  ! MODULE bgckmicrowave_mod (prefix='mwbg' category='1. High-level functionality')
  !
  ! :Purpose: Variables for microwave background check and quality control.
  !
  use burp_module

  implicit none
  save
  private

  ! public variables
  public :: scanpos, ztb, biasCorr, zlat,zlon,zenith, ilq,itt, ican, qcflag2, qcflag1
  public :: mxele, mxval, mxnt, nchan, mxsat, mxan, mxor, mxscan, zmisg
  public :: reportIndex
  public :: nval
  public :: satqc_debug, satqc_modlsqtt, satqc_useUnbiasedObsForClw 

  ! Public functions
  public :: satqc_getData, satqc_landIceMaskAtms, satqc_grossValueCheck
  public :: satqc_firstQcCheckAtms, satqc_nrlFilterAtms, satqc_writeBlocks

  ! Set array limits (ATMS: 22 chan, 96 FOV):
  integer, parameter :: mxele=24,mxval=30,mxnt=2800
  integer, parameter :: nchan=22
  integer, parameter :: mxsat=3
  integer, parameter :: mxan=22
  integer, parameter :: mxor=20
  integer, parameter :: mxscan=96

  ! Other variables:
  real, parameter    :: zmisg=9.9e09

  integer  :: error, reportIndex, nval

  integer, dimension(mxnt)             :: scanpos
  real, dimension(mxval*mxnt)          :: ztb, biasCorr
  real, dimension(mxnt)                :: zlat,zlon,zenith
  integer, dimension(mxnt)             :: ilq,itt
  integer, dimension(mxval*mxnt)       :: ican, qcflag2
  integer, dimension(mxnt,3)           :: qcflag1

  logical :: satqc_debug, satqc_modlsqtt, satqc_useUnbiasedObsForClw 

contains

  subroutine satqc_getData(rpt)
    !--------------------------------------------------------------------------------------
    ! Object:   This routine extracts the needed data from the blocks in the report:
    !             zenith(nt)       = satellite zenith angle (btyp=3072,ele=7024)
    !             ilq(nt)          = land/sea qualifier     (btyp=3072,ele=8012)
    !             itt(nt)          = terrain-type (ice)     (btyp=3072,ele=13039)
    !             zlat(nt)                                  (btyp=5120,ele=5002)
    !             zlon(nt)                                  (btyp=5120,ele=6002)
    !             ztb(nt*nval),    = brightness temperature (btyp=9248/9264,ele=12163) 
    !             scanpos(nt),     = scan position (fov)    (btyp=3072,ele=5043) 
    !             nval,            = number of channels     (btyp=9248/9264)
    !             nt,              = number of locations    (btyp=5120,etc.)
    !             qcflag1(nt,3)    = flag values for btyp=3072 block ele 033078, 033079, 033080
    !             qcflag2(nt*nval) = flag values for btyp=9248 block ele 033081
    !             ican(nt*nval)    = channel numbers btyp=9248 block ele 5042 (= 1-22)

    ! NOTE:  reportIndex = report number (from MAIN program) **** DO NOT MODIFY ****
    !       kk = variable for loops over locations (nt)
    !        j = variable for loops over nval (nval = 1 or nchan)

    type(BURP_RPT)         :: rpt
    type(BURP_BLOCK)       :: blk

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: indice, indice1, indice2, indice3, indice4, kk, j
    integer :: ind_lat,ind_lon

    integer :: ipos

    ! 1) Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 5120, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR -  Location/time (3D) block (btyp=5120) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()
   
    ind_lat = BURP_Find_Element(blk,5002,IOSTAT=error)
    ind_lon = BURP_Find_Element(blk,6002,IOSTAT=error)
    
    if ( (ind_lat > 0) .AND. (ind_lon > 0) ) then
      j = 1
      do kk =1, my_nt
        zlat(kk) = BURP_Get_Rval(blk,ind_lat,j,kk,error)
        zlon(kk) = BURP_Get_Rval(blk,ind_lon,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - lat, lon elements (5002,6002) missing in 3D block (btyp=5120). Report = ', reportIndex
      call abort()
    endif

    ! 2) Get info elements from the INFO block   BTYP = 3072
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 3072, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR - INFO block (btyp=3072) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,7024,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        zenith(kk) = BURP_Get_Rval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - satellite zenith angle missing in INFO block (ele=7024). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,8012,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        ilq(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - land/sea qualifier missing in INFO block (ele=8012). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,13039,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        itt(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - terrain-type missing in INFO block (ele=13039). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,5043,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        scanpos(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - scan position missing in INFO block (ele=5043). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,33078,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,1) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Geolocation quality code missing in INFO block (ele=33078). Report = ', reportIndex
      call abort()
    endif  

    indice = BURP_Find_Element(blk,33079,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,2) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Granule level QC flag missing in INFO block (ele=33079). Report = ', reportIndex
      call abort()
    endif  

    indice = BURP_Find_Element(blk,33080,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,3) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Scan level QC flag missing in INFO block (ele=33080). Report = ', reportIndex
      call abort()
    endif  

    ! 3) Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = nchan)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9248, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9264, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - DATA block (btyp 9248 or 9264) not found in report number ', reportIndex
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    nval = my_nval    ! set nval (#channels) for MAIN program

    indice1 = BURP_Find_Element(blk,12163,IOSTAT=error)
    indice2 = BURP_Find_Element(blk,33081,IOSTAT=error)
    indice3 = BURP_Find_Element(blk,5042,IOSTAT=error)
    indice4 = BURP_Find_Element(blk,12233,IOSTAT=error)
   
    if ( indice1 > 0 .and. indice2 > 0 .and. indice3 > 0 .and. indice4 > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          ztb(ipos)     = BURP_Get_Rval(blk,indice1,j,kk,error)
          qcflag2(ipos) = BURP_Get_Tblval(blk,indice2,j,kk,error)
          ican(ipos)    = BURP_Get_Tblval(blk,indice3,j,kk,error)
          biasCorr(ipos)= BURP_Get_Rval(blk,indice4,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Elements are missing in DATA block. Report = ', reportIndex
      if ( indice1 <= 0 )  write(*,*) '       Tb data             (012163) are missing!'
      if ( indice2 <= 0 )  write(*,*) '       Data level QC flags (033081) are missing!'
      if ( indice3 <= 0 )  write(*,*) '       Channel numbers     (005042) are missing!'
      if ( indice4 <= 0 )  write(*,*) '       Bias correction     (012233) are missing!'
      call abort()
    endif 

    return

  end subroutine satqc_getData


  subroutine satqc_writeBlocks(lsq,trn,riwv,rclw,ident,logicalFlags,lutb,rpt,rpt_out)
    ! Object:   This routine modifies the blocks in the input Report (rpt) 
    integer, intent(in), dimension(:)   :: lsq
    integer, intent(in), dimension(:)   :: trn
    real,    intent(in), dimension(:)   :: riwv
    real,    intent(in), dimension(:)   :: rclw
    integer, intent(in), dimension(:)   :: ident
    logical, intent(in), dimension(:,:) :: logicalFlags
    logical :: lutb
    type(BURP_RPT)         :: rpt
    type(BURP_RPT)         :: rpt_out

    type(BURP_BLOCK)       :: blk, blk_copy

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_btyp, my_bfam, iidata
    integer :: indice, indice1, indice2, j, kk, ipos

    Call BURP_Init(blk, B2=blk_copy, IOSTAT=error)

    ! 1) Read and modify the blocks in rpt and add them to rpt_out
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(rpt,BLOCK= blk,SEARCH_FROM= ref_blk,IOSTAT= error)
      if (error /= burp_noerr) call abort()
      
      if (ref_blk < 0) Exit

      Call BURP_Get_Property(blk, &
                  NELE   = my_nele, &
                  NVAL   = my_nval, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                  NT     = my_nt, &         ! 1 or number of locations in block
                  BTYP   = my_btyp, &
                  BFAM   = my_bfam, &
                  IOSTAT = error)
      if (error /= burp_noerr) call abort()

      ! 3D Block
      if (my_btyp == 5120) then     
      ! Set bit 6 in 24-bit global flags if any data rejected

      ! Extract the global flags, element 55200
      indice = BURP_Find_Element(blk,55200,IOSTAT=error)
        
      if ( indice > 0 ) then
        j = 1
        do kk =1, my_nt
          iidata = BURP_Get_Tblval(blk,indice,j,kk,error)
          if ( ANY(logicalFlags(kk,:)) ) iidata = IBSET(iidata,6)
          Call BURP_Set_Tblval(blk,indice,j,kk,iidata)
        end do
      else
        write(*,*) 'ERREUR - Global flag missing in 3D block (ele=55200). Report = ', reportIndex
        call abort()
      endif 
        
      blk_copy = blk
      Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
      if (error /= burp_noerr)  call abort()

      ! INFO block (mix of integer and real data)
      elseif (my_btyp == 3072) then


        ! Add new elements QC indent flag (ident), CLW (rclw) and ECMWF_SI (riwv)
        !   OPTION: replace land-sea qualifier (lsq) and terrain type (trn) with internal values
        !   NOTE: after setting real values (Rval), call BURP_Convert_Block()!!

        Call BURP_Resize_Block(blk, ADD_NELE = 3, IOSTAT = error)
        if (error /= burp_noerr)  call abort()
        
        Call BURP_Set_Element(blk, NELE_IND = my_nele+1, ELEMENT = 25174, IOSTAT = error)
        Call BURP_Set_Element(blk, NELE_IND = my_nele+2, ELEMENT = 13209, IOSTAT = error)
        Call BURP_Set_Element(blk, NELE_IND = my_nele+3, ELEMENT = 13208, IOSTAT = error)
        Call BURP_Encode_Block(blk)   ! encode the element numbers in the block

        j = 1
        do kk =1, my_nt
          iidata = ident(kk)
          Call BURP_Set_Rval  (blk, NELE_IND=my_nele+2,NVAL_IND=j,NT_IND=kk, RVAL=rclw(kk),IOSTAT=error)
          Call BURP_Set_Rval  (blk, NELE_IND=my_nele+3,NVAL_IND=j,NT_IND=kk, RVAL=riwv(kk),IOSTAT=error)
        end do
      
        Call BURP_Convert_Block(blk)

        j = 1
        do kk =1, my_nt
          iidata = ident(kk)
          Call BURP_Set_Tblval(blk, NELE_IND=my_nele+1,NVAL_IND=j,NT_IND=kk, TBLVAL=iidata,IOSTAT= error)
        end do

        if (satqc_modlsqtt) then
          indice1 = BURP_Find_Element(blk,  8012, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          indice2 = BURP_Find_Element(blk, 13039, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          if ( indice1 > 0 .and. indice2 > 0 ) then
            j = 1
            do kk =1, my_nt
              iidata = lsq(kk)
              Call BURP_Set_Tblval(blk,indice1,j,kk,iidata,error)
              if (error /= burp_noerr)  call abort()
              iidata = trn(kk)
              Call BURP_Set_Tblval(blk,indice2,j,kk,iidata,error)
              if (error /= burp_noerr)  call abort()
            end do
          else
            write(*,*) 'ERREUR - land/sea qualifier (ele=8012) and/or terrain type (ele=13039) not found in INFO block. Report = ', reportIndex
            call abort()
          endif                     
        endif

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()      

      !  DATA block
      elseif (my_btyp == 9248 .or. my_btyp ==9264) then 
        ! Modify Tb data if any data (ztb) were set to zmisg (lutb=.true.)

        if (lutb) then
          indice = BURP_Find_Element(blk, 12163, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          if ( indice > 0 ) then
            ipos = 0
            do kk =1, my_nt
              do j = 1, my_nval
                ipos = ipos + 1
                Call BURP_Set_Rval(blk,NELE_IND=indice,NVAL_IND=j,NT_IND=kk,RVAL=ztb(ipos),IOSTAT=error)  
                if (error /= burp_noerr)  call abort()
              enddo
            enddo
          else
            write(*,*) 'ERREUR - Cannot find Tb (ele=12163) in DATA block!. Report = ', reportIndex
            call abort()
          endif
!          if ( satqc_debug ) Call BURP_TO_STDOUT(blk, CONVERT =.FALSE.)
          blk_copy = blk
          Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
          if (error /= burp_noerr)  call abort() 
        else
          Call BURP_Write_Block(rpt_out,BLOCK=blk,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
          if (error /= burp_noerr)  call abort() 
        endif

      ! FLAG block
      elseif (my_btyp == 15392 .or. my_btyp == 15408) then 
        ! Modify data flag values (set bit 7) for rejected data    

        indice = BURP_Find_Element(blk, 212163, IOSTAT=error)
        if (error /= burp_noerr)  call abort()
        if ( indice > 0 ) then 
          do kk =1, my_nt
            do j = 1, my_nval
              iidata = BURP_Get_Tblval(blk,indice,j,kk,error)
              if (logicalFlags(kk,j)) iidata = IBSET(iidata,7)
              Call BURP_Set_Tblval(blk,indice,j,kk,iidata)
            enddo
          enddo
        else
          write(*,*) 'ERREUR - Data QC flags (ele=212163) not found in FLAG block. Report = ', reportIndex
          call abort()      
        endif

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! OTHER BLOCK 
      else 
        Call BURP_Write_Block(rpt_out,BLOCK=blk,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort() 
        
      endif      

    enddo BLOCKS

    return

  end subroutine satqc_writeBlocks


  subroutine satqc_landIceMaskAtms(mglg_file,npts,zlat,zlon,zlq,ztt,waterobs)
    ! Adapted from: land_ice_mask_ssmis.ftn90 of satqc_ssmis (D. Anselmo, S. Macpherson)
    !
    ! Object:   This routine sets waterobs array by performing a land/ice proximity check using
    !           using analysis MG and LG (or GL) fields used by the model which produces the trial field.
    !           The purpose of this check is to remove obs that reside close to coasts or ice,
    !           and so whose TBs may be contaminated.
    !           The GEM Global (glbhyb2) analysis contains MG and LG fields (on different grids).
    !
    !           NOTE: The 0.1 deg binary ice field check from land_ice_mask_ssmis.ftn90
    !           was removed. The land/sea qualifier (zlq) and terrain type (ztt) are modified
    !           to indicate proximity to land and sea-ice but are NOT changed in output BURP file.
    !
    !           In the application of this check, a 5x5 mesh, with spacing defined by rlat_km and
    !           rlon_km, is positioned with its center over an obs pt (2 grid pts on either side
    !           of the obs pt; size of mesh is equal to 4*rlat_km x 4*rlon_km). The values of MG
    !           and LG are evaluated at the grid points of this mesh. The maximum value of each
    !           determines whether the obs pt is too close to ice or land to be retained.
    !           **NOTE: the threshold value for MG has a very strong effect on the distance
    !                   from land that is permitted for an obs to be retained
    !
    !
    !      Maximum FOV             x---x---x---x---x     ^
    !         = 75km x 75km        |   |   |   |   |     |
    !         for Meso-sphere CHs  x---x---x---x---x     |
    !         = 74km x 47km        |   |   |   |   |     |
    !         for 19 GHz           x---x---o---x---x     | = 4*rlat_km
    !                              |   |   |   |   |     | = 4*40 km
    !                           ^  x---x---x---x---x     | = 160 km = 80 km north & south
    !                   rlat_km |  |   |   |   |   |     |
    !                           v  x---x---x---x---x     v
    !                                          <--->
    !                                         rlon_km
    !
    !                              <--------------->
    !                                 = 4*rlon_km
    !                                 = 4*40 km
    !                                 = 160 km = 80 km east & west
    !
    !
    !               MG value = 1.0  ==>  LAND       MG value = 0.0  ==>  OCEAN
    !               LG value = 1.0  ==>  ICE        LG value = 0.0  ==>  NO ICE
    !
    !
    ! Version:      Date:      Comment:
    ! --------      -----      --------
    !   0.1       16/08/12     Original adapted code.      S. Macpherson  
    !   0.2       01/03/14     Open mglg_file in R/O mode  S. Macpherson
    !
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! mglg_file  - input  -  name of file holding model MG and LG (or GL) fields
    ! npts       - input  -  number of input obs pts in report
    ! zlat       - input  -  array holding lat values for all obs pts in report
    ! zlon       - input  -  array holding lon values for all obs pts in report
    ! zlq        - in/out -  array holding land/sea qualifier values for all obs
    !                        pts of report (0 = land, 1 = sea)
    ! ztt        - in/out -  array holding terrain-type values for all obs pts
    !                        of current report (-1 land/open water, 0 = ice)
    ! waterobs   - output -  logical array identifying for each obs in current report
    !                        whether it is over open water, far from coast/ice
    ! mxlat      -internal-  number of grid pts in lat. direction for mesh
    ! mxlon      -internal-  number of grid pts in lon. direction for mesh
    ! rlat_km    -internal-  spacing desired between mesh grid points in km
    !                        along lat. direction
    ! rlon_km    -internal-  spacing desired between mesh grid points in km
    !                        along lon. direction
    ! dlat       -internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlat_km
    ! dlon       -internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlon_km
    ! rkm_per_deg -internal- distance in km per degree
    !                           = Earth radius * PI/180.0
    !                           = 6371.01 km * PI/180.0
    !                           = 111.195 km
    ! nlat,nlon  -internal-  used to define the lat/lon of the grid pts of mesh
    ! zlatbox    -internal-  lat values at all grid pts of mesh for all obs pts
    ! zlonbox    -internal-  lon values at all grid pts of mesh for all obs pts
    ! latmesh    -internal-  lat values at all grid pts of mesh for 1 obs pt
    ! lonmesh    -internal-  lon values at all grid pts of mesh for 1 obs pt
    ! mgintob    -internal-  interpolated MG values at all grid pts of mesh for 1 obs pt
    ! lgintob    -internal-  interpolated LG values at all grid pts of mesh for 1 obs pt
    ! mgintrp    -internal-  max. interpolated MG value on mesh for all obs pts
    ! lgintrp    -internal-  max. interpolated LG value on mesh for all obs pts
    ! MGthresh   -internal-  maximum allowable land fraction for obs to be kept
    ! LGthresh   -internal-  maximum allowable ice  fraction for obs to be kept
    !--------------------------------------------------------------------
    !  use var_declare
    implicit none

    ! Arguments:
    character(len=128), intent(in) :: mglg_file

    integer, intent(in)                   :: npts
    real,    intent(in),     dimension(:) :: zlat,zlon
    integer, intent(inout),  dimension(:) :: zlq, ztt

    logical, intent(out), dimension(:) :: waterobs

    ! Locals:
    integer, parameter :: mxlat=5,mxlon=5
    integer, parameter :: iungeo=50

    integer :: ier,key,istat
    integer :: ni,nj,nk,nilg,njlg
    integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    integer :: indx,ii,jj,kk
    integer :: nlat,nlon

    integer, dimension(10) :: alloc_status
  
    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195

    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1)  :: grtyp,grtyplg
  
    logical  :: llg

    ! F90 allocatable arrays:
    real, allocatable, dimension(:)   :: mg,lg
    real, allocatable, dimension(:)   :: latmesh,lonmesh
    real, allocatable, dimension(:)   :: mgintob,lgintob
    real, allocatable, dimension(:,:) :: zlatbox,zlonbox
    real, allocatable, dimension(:)   :: mgintrp,lgintrp
  
    ! RMNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid,gdidlg

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    integer :: idum1,idum2,idum3

    ! Allocate space for arrays holding values on mesh grid pts.
    alloc_status(:) = 0
    allocate ( latmesh(mxlat*mxlon), stat=alloc_status(1) )
    allocate ( lonmesh(mxlat*mxlon), stat=alloc_status(2) )
    allocate ( mgintob(mxlat*mxlon), stat=alloc_status(3) )
    allocate ( lgintob(mxlat*mxlon), stat=alloc_status(4) )
    allocate ( zlatbox(mxlat*mxlon,npts), stat=alloc_status(5) )
    allocate ( zlonbox(mxlat*mxlon,npts), stat=alloc_status(6) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'satqc_landIceMaskAtms: Memory allocation error '
      call abort()
    endif

    ! Open FST file.
    ier = fnom( iungeo,mglg_file,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    ! Read MG field.
    key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
    if ( key <  0 ) then
      write(*,*) 'satqc_landIceMaskAtms: The MG field is MISSING '
      call abort()
    end if

    allocate ( mg(ni*nj), stat=alloc_status(7) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'satqc_landIceMaskAtms: Memory allocation error '
      call abort()
    endif
    ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                idum18)


    ! Read LG field. Use GL field as backup.
    ! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.
    llg=.false.
    key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
    if ( key <  0 ) then
      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
      if ( key <  0 ) then
        write(*,*) 'satqc_landIceMaskAtms: No ice (LG or GL) fields found. Aborting! '
        call abort()
      else
        !write(*,*) 'satqc_landIceMaskAtms: The GL field was found and will be used.'
      endif
    else
      llg=.true.
    endif

    allocate ( lg(nilg*njlg), stat=alloc_status(8) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'satqc_landIceMaskAtms: Memory allocation error '
      call abort()
    endif
    
    if ( llg ) then
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
    else
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
    endif

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
                ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
                idum18)

    ! For each obs pt, define a grid of artificial pts surrounding it.
    nlat = ( mxlat - 1 ) / 2
    nlon = ( mxlon - 1 ) / 2

    dlat = rlat_km / rkm_per_deg
    do kk = 1, npts
      indx = 0

      do ii = -nlat, nlat
        rii = float(ii)
        xlat = zlat(kk) + rii*dlat
        xlat = max( -90.0, min(90.0,xlat) )
        xlatrad = xlat*pi/180.0

        do jj = -nlon, nlon
          dlon = rlon_km / ( rkm_per_deg*cos(xlatrad) )
          rjj = float(jj)
          indx = indx + 1
          xlon = zlon(kk) + rjj*dlon
          if ( xlon < -180. ) xlon = xlon + 360.
          if ( xlon >  180. ) xlon = xlon - 360.
          if ( xlon <    0. ) xlon = xlon + 360.
          zlatbox(indx,kk) = xlat
          zlonbox(indx,kk) = xlon
        end do

      end do
    end do


    ! Interpolate values from MG and LG field to grid pts of mesh centred over each obs pt.
    ! Determine for each obs pt, the max interpolated MG and LG value within the box
    ! surrounding it.
    ier    = ezsetopt('INTERP_DEGREE','LINEAR')
    gdid   = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
    gdidlg = ezqkdef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iungeo)

    allocate ( mgintrp(npts), stat=alloc_status(9) )
    allocate ( lgintrp(npts), stat=alloc_status(10) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'satqc_landIceMaskAtms: Memory allocation error '
      call abort()
    endif

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, npts

      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      mgintrp(kk) = maxval(mgintob(:))
      lgintrp(kk) = maxval(lgintob(:))

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.
    waterobs(:) = .false.   ! not over open water
    ztt(:) = -1             ! no ice (reset terain type)
    zlq(:) = 0              ! land   (reset land/sea qualifier)

    do kk = 1, npts
      if ( mgintrp(kk) < MGthresh ) zlq(kk) = 1  ! ocean point away from coast
      if ( lgintrp(kk) >= LGthresh .and. zlq(kk) == 1 ) ztt(kk) = 0  ! sea-ice affected point
      if ( lgintrp(kk)  < LGthresh .and. zlq(kk) == 1 ) then
        waterobs(kk) = .true.  ! water point not in close proximity to land or sea-ice
      end if
    end do

    ! Deallocate arrays and close FST file.
    alloc_status(:) = 0
    deallocate ( mgintrp, stat=alloc_status(1) )
    deallocate ( lgintrp, stat=alloc_status(2) )
    deallocate ( mg,      stat=alloc_status(3) )
    deallocate ( lg,      stat=alloc_status(4) )
    deallocate ( latmesh, stat=alloc_status(5) )
    deallocate ( lonmesh, stat=alloc_status(6) )
    deallocate ( mgintob, stat=alloc_status(7) )
    deallocate ( lgintob, stat=alloc_status(8) )
    deallocate ( zlatbox, stat=alloc_status(9) )
    deallocate ( zlonbox, stat=alloc_status(10) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'satqc_landIceMaskAtms: Memory deallocation error '
      call abort()
    endif
    ier = fstfrm(iungeo)
    ier = fclos(iungeo)

    return
  end subroutine satqc_landIceMaskAtms


  subroutine satqc_grossValueCheck(npts,ztb,grossrej)
    !  Object: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    !
    ! Variable Definitions:
    ! npts            - input  -  number of obs pts to process
    ! ztb             - input  -  Tbs from input BURP file
    ! grossrej        - output -  logical array defining which obs are to be rejected
    implicit none

    ! Arguments
    integer, intent(in) :: npts

    real,    intent(in),  dimension(:) :: ztb
    logical, intent(out), dimension(:) :: grossrej

    ! Locals
    integer :: ii, indx1, indx2

    grossrej(1:npts) = .true.
    indx1 = 1
    do ii = 1, npts

      indx2 = ii*nchan
      if ( all( ztb(indx1:indx2) > 50.0 ) .and. all( ztb(indx1:indx2) < 380.0 ) ) then
        grossrej(ii) = .false.
      end if
      indx1 = indx2 + 1

    end do

    return
  end subroutine satqc_grossValueCheck


  subroutine satqc_firstQcCheckAtms(zenith,ilq,itt,zlat,zlon,ztb,scanpos,stnid,nval,nt,lqc, &
               grossrej,lsq,trn,qcflag1,qcflag2,ican,blat,blon,lutb)
    !  This routine performs basic quality control checks on the data. It sets array
    !  lqc(nt,nchan) elements to .true. to flag data with failed checks.
    !
    !  The 7 QC checks are:
    !                 1) Invalid land/sea qualifier or terrain type,
    !                 2) Invalid field of view number,
    !                 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 5) Change in (computed) lsq,trn from (input) ilq,itt (from MG,LG fields)
    !                      ilq= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      itt=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      lsq= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      trn=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)
    !                 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !
    !  In most cases, lqc(ii,nchan) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (ztb) is set to missing_value 
    !  for checks 3 and 4 fails.
    !  use var_declare      ! nchan=22, zmisg, mxscan
    implicit none

    integer, intent(in), dimension(:)    :: ilq, itt, scanpos
    integer, intent(in), dimension(:)    :: ican, qcflag2
    integer, intent(in), dimension(:,:)  :: qcflag1
    integer, intent(in)                  :: nval, nt
    integer, intent(in), dimension(:)    :: lsq, trn

    logical, intent(in), dimension(:)    :: grossrej     ! dim(nt), true if 1 or more Tb fail gross error check
    logical, intent(out)                 :: lutb         ! true if Tb(ztb) are set to missing_value

    real, intent(in), dimension(:)       :: zlat, zlon
    real, intent(inout), dimension(:)    :: ztb, zenith
    
    integer, intent(in)                  :: blat, blon   ! NT box lat,lon (header)
     
    logical, intent(inout), dimension(:,:) :: lqc        ! dim(nt,nchan), lqc = .false. on input
     
    character(len=9), intent(in)         :: stnid
     
    !  Locals
    integer :: ii, jj, indx1, icount
    logical :: fail, fail1, fail2

    write(*,*) '============================================================================'
    write(*,*) 'satqc_firstQcCheckAtms: Processing data box: Stnid, lat, lon = ', stnid, blat, blon
    write(*,*) ' '

    lutb = .false.

    ! Global rejection checks

    ! Check if number of channels is correct
    if ( nval /= nchan ) then
      write(*,*) 'WARNING: Number of channels (',nval, ') is not equal to nchan (', nchan,')'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    endif

    ! Check for errors in channel numbers (should be 1-22 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,nt
      do jj = 1,nchan
        if ( ican(indx1+jj-1) /= jj ) fail = .true.
      enddo
      indx1 = indx1 + nchan
    enddo
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  ican(nt*nchan) array = ', ican(:)
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    endif

    ! 1) invalid land/sea qualifier or terrain type
    !  ilq = 0 (land),     1 (sea)
    !  itt = 0 (sea-ice), -1 otherwise
    !  lsq = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  trn = 0 (over or near analyzed sea-ice [LG]), -1 otherwise
    do ii = 1,nt
      fail = .false.
      if ( ilq(ii) < 0  .or. ilq(ii) > 2 ) fail = .true.
      if ( itt(ii) < -1 .or. itt(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
        write(*,*) '  ilq, itt, (lat, lon) = ', ilq(ii), itt(ii), '(',zlat(ii), zlon(ii),')'
      endif
      if ( ilq(ii) == 0 .and. itt(ii) == 0 ) then
        fail = .true.
        write(*,*) 'WARNING: Sea ice point (itt=0) at land point (ilq=0)!'
        write(*,*) ' lat, lon =  ', zlat(ii), zlon(ii)
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo

    do ii = 1,nt
      fail = .false.
      if ( lsq(ii) < 0  .or. lsq(ii) > 2 ) fail = .true.
      if ( trn(ii) < -1 .or. trn(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
        write(*,*) '  lsq, trn, (lat, lon) = ', lsq(ii), trn(ii), '(',zlat(ii), zlon(ii),')'
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo
 
    ! 2) invalid field of view number
    do ii = 1,nt
      fail = .false.
      if ( scanpos(ii) < 1  .or. scanpos(ii) > mxscan ) then
        fail = .true.
        write(*,*) 'WARNING: Invalid field of view! scanpos, lat, lon = ', scanpos(ii), zlat(ii), zlon(ii)
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo

    ! 3) satellite zenith angle missing or out of range (> 75 deg)
    !  If bad zenith, then set Tb (and zenith) = missing value
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zenith(ii) > 75.0 .or. zenith(ii) < 0. ) then
        fail = .true.
        write(*,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', zenith(ii), zlat(ii), zlon(ii)
        zenith(ii) = zmisg
        lutb = .true.
      endif
      do jj = 1,nchan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchan
    enddo

    ! 4) Lat,lon check
    ! Check for undecoded BURP file integer values of lat,lon = 0,0
    ! (usually associated with missing zenith angle and erroneous Tb=330K)

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zlat(ii) == -90.0  .and. zlon(ii) == -180.0 ) then
        fail = .true.
        icount =  icount + 1
        lutb = .true.
      endif
      do jj = 1,nchan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchan
    enddo
    if ( icount > 0 ) write(*,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( abs(zlat(ii)) > 90.0  .or. abs(zlon(ii)) > 180.0 ) then
        fail = .true.
        icount =  icount + 1
        lutb = .true.
      endif
      do jj = 1,nchan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchan
    enddo
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    icount = 0
    do ii = 1,nt
      fail = .false.
      if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) fail = .true.
      if ( fail ) then
        icount =  icount + 1
      endif
    enddo
    if ( icount > 0 ) write(*,*) 'INFO: Num. pts with land/sea qualifier or terrain type changed (MG,LG) = ', icount

    ! 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !  33078 Geolocation quality code     qcflag1(ii,1)  code value = 0-15 (0= OK, 15=misg)
    !  33079 Granule level quality flags  qcflag1(ii,2)  16 bit flag  (start bit 6(2^5)=32) (misg=2^16-1 = 65535)
    !  33080 Scan level quality flags     qcflag1(ii,3)  20 bit flag  (start bit 7(2^6)=64) (misg=2^20-1) 
    !  33081 Channel data quality flags   qcflag2        12 bit flag  (start bit 3(2^2)=4)  (misg=2^12-1)
    !
    !  See http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/2010edition/BUFRver16/BUFR_16_0_0_TableD.pdf

    indx1 = 1
    do ii = 1,nt 
      fail1 = .false.
      fail = .false.
      if ( (qcflag1(ii,1) > 0) .or. (qcflag1(ii,2) >= 32) .or. (qcflag1(ii,3) >= 64) ) then
         write(*,*) 'WARNING: INFO BLOCK QC flag(s) indicate problem with data'
         write(*,*) ' ele33078 = ',qcflag1(ii,1),' ele33079 = ',qcflag1(ii,2),' ele33080 = ', qcflag1(ii,3)
         write(*,*) ' lat, lon = ', zlat(ii), zlon(ii)
         fail1 = .true.
         if ( grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
      endif
      do jj = 1,nchan
        fail2 = .false.
        if ( qcflag2(indx1+jj-1) >= 4 ) then
          !write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 = ', qcflag2(indx1+jj-1)
          !write(*,*) '    Lat, lon, channel = ', zlat(ii), zlon(ii), ican(indx1+jj-1)
          fail2 = .true.
          fail = .true.
          !if ( (.not. fail1) .and. grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
        endif
        if ( fail2 .or. fail1 ) lqc(ii,jj) = .true.
      enddo
      if ( fail ) write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 >= 4 for one or more channels! lat, lon = ', zlat(ii), zlon(ii)
      indx1 = indx1 + nchan
    enddo
     
    write(*,*) 'satqc_firstQcCheckAtms: Total number of data processed in this box = ', nt*nchan
    write(*,*) '         Total number of data flagged in this box   = ', COUNT(lqc)
    write(*,*) ' '

    return
  end subroutine satqc_firstQcCheckAtms


  subroutine satqc_nrlFilterAtms(ier, ni, tb23, bcor23, tb31, bcor31, tb50, bcor50, &
                   tb89, bcor89, tb165, bcor165, pangl, plat, ilansea, iglace, &
                   waterobs, grossrej, clw, si_ecmwf, si_bg, iNumSeaIce, iRej, SeaIce)
    !OBJET          Compute the following parameters using 5 ATMS channels:
    !                  - sea ice, 
    !                  - cloud liquid water (clw), 
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                o  open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and iglace (itt or terrain type) is changed accordingly
    !                o  clw are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                o  clw and si only computed over open water away from coasts and sea-ice
    !                o  clw and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      ier         - output - error return code for each location:
    !                                        0, ok,  
    !                                        1, input parameter out of range or grossrej=.true. 
    !               ni          - input  -  number of points to process (= NT)
    !               tb23        - input  -  23Ghz brightness temperature (K) -- ch. 1
    !               tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2
    !               tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3
    !               tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               pangl       - input  -  satellite zenith angle (deg.)
    !               plat        - input  -  latitude (deg.)
    !               ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
    !               iglace      - in/out -  terrain type (0=ice, -1 otherwise)
    !               waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               grossrej    - input  -  .true. if any channel had a gross error from satqc_grossValueCheck
    !               clw         - output -  cloud liquid water (kg/m**2) from tb23 & tb31
    !               si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
    !               si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
    !               iNumSeaIce  - in/out -  running counter for number of open water points
    !                                       with sea-ice detected (from algorithm)
    !               iRej        - in/out -  running counter for number of locations with bad
    !                                       pangl, plat, ilansea, or with grossrej=true
    !               SeaIce      - output -  computed sea-ice fraction from tb23 & tb50 
    !
    !               ice         - internal -  sea ice
    !             
    !
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is set to the missing value, i.e. -99.
    !
    implicit none

    integer    ::  i

    integer, intent(in)                   ::  ni
    integer, intent(inout)                ::  iNumSeaIce
    integer, intent(in), dimension(:)     ::  ilansea
    integer, intent(out), dimension(:)    ::  ier
    integer, intent(inout), dimension(:)  ::  iglace
    integer, intent(inout)                ::  iRej
    
    
    logical, intent(in), dimension(:)     ::  grossrej
    logical, intent(inout), dimension(:)  ::  waterobs

    real, intent(in),  dimension(:)  ::  tb23, tb31, tb50, tb89, tb165, pangl, plat
    real, intent(in),  dimension(:)  ::  bcor23, bcor31, bcor50, bcor89, bcor165
    real, intent(out), dimension(:)  ::  clw, si_ecmwf, si_bg, SeaIce 

    real, dimension(ni)              ::  ice

    real       ::  aa, deltb, abslat, cosz
    real       ::  t23, t31, t50, t89, t165

    real, parameter  ::  rmisg = -99.

    ier = 0

    ! 1) Initialise parameters:
    do i = 1, ni
      ice(i)      = rmisg
      clw(i)      = rmisg
      si_ecmwf(i) = rmisg
      si_bg(i)    = rmisg
      SeaIce(i)   = 0.0
    enddo

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( pangl(i)   .lt.   0.  .or. &
           pangl(i)   .gt.  70.  .or. &
           plat(i)    .lt. -90.  .or. & 
           plat(i)    .gt.  90.  .or. &  
           ilansea(i) .lt.   0   .or. & 
           ilansea(i) .gt.   1        ) then
         ier(i) = 1
      endif

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)       
      if ( grossrej(i) ) ier(i) = 1 

    enddo

    ! 3) Compute parameters:
    do i = 1, ni

      if ( ier(i) .eq. 0 ) then

        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))

        if ( satqc_useUnbiasedObsForClw ) then
          t23 = tb23(i)
          t31 = tb31(i)
          t50 = tb50(i)
          t89 = tb89(i)
          t165 = tb165(i)
        else
          t23 = tb23(i) - bcor23(i)
          t31 = tb31(i) - bcor31(i)
          t50 = tb50(i) - bcor50(i)
          t89 = tb89(i) - bcor89(i)
          t165 = tb165(i) - bcor165(i)
        end if
        deltb = t89 - t165

        ! Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
        if ( ilansea(i) .eq. 1 ) then  ! water point

          if ( abslat .lt. 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          endif
          
          SeaIce(i) = ice(i)
          
          if ( ice(i) .ge. 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            iglace(i) = 0
          endif
          
        endif

        ! Compute CLW and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 .lt. 284. .and. t31 .lt. 284. ) then
            aa = 8.24 - (2.622 - 1.846*cosz)*cosz
            clw(i) = aa + 0.754*alog(285.0-t23) - 2.265*alog(285.0-t31)
            clw(i) = clw(i)*cosz
            if ( clw(i) .lt. 0.0 ) clw(i) = 0.0
          endif
          si_ecmwf(i) = deltb - (-46.94 + 0.248*pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104*pangl(i))
        endif

      else  ! ier(i) .eq. 1 case
         iRej = iRej + 1

      endif ! if ( ier(i) .eq. 0 )

      if ( satqc_debug .and. (i .le. 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(*,*) ' ier(i),ice(i),clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      endif

    enddo   ! i loop over ni points

    return

  end subroutine satqc_nrlFilterAtms


  subroutine read_coeff(sats,chans,fovbias,coeff,nsat,nchan,nfov,npred,cinstrum,maxpred,iun,coeff_file,ptypes)
    ! max # of satellites, # of channels (24), # of FOV (60)
    !use var_declare,  only : mxsat, mxan, mxscan, fnom, fclos
    ! RETURNS:

    !   sats(nsat)            = satellite names
    !   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
    !   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i,j,1)          = regression constant
    !   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients

    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty

    !   maxpred (input) is max number of predictors
    implicit none

    ! IN
    integer                        :: maxpred, iun
    character(len=90)              :: coeff_file

    ! OUT
    character(len=9), dimension(mxsat) :: sats        ! dim(maxsat), satellite names
    integer*4, dimension(mxsat,mxan)   :: chans       ! dim(maxsat, maxchan), channel numbers
    real, dimension(mxsat,mxan,mxscan) :: fovbias     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
    real, dimension(mxsat,mxan,maxpred+1) :: coeff    ! dim(maxsat,maxchan,maxpred+1)
    integer                            :: nsat, nfov, nbscan
    integer, dimension(mxsat)          :: nchan       ! dim(maxsat), number of channels
    integer, dimension(mxsat,mxan)     :: npred       ! dim(maxsat, maxchan), number of predictors
    character(len=5)                   :: cinstrum    ! string: instrument (e.g. SSMIS)
    character(len=2), dimension(mxsat,mxan,maxpred)  :: ptypes ! dim(maxsat,maxchan,maxpred)

    ! LOCAL
    character(len=8)               :: sat
    character(len=120)             :: line
    integer*4                      :: chan
    integer                        :: ndata, nbfov, nbpred, i, j, k, ier, istat, ii
    logical                        :: newsat
    real                           :: dummy

    coeff    = 0.0
    fovbias  = 0.0
    sats     = 'XXXXXXXXX'
    cinstrum = 'XXXXX'
    chans    = 0
    npred    = 0
    nsat     = 0
    nchan    = 0
    nfov     = 0
    ptypes   = 'XX'

    nbscan = mxscan

    ier = FNOM(iun,coeff_file,'FMT',0)

    IF (ier == 0) THEN

      WRITE(*,*)
      WRITE(*,*) 'Bias correction coefficient file open = ', coeff_file

      READ(iun,*,IOSTAT=istat)
      IF ( istat < 0 ) THEN
        WRITE(*,*) '  ERROR- File appears empty.'
        RETURN
      END IF
      REWIND(iun)

      ii = 0

      ! Loop over the satellites/channels in the file
      do
        read(iun,'(A)',IOSTAT=istat) line
        if ( istat < 0 ) EXIT
        if ( line(1:3) == 'SAT' ) then
          newsat = .true.
          read(line,'(T53,A8,1X,A5,1X,I6,1X,I8,1X,I2,1X,I3)',IOSTAT=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading data from SATELLITE line in coeff file!'
            return
          endif
          do i = 1, mxsat
            if ( trim(sats(i)) == trim(sat) ) then
              newsat = .false.
              ii = i
            endif
          end do
          if ( newsat ) then
            ii = ii + 1
            if ( ii > mxsat ) then
              write(*,*) ' ERROR - max number of satellites exceeded in coeff file!'
              return
            endif
            sats(ii) = sat
            if (ii > 1) nchan(ii-1) = j
            j = 1
          else
            j = j + 1
          endif
          chans(ii, j) = chan
          npred(ii, j) = nbpred
          if ( nbpred > maxpred ) then
            write(*,*) ' ERROR - max number of predictors exceeded in coeff file!'
            return
          endif
          read(iun,'(A)',IOSTAT=istat) line
          if ( line(1:3) /= 'PTY' ) then
            write(*,*) ' ERROR - list of predictors is missing in coeff file!'
            return
          endif
          if ( nbpred > 0 ) then
            read(line,'(T8,6(1X,A2))',IOSTAT=istat) (ptypes(ii,j,k),k=1,nbpred)
            if ( istat /= 0 ) then
              write(*,*) ' ERROR - reading predictor types from PTYPES line in coeff file!'
              return
            endif
          endif
          read(iun,*,IOSTAT=istat) (fovbias(ii,j,k),k=1,nbfov)
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading fovbias in coeff file!'
            return
          endif
          if ( nbpred > 0 ) then
            read(iun,*,IOSTAT=istat) (coeff(ii,j,k),k=1,nbpred+1)
          else
            read(iun,*,IOSTAT=istat) dummy
          endif
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading coeff in coeff file!'
            return
          endif

        else
          EXIT
        endif

      end do

      if ( ii == 0 ) then
        write(*,*) ' ERROR - No data read from coeff file!'
        return
      endif

      nsat      = ii
      nfov      = nbfov
      nchan(ii) = j
      if ( nbscan /= 0 ) then
        if ( nfov /= mxscan ) then
          write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (mxscan).'
          write(*,*) '         nfov = ', nfov
          write(*,*) '       mxscan = ', mxscan
          !call abort()
        endif
      else ! nbscan = 0 case
        if ( nfov /= 1 ) then
          write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (1).'
          write(*,*) '         nfov = ', nfov
          !call abort()
        endif
      endif

      write(*,*) ' '
      write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
      write(*,*) ' '
      write(*,*) ' Number of satellites =     ', nsat
      write(*,*) ' Number of FOV =            ', nfov
      write(*,*) ' Max number of predictors = ', maxval(npred)
      write(*,*) ' '
      do i = 1, nsat
        write(*,*) '  Satellite = ' // sats(i)
        write(*,*) '     Number of channels = ', nchan(i)
        write(*,*) '     predictors, fovbias, coeff for each channel: '
        do j = 1, nchan(i)
          write(*,*) i, chans(i,j)
          if ( npred(i,j) > 0 ) then 
            write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
          else
            write(*,'(A)') 'No predictors'
          endif
          write(*,*) (fovbias(i,j,k),k=1,nfov)
          write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
        end do
      end do
      write(*,*) ' '

      ier = FCLOS(iun)

    ELSE

      write(*,*) 'READ_COEFF: ERROR - Problem opening the coeff file!'

    ENDIF

    RETURN

  end subroutine read_coeff


end module bgcksatqc_mod
