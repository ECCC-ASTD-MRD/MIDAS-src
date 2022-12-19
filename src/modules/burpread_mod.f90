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

module  burpread_mod
  ! MODULE burpread_mod (prefix='brpr' category='3. Observation input/output')
  !
  ! :Purpose: To read and update BURP observation files. Data is stored in 
  !           obsSpaceData object.
  !

use codePrecision_mod
use bufr_mod
use burp_module
use ObsSpaceData_mod
use MathPhysConstants_mod
use earthConstants_mod
use utilities_mod
use obsUtil_mod
use obsVariableTransforms_mod
use obsFilter_mod
use tovs_nl_mod
use kdtree2_mod
use codtyp_mod

implicit none
save

private

! public procedures
public :: brpr_readBurp, brpr_updateBurp, brpr_getTypeResume, brpr_addCloudParametersandEmissivity
public :: brpr_addElementsToBurp, brpr_updateMissingObsFlags, brpr_burpClean

! Namelist variables
INTEGER          :: NELEMS
INTEGER          :: NELEMS_SFC
INTEGER          :: BNBITSOFF
INTEGER          :: BNBITSON
INTEGER          :: BBITOFF(15)
INTEGER          :: BBITON(15)
INTEGER          :: NELEMS_GPS
INTEGER          :: LISTE_ELE_GPS(20)
INTEGER          :: BLISTELEMENTS(20)
INTEGER          :: BLISTELEMENTS_SFC(20)
INTEGER          :: BN_ITEMS
CHARACTER(len=3) :: BITEMLIST(20)
CHARACTER(len=7) :: TYPE_RESUME = 'UNKNOWN'
LOGICAL          :: ENFORCE_CLASSIC_SONDES
LOGICAL          :: UA_HIGH_PRECISION_TT_ES
LOGICAL          :: UA_FLAG_HIGH_PRECISION_TT_ES
LOGICAL          :: READ_QI_GA_MT_SW
logical          :: addBtClearToBurp
integer          :: clwFgElementId
integer          :: btClearElementId

CONTAINS

  character(len=7) function brpr_getTypeResume
    brpr_getTypeResume=TYPE_RESUME
  end function brpr_getTypeResume

  subroutine brpr_updateBurp(obsdat,familytype,brp_file,filenumb)
    !
    !:Purpose: To update variables relative to assimilation in burp files

    !***************************************************************************
    !
    !          WHEN SEARCHING FOR A SPECIFIC BLOCK BY ITS BTYP, VALUES OF
    !          BIT 0 TO 3 ARE IRRELEVANT WHILE BIT 4 IS 0 FOR GLOBAL AND 1
    !          FOR REGIONAL MODEL. HERE, WE SEARCH BLOCK BY THEIR FIRST
    !          10 BITS (BIT 5 TO 14).
    !
    !***************************************************************************

    IMPLICIT NONE

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat ! obsSpaceData object
    CHARACTER(len=*)       :: FAMILYTYPE ! type of family('UA','SF','AI','SW','TO', ...)
    CHARACTER(len=*)       :: BRP_FILE   ! name of burp file 
    integer                :: FILENUMB

    ! Locals:
    type(kdtree2), pointer            :: tree
    integer, parameter                :: maxNumSearch = 100
    integer                           :: numFoundSearch, bodyCount, resultIndex
    type(kdtree2_result)              :: searchResults(maxNumSearch)
    real(kdkind)                      :: maxRadius = 0.000001d0
    real(kdkind)                      :: refPosition(2)
    real(kdkind), allocatable         :: PPPandVNM(:,:)
    integer, allocatable              :: bodyIndexList(:)

    INTEGER,  PARAMETER    :: NBLOC_LIST = 9
    INTEGER                :: LNMX

    TYPE(BURP_FILE)        :: FILE_IN
    TYPE(BURP_RPT)         :: RPT_IN,CP_RPT
    TYPE(BURP_BLOCK)       :: BLOCK_IN,BLOCK_OMA,BLOCK_OMP,BLOCK_OER,BLOCK_FGE,BLOCK_FLG,BLOCK_FSO
    TYPE(BURP_BLOCK)       :: BLOCK_OMA_SFC,BLOCK_OMP_SFC,BLOCK_OER_SFC,BLOCK_FGE_SFC,BLOCK_FLG_SFC,BLOCK_FSO_SFC
    TYPE(BURP_BLOCK)       :: Block_FLG_CP,BLOCK_OBS_MUL_CP,BLOCK_MAR_MUL_CP,BLOCK_OBS_SFC_CP,BLOCK_MAR_SFC_CP
    TYPE(BURP_BLOCK)       :: BLOCK_GEN, BLOCK_OBS_BND,BLOCK_MAR_BND,BLOCK_ORB

    CHARACTER(LEN=5)       :: FAMILYTYPE2
    CHARACTER(LEN=9)       :: OPT_MISSING
    integer                :: BTYP,BFAM,BTYP10,BTYP10FLG_uni,BTYP10obs_uni 
    integer                :: BTYP10DES,BTYP10INF,BTYP10OBS,BTYP10FLG

    integer                :: NB_RPTS,REF_RPT,REF_BLK,COUNT
    INTEGER, ALLOCATABLE   :: ADDRESS(:)
    real                   :: VCOORD

    integer                :: NBELE,NVALE,NTE, numHeader
    integer                :: J,JJ,K,KK,KI,IL,Jo,ERROR,OBSN,KOBSN,ITEM
    integer                :: IND_ELE,IND_VCOORD
    integer                :: IND_ELE_MAR,IND_ELEU,IND_ELEF,IND_ELE_stat,IND_ELE_tth,IND_ELE_esh
    integer                :: IND_LAT,IND_LON,IND_TIME,IND_obsClear

    integer                :: vcord_type(10),SUM
    real                   :: ELEVFACT
    integer                :: status ,idtyp,lati,long,dx,dy,elev, &
                              drnd,date_h,hhmm_h,oars,runn
    integer                :: IND055200

    integer                :: iele,NELE,NELE_SFC,NVAL,NT,NELE_INFO
    integer                :: bit_alt,btyp_offset,btyp_offset_uni
    integer                :: BKNAT,BKTYP,BKSTP
    character(len = 5)     :: BURP_TYP
    CHARACTER(LEN=9)       :: STNID,STN_RESUME,STID
    LOGICAL                :: HIRES,HIPCS
    integer                :: NDATA_SF
    integer                :: IFLAG,BITSflagoff

    integer                :: OBS_START,SAVE_OBS
    integer                :: IL_INDEX,IRLN,INLV,LK,VNM
    real                   :: OBS,OMA,OMP,OER,FSO,FGE,OBSVA,CONVFACT, BCOR, obsClear
    integer                :: FLG,TIME,ILEMU,ILEMV,ILEMD,VCOORD_POS,ILEMZBCOR,ILEMTBCOR,ILEMHBCOR

    integer                :: BLOCK_LIST(NBLOC_LIST),bl

    integer                :: new_bktyp,post_bit,STATUS_HIRES,BIT_STATUS,FILEN
    LOGICAL                :: REGRUP,WINDS,OMA_SFC_EXIST,OMA_ALT_EXIST

    integer                :: LISTE_ELE(20),LISTE_ELE_SFC(20),is_in_list
    integer                :: ADDSIZE
    
    LOGICAL                :: LBLOCK_OER_CP, LBLOCK_FGE_CP
    TYPE(BURP_BLOCK)       :: BLOCK_OER_CP, BLOCK_FGE_CP
    logical                :: FSOFound
    logical                :: btClearElementFound 

    ! ensure kdtrees object is null
    nullify(tree)

    STATUS_HIRES = 0
    FAMILYTYPE2= 'SCRAP'
    vcord_type(:)=-1
    vcord_type(1)=0
    NELE_INFO=1
    NELE_SFC=0
    NELE=0
    ILEMU=11003
    ILEMV=11004
    ILEMD=11001
    ILEMZBCOR=15234 ! bcor element for GP  ZTD observations
    ILEMTBCOR=12204 ! bcor element for altitude TT observations
    ILEMHBCOR=12243 ! bcor element for altitude ES observations
    ELEVFACT=0.
    BNBITSOFF=0
    BNBITSON=0
    ENFORCE_CLASSIC_SONDES=.false.
    UA_HIGH_PRECISION_TT_ES=.false.
    UA_FLAG_HIGH_PRECISION_TT_ES=.false.
    ADDSIZE=100000
    LNMX=100000
    LISTE_ELE_SFC(:)=-1
    LISTE_ELE(:)=-1
    SELECT CASE(trim(FAMILYTYPE))
      CASE('UA')
        BURP_TYP='multi'
        vcord_type(1)=7004

        LISTE_ELE_SFC(1:6) = (/12004,11011,11012,10051,10004,12203/)
        NELE_SFC=6
        call BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'UA'
        LISTE_ELE(1:5) = (/12001,11001,11002,12192,10194/)
        NELE=5
        ENFORCE_CLASSIC_SONDES=.false.
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=10000
      CASE('AI')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:4) = (/12001,12192,11001,11002/)
        NELE=4
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=5000
      CASE('AL')
        BURP_TYP='uni'
        vcord_type(1)=7071

        LISTE_ELE(1:4) = (/10004,40030,12001,5021/)
        NELE=4
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=5000
      CASE('SW')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        vcord_type(2)=-1
!       ADDSIZE=5000
!       ' AMVS-> ADDSIZE 600000 oct 2014 pik'
        ADDSIZE=600000
      CASE('SF')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELEMS_SFC=6
        LISTE_ELE_SFC(1:NELEMS_SFC) = (/10004,12004,10051,12203,11011,11012/)
        BLISTELEMENTS_SFC(1:NELEMS_SFC) = LISTE_ELE_SFC(1:NELEMS_SFC)

        call BRPACMA_NML('namburp_sfc') ! read NELEMS_SFC, BLISTELEMENTS_SFC(1:NELEMS_SFC)
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'SFC'
        WINDS=.TRUE.
        IF (trim(FAMILYTYPE) == 'GP') WINDS=.FALSE.
        ILEMU=11215
        ILEMV=11216
        ILEMD=11011
        ADDSIZE=5000
      CASE('GP')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELEMS_GPS=6
        LISTE_ELE_GPS(1:NELEMS_GPS) = (/10004,12004,12203,15031,15032,15035/)

        call BRPACMA_NML('namburp_sfc') ! read NELEMS_GPS, LISTE_ELE_GPS(1:NELEMS_GPS)
        NELE_SFC=NELEMS_GPS             !   -- ignore NELEMS_SFC, BLISTELEMENTS_SFC(1:NELEMS_SFC)
        BLISTELEMENTS_SFC(1:NELEMS_GPS) = LISTE_ELE_GPS(1:NELEMS_GPS)

        FAMILYTYPE2= 'SFC'
        WINDS=.FALSE.
        ADDSIZE=5000      
      CASE('SC')
        BURP_TYP='uni'
        LISTE_ELE_SFC(1:2) = (/11012,11011/)
        NELE=2
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS

        FAMILYTYPE2='SCAT'
        WINDS=.TRUE.
        ILEMU=11215
        ILEMV=11216
        ILEMD=11011
        ADDSIZE=5000
      CASE('PR')
        BURP_TYP='multi'
        vcord_type(1)=7006
        ELEVFACT=1.

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2

        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=10000
      CASE('RO')
        BURP_TYP='multi'
        vcord_type(1)=7007
        vcord_type(2)=7040

        LISTE_ELE(1:2) = (/15036,15037/)
        NELE=2

        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.FALSE.
        !================GPS-RO CANNOT BE FILTERED=======
        BNBITSOFF=0
        BNBITSON=0
        !================GPS-RO CANNOT BE FILTERED=======
        NELE_INFO=16
      CASE('TO')
        BURP_TYP='multi'
        vcord_type(1)=5042
        vcord_type(2)=2150
        
        LISTE_ELE(1:1) = (/12163/)
        NELE=1

        CALL BRPACMA_NML('namburp_tovs')
        NELE=NELEMS

        if ( addBtClearToBurp ) then
          btClearElementFound = .false.

          elementLoop: do iele = 1, NELE
            if ( BLISTELEMENTS(iele) == btClearElementId ) then
              btClearElementFound = .true.
              exit elementLoop
            end if
          end do elementLoop

          if ( .not. btClearElementFound ) &
            call utl_abort('brpr_updateBurp: btClearElement element should be in namelist.')
        end if

        NELE_INFO=23
        WINDS=.FALSE.
        ADDSIZE=600000
      CASE('CH')
        BURP_TYP='multi'   ! Both 'multi' and 'uni' are possible for this family.
                           ! 'uni' level data are assumed not to have any accompanynig vertical 
                           ! coordinate element in addition to having only one level.
        vcord_type(1:8) = (/7004,7204,7006,7007,5042,2150,2071,0/)  ! 0 must be at end.

        LISTE_ELE_SFC(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                                15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE_SFC=19
        call BRPACMA_NML('namburp_chm_sfc')
        NELE_SFC=NELEMS_SFC
        WINDS=.FALSE.
        ADDSIZE=10000

        FAMILYTYPE2='CH'
        LISTE_ELE(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                                15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE=19
        call BRPACMA_NML('namburp_chm')
        NELE=NELEMS
      CASE('GO','MI')
        call utl_abort('brpr_updateBurp: unknown familyType : ' // trim(familyType))
    END SELECT
    LISTE_ELE    (1:NELE    )=BLISTELEMENTS(1:NELE)
    LISTE_ELE_SFC(1:NELE_SFC)=BLISTELEMENTS_SFC(1:NELE_SFC)
    if(NELE     > 0)write(*,*)  ' LISTE_ELE =',LISTE_ELE
    if(NELE_SFC > 0)write(*,*)  ' LISTE_ELE_SFC =',LISTE_ELE_SFC
    if(BNBITSOFF > 0 .or. BNBITSON > 0)write(*,*)  ' BNBITSON BNBITSOFF SIZE OF CP_RPT   =',BNBITSON,BNBITSOFF,LNMX*8


    TYPE_RESUME='POSTALT'
    BN_ITEMS=1
    BITEMLIST(1)='OMA'
    call BRPACMA_NML('namburp_update')
    write(*,*) ' BN_ITEMS   =',BN_ITEMS
    write(*,*) ' ITEMS TO ADD IN BURP FILE REPORTS =', BITEMLIST(1:BN_ITEMS)
    write(*,*) ' BTYP OF UPDATED BURP FILE=', TYPE_RESUME

    ! check if there is FSO calculation
    FSOFound = .false.
    do item = 1, BN_ITEMS
      if ( BITEMLIST(item) == 'FSO' ) FSOFound = .true.
    end do

    SELECT CASE( trim(TYPE_RESUME))
      CASE("BGCKALT", "POSTALT")
        BIT_STATUS = 12 
      CASE("DERIALT")
        BIT_STATUS = 11 
    END SELECT 

    if (trim(BURP_TYP) == 'uni') then
      btyp_offset=256
    else
      btyp_offset=0
    end if

    if (trim(familytype) == 'AL') then
      btyp_offset=255
    end if

    if ( TRIM(FAMILYTYPE2) == 'SCAT') then
      btyp_offset= 0
      btyp_offset_uni= 256 +0
    elseif ( TRIM(FAMILYTYPE2) == 'SFC') then
      btyp_offset= 0
      btyp_offset_uni= 256 +32
    elseif ( TRIM(FAMILYTYPE2) == 'UA') then
      btyp_offset_uni= 256 +32
    elseif ( TRIM(FAMILYTYPE2) == 'CH') then
      btyp_offset_uni= 256
    else
      btyp_offset_uni= -999 !set to -999 when not used
    end if


    write(*,*) '----------------------------------------------------'
    write(*,*) '-----------     BEGIN brpr_updateBurp   ------------'
    write(*,*) 'FAMILYTYPE   =',FAMILYTYPE
    write(*,*) 'BURP_TYP btyp_offset    =',BURP_TYP, btyp_offset
    write(*,*) 'BURP_TYP btyp_offset_uni=',BURP_TYP, btyp_offset_uni
    write(*,*) '----------------------------------------------------'


    ! initialisation

    SUM=0
    opt_missing = 'MISSING'


    call BURP_Set_Options( &
       & REAL_OPTNAME       = opt_missing, &
       & REAL_OPTNAME_VALUE = MPC_missingValue_R4, &
       & CHAR_OPTNAME       = 'MSGLVL', &
       & CHAR_OPTNAME_VALUE = 'FATAL', &
       & IOSTAT          = error )
    call handle_error(error, "brpr_updateBurp: BURP_Set_Options")

    call BURP_Init(File_in      ,IOSTAT=error)
    call handle_error(error, "brpr_updateBurp: BURP_Init FILE_in")

    call BURP_Init(Rpt_in,CP_RPT)
    call BURP_Init(Block_in)
    call BURP_Init(BLOCK_OMA)
    call BURP_Init(BLOCK_OMP)
    call BURP_Init(BLOCK_OER)
    call BURP_Init(BLOCK_FGE)
    call BURP_Init(BLOCK_FSO)
    call BURP_Init(BLOCK_OMA_SFC)
    call BURP_Init(BLOCK_OMP_SFC)
    call BURP_Init(BLOCK_OER_SFC)
    call BURP_Init(BLOCK_FGE_SFC)
    call BURP_Init(BLOCK_FSO_SFC)
    call BURP_Init(BLOCK_FLG_SFC)
    
    call BURP_Init(BLOCK_FLG)
    call BURP_Init(Block_FLG_CP)

    call BURP_Init(BLOCK_OBS_MUL_CP)
    call BURP_Init(BLOCK_MAR_MUL_CP)

    call BURP_Init(BLOCK_OBS_SFC_CP)
    call BURP_Init(BLOCK_MAR_SFC_CP)

    Call BURP_Init(BLOCK_GEN)
    Call BURP_Init(BLOCK_OBS_BND)
    Call BURP_Init(BLOCK_MAR_BND)
    Call BURP_Init(BLOCK_ORB)

    ! opening file
    ! ------------
    write(*,*) 'OPENING BURP FILE FOR UPDATE = ', trim(brp_file)

    call BURP_New(File_in, FILENAME = brp_file, &
                   & MODE    = FILE_ACC_APPEND, &
                   & IOSTAT  = error )
    call handle_error(error, "brpr_updateBurp: BURP_new error while opening burp file " // trim(brp_file))

    ! obtain input burp file number of reports
    ! ----------------------------------------
    call BURP_Get_Property(File_in, NRPTS=nb_rpts)
    call BURP_Init(Rpt_in)

    write(*,*) '-----------------------------------------'
    write(*,*) 'IOSTAT    =',error
    write(*,*) 'NUMBER OF REPORTS IN FILE  = ',nb_rpts
    write(*,*) '-----------------------------------------'

    ! scan input burp file to get all reports address
    ! -----------------------------------------------

    Allocate(address(nb_rpts))
    address(:) = 0
    count      = 0
    ref_rpt    = 0
    bit_alt    = 0
    stn_resume='NOT_FOUND'

    do
      ref_rpt = BURP_Find_Report(File_in, &
                 & REPORT      = Rpt_in, &
                 & SEARCH_FROM = ref_rpt, &
                 & IOSTAT      = error)
      call handle_error(error, "brpr_updateBurp: BURP_Find_Report #1")
      call burp_get_property(Rpt_in, STNID = stnid )
      IF ( stnid(1:2) == ">>" ) then
        STN_RESUME=stnid
        SELECT CASE(stnid)
          CASE(">>BGCKALT", ">>POSTALT")
            bit_alt=1
          CASE(">>DERIALT")
            bit_alt=2
        END SELECT 
      END IF

      if (ref_rpt < 0) Exit
      if (count == nb_rpts) then
        write(*,*) 'brpr_updateBurp: ERROR: count = nb_rpts:',count,nb_rpts
        exit
      end if
      count = count + 1
      address(count) = ref_rpt
    end do

    if (stn_resume == 'NOT_FOUND') then
      write(*,*) 'brpr_updateBurp: WARNING: No RESUME record found in this file, ' //  &
                 'check if found during reading of all files'
      ! try to get value from previously read file
      if ( type_resume /= 'UNKNOWN' ) then
        stn_resume = '>>' // type_resume
        SELECT CASE(stn_resume)
          CASE(">>BGCKALT", ">>POSTALT")
            bit_alt=1
          CASE(">>DERIALT")
            bit_alt=2
          CASE DEFAULT
            write(*,*) 'brpr_updateBurp: WARNING: Unknown RESUME record found, assume BGCKALT'
            stn_resume = '>>BGCKALT'
            bit_alt=1
        END SELECT 
      else
        write(*,*) 'brpr_updateBurp: WARNING: No file read has RESUME record, assume BGCKALT'
        stn_resume = '>>BGCKALT'
        bit_alt=1
      end if
    end if

    write(*,'(a9,1x,a16,1x,i2)' )STN_RESUME,' bit_alt==== >  ',bit_alt

    BTYP10obs     = 291 -btyp_offset
    BTYP10obs_uni = 291 -btyp_offset_uni
    if (bit_alt == 2) btyp10obs     =  BTYP10obs - 2
    if (bit_alt == 2) btyp10obs_uni =  BTYP10obs_uni - 2

    BTYP10flg = 483 -btyp_offset
    BTYP10flg_uni = 483 -btyp_offset_uni
    if (bit_alt == 2) BTYP10flg =  BTYP10flg  - 2
    if (bit_alt == 2) BTYP10flg_uni =  BTYP10flg_uni  - 2
    BTYP10des = 160

    BTYP10inf = 96

    write(*, *)  ' NUMBER OF VALID REPORTS IN FILE = ',count
    write(*, *)  ' BTYP10obs BTYP10obs_uni         = ',BTYP10obs,BTYP10obs_uni
    BITSflagoff=0
    DO J = 1, Bnbitsoff
      BITSflagoff = IBSET ( BITSflagoff, 13-BBITOFF(J) )
    END DO

    if ( count > 0 ) then

      OBS_START=1
      SAVE_OBS=1
      DO Jo=1,obs_numHeader(obsdat)
        filen= obs_headElem_i(obsdat,OBS_OTP,Jo)
        if ( filen == filenumb) then
          OBS_START=Jo
          SAVE_OBS=Jo
          exit
        end if
      END DO
      write(*, *)  '  FILE = ',trim(brp_file),' OBS_START= ',OBS_START

      ! Create a new report.
      ! The factor 12 before 'LNMX' is arbitrary.
      ! We increase it from time to time as we encounter some
      ! problems.
      call BURP_New(Cp_rpt, ALLOC_SPACE=12*LNMX, IOSTAT=error)
      call handle_error(error, "brpr_updateBurp: error while allocating Cp_rpt")

      ! LOOP ON REPORTS
      REPORTS: do kk = 1, count

        call BURP_Get_Report(File_in, &
                & REPORT    = Rpt_in, &
                & REF       = address(kk), &
                & IOSTAT    = error)
        call handle_error(error, "brpr_updateBurp: BURP_Get_Report #1") 
        call burp_get_property(Rpt_in, &
               STNID = stnid ,TEMPS =hhmm_h,FLGS = status ,IDTYP =idtyp,LATI = lati &
               ,LONG = long ,DX = dx ,DY = dy,ELEV=elev,DRND =drnd,DATE =date_h &
               ,OARS =oars,RUNN=runn ,IOSTAT=error)
        call handle_error(error, "brpr_updateBurp: burp_get_property #1")

        IF ( stnid(1:2) == ">>" ) THEN
          write(*,*)  ' RESUME RECORD POSITION IN BURP FILE =',stnid,kk
          call BURP_Copy_Header(TO=Cp_rpt,FROM=Rpt_in)
          call BURP_Init_Report_Write(File_in,Cp_Rpt, IOSTAT=error)
          call handle_error(error, "brpr_updateBurp: Burp_Init_Report_Write #1")
          call BURP_Set_Property(Cp_Rpt,STNID =">>"//TYPE_RESUME)
          call BURP_Delete_Report(File_in,Rpt_in, IOSTAT=error)
          call handle_error(error, "brpr_updateBurp: BURP_Delete_Report #1")
          call BURP_Write_Report(File_in,Cp_rpt, IOSTAT=error)
          call handle_error(error, "brpr_updateBurp: BURP_Write_report #1")
          cycle REPORTS
        ELSE
          !write(*,*) ' UPDATING STN IN BURP FILE =',  TRIM(FAMILYTYPE),KK,stnid,lati,LONG,dx,DY,elev,idtyp
        END IF
        call BURP_Copy_Header(TO=Cp_rpt,FROM=Rpt_in)
        call BURP_Init_Report_Write(File_in,Cp_Rpt, IOSTAT=error)
        call handle_error(error, "brpr_updateBurp: Burp_Init_Report_Write #2")

        ! FIRST LOOP ON BLOCKS

        ref_blk = 0

        LBLOCK_OER_CP=.false.
        LBLOCK_FGE_CP=.false.
        HIRES=.FALSE.
        HIPCS=.FALSE.
        REGRUP=.false.
        NDATA_SF=-1
        !WRITE(*,*)'  record number =',kk,' obs_start =',obs_start
        BLOCK_LIST(:)=-1
        BLOCKS0: do
          ref_blk = BURP_Find_Block(Rpt_in, &
                     & BLOCK       = Block_in, &
                     & SEARCH_FROM = ref_blk, &
                     & IOSTAT      = error)
          call handle_error(error, "brpr_updateBurp: BURP_Find_Block #1")

          if (ref_blk < 0) EXIT BLOCKS0

          call BURP_Get_Property(Block_in, &
                                 & NELE   = nbele, &
                                 & NVAL   = nvale, &
                                 & NT     = nte, &
                                 & BFAM   = bfam, &
                                 & BTYP   = btyp, &
                                 & BKTYP  = bktyp, &
                                 & BKNAT  = BKNAT, &
                                 & BKSTP  = BKSTP, &
                                 & IOSTAT = error)
          call handle_error(error, "brpr_updateBurp: BURP_Get_Property #2")
          if(trim(familytype) == 'AL')then

            ! Fudge the block type, because the data are simulated
            if(btyp == 1024) then
              btyp=1152
            else if(btyp == 7168)then
              btyp=7296
            end if
          end if

          btyp10    = ishft(btyp,-5)
          if ( btyp10 ==  BTYP10des ) then
            Block_FLG_CP=BLOCK_IN
            BLOCK_LIST(1)=BTYP
            REGRUP=.TRUE.
          elseif ( btyp10 == btyp10obs_uni .and. bkstp <= 4 ) then
            BLOCK_OBS_SFC_CP=BLOCK_IN
            BLOCK_LIST(2)=BTYP
            NDATA_SF=0
          elseif ( btyp10 == btyp10flg_uni .and. bkstp <= 4) then
            BLOCK_MAR_SFC_CP=BLOCK_IN
            BLOCK_LIST(3)=BTYP
          elseif ( btyp10 == btyp10obs .and. bfam == 0 ) then
            BLOCK_LIST(4)=BTYP
            BLOCK_OBS_MUL_CP=BLOCK_IN
          elseif ( btyp10 == btyp10flg ) then
            BLOCK_LIST(5)=BTYP
            BLOCK_MAR_MUL_CP=BLOCK_IN
          elseif ( (btyp10 == btyp10inf) .or. (btyp10 - btyp10inf == 1) ) then
            BLOCK_LIST(6) = BTYP
            BLOCK_GEN     = BLOCK_IN
          else if (trim(familytype) == 'RO' .and. bfam == 0 .and. btyp ==  9217) then
            BLOCK_LIST(7) = BTYP
            BLOCK_OBS_BND = BLOCK_IN
          else if (trim(familytype) == 'RO' .and. bfam == 0 .and. btyp == 15361) then
            BLOCK_LIST(8) = BTYP
            BLOCK_MAR_BND = BLOCK_IN
          else if (trim(familytype) == 'RO' .and. bfam == 0 .and. btyp ==  9220) then
            BLOCK_LIST(9) = BTYP
            BLOCK_ORB     = BLOCK_IN
          else
            !WRITE(*, *)' POUR STATION bloc NON CONNU: ',STNID,ref_blk,bfam,familytype
          end if
          if (bfam == 10.and.bkstp == 14) then
              LBLOCK_OER_CP=.true.
              BLOCK_OER_CP=BLOCK_IN
          else  if (bfam == 10.and.bkstp == 15) then
              LBLOCK_FGE_CP=.true.
              BLOCK_FGE_CP=BLOCK_IN
          end if
        end do BLOCKS0
        if ( TYPE_RESUME == 'POSTALT' .or. TYPE_RESUME == 'BGCKALT') THEN
          post_bit=2
        else
          post_bit=0
        end if
        
        BLOCKS1: do bl=1,NBLOC_LIST

          if( BLOCK_LIST(bl) < 0 ) cycle
          ref_blk = BURP_Find_Block(Rpt_in, &
                               & BLOCK       = Block_in, &
                               & BTYP = BLOCK_LIST(bl), &
                               & IOSTAT      = error)
          call handle_error(error, "brpr_updateBurp: BURP_Find_Block #2")

          if (ref_blk < 0) cycle BLOCKS1

          call BURP_Get_Property(Block_in, &
                           & NELE   = nbele, &
                           & NVAL   = nvale, &
                           & NT     = nte, &
                           & BFAM   = bfam, &
                           & BTYP   = btyp, &
                           & BKTYP  = bktyp, &
                           & BKNAT  = BKNAT, &
                           & IOSTAT = error)
          call handle_error(error, "brpr_updateBurp: BURP_Get_Property #3")
          if(trim(familytype) == 'AL')then

            ! Fudge the block type, because the data are simulated
            if(btyp == 1024) then
              btyp=1152
            else if(btyp == 7168)then
              btyp=7296
            end if
          end if

          ! observation block (btyp = 0100 100011X XXXX)
          !======================================================
          btyp10    = ishft(btyp,-5)
          !======================================================

          OBS_START=SAVE_OBS
          !if ( btyp10 - btyp10obs_uni == 0 .and. bfam == 0 ) then
          if ( bl == 2 ) then
            OBSN=OBS_START

            NDATA_SF=0
            new_bktyp=bktyp
            if ( post_bit > 0 ) then
              new_bktyp=IBSET(bktyp,post_bit)

              if ( FSOFound ) then
                ! to correct the btyp of SC and UA surface observation for the block
                ! of observation value and flag
                call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BKTYP =new_bktyp, BKSTP=0)
                call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BKTYP =new_bktyp, BKSTP=0)
              else
                call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BKTYP =new_bktyp)
                call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BKTYP =new_bktyp)
              end if

            end if

            il_index=0
             
            !call BURP_Delete_BLOCK(Rpt_in,BLOCK=Block_in)

            OMA_SFC_EXIST=.true.

            if (.not.WINDS) then

               call BURP_New(BLOCK_OMA_SFC,NELE =NBELE,NVAL=nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMA #1") 
               call BURP_New(BLOCK_OMP_SFC,NELE =NBELE,NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMP_SFC #1") 
               call BURP_New(BLOCK_OER_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OER_SFC #1") 
               call BURP_New(BLOCK_FGE_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                    ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FGE_SFC #1") 
               call BURP_New(BLOCK_FSO_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FSO_SFC #1")
            else

               IND_eleu  = BURP_Find_Element(Block_in, ELEMENT=11215)
               IND_elef  = BURP_Find_Element(Block_in, ELEMENT=11011)

               call BURP_New(BLOCK_OMA_SFC,NELE =NBELE+2,NVAL=nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMA_SFC #2") 
               call BURP_New(BLOCK_OMP_SFC,NELE =NBELE+2,NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMP_SFC #2") 
               call BURP_New(BLOCK_OER_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OER_SFC #2") 
               call BURP_New(BLOCK_FGE_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FGE_SFC #2") 
               call BURP_New(BLOCK_FSO_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2 &
                             ,IOSTAT = error)
               call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FSO_SFC #2")
               ILEMU = 11215
               ILEMV = 11216
               call BURP_Set_Element( BLOCK_OMA_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA_SFC ILEMU") 
               call BURP_Set_Element( BLOCK_OMA_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA_SFC ILEMV") 

               call BURP_Set_Element( BLOCK_OMP_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP_SFC ILEMU") 
               call BURP_Set_Element( BLOCK_OMP_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP_SFC ILEMV") 
   
               call BURP_Set_Element( BLOCK_OER_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER_SFC ILEMU") 
               call BURP_Set_Element( BLOCK_OER_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER_SFC ILEMV") 

               call BURP_Set_Element( BLOCK_FGE_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE_SFC ILEMU") 
               call BURP_Set_Element( BLOCK_FGE_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE_SFC ILEMV") 

               call BURP_Set_Element( BLOCK_FSO_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO_SFC ILEMU")
               call BURP_Set_Element( BLOCK_FSO_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)
               call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO_SFC ILEMV")

               do k =1,nte
                 call BURP_Set_Rval(Block_OMA_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_OMA_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_OMP_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_OMP_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_OER_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_OER_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_FGE_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_FGE_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 call BURP_Set_Rval(Block_FSO_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 )
                 call BURP_Set_Rval(Block_FSO_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 )
               end do
               il_index=2

               IF (IND_eleu < 0 .and. IND_elef > 0) THEN
                 call BURP_Resize_Block(BLOCK_OBS_SFC_CP,ADD_NELE = 2 ,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp: BURP_Resize_Block #1") 
                 call BURP_Set_Element( BLOCK_OBS_SFC_CP,NELE_IND = nbele+1,ElEMENT=ILEMU,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OBS_SFC_CP ILEMU") 
                 call BURP_Set_Element( BLOCK_OBS_SFC_CP,NELE_IND = nbele+2,ElEMENT=ILEMV,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OBS_SFC_CP ILEMV") 
                 do k =1,nte
                   call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =nbele+1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                   call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =nbele+2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 end do

                 call BURP_Resize_Block(BLOCK_MAR_SFC_CP,ADD_NELE = 2 ,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp:  BURP_Resize_Block #2") 
                 call BURP_Set_Element( BLOCK_MAR_SFC_CP,NELE_IND = nbele+1,ElEMENT=ILEMU+200000,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_MAR_SFC_CP ILEMU") 
                 call BURP_Set_Element( BLOCK_MAR_SFC_CP,NELE_IND = nbele+2,ElEMENT=ILEMV+200000,IOSTAT=error)
                 call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_MAR_SFC_CP ILEMV") 
                 do k =1,nte
                   call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =nbele+1,NVAL_IND =1,NT_IND = k,tblval = 0 ) 
                   call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =nbele+2,NVAL_IND =1,NT_IND = k,tblval = 0 ) 
                 end do
                 il_index=2
                 nbele=nbele +2
               end if

             end if

!=========================================
             REGRUP_SFC: do k=1,nte
!=========================================
               KOBSN=0

!-------------------------------------------
               elems_sfc: do IL = 1, NBELE  
!-------------------------------------------

                 iele=-1
                 iele=BURP_Get_Element(BLOCK_OBS_SFC_CP,INDEX =il,IOSTAT= error)
                 call handle_error(error, "brpr_updateBurp: BURP_GET_ELEMENT #1") 

                 IND_ELE_MAR= BURP_Find_Element(Block_MAR_SFC_CP, ELEMENT=iele+200000)
                 if (IND_ele_mar <= 0 ) cycle

                 IND_ele  = BURP_Find_Element(BLOCK_OBS_SFC_CP, ELEMENT=iele)
                 if (ind_ele==-1) call handle_error(IND_ELE, "brpr_updateBurp: element not found in BLOCK_OBS_SFC_CP")

                 if ( k == 1 ) then
                   if( OMA_SFC_EXIST ) then
                     if (iele /= ILEMU .and. iele /= ILEMV) then
                       il_index=il_index +1
                       call BURP_Set_Element (BLOCK_OMA_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                       call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA_SFC")
                       call BURP_Set_Element (BLOCK_OMP_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                       call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP_SFC")
                       call BURP_Set_Element (BLOCK_OER_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                       call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER_SFC")
                       call BURP_Set_Element (BLOCK_FGE_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                       call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE_SFC")
                       call BURP_Set_Element (BLOCK_FSO_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                       call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO_SFC")
                     end if
                   end if
                 end if

                 is_in_list=-1
                 is_in_list=FIND_INDEX(LISTE_ELE_SFC,iele)
                 if (is_in_list < 0   .and. iele  /= ILEMU .and. iele /= ILEMV) cycle ELEMS_SFC
                 IND_ele_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele)
                 if (IND_ele_stat==-1) call handle_error(IND_ele_stat, "brpr_updateBurp: element not found in BLOCK_OMA_SFC #1")
                 call BURP_Set_Rval(Block_OMA_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                 call BURP_Set_Rval(Block_OMP_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                 call BURP_Set_Rval(Block_OER_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                 call BURP_Set_Rval(Block_FGE_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                 call BURP_Set_Rval(Block_FSO_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4) 

                 IFLAG =  BURP_Get_Tblval(Block_MAR_SFC_CP,NELE_IND = IND_ele_mar,NVAL_IND = 1, NT_IND = k)
                 OBSVA =  BURP_Get_Rval  (Block_OBS_SFC_CP,NELE_IND = IND_ele    ,NVAL_IND = 1, NT_IND = k)
                 if (OBSVA == MPC_missingValue_R4 .and. iele /= ILEMU  .and. iele /= ILEMV ) cycle
                 if(iand(IFLAG,BITSflagoff) /= 0) cycle

                 if (OBSN > obs_numHeader(obsdat) ) write(*,*)  ' debordement  surface OBS_START=',OBS_START
                 if (OBSN > obs_numHeader(obsdat))  cycle

                 IRLN=obs_headElem_i(obsdat,OBS_RLN,OBSN )
                 INLV=obs_headElem_i(obsdat,OBS_NLV,OBSN )

                 IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele, IOSTAT=error)
                 if (IND_ELE_stat==-1) call handle_error(error, "brpr_updateBurp: element not found in BLOCK_OMA_SFC #2")
                 STID=obs_elem_c(obsdat,'STID',obs_start)
                 if ( STID /= stnid ) cycle

                 OBSDATA: do LK=IRLN,IRLN+INLV-1

                   VNM=obs_bodyElem_i(obsdat,OBS_VNM ,LK)
                   if( VNM == iele ) then 
                     OBS=obs_bodyElem_r(obsdat,OBS_VAR,LK)
                     OMA=obs_bodyElem_r(obsdat,OBS_OMA ,LK)
                     OMP=obs_bodyElem_r(obsdat,OBS_OMP ,LK)
                     OER=obs_bodyElem_r(obsdat,OBS_OER ,LK)
                     FGE=obs_bodyElem_r(obsdat,OBS_HPHT,LK)
                     if ( obs_columnActive_RB(obsdat,OBS_FSO) ) then
                       FSO=obs_bodyElem_r(obsdat,OBS_FSO,LK)
                     else
                       FSO = MPC_missingValue_R4
                     end if
                     if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) then
                       BCOR = obs_bodyElem_r(obsdat,OBS_BCOR,LK)
                     else
                       BCOR = MPC_missingValue_R4
                     end if
                     FLG=obs_bodyElem_i(obsdat,OBS_FLG ,LK)
                     KOBSN= KOBSN + 1
                     SUM=SUM +1
                     call BURP_Set_Rval( Block_OER_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OER   ) 
                     call BURP_Set_Rval( Block_FGE_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FGE  ) 
                     call BURP_Set_Rval( Block_FSO_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FSO  )

                     IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele)
                     if (IND_ELE_stat==-1) call handle_error(IND_ELE_stat, "brpr_updateBurp: element not found in BLOCK_OMA_SFC #3")
                     call BURP_Set_Rval( Block_OMA_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OMA)

                     IND_ELE_stat  = BURP_Find_Element(BLOCK_OMP_SFC, ELEMENT=iele)
                     if (IND_ELE_stat==-1)call handle_error(IND_ELE_stat, "brpr_updateBurp: element not found in BLOCK_OMP_SFC #1")
                     call BURP_Set_Rval( Block_OMP_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OMP)

                     IND_ELE_stat  = BURP_Find_Element(BLOCK_FSO_SFC, ELEMENT=iele)
                     if (IND_ELE_stat==-1)call handle_error(IND_ELE_stat, "brpr_updateBurp: element not found in BLOCK_FSO_SFC #1")
                     call BURP_Set_Rval( Block_FSO_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FSO)
   
                     call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =IND_ele_mar,NVAL_IND =1,NT_IND = k ,TBLVAL= FLG) 
   
                     !OBS=obs_bodyElem_r(obsdat,OBS_VAR,LK)
                     IND_ele  = BURP_Find_Element(BLOCK_OBS_SFC_CP, ELEMENT=iele)
                     if (IND_ELE==-1) call handle_error(IND_ele, "brpr_updateBurp: element not found in BLOCK_OBS_SFC #1")
                     call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =IND_ele,NVAL_IND =1,NT_IND = k,RVAL = OBS ) 

                     if (iele == BUFR_NEZD) then
                       IND_ele  = BURP_Find_Element(BLOCK_OBS_SFC_CP, ELEMENT=ILEMZBCOR)
                       if (IND_ele > 0 .and. obs_columnActive_RB(obsdat,OBS_BCOR)) &
                            call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =IND_ele,NVAL_IND =1,NT_IND = k,RVAL = BCOR ) 
                     end if
                        
                     exit
                   end if

                 end do  OBSDATA
!-------------------------------------------
               end do  elems_sfc
!-------------------------------------------

               if ( REGRUP .and. KOBSN > 0   )  then
                 STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START )
                 STATUS=IBSET(STATUS,BIT_STATUS)
                 ind055200 = BURP_Find_Element(Block_FLG_CP, ELEMENT=055200)
                 if (ind055200==-1)call handle_error(ind055200, "brpr_updateBurp: element 055200 not found in Block_FLG_CP")
                 call BURP_Set_tblval( Block_FLG_CP, NELE_IND =ind055200,NVAL_IND =1,NT_IND = k ,TBLVAL = STATUS ) 
                 OBSN=OBSN +1
                 OBS_START=OBS_START +1
               end if
!==============================================
             end do  REGRUP_SFC
!=============================================

             do item=1,BN_ITEMS

               if ( BITEMLIST(item) == 'OMA') then
                 call BURP_Reduce_Block(BLOCK_OMA_SFC, NEW_NELE =il_index )
                 call BURP_Write_Block( CP_RPT, BLOCK_OMA_SFC,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                 call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OMA_SFC") 
                 cycle
               end if
               if ( BITEMLIST(item) == 'OMP') then
                 call BURP_Reduce_Block(BLOCK_OMP_SFC, NEW_NELE =il_index )
                 call BURP_Write_Block( CP_RPT, BLOCK_OMP_SFC,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                 call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OMP_SFC") 
                 cycle
               end if
              if ( BITEMLIST(item) == 'OER') then
                if (.not.LBLOCK_OER_CP) then
                  call BURP_Reduce_Block(BLOCK_OER_SFC, NEW_NELE =il_index )
                  call BURP_Write_Block( CP_RPT, BLOCK_OER_SFC,&
                       ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OER_SFC") 
                else
                  call BURP_Set_Property(BLOCK_OER_CP ,BKTYP =new_bktyp)
                  call BURP_Write_Block( CP_RPT, BLOCK_OER_CP,&
                       ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OER_CP")                 
                end if
                cycle
              end if
              if ( BITEMLIST(item) == 'FGE') then
                if (.not.LBLOCK_FGE_CP) then
                  call BURP_Reduce_Block(BLOCK_FGE_SFC, NEW_NELE =il_index )
                  call BURP_Write_Block( CP_RPT, BLOCK_FGE_SFC,&
                       ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FGE_SFC") 
                else
                  call BURP_Set_Property(BLOCK_FGE_CP ,BKTYP =new_bktyp)
                  call BURP_Write_Block( CP_RPT, BLOCK_FGE_CP,&
                       ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FGE_CP")                 
                end if
                cycle
              end if

              if ( BITEMLIST(item) == 'FSO') then
                call BURP_Reduce_Block(BLOCK_FSO_SFC, NEW_NELE =il_index )
                call BURP_Write_Block( CP_RPT, BLOCK_FSO_SFC,&
                     ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FSO_SFC")
                cycle
              end if

            end do

            call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BFAM =0)
            call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BFAM =0)

            call BURP_Write_Block( CP_RPT, BLOCK_OBS_SFC_CP,&
                                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OBS_SFC_CP") 
            call BURP_Write_Block( CP_RPT, BLOCK_MAR_SFC_CP,&
                                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_MAR_SFC_CP") 
            IF ( KOBSN > 0 .and. .not. REGRUP ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
              STATUS=IBSET(STATUS,BIT_STATUS)
              call BURP_Set_Property(CP_RPT ,FLGS =STATUS)
            END IF

            SAVE_OBS=OBS_START
            NDATA_SF=KOBSN
            if (BLOCK_LIST(4) == -1 .and.  KOBSN > 0 .and. .not. regrup ) THEN
              SAVE_OBS=SAVE_OBS+1
              OBS_START=OBS_START+1
            end if

          end if  ! bl == 2


          !if ( btyp10 - btyp10obs == 0 .and. bfam == 0 ) then
          if ( bl == 4 ) then
            ILEMU=11003
            ILEMV=11004
            new_bktyp=bktyp
            if ( post_bit  > 0 ) then
              new_bktyp=IBSET(bktyp,post_bit)
              call BURP_Set_Property(BLOCK_OBS_MUL_CP ,BKTYP =new_bktyp)
              call BURP_Set_Property(BLOCK_MAR_MUL_CP ,BKTYP =new_bktyp)
            end if
            OBSN=OBS_START
            NVAL=NVALE ;  NT=NTE

            il_index=1
            call BURP_New(BLOCK_OMA, NELE =1, NVAL =nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                  ,IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMA") 
            call BURP_New(BLOCK_OMP, NELE =1, NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                  ,IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OMP") 
            call BURP_New(BLOCK_OER, NELE =1, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                  ,IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_OER") 
            call BURP_New(BLOCK_FGE, NELE =1, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                  ,IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FGE") 
            call BURP_New(BLOCK_FSO, NELE =1, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2  &
                  ,IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_New BLOCK_FSO")

            VCOORD_POS=0
            k=0
            IND_VCOORD=-1
            do while (vcord_type(k+1) /= -1 .and. IND_VCOORD == -1)
               k=k+1 
               IND_VCOORD  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=vcord_type(k))
            end do
            IF ( IND_VCOORD > 0 ) then
              IF (trim(FAMILYTYPE) == trim('CH')) THEN
                 ELEVFACT=0.0
                 IF (vcord_type(k) == 7006) ELEVFACT=1.0
              END IF
              call BURP_Set_Element(BLOCK_OMA,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA") 
              call BURP_Set_Element(BLOCK_OMP,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP") 
              call BURP_Set_Element(BLOCK_OER,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER") 
              call BURP_Set_Element(BLOCK_FGE,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE") 
              call BURP_Set_Element(BLOCK_FSO,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO #1")
              VCOORD_POS=1
            ELSE IF (IND_VCOORD == -1) then
              !write(*,*) '  PAS DE COORDONNEE VERTICALE famille ',trim(FAMILYTYPE)
              il_index=0
            end if
            VCOORD = -999.

            IND_eleu  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=ILEMU)
            IND_elef  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=ILEMD)

            OMA_ALT_EXIST=.false.
            if(WINDS .and. IND_eleu < 0 .and. IND_elef > 0) then
              call BURP_Resize_Block(BLOCK_OMA,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OMA #1") 
              call BURP_Set_Element( BLOCK_OMA,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA ILEMU") 
              call BURP_Set_Element( BLOCK_OMA,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA ILEMV") 
              OMA_ALT_EXIST=.true.

              call BURP_Resize_Block(BLOCK_OMP,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OMP #1") 
              call BURP_Set_Element( BLOCK_OMP,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP ILEMU") 
              call BURP_Set_Element( BLOCK_OMP,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP ILEMV") 

              call BURP_Resize_Block(BLOCK_OER,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OER #1") 
              call BURP_Set_Element( BLOCK_OER,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER ILEMU") 
              call BURP_Set_Element( BLOCK_OER,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER ILEMV") 

              call BURP_Resize_Block(BLOCK_FGE,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_FGE #1") 
              call BURP_Set_Element( BLOCK_FGE,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE ILEMU") 
              call BURP_Set_Element( BLOCK_FGE,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE ILEMV") 

              call BURP_Resize_Block(BLOCK_FSO,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_FSO #1")
              call BURP_Set_Element( BLOCK_FSO,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO ILEMU")
              call BURP_Set_Element( BLOCK_FSO,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO ILEMV")

              call BURP_Resize_Block(BLOCK_OBS_MUL_CP,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OBS_MUL_CP") 
              call BURP_Set_Element( BLOCK_OBS_MUL_CP,NELE_IND = nbele+1,ElEMENT=ILEMU,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OBS_MUL_CP ILEMU") 
              call BURP_Set_Element( BLOCK_OBS_MUL_CP,NELE_IND = nbele+2,ElEMENT=ILEMV,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OBS_MUL_CP ILEMV") 

              call BURP_Resize_Block(BLOCK_MAR_MUL_CP,ADD_NELE = 2 ,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_MAR_MUL_CP") 
              call BURP_Set_Element( BLOCK_MAR_MUL_CP,NELE_IND = nbele+1,ElEMENT=ILEMU+200000,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_MAR_MUL_CP ILEMU") 
              call BURP_Set_Element( BLOCK_MAR_MUL_CP,NELE_IND = nbele+2,ElEMENT=ILEMV+200000,IOSTAT=error)
              call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_MAR_MUL_CP ILEMV") 

              do k=1,nte
                do jj=1,nvale
                  call BURP_Set_Rval( Block_OMA,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_OMA,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_OMP,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_OMP,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_OER,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_OER,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_FGE,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_FGE,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  call BURP_Set_Rval( Block_FSO,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  )
                  call BURP_Set_Rval( Block_FSO,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  )

                  call BURP_Set_Rval(  BLOCK_OBS_MUL_CP, NELE_IND =nbele+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                  call BURP_Set_Rval(  BLOCK_OBS_MUL_CP, NELE_IND =nbele+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                  call BURP_Set_tblval(BLOCK_MAR_MUL_CP, NELE_IND =nbele+1 ,NVAL_IND =jj , NT_IND  = k , TBLVAL = 0  ) 
                  call BURP_Set_tblval(BLOCK_MAR_MUL_CP, NELE_IND =nbele+2 ,NVAL_IND =jj , NT_IND  = k , TBLVAL = 0  ) 
                end do
              end do
              nbele=nbele+2

              il_index=il_index+2
            end if
            !call BURP_Delete_BLOCK(Rpt_in,BLOCK=Block_in)

            ! LAT LON TIME IN DATA BLOCK
            IND_LAT   = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=5001)
            IND_LON   = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=6001)
            IND_TIME  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=4015)
            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES=.true.
            if(ENFORCE_CLASSIC_SONDES) hires=.false.

            if( (TRIM(FAMILYTYPE2) == 'UA') .and. UA_HIGH_PRECISION_TT_ES ) HIPCS=.true.

            !print * , ' hires =true ? ndata_sf ',stnid,hires,NDATA_SF

            if ( HIRES .AND. NDATA_SF  > 0 ) OBS_START =OBS_START +1
            OBSN=OBS_START
            STATUS_HIRES=0

            regrup_LOOP: do k=1,nte
              KOBSN=0

              levels: do j=1,nvale

                if (OBSN > obs_numHeader(obsdat)) then 
                  write(*,*) ' debordement  altitude OBSN=',OBSN
                else
                  IRLN=obs_headElem_i(obsdat,OBS_RLN,OBSN)
                  INLV=obs_headElem_i(obsdat,OBS_NLV,OBSN)
                  if ((j == 1 .or. HIRES) .and. INLV > 0) then
                    if (allocated(PPPandVNM))     deallocate(PPPandVNM)
                    if (allocated(bodyIndexList)) deallocate(bodyIndexList)
                    allocate(PPPandVNM(2,INLV))
                    allocate(bodyIndexList(INLV))
                    bodyCount = 0
                    do LK = IRLN, IRLN+INLV-1
                      bodyCount = bodyCount + 1
                      PPPandVNM(1,bodyCount) = obs_bodyElem_r(obsdat,OBS_PPP,LK) - (ELEV-400.)*ELEVFACT
                      PPPandVNM(2,bodyCount) = real(obs_bodyElem_i(obsdat,OBS_VNM,LK),8)
                      bodyIndexList(bodyCount) = lk
                    end do
                    if (associated(tree)) call kdtree2_destroy(tree)
                      tree => kdtree2_create(PPPandVNM, sort=.true., rearrange=.true.)
                  end if 
                end if

                !pikpik
                if(HIRES) KOBSN=0
                !pikpik

                ! Loop over elements to add all to the new blocks and the initialize values 
                elems0: do IL = 1, NBELE

                  iele=BURP_Get_Element(BLOCK_OBS_MUL_CP,INDEX =il,IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Get_Element BLOCK_OBS_MUL_CP") 

                  IND_ELE_MAR= BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=iele+200000)
                  if (IND_ele_mar <  0 ) cycle

                  IND_ele       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele)
                  if (IND_ele == IND_LAT  .and. hires ) cycle
                  if (IND_ele == IND_LON  .and. hires ) cycle
                  if (IND_ele == IND_TIME .and. hires ) cycle
                  IND_ELE_STAT=-1
                  IND_ele_STAT  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele)

                  if(j == 1 .and. il /= ind_vcoord .and. IND_ELE_STAT < 1 ) then

                    il_index=il_index +1
                    call BURP_Resize_Block(BLOCK_OMA,ADD_NELE = 1 ,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OMA #2") 
                    call BURP_Set_Element (BLOCK_OMA,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMA")

                    call BURP_Resize_Block(BLOCK_OMP,ADD_NELE = 1 ,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OMP #2") 
                    call BURP_Set_Element (BLOCK_OMP,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OMP")

                    call BURP_Resize_Block(BLOCK_OER,ADD_NELE = 1 ,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_OER #2") 
                    call BURP_Set_Element (BLOCK_OER,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_OER")

                    call BURP_Resize_Block(BLOCK_FGE,ADD_NELE = 1 ,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_FGE #2") 
                    call BURP_Set_Element (BLOCK_FGE,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FGE")

                    call BURP_Resize_Block(BLOCK_FSO,ADD_NELE = 1 ,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Resize_Block BLOCK_FSO #2")
                    call BURP_Set_Element (BLOCK_FSO,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)
                    call handle_error(error, "brpr_updateBurp: BURP_Set_Element BLOCK_FSO #2")

                    do ki=1,nte
                      do jj=1,nvale
                        call BURP_Set_Rval( Block_OMA, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                        call BURP_Set_Rval( Block_OMP, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                        call BURP_Set_Rval( Block_OER, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                        call BURP_Set_Rval( Block_FGE, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                        call BURP_Set_Rval( Block_FSO, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 )
                      end do
                    end do

                  end if

                end do elems0

                ! Loop over elements to update the values
                elems: do IL = 1, NBELE  

                  iele=-1
                  iele=BURP_Get_Element(BLOCK_OBS_MUL_CP,INDEX =il,IOSTAT= error)
                  call handle_error(error, "brpr_updateBurp: BURP_Get_Element BLOCK_OBS_MUL_CP") 

                  IND_ELE_MAR= BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=iele+200000)
                  if (IND_ele_mar <  0 ) cycle

                  IND_ele       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele)
                  if (IND_ele == IND_LAT  .and. hires ) cycle
                  if (IND_ele == IND_LON  .and. hires ) cycle
                  if (IND_ele == IND_TIME .and. hires ) cycle
                  IND_ELE_STAT=-1
                  IND_ele_STAT  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele)
                  VCOORD =  BURP_Get_Rval(BLOCK_OBS_MUL_CP, &
                                      &   NELE_IND = IND_VCOORD, &
                                      &   NVAL_IND = j, &
                                      &   NT_IND   = k)
                  IF (il == IND_VCOORD) THEN
                    call BURP_Set_Rval( Block_OMA, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                    call BURP_Set_Rval( Block_OMP, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                    call BURP_Set_Rval( Block_OER, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                    call BURP_Set_Rval( Block_FGE, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                    call BURP_Set_Rval( Block_FSO, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  )
                  END IF

                  IND_ele       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele)
                  if (IND_ELE==-1)call handle_error(IND_ele, "brpr_updateBurp: element not found in BLOCK_OBS_MUL_CP")

                  is_in_list=-1
                  is_in_list=FIND_INDEX(LISTE_ELE,iele)
                  if (is_in_list < 0   .and. iele  /= ILEMU .and. iele /= ILEMV) cycle

                  IFLAG =  BURP_Get_Tblval(Block_MAR_MUL_CP,NELE_IND = IND_ELE_MAR,NVAL_IND = J, NT_IND = k)
                  OBSVA =  BURP_Get_Rval  (Block_OBS_MUL_CP,NELE_IND = IND_ele    ,NVAL_IND = J, NT_IND = k)
                  if(iand(iflag,BITSflagoff) /= 0) CYCLE ELEMS     

                  OBSVA =  BURP_Get_Rval  (Block_OBS_MUL_CP,NELE_IND = IND_ele    ,NVAL_IND = J, NT_IND = k)
                  !if( idtyp == 168 ) print * , ' bobossmi avant vcoord obsva    stnid =',  VCOORD,OBSVA,stnid

                  if (VCOORD == MPC_missingValue_R4 ) CYCLE ELEMS
                  if (OBSVA == MPC_missingValue_R4 .and. iele /= ILEMU  .and. iele /= ILEMV ) CYCLE ELEMS

                  if (OBSN > obs_numHeader(obsdat)) write(*,*) ' debordement  altitude OBSN=',OBSN
                  if (OBSN > obs_numHeader(obsdat))  cycle
                  TIME=obs_headElem_i(obsdat,OBS_ETM,OBSN)
                  STID=obs_elem_c(obsdat,'STID',OBSN)
                  if ( STID /= stnid ) cycle

                  convfact=1.
                  if (iele == 10194) convfact=1./ec_rg

                  if (INLV > 0) then
                    refPosition(1) = vcoord
                    refPosition(2) = iele
                    call kdtree2_r_nearest(tp=tree,  &
                                           qv=refPosition, r2=maxRadius, &
                                           nfound=numFoundSearch,        &
                                           nalloc=maxNumSearch,          &
                                           results=searchResults)
                    if (numFoundSearch == 0) cycle ELEMS
                    if (numFoundSearch > 1) then
                      write(*,*) 'vcoord, iele = ', vcoord, iele
                      do resultIndex = 1, numFoundSearch
                        write(*,*) 'ppp = ', PPPandVNM(1,searchResults(resultIndex)%idx)
                        write(*,*) 'vnm = ', PPPandVNM(2,searchResults(resultIndex)%idx)
                      end do
                      write(*,*) 'brpr_updateBurp: multiple obs matches found, taking closest'
                    end if
                    lk = bodyIndexList(searchResults(1)%idx)
                  else
                    cycle ELEMS
                  end if

                 
                  OBS=obs_bodyElem_r(obsdat,OBS_VAR,LK)*convfact
                  OMA=obs_bodyElem_r(obsdat,OBS_OMA,LK)
                  OMP=obs_bodyElem_r(obsdat,OBS_OMP,LK)
                  OER=obs_bodyElem_r(obsdat,OBS_OER,LK)
                  FGE=obs_bodyElem_r(obsdat,OBS_HPHT,LK)
                  if ( obs_columnActive_RB(obsdat,OBS_FSO) ) then
                    FSO=obs_bodyElem_r(obsdat,OBS_FSO,LK)
                  else
                    FSO = MPC_missingValue_R4
                  end if
                  if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) then
                    BCOR = obs_bodyElem_r(obsdat,OBS_BCOR,LK)
                  else
                    BCOR = MPC_missingValue_R4
                  end if
                  if ( obs_columnActive_RB(obsdat,OBS_BTCL) ) then
                    obsClear = obs_bodyElem_r(obsdat,OBS_BTCL,LK)
                  else
                    obsClear = MPC_missingValue_R4
                  end if
                  FLG=obs_bodyElem_i(obsdat,OBS_FLG,LK)
                  KOBSN= KOBSN + 1
                  IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele)
                  if (IND_ELE_stat==-1)call handle_error(IND_ELE_stat, "brpr_updateBurp: element not found in BLOCK_OMA")
                  if  ( OMA /= MPC_missingValue_R4 ) then
                    OMA=OMA*convfact
                  end if

                  call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )

                  if(HIPCS) then  
                    IND_ELE_tth   = BURP_Find_Element(BLOCK_OMA, ELEMENT=12101)
                    if (IND_ELE_tth==-1) call handle_error(IND_ELE_tth, "brpr_updateBurp: element 12101 not found in BLOCK_OMA")
                    IND_ELE_esh   = BURP_Find_Element(BLOCK_OMA, ELEMENT=12239)
                    if (IND_ELE_esh==-1) call handle_error(IND_ELE_esh, "brpr_updateBurp: element 12239 not found in BLOCK_OMA")
                    if(iele == 12001) then
                      if (IND_ele_tth==-1) write(*,*) 'PB1 SYLVAIN TTH' 
                      call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )
                    end if
                    if(iele == 12192) then
                      if (IND_ele_esh==-1) write(*,*) 'PB1 SYLVAIN ESH'
                      call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )
                    end if
                  endif

                  !if(trim(familytype) == 'TO' )print *,' bingo  stnid kk vnm ppp flg omp ',stnid,kk,vnm,ppp,flg,omp,oma
                  SUM=SUM +1
                  IND_ELE_stat  = BURP_Find_Element(BLOCK_OMP, ELEMENT=iele)
                  if (IND_ELE_stat==-1)call handle_error(IND_ELE_stat, "element not found in BLOCK_OMP")
                  if  ( OMP /= MPC_missingValue_R4 ) then
                    OMP=OMP*convfact
                  end if
                  call BURP_Set_Rval(  Block_OMP,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OMP)

                  if(HIPCS) then
                    if(iele == 12001) call BURP_Set_Rval(Block_OMP,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OMP )
                    if(iele == 12192) call BURP_Set_Rval(Block_OMP,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OMP )
                  endif

                  call BURP_Set_Rval(  Block_OER,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OER  ) 

                  if(HIPCS) then
                    if(iele == 12001) call BURP_Set_Rval(Block_OER,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OER )
                    if(iele == 12192) call BURP_Set_Rval(Block_OER,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OER )
                  endif

                  call BURP_Set_Rval(  Block_FGE,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = FGE ) 

                  if(HIPCS) then
                    if(iele == 12001) call BURP_Set_Rval(Block_FGE,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = FGE )
                    if(iele == 12192) call BURP_Set_Rval(Block_FGE,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = FGE )
                  endif

                  call BURP_Set_Rval(  Block_FSO,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )

                  if(HIPCS) then
                    if(iele == 12001) call BURP_Set_Rval(Block_FSO,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )
                    if(iele == 12192) call BURP_Set_Rval(Block_FSO,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )
                  endif

                  IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=iele+200000, IOSTAT=error)
                  if (IND_ELE_mar==-1)call handle_error(IND_ele_mar, "brpr_updateBurp: element not found in Block_MAR_MUL_CP")

                  call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )

                  if(HIPCS) then
                    if(iele == 12001) then
                      IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=212101)
                      if (IND_ELE_mar==-1)call handle_error(IND_ele_mar, "brpr_updateBurp: element 212001 not found in Block_MAR_MUL_CP")
                      call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )
                    endif
                    if(iele == 12192) then
                      IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=212239)
                      if (IND_ELE_mar==-1)call handle_error(IND_ele_mar, "brpr_updateBurp: element 212239 not found in Block_MAR_MUL_CP")
                      call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )
                    endif
                  endif

                  IND_ele = -1
                  if (iele == BUFR_NBT3) then
                    IND_ele  = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=12233)
                  elseif (iele == BUFR_NETT) then
                    IND_ele  = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=ILEMTBCOR)
                  elseif (iele == BUFR_NEES) then
                    IND_ele  = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=ILEMHBCOR)
                  end if
                      
                  if (IND_ele > 0 .and. obs_columnActive_RB(obsdat,OBS_BCOR)) then
                       call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_ele,NVAL_IND =j,NT_IND = k,RVAL = BCOR)
                  end if

                  IND_obsClear =  BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=btClearElementId)
                  if ( IND_obsClear > 0 .and. obs_columnActive_RB(obsdat,OBS_BTCL) ) then
                    Call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_obsClear,NVAL_IND =j,NT_IND = k,RVAL = obsClear) 
                  end if

                  IND_ele  = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=iele)
                  if (IND_ele==-1)call handle_error(IND_ele, "brpr_updateBurp: element not found in Block_OBS_MUL_CP")

                  call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_ele,NVAL_IND =j,NT_IND = k,RVAL = OBS)
                  
                  if (HIPCS) then
                     if (iele == 12001) then
                        IND_ele = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=12101, IOSTAT=error)
                        if (IND_ele==-1)call handle_error(IND_ele, "brpr_updateBurp: element 12101 not found in Block_OBS_MUL_CP")
                        Call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_ele,NVAL_IND =j,NT_IND = k,RVAL = OBS)
                     end if
                     if (iele == 12192) then
                        IND_ele = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=12239, IOSTAT=error)
                        if (IND_ele==-1)call handle_error(IND_ele, "brpr_updateBurp: element 12239 not found in Block_OBS_MUL_CP")
                        Call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_ele,NVAL_IND =j,NT_IND = k,RVAL = OBS)
                     end if
                  end if
                  
                  IF (HIRES .and. KOBSN > 0 ) THEN
                    STATUS=obs_headElem_i(obsdat,OBS_ST1,OBSN )
                    STATUS_HIRES=ior(STATUS_HIRES,STATUS)
                  END IF

                end do ELEMS

                IF (HIRES .and. KOBSN > 0 ) OBSN=OBSN +1

              end do LEVELS

              if ( REGRUP .and. KOBSN > 0   )  then
                STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START )
                STATUS=IBSET(STATUS,BIT_STATUS)
                ind055200 = BURP_Find_Element(Block_FLG_CP, ELEMENT=055200)
                if (ind055200==-1)call handle_error(ind055200, "brpr_updateBurp: element 55200 not found in Block_FLG_CP")
                call BURP_Set_tblval( Block_FLG_CP, NELE_IND =ind055200,NVAL_IND =1,NT_IND = k ,TBLVAL = STATUS ) 
                OBSN=OBSN +1
                OBS_START=OBS_START +1
              end if

            end do regrup_LOOP

            IF (HIRES .and. KOBSN > 0 .and. .not. regrup ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
!pik 8-2014   STATUS=IBSET(STATUS,BIT_STATUS)
              STATUS_HIRES=IBSET(STATUS_HIRES,BIT_STATUS)
              call BURP_Set_Property(CP_RPT ,FLGS =STATUS_HIRES)
            END IF
            IF (HIRES )OBS_START=OBSN
            IF (HIRES )SAVE_OBS=OBS_START
            IF (REGRUP)OBS_START=OBSN

            IF ( .not. HIRES .and. .not. regrup .and. KOBSN > 0 ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
              STATUS=IBSET(STATUS,BIT_STATUS)
              OBS_START=OBSN +1
              OBSN=OBSN +1
              call BURP_Set_Property(CP_RPT ,FLGS =STATUS)
            END IF
            IF ( .not. HIRES .and. .not. regrup .and. KOBSN == 0 ) THEN
              write(*,*) ' KOBSN=0  stnid=',stnid
            END IF
            IF ( .not. HIRES .and.  regrup .and. KOBSN == 0 ) THEN
              numHeader = obs_numHeader(obsdat)
              write(*,*)' KOBSN=0 regrup stnid=',kk,stnid,obsn,numHeader
            END IF

            IF (REGRUP) SAVE_OBS=OBS_START

            call BURP_Reduce_Block(BLOCK_OMA, NEW_NELE =il_index )
            call BURP_Reduce_Block(BLOCK_OMP, NEW_NELE =il_index )
            call BURP_Reduce_Block(BLOCK_OER, NEW_NELE =il_index )
            call BURP_Reduce_Block(BLOCK_FGE, NEW_NELE =il_index )
            call BURP_Reduce_Block(BLOCK_FSO, NEW_NELE =il_index )
            call BURP_Set_Property(BLOCK_OBS_MUL_CP ,BFAM =0)
            call BURP_Set_Property(BLOCK_MAR_MUL_CP ,BFAM =0)

            call BURP_Write_Block( CP_RPT, BLOCK_OBS_MUL_CP,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OBS_MUL_CP") 
            call BURP_Write_Block( CP_RPT, BLOCK_MAR_MUL_CP,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_MAR_MUL_CP") 

            do item=1,BN_ITEMS

              if ( BITEMLIST(item) == 'OMA') then
                call BURP_Write_Block( CP_RPT, BLOCK_OMA,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OMA") 
              end if
              if ( BITEMLIST(item) == 'OMP') then
                call BURP_Write_Block( CP_RPT, BLOCK_OMP,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OMP") 
              end if
              if ( BITEMLIST(item) == 'OER') then
                if (.not.LBLOCK_OER_CP) then
                   call BURP_Write_Block( CP_RPT, BLOCK_OER,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                   call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OER") 
                else
                   call BURP_Set_Property(BLOCK_OER_CP ,BKTYP =new_bktyp)
                   call BURP_Write_Block( CP_RPT, BLOCK_OER_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
                   call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OER_CP") 
                  end if
              end if
              if ( BITEMLIST(item) == 'FGE') then
                if (.not.LBLOCK_FGE_CP) then
                   call BURP_Write_Block( CP_RPT, BLOCK_FGE,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                   call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FGE") 
                else
                   call BURP_Set_Property(BLOCK_FGE_CP ,BKTYP =new_bktyp)
                   call BURP_Write_Block( CP_RPT, BLOCK_FGE_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)
                   call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FGE_CP") 
                  end if
              end if
              if ( BITEMLIST(item) == 'FSO') then
                 call BURP_Write_Block( CP_RPT, BLOCK_FSO,&
                        ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                 call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_FSO")
              end if

            end do

            if (regrup ) OBS_START = OBSN
            if( .not. hires ) SAVE_OBS = OBS_START

          end if  ! bl == 4

          if ( bl == 6 ) then
            call BURP_Write_Block( CP_RPT, BLOCK_GEN, ENCODE_BLOCK = .FALSE., &
                                   CONVERT_BLOCK = .FALSE., IOSTAT = error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_GEN")
          end if

          ! descriptor block (btyp = 0010 100000X XXXX) 

          BTYP10des = 160

          !if ( BTYP10 - BTYP10des == 0 ) then
          if (  bl == 1 ) then
            OBS_START=SAVE_OBS
          end if

          !==================== IASI  SPECIAL BLOCK==================
          if ( (BTYP == 9217 .or. BTYP == 15361) .and.  IDTYP == 186   ) then
            call BURP_Write_Block( CP_RPT, BLOCK_in, ENCODE_BLOCK = .FALSE., &
                                   CONVERT_BLOCK = .FALSE., IOSTAT= error)
            call handle_error(error, "brpr_updateBurp: BURP_Write_Block bloc special IASI")
          end if
          !==================== IASI  SPECIAL BLOCK==================

          !==================== GPSRO BLOCKS TO KEEP IF THEY EXIST===
          if ( IDTYP == codtyp_get_codtyp('ro') ) then
            if (BTYP ==  9217) then
              Call BURP_Write_Block( CP_RPT, BLOCK_OBS_BND, &
                   ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., &
                   IOSTAT= error)
              call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_OBS_BND")
            end if
            if (BTYP == 15361) then
              Call BURP_Write_Block( CP_RPT, BLOCK_MAR_BND, &
                   ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., &
                   IOSTAT= error)
              call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_MAR_BND")
            end if
            if (BTYP ==  9220) then
              Call BURP_Write_Block( CP_RPT, BLOCK_ORB, &
                   ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., &
                   IOSTAT= error)
              call handle_error(error, "brpr_updateBurp: BURP_Write_Block CP_RPT BLOCK_ORB")
            end if
          end if
          !==================== GPSRO BLOCKS=========================
        end do BLOCKS1

        if (  REGRUP )  then
          call BURP_Write_Block( CP_RPT, Block_FLG_CP, ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE.)
        end if

        call BURP_Delete_Report(File_in,Rpt_in, IOSTAT=error)
        call handle_error(error, "brpr_updateBurp: BURP_Delete_Report")
        call BURP_Write_Report(File_in,CP_RPT,IOSTAT= error)
        call handle_error(error, "brpr_updateBurp: BURP_Write_Report") 

      end do REPORTS

    end if

    call cleanup()

    write(*,*) ' BURPFILE  UPDATED SUM = ',trim(brp_file),SUM

  contains

    !--------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(26)

      errors(:) = 0
      deallocate(address, stat=errors(1))
      call BURP_Free(Rpt_in,CP_RPT, iostat=errors(2))
      call BURP_Free(Block_in, iostat=errors(3))
      call BURP_Free(Block_OMA, iostat=errors(4))
      call BURP_Free(Block_OMP, iostat=errors(5))
      call BURP_Free(Block_OER, iostat=errors(6))
      call BURP_Free(Block_FGE, iostat=errors(7))
      call BURP_Free(Block_FSO, iostat=errors(8))
      call BURP_Free(Block_OMA_SFC, iostat=errors(9))
      call BURP_Free(Block_OMP_SFC, iostat=errors(10))
      call BURP_Free(Block_OER_SFC, iostat=errors(11))
      call BURP_Free(Block_FGE_SFC, iostat=errors(12))
      call BURP_Free(Block_FSO_SFC, iostat=errors(13))
      call BURP_Free(Block_FLG_SFC, iostat=errors(14))
      call BURP_Free(Block_FLG, iostat=errors(15))
      call BURP_Free(Block_FLG_CP, iostat=errors(16))
      call BURP_Free(Block_MAR_MUL_CP, iostat=errors(17))
      call BURP_Free(Block_MAR_SFC_CP, iostat=errors(18))
      call BURP_Free(Block_OBS_MUL_CP, iostat=errors(19))
      call BURP_Free(Block_OBS_SFC_CP, iostat=errors(20))
      call BURP_Free(Block_GEN, iostat=errors(21))
      call BURP_Free(Block_OBS_BND, iostat=errors(23))
      call BURP_Free(Block_MAR_BND, iostat=errors(24))
      call BURP_Free(Block_ORB, iostat=errors(25))
      call BURP_Free(File_in, iostat=errors(26))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_updateBurp: error while deallocating memory: ", errors(:)
      if (associated(tree)) call kdtree2_destroy(tree)
      if (allocated(PPPandVNM)) deallocate(PPPandVNM)
      if (allocated(bodyIndexList)) deallocate(bodyIndexList)
      
    end subroutine cleanup

    !--------- handle_error -----
    subroutine handle_error(icode, errorMessage)
      implicit none
      character (len=*) :: errorMessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errorMessage))
      end if
    end subroutine handle_error

  end subroutine brpr_updateBurp


  SUBROUTINE BRPACMA_NML(NML_SECTION, beSilent_opt)

    IMPLICIT NONE

    logical, optional  :: beSilent_opt
    logical            :: beSilent
    INTEGER            :: NULNAM,IER,FNOM,FCLOS
    CHARACTER *256     :: NAMFILE
    CHARACTER(len = *) :: NML_SECTION

    NAMELIST /NAMBURP_FILTER_CONV/NELEMS,    BLISTELEMENTS,    BNBITSOFF,BBITOFF,BNBITSON,BBITON, &
         ENFORCE_CLASSIC_SONDES,UA_HIGH_PRECISION_TT_ES,UA_FLAG_HIGH_PRECISION_TT_ES,READ_QI_GA_MT_SW
    NAMELIST /NAMBURP_FILTER_SFC/ NELEMS_SFC,BLISTELEMENTS_SFC,BNBITSOFF,BBITOFF,BNBITSON,BBITON,NELEMS_GPS,LISTE_ELE_GPS
    NAMELIST /NAMBURP_FILTER_TOVS/NELEMS,BLISTELEMENTS,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_FILTER_CHM_SFC/NELEMS_SFC,BLISTELEMENTS_SFC,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_FILTER_CHM/NELEMS,BLISTELEMENTS,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_UPDATE/BN_ITEMS, BITEMLIST,TYPE_RESUME

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    NAMFILE=trim("flnml")
    nulnam=0
    IER=FNOM(NULNAM,NAMFILE,'R/O',0)
    write(*,*) ' READ NML_SECTION =',trim(NML_SECTION)

    SELECT CASE(trim(NML_SECTION))
      CASE( 'namburp_sfc')
        READ(NULNAM,NML=NAMBURP_FILTER_SFC)
        if (.not.beSilent) write(*,nml=NAMBURP_FILTER_SFC)
      CASE( 'namburp_conv')
        READ(NULNAM,NML=NAMBURP_FILTER_CONV)
        if (.not.beSilent) write(*,nml=NAMBURP_FILTER_CONV)
      CASE( 'namburp_tovs')
        READ(NULNAM,NML=NAMBURP_FILTER_TOVS)
        if (.not.beSilent) write(*,nml=NAMBURP_FILTER_TOVS)
      CASE( 'namburp_chm_sfc')
        READ(NULNAM,NML=NAMBURP_FILTER_CHM_SFC)
        if (.not.beSilent) write(*,nml=NAMBURP_FILTER_CHM_SFC)
      CASE( 'namburp_chm')
        READ(NULNAM,NML=NAMBURP_FILTER_CHM)
        if (.not.beSilent) write(*,nml=NAMBURP_FILTER_CHM)
      CASE( 'namburp_update')
        READ(NULNAM,NML=NAMBURP_UPDATE)
        if (.not.beSilent) write(*,nml=NAMBURP_UPDATE)
    END SELECT

    ier=FCLOS(NULNAM)
  END SUBROUTINE BRPACMA_NML


  subroutine brpr_readBurp(obsdat,familytype,brp_file,filenumb)
    !
    !:Purpose: Select variables relative to airs in burp file. Read burp file.

    !***********************************************************************
    !
    !          WHEN SEARCHING FOR A SPECIFIC BLOCK BY ITS BTYP, VALUES OF
    !          BIT 0 TO 3 ARE IRRELEVANT WHILE BIT 4 IS 0 FOR GLOBAL AND 1
    !          FOR REGIONAL MODEL. HERE, WE SEARCH BLOCK BY THEIR FIRST
    !          10 BITS (BIT 5 TO 14).
    !
    !***********************************************************************
    IMPLICIT NONE

    ! Arguments
    type (struct_obs), intent(inout) :: obsdat
    CHARACTER *2           :: FAMILYTYPE
    CHARACTER(LEN=128)     :: BRP_FILE ! name of burp file
    integer                :: FILENUMB

    ! Locals:
    CHARACTER *2           :: UNI_FAMILYTYPE

    TYPE(BURP_FILE)        :: FILE_IN
    TYPE(BURP_RPT)         :: RPT_IN
    TYPE(BURP_BLOCK)       :: BLOCK_IN
    
    CHARACTER(LEN=5)       :: FAMILYTYPE2
    CHARACTER(LEN=9)       :: OPT_MISSING
    integer                :: BTYP,BFAM,BKSTP,BTYP10,BTYP10FLG_uni,BTYP10obs_uni 
    integer                :: BTYP10DES,BTYP10INF,BTYP10OBS,BTYP10FLG

    integer                :: NB_RPTS,REF_RPT,REF_BLK,COUNT, numHeader
    INTEGER, ALLOCATABLE   :: ADDRESS(:)

    real   , ALLOCATABLE   :: OBSVALUE(:,:,:),OBSVALUE_SFC(:,:,:)
    real   , ALLOCATABLE   :: OBSERV  (:,:),    OBSERV_SFC(:,:)
    real   , ALLOCATABLE   :: BiasCorrection_sfc(:,:,:)
    real   , ALLOCATABLE   :: BCOR_SFC(:,:)
    INTEGER, PARAMETER     :: MAXRONVAL=500
    real                   :: ROLAT(MAXRONVAL), ROLON(MAXRONVAL)

    INTEGER, ALLOCATABLE   :: MTVAL(:)
    INTEGER, ALLOCATABLE   :: HAVAL(:), GAVAL(:), QI1VAL(:) ,QI2VAL(:), LSVAL(:)
    real(pre_obsReal) , ALLOCATABLE  :: azimuth(:)
    INTEGER, ALLOCATABLE   :: QCFLAG  (:,:,:),  QCFLAG_SFC(:,:,:)
    INTEGER, ALLOCATABLE   :: QCFLAGS (:,:),   QCFLAGS_SFC(:,:)
    integer, allocatable   :: hiresTimeFlag(:,:), hiresLatFlag(:,:)

    real   , ALLOCATABLE   :: VCOORD  (:,:),  VCOORD_SFC(:)
    real   , ALLOCATABLE   :: VCORD   (:)

    INTEGER, ALLOCATABLE   :: LAT(:),LON(:),HHMM(:),DATE(:),GLBFLAG(:)
    real   , ALLOCATABLE   :: HLAT(:,:),  HLON(:,:),  HTIME(:,:)
    INTEGER, ALLOCATABLE   :: PHASE(:,:)
    INTEGER, ALLOCATABLE   :: dataQcFlagLEV(:),dataQcFlag2(:,:)
    integer                :: dataQcFlagLEV_sfc(1)
    INTEGER, ALLOCATABLE   :: dataCloudFracLEV(:),dataCloudFrac(:,:)
    integer                :: dataCloudFracLEV_sfc(1)
    real   , ALLOCATABLE   :: HLAT_SFC(:),HLON_SFC(:),HTIME_SFC(:)

    real   , ALLOCATABLE   :: RINFO(:,:)
    real   , ALLOCATABLE   :: TRINFO(:)

    real   , ALLOCATABLE   :: EMIS(:,:), SURF_EMIS(:)
    real   , ALLOCATABLE   :: BCOR(:,:), BiasCorrection(:,:)
    real   , ALLOCATABLE   :: BCOR2(:,:,:), BiasCorrection2(:,:)

    REAL(pre_obsReal), ALLOCATABLE  :: CFRAC(:,:)

    REAL(pre_obsReal), ALLOCATABLE :: RADMOY(:,:,:)
    REAL(pre_obsReal), ALLOCATABLE :: radstd(:,:,:)

    integer                :: LISTE_INFO(29),LISTE_ELE(20),LISTE_ELE_SFC(20)
    
    integer                :: NBELE,NVALE,NTE
    integer                :: J,JJ,K,KK,KL,IL,ERROR,OBSN
    integer                :: info_elepos,IND_ELE,IND_VCOORD,IND_QCFLAG,IND_SW
    integer                :: IND055200,IND4208,ind4197,IND5002,IND6002,ind_al,IND5001,IND6001
    integer                :: IND_LAT,IND_LON,IND_TIME,IND_EMIS,IND_BCOR,IND_PHASE,IND_BCOR_TT,IND_BCOR_HU
    integer                :: FLAG_PASSAGE1,FLAG_PASSAGE2,FLAG_PASSAGE3,FLAG_PASSAGE4
    integer                :: IND_dataQcFlag0, IND_dataQcFlag1, IND_dataQcFlag2, IND_dataCloudFrac 

    integer                :: vcord_type(10),SUM,vcoord_type
    REAL(pre_obsReal)      :: RELEV,XLAT,XLON,RELEV2
    real                   :: XTIME,ROLAT0,ROLON0,ROLAT1,ROLON1
    integer                :: status ,idtyp,lati,long,dx,dy,elev
    integer                :: drnd,date_h,hhmm_h,oars,runn,YMD_DATE,HM,kstamp,kstamp2,HM_SFC,YMD_DATE_SFC

    integer                :: iele,NELE,NELE_SFC,NVAL,NT,NELE_INFO,LN
    integer                :: bit_alt,btyp_offset,btyp_offset_uni
    character(len = 5)     :: BURP_TYP
    CHARACTER(LEN=9)       :: STNID,STN_RESUME
    LOGICAL                :: HIRES,HIRES_SFC,HIPCS,phasePresent,LOK,LROK
    integer                :: NDATA,NDATA_SF
    integer                :: IER,date2,time2,time_sonde,NEWDATE
    real                   :: RAD_MOY,RAD_STD
    integer                :: iclass,NCHANAVHRR,NCLASSAVHRR,ichan,iobs,inorm
    integer                :: infot
    integer                :: ILEMZBCOR, ILEMTBCOR, ILEMHBCOR
    
    
    LISTE_INFO(1:29) = (/ 1007,002019,007024,007025 ,005021, 005022, 008012, &
        013039,020010,2048,2022,33060,33062,33039,10035,10036,08046,5043, &
        013209,clwFgElementId,1033,2011,4197,5040,33078,33079,33080,020029,25174 /)

    RELEV2=0.0
    FAMILYTYPE2= 'SCRAP'
    vcord_type(:)=-1
    vcord_type(1)=0
    NELE_INFO=1
    NELE_SFC=0
    NELE=0
    BNBITSOFF=0
    BNBITSON=0
    ILEMZBCOR=15234 ! bcor element for GP  ZTD observations
    ILEMTBCOR=12204 ! bcor element for altitude TT observations
    ILEMHBCOR=12243 ! bcor element for altitude ES observations
    ENFORCE_CLASSIC_SONDES=.false.
    UA_HIGH_PRECISION_TT_ES=.false.
    UA_FLAG_HIGH_PRECISION_TT_ES=.false.
    READ_QI_GA_MT_SW=.false.
    UNI_FAMILYTYPE = 'SF'
    LISTE_ELE_SFC(:)=-1
    LISTE_ELE(:)=-1
    SELECT CASE(trim(FAMILYTYPE))
      CASE('UA')
        BURP_TYP='multi'
        vcord_type(1)=7004

        LISTE_ELE_SFC(1:6) = (/12004,11011,11012,10051,10004,12203/)
        NELE_SFC=6
        call BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'UA'
        LISTE_ELE(1:5) = (/12001,11001,11002,12192,10194/)
        NELE=5
        ENFORCE_CLASSIC_SONDES=.false.
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        NELE_INFO=23
      CASE('AI')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:4) = (/12001,12192,11001,11002/)
        NELE=4
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('AL')
        BURP_TYP='uni'
        vcord_type(1)=7071

        LISTE_ELE(1) = 40030
        NELE=1
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('SW')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('SF')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELEMS_SFC=6
        LISTE_ELE_SFC(1:NELEMS_SFC) = (/10004,12004,10051,12203,11011,11012/)
        BLISTELEMENTS_SFC(1:NELEMS_SFC) = LISTE_ELE_SFC(1:NELEMS_SFC) ! default list

        call BRPACMA_NML('namburp_sfc') ! read NELEMS_SFC, BLISTELEMENTS_SFC(1:NELEMS_SFC)
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'SFC'
      CASE('GP')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELEMS_GPS=6
        LISTE_ELE_GPS(1:NELEMS_GPS) = (/10004,12004,12203,15031,15032,15035/) ! default list

        call BRPACMA_NML('namburp_sfc') ! read NELEMS_GPS, LISTE_ELE_GPS(1:NELEMS_GPS)
        NELE_SFC=NELEMS_GPS             !   -- ignore NELEMS_SFC, BLISTELEMENTS_SFC(1:NELEMS_SFC)
        BLISTELEMENTS_SFC(1:NELEMS_GPS) = LISTE_ELE_GPS(1:NELEMS_GPS)

        FAMILYTYPE2= 'SFC'
        UNI_FAMILYTYPE = 'GP'
      CASE('SC')
        vcord_type(1)=0
        BURP_TYP='uni'
        LISTE_ELE_SFC(1:2) = (/11012,11011/)
        NELE=2
        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS

        FAMILYTYPE2= 'UASFC2'
      CASE('PR')
        BURP_TYP='multi'
        vcord_type(1)=7006

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2

        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('RO')
        BURP_TYP='multi'
        vcord_type(1)=7007
        vcord_type(2)=7040

        LISTE_ELE(1:2) = (/15036,15037/)
        NELE=2

        call BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        !================GPS-RO CANNOT BE FILTERED=======
        BNBITSOFF=0
        BNBITSON=0
        !================GPS-RO CANNOT BE FILTERED=======
        NELE_INFO=18
      CASE('TO')
        BURP_TYP='multi'
        vcord_type(1)=5042
        vcord_type(2)=2150

        LISTE_ELE(1:1) = (/12163/)
        NELE=1

        call BRPACMA_NML('namburp_tovs')
        NELE=NELEMS

        NELE_INFO=29
     CASE('CH')

        BURP_TYP='multi'  ! Both 'multi' and 'uni' are possible for this family.
                          ! 'uni' level data are assumed not to have any accompanynig vertical 
                          ! coordinate element in addition to having only one level.
        vcord_type(1:8) = (/7004,7204,7006,7007,5042,2150,2071,0/) ! 0 must be at end.
        NELE_INFO=18

        UNI_FAMILYTYPE = 'CH'
        LISTE_ELE_SFC(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                          15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE_SFC=19
        call BRPACMA_NML('namburp_chm_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2='CH'
        LISTE_ELE(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                          15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE=19
        call BRPACMA_NML('namburp_chm')
        NELE=NELEMS 
      CASE('GO','MI')
        call utl_abort('brpr_readBurp: unknown familyType : ' // trim(familyType))
    END SELECT

    LISTE_ELE    (1:NELE    )=BLISTELEMENTS(1:NELE)
    LISTE_ELE_SFC(1:NELE_SFC)=BLISTELEMENTS_SFC(1:NELE_SFC)
    if (NELE     > 0) then
      write(*,*)  ' LISTE_ELE =',LISTE_ELE
      call ovt_setup(LISTE_ELE(1:NELE))
    end if
    if (NELE_SFC > 0) then
      write(*,*)  ' LISTE_ELE_SFC =',LISTE_ELE_SFC(1:NELE_SFC)
      call ovt_setup(LISTE_ELE_SFC(1:NELE_SFC))
    end if

    write(*,*) ' BNBITSON BNBITSOFF     =',BNBITSON,BNBITSOFF

    btyp_offset_uni=-999
    btyp_offset=-999
    if (trim(BURP_TYP) == 'uni') then
      btyp_offset=256
    else
      btyp_offset=0
    end if

    if (trim(familytype) == 'AL') then
      btyp_offset=255
    end if

    if (TRIM(FAMILYTYPE2) == 'SFC') then
      btyp_offset= btyp_offset+32
      btyp_offset_uni= 256 +32
    elseif ( TRIM(FAMILYTYPE2) == 'UA') then
      btyp_offset_uni= 256 +32
    elseif (TRIM(FAMILYTYPE2) == 'CH') then
      btyp_offset_uni= 256
    else
      btyp_offset_uni= -999 ! set to -999 when not used
    end if

    write(*,*) '-----------------------------------------------'
    write(*,*) '-----------  BEGIN brpr_readBurp   ------------'
    write(*,*) 'FAMILYTYPE vcord_type   =',FAMILYTYPE,vcord_type
    write(*,*) 'BURP_TYP btyp_offset    =',BURP_TYP, btyp_offset
    write(*,*) '-----------------------------------------------'


    ! initialisation
    ! --------------
    SUM=0
    flag_passage1 = 0
    flag_passage2 = 0
    flag_passage3 = 0
    flag_passage4 = 0
    ! initialiase the qc flag2 indice

    opt_missing = 'MISSING'


    call BURP_Set_Options( &
       & REAL_OPTNAME       = opt_missing, &
       & REAL_OPTNAME_VALUE = MPC_missingValue_R4, &
       & CHAR_OPTNAME       = 'MSGLVL', &
       & CHAR_OPTNAME_VALUE = 'FATAL', &
       & IOSTAT          = error )
    call handle_error(error, "brpr_readBurp: BURP_Set_Options")

    call BURP_Init(File_in      ,IOSTAT=error)
    call handle_error(error, "brpr_readBurp: BURP_init File_in")
    call BURP_Init(Rpt_in)
    call BURP_Init(Block_in)

    ! opening file
    write(*,*) 'OPENING BURP FILE FOR READING = ', trim(brp_file)

    call BURP_New(File_in, FILENAME = brp_file, &
                   & MODE    = FILE_ACC_READ, &
                   & IOSTAT  = error )
    call handle_error(error, "brpr_readBurp: Burp_New. Problem opening " // trim(brp_file))


    ! obtain input burp file number of reports

    call BURP_Get_Property(File_in, NRPTS=nb_rpts)

    write(*,*) '-----------------------------------------'
    write(*,*) 'IOSTAT    =',error
    write(*,*) 'NUMBER OF REPORTS IN FILE  = ',nb_rpts
    write(*,*) '-----------------------------------------'


    ! scan input burp file to get all reports address

    Allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0
    bit_alt = 0
    stn_resume='NOT_FOUND'

    do
      ref_rpt = BURP_Find_Report(File_in, &
                    & REPORT      = Rpt_in,  &
                    & SEARCH_FROM = ref_rpt, &
                    & IOSTAT      = error)
      call handle_error(error, "brpr_readBurp: BURP_Find_Report")
      call burp_get_property(Rpt_in, STNID = stnid )
      IF ( stnid(1:2) == ">>" ) then
        STN_RESUME=stnid
        TYPE_RESUME=STN_RESUME(3:9)
        SELECT CASE(stnid)
          CASE(">>BGCKALT", ">>POSTALT")
            bit_alt=1
          CASE(">>DERIALT")
            bit_alt=2
          CASE DEFAULT
            write(*,*) 'brpr_readBurp: WARNING: Unknown RESUME record found, assume BGCKALT'
            bit_alt=1
        END SELECT 
      END IF

      if (ref_rpt < 0) Exit

      if (count == nb_rpts) then
        write(*,*) 'brpr_readBurp: ERROR: count = nb_rpts:',count,nb_rpts
        exit
      end if

      count = count + 1
      address(count) = ref_rpt
    end do

    if (stn_resume == 'NOT_FOUND') then
      write(*,*) 'brpr_readBurp: WARNING: No RESUME record found in this file, ' //  &
                 'check if already read in another file'
      ! try to get value from previously read file
      if ( type_resume /= 'UNKNOWN' ) then
        stn_resume = '>>' // type_resume
        SELECT CASE(stn_resume)
          CASE(">>BGCKALT", ">>POSTALT")
            bit_alt=1
          CASE(">>DERIALT")
            bit_alt=2
          CASE DEFAULT
            write(*,*) 'brpr_readBurp: WARNING: Unknown RESUME record found, assume BGCKALT'
            stn_resume = '>>BGCKALT'
            bit_alt=1
        END SELECT 
      else
        write(*,*) 'brpr_readBurp: WARNING: No file read has RESUME record, assume BGCKALT'
        stn_resume = '>>BGCKALT'
        bit_alt=1
      end if
    end if

    write(*,*) STN_RESUME,' bit_alt==== >  ',bit_alt


    BTYP10obs     = 291 -btyp_offset
    BTYP10obs_uni = 291 -btyp_offset_uni
    if (bit_alt == 2) btyp10obs =  BTYP10obs - 2
    if (bit_alt == 2) btyp10obs_uni =  BTYP10obs_uni - 2

    BTYP10flg     = 483 -btyp_offset
    BTYP10flg_uni = 483 -btyp_offset_uni
    if (bit_alt == 2) BTYP10flg =  BTYP10flg  - 2
    if (bit_alt == 2) BTYP10flg_uni =  BTYP10flg_uni  - 2

    write(*, *)  ' NUMBER OF VALID REPORTS IN FILE = ',count
    write(*, *)  ' BTYP10obs BTYP10obs_uni         = ',BTYP10obs,BTYP10obs_uni

    if ( count > 0 ) then
      ! LOOP ON REPORTS
      REPORTS: do kk = 1, count
            
        call BURP_Get_Report(File_in, &
                      & REPORT    = Rpt_in, &
                      & REF       = address(kk), &
                      & IOSTAT    = error)
        call handle_error(error, "brpr_readBurp: BURP_Get_Report") 
        call burp_get_property(Rpt_in, &
               STNID = stnid ,TEMPS =hhmm_h,FLGS = status ,IDTYP =idtyp,LATI = lati &
               ,LONG = long ,DX = dx ,DY = dy,ELEV=elev,DRND =drnd,DATE =date_h &
               ,OARS =oars,RUNN=runn ,IOSTAT=error)
        call handle_error(error, "brpr_readBurp: burp_get_property")
        IF ( stnid(1:2) == ">>" ) cycle
        !  LOOP ON BLOCKS

        ! Ensure the dataQcFlag arrays are not allocated before looping over blocks
        if (allocated(dataQcFlag2)) deallocate(dataQcFlag2)
        if (allocated(dataQcFlagLEV)) deallocate(dataQcFlagLEV)
        
        ! Ensure the dataCloudFrac arrays are not allocated before looping over blocks
        if (allocated(dataCloudFrac)) deallocate(dataCloudFrac)
        if (allocated(dataCloudFracLEV)) deallocate(dataCloudFracLEV)
        
        ref_blk = 0

        HIRES=.FALSE.
        HIPCS=.FALSE.
        HIRES_SFC=.FALSE.
        phasePresent = .false.
        LROK = .false.

        BLOCKS1: do

          ref_blk = BURP_Find_Block(Rpt_in, &
                    & BLOCK       = Block_in, &
                    & SEARCH_FROM = ref_blk, &
                    & IOSTAT      = error)
          call handle_error(error, "brpr_readBurp: BURP_Find_Block")

          if (ref_blk < 0) EXIT BLOCKS1

          call BURP_Get_Property(Block_in, &
                      & NELE   = nbele, &
                      & NVAL   = nvale, &
                      & NT     = nte, &
                      & BFAM   = bfam, &
                      & BTYP   = btyp, &
                      & BKSTP  = BKSTP, &
                      & IOSTAT = error)
          call handle_error(error, "brpr_readBurp: BURP_Get_Property")

          ! Read slant latlon if type is RO
          if (trim(familytype) == 'RO' .and. bfam == 0 .and. LROK == .FALSE.) then
            ROLAT0 = 0.01*lati- 90.
            ROLON0 = 0.01*long
            if (ROLON0 > 180.) ROLON0 = ROLON0-360.
            ROLAT(:) = ROLAT0
            ROLON(:) = ROLON0
            IND5001 = BURP_Find_Element(Block_in, ELEMENT = 5001)
            IND6001 = BURP_Find_Element(Block_in, ELEMENT = 6001)
            if (IND5001 > 0 .and. IND6001 > 0) then
              do j = 1, nvale
                ROLAT1 = BURP_Get_Rval(Block_in, &
                                       NELE_IND = IND5001, &
                                       NVAL_IND = j, &
                                       NT_IND = 1, IOSTAT = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Rval Block_in 5001")
                ROLON1 = BURP_Get_Rval(Block_in, &
                                       NELE_IND = IND6001, &
                                       NVAL_IND = j, &
                                       NT_IND = 1, IOSTAT = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Rval Block_in 6001")
                lok = ( -90.1 < ROLAT1 .and. ROLAT1 <  90.1) .and. &
                      (-180.1 < ROLON1 .and. ROLON1 < 360.1)
                if (lok .and. j<=MAXRONVAL) then
                  ROLAT(j) = ROLAT1
                  ROLON(j) = ROLON1
                  LROK = .TRUE.
                end if
              end do
            end if
          end if

          ! observation block (btyp = 0100 100011X XXXX)
          if(trim(familytype) == 'AL')then

            ! Fudge the block type, because the data are simulated
            if(btyp == 1024) then
              btyp=1152
            else if(btyp == 7168)then
              btyp=7296
            end if
          end if
          btyp10    = ishft(btyp,-5)
          if ( btyp10 == btyp10obs_uni .and. bkstp <= 4 ) then   ! FAMILYTYPE = SF, GP data blocks

            flag_passage3 = 1

            allocate(obsvalue_sfc(NELE_SFC,1,nte))
            allocate(BiasCorrection_sfc(NELE_SFC,1,nte))
            allocate( OBSERV_SFC(NELE_SFC,1) )
            allocate( BCOR_SFC(NELE_SFC,1) )

            allocate( vcoord_sfc(1))

            vcoord_type=0
            vcoord_SFC(:)       = 0
            obsvalue_sfc(:,:,:) = MPC_missingValue_R4
            BiasCorrection_sfc(:,:,:) = MPC_missingValue_R4
            IND_LAT   = BURP_Find_Element(Block_in, ELEMENT=5001)
            IND_LON   = BURP_Find_Element(Block_in, ELEMENT=6001)
            IND_TIME  = BURP_Find_Element(Block_in, ELEMENT=4015)
            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES_SFC=.true.
            if (HIRES_SFC) allocate(HLAT_SFC(nte),HLON_SFC(nte),HTIME_SFC(nte) )
            IF (HIRES_SFC) THEN
              do k=1,nte
                HLAT_SFC(k) =BURP_Get_Rval(Block_in, &
                                         & NELE_IND = IND_LAT, &
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
                HLON_SFC(k) =BURP_Get_Rval(Block_in, &
                                         & NELE_IND = IND_LON, &
        
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
                HTIME_SFC(k)=BURP_Get_Rval(Block_in, &
                                         & NELE_IND = IND_TIME, &
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
              end do
            END IF

            do IL = 1, NELE_SFC

              iele = LISTE_ELE_SFC(IL)
              IND_ele  = BURP_Find_Element(Block_in, ELEMENT=iele)
              if (IND_ele < 0 ) cycle

              do k=1,nte
                obsvalue_sfc(IL,1,k) =  BURP_Get_Rval(Block_in, &
                                         & NELE_IND = IND_ele, &
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
              end do
              
              IND_ele = -1
              if (iele == BUFR_NEZD) then
                IND_ele  = BURP_Find_Element(Block_in, ELEMENT=ILEMZBCOR)
              end if
              if (IND_ele > 0) then
                do k=1,nte
                  BiasCorrection_sfc(IL,1,k) = BURP_Get_Rval(Block_in, &
                                                & NELE_IND = IND_ele, &
                                                & NVAL_IND = 1, &
                                                & NT_IND   = k)
                end do
              end if

            end do
            
          end if

          if ( btyp10 == btyp10flg_uni .and. bkstp <= 4 ) then  ! FAMILYTYPE = SF, GP flag blocks

            flag_passage4 = 1

            allocate( qcflag_sfc (NELE_SFC,1,nte))
            allocate( qcflags_SFC(NELE_SFC,1) )
            QCFLAGS_SFC(:,:)=0

            do IL = 1, NELE_SFC
              iele=LISTE_ELE_SFC(IL) + 200000
              IND_QCFLAG  = BURP_Find_Element(Block_in, ELEMENT=iele)
              if (IND_QCFLAG < 0 ) cycle
              DO k=1,nte
                QCFLAG_sfc(IL,1,k) =  BURP_Get_Tblval(Block_in, &
                                         & NELE_IND = IND_QCFLAG, &
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
                SUM = SUM +1                
              END DO
            end do
          end if

          if ( btyp10 == btyp10obs .and. bfam == 0 ) then ! non sfc-type DATA block

            flag_passage3 = 1

            NVAL=NVALE ;  NT=NTE
            allocate(obsvalue(NELE,nvale,nte),VCOORD(nvale,nte))
            allocate(OBSERV(NELE,nvale))
            allocate(VCORD(nvale))
                           
            obsvalue(:,:,:) = MPC_missingValue_R4
            OBSERV  (:,:)   = MPC_missingValue_R4
            VCOORD  (:,:)   = 0.
            VCORD   (:)     = 0.

            k=0
            IND_VCOORD=-1
            do while (vcord_type(k+1) /= -1 .and. IND_VCOORD == -1)
               k=k+1 
               IND_VCOORD  = BURP_Find_Element(Block_in, ELEMENT=vcord_type(k))
            end do
            vcoord_type=0
            IF (IND_VCOORD > 0) vcoord_type=vcord_type(k)
            !if (IND_VCOORD == -1)write(*,*) 'PAS DE COORDONNEE VERTICALE STNID=',STNID,trim(FAMILYTYPE)

            ! LAT LON TIME IN DATA BLOCK
            IND_LAT   = BURP_Find_Element(Block_in, ELEMENT=5001 )
            IND_LON   = BURP_Find_Element(Block_in, ELEMENT=6001 )
            IND_TIME  = BURP_Find_Element(Block_in, ELEMENT=4015 )
            IND_EMIS  = BURP_Find_Element(Block_in, ELEMENT=55043)

            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES=.true.

            IND_BCOR    = -1
            IND_BCOR_TT = -1
            IND_BCOR_HU = -1
            IND_PHASE   = -1
            phasePresent = .false.
            if ( FAMILYTYPE == 'TO' ) IND_BCOR  = BURP_Find_Element(Block_in, ELEMENT=12233)
            if ( FAMILYTYPE == 'AI' .or. FAMILYTYPE2 == 'UA') then
               IND_BCOR_TT  = BURP_Find_Element(Block_in, ELEMENT=ILEMTBCOR)
               IND_BCOR_HU  = BURP_Find_Element(Block_in, ELEMENT=ILEMHBCOR)
            end if
            if ( FAMILYTYPE == 'AI' ) then
               IND_PHASE = BURP_Find_Element(Block_in, ELEMENT=8004)
               if (IND_PHASE > 0) then
                  allocate( phase(nvale,nte) )
                  phase(:,:) = MPC_missingValue_R4
                  phasePresent = .true.
               end if
            end if

            if( (TRIM(FAMILYTYPE2) == 'UA') .and. UA_HIGH_PRECISION_TT_ES ) HIPCS=.true.

            if(HIRES) allocate(HLAT(nvale,nte),HLON(nvale,nte),HTIME(nvale,nte) )

            ! If ATMS or AMSUA, or AMSUB read the element 33081 or 33082
            IND_dataQcFlag2 = -1
            if ( idtyp == 192 .or. idtyp == 164 .or. idtyp == 181 .or. idtyp == 182) then
              IND_dataQcFlag0 = BURP_Find_Element(Block_in, ELEMENT=33081)
              IND_dataQcFlag1 = BURP_Find_Element(Block_in, ELEMENT=33032)
              if ( IND_dataQcFlag0 > 0 .and. IND_dataQcFlag1 > 0 ) then 
                call  utl_abort('readBurp : Got two valid indices for IND_dataQcFlag2 in family' // trim(familyType))
              elseif ( IND_dataQcFlag0 > 0 .and. IND_dataQcFlag1 < 0 ) then 
                IND_dataQcFlag2 = IND_dataQcFlag0
              elseif ( IND_dataQcFlag0 < 0 .and. IND_dataQcFlag1 > 0 ) then 
                IND_dataQcFlag2 = IND_dataQcFlag1
              end if
            end if
            
            ! Allocate arrays for dataQcFlag if they are found in the file
            if (IND_dataQcFlag2 > 0) then
              allocate(dataQcFlag2(nvale,nte))
              allocate(dataQcFlagLEV(nvale))
              dataQcFlag2(:,:) = MPC_missingValue_INT
              dataQcFlagLEV(:) = MPC_missingValue_INT
            end if

            ! if CSR data idtyp = 185, then read ele 020081
            IND_dataCloudFrac = -1
            if ( idtyp == 185) then
              IND_dataCloudFrac = BURP_Find_Element(Block_in, ELEMENT=020081)
            end if 
            ! Allocate arrays for dataCloudFrac if they are found in the file
            if (IND_dataCloudFrac > 0) then
              allocate(dataCloudFrac(nvale,nte))
              allocate(dataCloudFracLEV(nvale))
              dataCloudFrac(:,:) = MPC_missingValue_INT
              dataCloudFracLEV(:) = MPC_missingValue_INT
            end if


            allocate(EMIS(nvale,nte))
            allocate(SURF_EMIS(nvale))
            EMIS(:,:)       = MPC_missingValue_R4
            
            OBSVALUE(:,:,:) = MPC_missingValue_R4
            
            if (IND_BCOR > 0) then                              ! TOVS
               allocate(BCOR(nvale,nte))
               allocate(BiasCorrection(nele,nvale))
               BiasCorrection(:,:) = 0.0
               BCOR(:,:) =  MPC_missingValue_R4
            elseif (IND_BCOR_TT > 0 .or. IND_BCOR_HU > 0) then  ! conventional (UA or AI)
               allocate(BCOR2(nele,nvale,nte))
               allocate(BiasCorrection2(nele,nvale))
               BCOR2(:,:,:) =  MPC_missingValue_R4
            end if
           

            ! Get the observations and conventional data bias corrections for each element in LISTE_ELE
              
            do IL = 1, NELE

              iele = LISTE_ELE(IL)
              IND_ele = BURP_Find_Element(Block_in, ELEMENT=iele)
              if (IND_ele < 0 ) cycle

              if (HIPCS .and. iele == 12001) then
                IND_ele = BURP_Find_Element(Block_in, ELEMENT=12101)
                if (IND_ele == -1) call handle_error(IND_ele, "brpr_readBurp: cannot find element 12101 in Block_in")
              end if
              if (HIPCS .and. iele == 12192) then
                IND_ele = BURP_Find_Element(Block_in, ELEMENT=12239)
                if (IND_ele == -1) call handle_error(IND_ele, "brpr_readBurp: cannot find element 12239 in Block_in")
              end if
              do k=1,nte
                do j=1,nvale
                  obsvalue(IL,j,k) =  BURP_Get_Rval(Block_in,NELE_IND=IND_ele,NVAL_IND=j,NT_IND=k)
                  if (iele == BUFR_NETT .and. IND_BCOR_TT > 0) &
                    & BCOR2(IL,j,k) =  BURP_Get_Rval(Block_in,NELE_IND=IND_BCOR_TT,NVAL_IND=j,NT_IND=k)
                  if (iele == BUFR_NEES .and. IND_BCOR_HU > 0) &
                    & BCOR2(IL,j,k) =  BURP_Get_Rval(Block_in,NELE_IND=IND_BCOR_HU,NVAL_IND=j,NT_IND=k)
                end do
              end do
              
            end do
            
            ! Get other needed elements including vccord, TOVS bias corrections and AI phase of flight
            do k=1,nte
              do j=1,nvale
                  IF (HIRES) THEN
                    HLAT(j,k) = BURP_Get_Rval(Block_in,NELE_IND = IND_LAT,NVAL_IND = j,NT_IND = k)
                    HLON(j,k) = BURP_Get_Rval(Block_in,NELE_IND = IND_LON,NVAL_IND = j,NT_IND = k)
                    HTIME(j,k)= BURP_Get_Rval(Block_in,NELE_IND = IND_TIME,NVAL_IND = j,NT_IND = k)
                  END IF
                  if ( phasePresent ) then
                    phase(j,k) = BURP_Get_Tblval(Block_in,NELE_IND = IND_phase,NVAL_IND = j,NT_IND = k)
                  end if
                  if ( IND_dataQcFlag2 > 0 ) then
                    dataQcFlag2(j,k) = BURP_Get_Tblval(Block_in,NELE_IND = IND_dataQcFlag2,NVAL_IND = j,NT_IND = k)
                  end if
                  if ( IND_dataCloudFrac > 0 ) then
                    dataCloudFrac(j,k) = BURP_Get_Tblval(Block_in,NELE_IND = IND_dataCloudFrac,NVAL_IND = j,NT_IND = k)
                  end if
                  IF (IND_EMIS > 0) THEN
                    EMIS(j,k) = BURP_Get_Rval(Block_in,NELE_IND = IND_EMIS,NVAL_IND = j,NT_IND = k)
                  END IF
                  IF (IND_BCOR > 0) THEN
                    BCOR(j,k) = BURP_Get_Rval(Block_in,NELE_IND = IND_BCOR,NVAL_IND = j,NT_IND = k)
                  END IF
                  if (IND_VCOORD > 0) then
                    VCOORD(j,k) = BURP_Get_Rval(Block_in,NELE_IND = IND_VCOORD,NVAL_IND = j,NT_IND = k)
                  end if
              end do
            end do

!==================================================================================
!
! Lire les divers elements de la famille SW
!
            allocate(qi1val(nte))
            allocate(qi2val(nte))
            allocate(mtval(nte))
            allocate(lsval(nte))
            allocate(haval(nte))
            allocate(gaval(nte))
            QI1VAL(:) = 0
            QI2VAL(:) = 0
            MTVAL (:) = 0
            LSVAL (:) = 0
            HAVAL (:) = 0
            GAVAL (:) = 0

            if (TRIM(FAMILYTYPE) == 'SW' .and. READ_QI_GA_MT_SW) then

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=33007)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                QI1VAL(k)= BURP_Get_Tblval(Block_in, &
                                           NELE_IND = IND_SW, &
                                           NVAL_IND = 1, &
                                           NT_IND   = k, &
                                           IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 33007")
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=33194)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                QI2VAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 33194")
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=2023)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                MTVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 2023")
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=8012)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                LSVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 8012")
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=13039)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                GAVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 13039")
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=2163)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                HAVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 2163")
              end do

            !====================================================================
            !
            ! Lire les divers elements de la famille AL
            !
            else if (trim(familytype) == 'AL') then
              if (.not. allocated(azimuth)) then
                allocate(azimuth(nte))
                azimuth(:) = 0.
              end if

              ! Read in the azimuth)
              ind_al=burp_find_element(block_in, element=BUFR_NEAZ)
              if (ind_al <= 0 ) cycle

              do k = 1, nte
                azimuth(k)= BURP_Get_Rval(Block_in, &
                                            nele_ind = ind_al, &
                                            nval_ind = 1, &
                                            nt_ind   = k, &
                                            iostat   = error)
                call handle_error(error, "brpr_readBurp: BURP_Get_Rval BUFR_NEAZ")
              end do
            end if ! AL
            !
            !====================================================================

          end if


          ! flag block (btyp = 0111 100011X XXXX)
          if ( btyp10 == btyp10flg ) then    ! non-sfc type flag block

            flag_passage4 = 1
            allocate(qcflag( NELE,nvale,nte))
            allocate(qcflags(NELE,nvale) )
            QCFLAG (:,:,:)  = 0
            QCFLAGS(:,:)    = 0

            allocate(hiresTimeFlag(nvale,nte))
            allocate(hiresLatFlag(nvale,nte))
            hiresTimeFlag(:,:) = 0
            hiresLatFlag(:,:) = 0

            do IL = 1, NELE

              iele=LISTE_ELE(IL)

              IND_QCFLAG  = BURP_Find_Element(Block_in, ELEMENT=200000+iele)
              if (IND_QCFLAG <= 0 ) cycle

              if (UA_FLAG_HIGH_PRECISION_TT_ES) then
                if (HIPCS .and. iele == 12001) then
                  IND_QCFLAG = BURP_Find_Element(Block_in, ELEMENT=212101)
                  if (IND_QCFLAG == -1) call handle_error(IND_QCFLAG, "brpr_readBurp: cannot find element 212101 in Block_in")
                end if
                if (HIPCS .and. iele == 12192) then
                  IND_QCFLAG = BURP_Find_Element(Block_in, ELEMENT=212239)
                  if (IND_QCFLAG == -1) call handle_error(IND_QCFLAG, "brpr_readBurp: cannot find element 212239 in Block_in")
                end if
              end if

              do k = 1, nte
                do j = 1, nvale
                  QCFLAG(IL,j,k)= BURP_Get_Tblval(Block_in, &
                                  &   NELE_IND = IND_QCFLAG, &
                                  &   NVAL_IND = j, &
                                  &   NT_IND   = k, &
                                  &   IOSTAT   = error)
                  call handle_error(error, "brpr_readBurp: BURP_Get_Tblval IND_QCFLAG")
                  SUM = SUM +1
                end do
              end do

            end do

            ! read the hires time and latitude flags, needed for UA thinning procedure
            IND_QCFLAG = BURP_Find_Element(Block_in, ELEMENT=204015)
            if (IND_QCFLAG > 0) then
              do k = 1, nte
                do j = 1, nvale
                  hiresTimeFlag(j,k)= BURP_Get_Tblval(Block_in, &
                                  &   NELE_IND = IND_QCFLAG, &
                                  &   NVAL_IND = j, &
                                  &   NT_IND   = k, &
                                  &   IOSTAT   = error)
                  call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 204015")
                  SUM = SUM +1
                end do
              end do
            end if
            IND_QCFLAG = BURP_Find_Element(Block_in, ELEMENT=205001)
            if (IND_QCFLAG > 0) then
              do k = 1, nte
                do j = 1, nvale
                  hiresLatFlag(j,k)= BURP_Get_Tblval(Block_in, &
                                  &   NELE_IND = IND_QCFLAG, &
                                  &   NVAL_IND = j, &
                                  &   NT_IND   = k, &
                                  &   IOSTAT   = error)
                  call handle_error(error, "brpr_readBurp: BURP_Get_Tblval 205001")
                  SUM = SUM +1
                end do
              end do
            end if

          end if

          ! info block (btyp = 0001 100000X XXXX) 
          BTYP10inf = 96

          if ( (btyp10 == btyp10inf) .or. (btyp10 - btyp10inf == 1) ) then

            allocate( RINFO(NELE_INFO,nte))
            allocate(TRINFO(NELE_INFO))

            flag_passage2 = 1

            do kl=1,NELE_INFO
              info_elepos = BURP_Find_Element(Block_in, &
                                 & ELEMENT  = LISTE_INFO(kl) )
              if (  info_elepos >= 0 )then
                
                do k =1 , nte
                  RINFO(kl,k)= BURP_Get_rval(Block_in, &
                              &   NELE_IND = info_elepos, &
                              &   NVAL_IND = 1, &
                              &   NT_IND   = k, &
                              &   IOSTAT   = error)
                  call handle_error(error, "brpr_readBurp: BURP_Get_rval info_elepos")
                  if  (RINFO(kl,k) == MPC_missingValue_R4)  THEN    
                    infot= BURP_Get_tblval(Block_in, &
                            &   NELE_IND = info_elepos, &
                            &   NVAL_IND = 1, &
                            &   NT_IND   = k, &
                            &   IOSTAT   = error)
                    call handle_error(error, "brpr_readBurp: BURP_Get_tblval info_elepos")
                    if (infot /= -1) RINFO(kl,k) =real(infot)
                  END IF

                end do

              else
                RINFO(kl,1:nte)=MPC_missingValue_R4
              end if

            end do

          end if


          ! descriptor block (btyp = 0010 100000X XXXX) 
          BTYP10des = 160

          if ( BTYP10 == BTYP10des ) then

            flag_passage1 = 1

            allocate(GLBFLAG(nte))
            allocate(    lat(nte))
            allocate(    lon(nte))
            allocate(   date(nte))
            allocate(   hhmm(nte))

            ! DATE  004208  HHMM 004197  STATUS 055200  LAT 005002  LON 006002  DELAY 004195
            ind055200 = BURP_Find_Element(Block_in, ELEMENT=055200)
            if (ind055200 == -1) call handle_error(ind055200, "brpr_readBurp: cannot find element 55200 in Block_in")
            ind5002   = BURP_Find_Element(Block_in, ELEMENT=5002)
            if (ind5002 == -1) call handle_error(ind5002, "brpr_readBurp: cannot find element 5002 in Block_in")
            ind6002   = BURP_Find_Element(Block_in, ELEMENT=6002)
            if (ind6002 == -1) call handle_error(ind6002, "brpr_readBurp: cannot find element 6002 in Block_in")
            ind4208   = BURP_Find_Element(Block_in, ELEMENT=4208)
            if (ind4208 == -1) call handle_error(ind4208, "brpr_readBurp: cannot find element 4208 in Block_in")
            ind4197   = BURP_Find_Element(Block_in, ELEMENT=4197)
            if (ind4197 == -1) call handle_error(ind4197, "brpr_readBurp: cannot find element 4197 in Block_in")

            do k = 1, nte
              LAT(k) =  BURP_Get_Tblval(Block_in, &
                         &   NELE_IND = ind5002, &
                         &   NVAL_IND = 1, &
                         &   NT_IND   = k)
              LON(k) =  BURP_Get_Tblval(Block_in, &
                         &   NELE_IND = ind6002, &
                         &   NVAL_IND = 1, &
                         &   NT_IND   = k)
              HHMM(k) =  BURP_Get_Tblval(Block_in, &
                          &   NELE_IND = ind4197, &
                          &   NVAL_IND = 1, &
                          &   NT_IND   = k)
              DATE(k) =  BURP_Get_Tblval(Block_in, &
                          &   NELE_IND = ind4208, &
                          &   NVAL_IND = 1, &
                          &   NT_IND   = k)
              GLBFLAG(k) =  BURP_Get_Tblval(Block_in, &
                             &   NELE_IND = ind055200, &
                             &   NVAL_IND = 1, &
                             &   NT_IND   = k)
            end do

          end if

          !==================== IASI  SPECIAL BLOCK==================
          if ( BTYP == 9217 .and.  IDTYP == 186   ) then
            NCLASSAVHRR=obs_getNclassAvhrr()
            NCHANAVHRR=obs_getNchanAvhrr()
            if (.not. allocated(CFRAC) ) allocate( CFRAC(NCLASSAVHRR,nte) )
            if (.not. allocated(RADMOY)) allocate(RADMOY(NCLASSAVHRR,NCHANAVHRR,nte))
            if (.not. allocated(radstd)) allocate(radstd(NCLASSAVHRR,NCHANAVHRR,nte))

            RADMOY(:,:,:)=MPC_missingValue_R4
            RADSTD(:,:,:)=MPC_missingValue_R4
            CFRAC(:,:)=MPC_missingValue_R4

            IASIQUAL: DO k = 1, nte
              iclass=1

              NVALS :do j=1,nvale
                DO  il=1,nbele
                  iele=BURP_Get_Element(Block_in,INDEX =il)
                  SELECT CASE(iele)
                    CASE(25085)
                      CFRAC(iclass,k)= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      call handle_error(error, "brpr_readBurp: BURP_Get_RVAL 25085")
                    CASE(5042)
                      !ICHAN= BURP_Get_TBLVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      ICHAN= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      call handle_error(error, "brpr_readBurp: BURP_Get_RVAL 5042")
                    CASE(25142)
                      !INORM= BURP_Get_TBLVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      INORM= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      call handle_error(error, "brpr_readBurp: BURP_Get_RVAL 25142")
                    CASE(14047)
                      RAD_MOY=BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      call handle_error(error, "brpr_readBurp: BURP_Get_RVAL 14047")
                      RADMOY(iclass,ICHAN,k)=RAD_MOY * 10.d0**(-1.d0 * INORM ) * 100000.d0
                    CASE(14048)
                      RAD_STD=BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      call handle_error(error, "brpr_readBurp: BURP_Get_RVAL 14048")
                      RADSTD(iclass,ICHAN,k)= RAD_STD * 10.d0**(-1.d0 * INORM ) * 100000.d0
                      IF (ICHAN==NCHANAVHRR) iclass=iclass+1
                      IF ( iclass==(NCLASSAVHRR+1) ) EXIT NVALS
                  END SELECT
                END DO
              end do NVALS

            END DO IASIQUAL

          end if

        end do BLOCKS1

        ! Loop over observations/locations (nte > 1 for families with grouped data like AI and TO)
        ! Fill obsSpaceData HEADER and BODY
        do k = 1, nte

          IF (  allocated(lat) ) then

            XLON    = (lon(k)*1.-18000.)*.01
            XLAT    = (lat(k)*1.- 9000.)*.01
            IF ( xlon  < 0. ) xlon  = 360. + xlon
            XLON    = XLON*MPC_RADIANS_PER_DEGREE_R8
            XLAT    = XLAT*MPC_RADIANS_PER_DEGREE_R8
            YMD_DATE=date(k)
            HM      =hhmm(k)
            STATUS  =GLBFLAG(K)
            RELEV   =REAL(ELEV) - 400.
          ELSE
            XLON    =.01*LONG
            XLAT    =LATI*.01 -90.
            XLON    = XLON*MPC_RADIANS_PER_DEGREE_R8
            XLAT    = XLAT*MPC_RADIANS_PER_DEGREE_R8
            YMD_DATE=date_h
            HM      =hhmm_h
            RELEV   =REAL(ELEV,pre_obsReal) - 400.
          END IF

          if (allocated(RINFO))      TRINFO(1:NELE_INFO)     =RINFO       (1:NELE_INFO,k)


          if(ENFORCE_CLASSIC_SONDES) hires=.false.

          IF  (HIRES ) THEN  ! LAT LON TIME IN DATA BLOCK (e.g. UA family)

            if (allocated(EMIS)) SURF_EMIS(1:NVAL)        = EMIS(1:NVAL,k)
            if (allocated(BCOR)) BiasCorrection(1,1:NVAL) = BCOR(1:NVAL,k)

            NDATA   =0
            NDATA_SF=0
            IF ( allocated(obsvalue_sfc) ) THEN
              OBSERV_SFC(1:NELE_SFC,1:1) = obsvalue_sfc(1:NELE_SFC,1:1,k)
              BCOR_SFC(1:NELE_SFC,1:1) = BiasCorrection_sfc(1:NELE_SFC,1:1,k)
              QCFLAGS_sfc(1:NELE_SFC,1:1) = qcflag_sfc(1:NELE_SFC,1:1,k)
              IF ( HIRES_SFC) THEN
                XLAT=HLAT_SFC(k);XLON=HLON_SFC(k);XTIME=HTIME_SFC(k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON
                ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                call INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)
                time2=time_sonde/10000
                YMD_DATE_SFC=date2
                HM_SFC=time2
              END IF

              ! dataQcFlagLev does not exist for surface data
              dataQcFlagLev_sfc(:) = MPC_missingValue_INT
              dataCloudFracLev_sfc(:) = MPC_missingValue_INT

              NDATA_SF= WRITE_BODY(obsdat,UNI_FAMILYTYPE,RELEV,vcoord_sfc,vcoord_type, &
                                   OBSERV_SFC,qcflags_sfc,NELE_SFC,1,LISTE_ELE_SFC, &
                                   dataQcFlagLEV_sfc,dataCloudFracLev_sfc, ROLAT,ROLON,& 
                                   BiasCorrection_opt=BCOR_SFC)

              IF ( NDATA_SF > 0) THEN
                call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE_SFC,HM_SFC,idtyp,STATUS,RELEV,FILENUMB)
                OBSN=obs_numHeader(obsdat)
                if (obs_columnActive_IH(obsdat,obs_prfl)) call obs_headSet_i(obsdat,OBS_PRFL,OBSN,kk)
                if (obs_columnActive_IH(obsdat,obs_hdd)) call obs_headSet_i(obsdat,OBS_HDD,OBSN,date_h)
                if (obs_columnActive_IH(obsdat,obs_hdt)) call obs_headSet_i(obsdat,OBS_HDT,OBSN,hhmm_h)
                call obs_setFamily(obsdat,trim(FAMILYTYPE),  OBSN )
                call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA_SF)
                IF (OBSN > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                END IF

                ! write info block to header (same values for all headers associated with report)
                if (allocated(TRINFO))  then
                  call writeInfo(obsdat,familytype, TRINFO,LISTE_INFO,NELE_INFO  )
                else
                  call setInfoToMissing(obsdat)
                end if

              END IF

            END IF

            IF ( allocated(obsvalue) ) THEN

              if (.not. allocated(dataQcFlag2)) then
                if (.not. allocated(dataQcFlagLEV)) allocate(dataQcFlagLEV(1))
                dataQcFlagLEV(:) = MPC_missingValue_INT
              end if

              if (.not. allocated(dataCloudFrac)) then
                if (.not. allocated(dataCloudFracLEV)) allocate(dataCloudFracLEV(1))
                dataCloudFracLEV(:) = MPC_missingValue_INT
              end if

              ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
              do  JJ =1,nval
                OBSERV(1:NELE,1:1) = obsvalue(1:NELE,jj:jj,k)
                if (allocated(BCOR2))  BiasCorrection2(1:NELE,1:1) = BCOR2(1:NELE,jj:jj,k)
                if (allocated(qcflag)) QCFLAGS(1:NELE,1:1) = qcflag(1:NELE,jj:jj,k)
                if (allocated(dataQcFlag2)) dataQcFlagLEV(1:NVAL) = dataQcFlag2(1:NVAL,k)
                if (allocated(dataCloudFrac)) dataCloudFracLEV(1:NVAL) = dataCloudFrac(1:NVAL,k)

                XLAT=HLAT(jj,k);XLON=HLON(jj,k);XTIME=HTIME(jj,k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON

                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                call INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)

                time2=time_sonde/10000

                VCORD(1)=VCOORD(jj,k)
                if (allocated(BCOR)) then
                  NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type, &
                                    OBSERV,qcflags,NELE,1,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV, &
                                    ROLAT,ROLON, SURF_EMIS_opt=SURF_EMIS, BiasCorrection_opt=BiasCorrection)
                elseif (allocated(BCOR2)) then
                  NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type, &
                                    OBSERV,qcflags,NELE,1,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV, &
                                    ROLAT,ROLON, SURF_EMIS_opt=SURF_EMIS, BiasCorrection_opt=BiasCorrection2)
                else
                  NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type, &
                                    OBSERV,qcflags,NELE,1,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV, ROLAT,ROLON, &
                                    SURF_EMIS_opt=SURF_EMIS)
                end if

                IF (NDATA > 0) THEN
                  if ( phasePresent ) then
                    call WRITE_HEADER(obsdat,STNID,XLAT,XLON,date2,time2,idtyp,STATUS,RELEV,FILENUMB,phase(jj,k))
                  else
                    call WRITE_HEADER(obsdat,STNID,XLAT,XLON,date2,time2,idtyp,STATUS,RELEV,FILENUMB)
                  end if
!==================================================================================
!
! Ajoute qivals dans les argument de WRITE_QI

                  if (TRIM(FAMILYTYPE) == 'SW' .and. READ_QI_GA_MT_SW) call WRITE_QI(obsdat,QI1VAL(k),QI2VAL(k),MTVAL(k),LSVAL(k),HAVAL(k),GAVAL(k))

                  if(trim(familytype) == 'AL')call write_al(obsdat, azimuth(k))

                  OBSN=obs_numHeader(obsdat)
                  if (obs_columnActive_IH(obsdat,obs_prfl)) call obs_headSet_i(obsdat,OBS_PRFL,OBSN,kk)
                  if (obs_columnActive_IH(obsdat,obs_hdd))  call obs_headSet_i(obsdat,OBS_HDD,OBSN,date_h)
                  if (obs_columnActive_IH(obsdat,obs_hdt))  call obs_headSet_i(obsdat,OBS_HDT,OBSN,hhmm_h)
                  if (obs_columnActive_IH(obsdat,obs_tflg)) call obs_headSet_i(obsdat,OBS_TFLG,OBSN,hiresTimeFlag(jj,k))
                  if (obs_columnActive_IH(obsdat,obs_lflg)) call obs_headSet_i(obsdat,OBS_LFLG,OBSN,hiresLatFlag(jj,k))
                  call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                  call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA)
                  IF (OBSN > 1 ) THEN
                    LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                    call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                  ELSE
                    call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                  END IF

                  ! write info block to header (same values for all headers associated with report)
                  if (allocated(TRINFO))  then
                    call writeInfo(obsdat,familytype, TRINFO,LISTE_INFO,NELE_INFO  )
                  else
                    call setInfoToMissing(obsdat)
                  end if

                END IF

              end do

            END IF

          ELSE  ! not HIRES

            if (allocated(EMIS)) SURF_EMIS(1:NVAL)        = EMIS(1:NVAL,k)
            if (allocated(BCOR)) BiasCorrection(1,1:NVAL) = BCOR(1:NVAL,k)
            NDATA   =0
            NDATA_SF=0
            IF ( allocated(obsvalue_sfc) ) THEN
              IF ( HIRES_SFC) THEN
                XLAT=HLAT_SFC(k);XLON=HLON_SFC(k);XTIME=HTIME_SFC(k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON
                ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                call INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)
                time2=time_sonde/10000
                YMD_DATE=date2
                HM=time2
              END IF

              OBSERV_SFC (1:NELE_SFC,1:1) = obsvalue_sfc(1:NELE_SFC,1:1,k)
              QCFLAGS_sfc(1:NELE_SFC,1:1) = qcflag_sfc  (1:NELE_SFC,1:1,k)
              BCOR_SFC(1:NELE_SFC,1:1) = BiasCorrection_sfc (1:NELE_SFC,1:1,k)

              ! dataQcFlagLev does not exist for surface data
              dataQcFlagLev_sfc(:) = MPC_missingValue_INT
              dataCloudFracLev_sfc(:) = MPC_missingValue_INT

              NDATA_SF= WRITE_BODY(obsdat,UNI_FAMILYTYPE,RELEV,vcoord_sfc,vcoord_type, &
                                   OBSERV_sfc,qcflags_sfc,NELE_SFC,1,LISTE_ELE_SFC, &
                                   dataQcFlagLEV_sfc,dataCloudFracLev_sfc, ROLAT,ROLON, &
                                   BiasCorrection_opt=BCOR_SFC)

              IF ( NDATA_SF > 0) THEN
                call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE,HM,idtyp,STATUS,RELEV,FILENUMB)
                OBSN=obs_numHeader(obsdat) 
                if (obs_columnActive_IH(obsdat,obs_prfl)) call obs_headSet_i(obsdat,OBS_PRFL,OBSN,kk)
                if (obs_columnActive_IH(obsdat,obs_hdd)) call obs_headSet_i(obsdat,OBS_HDD,OBSN,date_h)
                if (obs_columnActive_IH(obsdat,obs_hdt)) call obs_headSet_i(obsdat,OBS_HDT,OBSN,hhmm_h)
                call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                call obs_headSet_i(obsdat,OBS_NLV ,OBSN,NDATA_SF)
                IF (OBSN  > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                END IF
              END IF
            END IF
          
            IF ( allocated(obsvalue) ) THEN
              OBSERV(1:NELE,1:NVAL) = obsvalue(1:NELE,1:NVAL,k)
              if (allocated(BCOR2))  BiasCorrection2(1:NELE,1:NVAL) = BCOR2(1:NELE,1:NVAL,k)
              if (allocated(qcflag)) QCFLAGS(1:NELE,1:NVAL) = qcflag(1:NELE,1:NVAL,k)
              if (allocated(dataQcFlag2)) then
                dataQcFlagLEV(1:NVAL) = dataQcFlag2(1:NVAL,k)
              else
                if (.not. allocated(dataQcFlagLev)) allocate(dataQcFlagLEV(1))
                dataQcFlagLEV(:) = MPC_missingValue_INT
              end if
              if (allocated(dataCloudFrac)) then
                dataCloudFracLEV(1:NVAL) = dataCloudFrac(1:NVAL,k)
              else
                if (.not. allocated(dataCloudFracLev)) allocate(dataCloudFracLEV(1))
                dataCloudFracLEV(:) = MPC_missingValue_INT
              end if

              VCORD(1:NVAL) = VCOORD(1:NVAL,k)

              !CASES DEPENDING ON WETHER ON NOT WE HAVE MW DATA
              if (allocated(BCOR)) then
                NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type,OBSERV, &
                                  qcflags,NELE,NVAL,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV,ROLAT,ROLON, &
                                  SURF_EMIS_opt=SURF_EMIS,BiasCorrection_opt=BiasCorrection)
              elseif (allocated(BCOR2)) then
                NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type,OBSERV, &
                                  qcflags,NELE,NVAL,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV, ROLAT,ROLON, &
                                  SURF_EMIS_opt=SURF_EMIS,BiasCorrection_opt=BiasCorrection2)
              else
                NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type,OBSERV, &
                                  qcflags,NELE,NVAL,LISTE_ELE,dataQcFlagLEV,dataCloudFracLEV, &
                                  ROLAT,ROLON, SURF_EMIS_opt=SURF_EMIS)
              end if
              
              IF (NDATA > 0) THEN

                IF (NDATA_SF == 0) THEN
                  if ( phasePresent ) then
                    call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE,HM,idtyp,STATUS,RELEV,FILENUMB,phase(1,k))
                  else
                    call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE,HM,idtyp,STATUS,RELEV,FILENUMB)
                  end if
                 
!==================================================================================
!
! Ajoute qivals dans les argument de WRITE_QI

                  if (TRIM(FAMILYTYPE) == 'SW') call WRITE_QI(obsdat,QI1VAL(k),QI2VAL(k),MTVAL(k),LSVAL(k),HAVAL(k),GAVAL(k))

                  if(trim(familytype) == 'AL')call write_al(obsdat, azimuth(k))

                  OBSN=obs_numHeader(obsdat) 
                  call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                END IF
                OBSN=obs_numHeader(obsdat) 
                call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA+NDATA_SF)
                IF (OBSN   > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                END IF

              END IF

            END IF

            !============ IASI =====================================
            if ( allocated(RADMOY) .and. NDATA > 0 ) then
              OBSN=obs_numHeader(obsdat)

              iclass=1
              do iobs=OBS_CF1,OBS_CF7
                if(obs_columnActive_RH(obsdat,iobs)) then
                  call obs_headSet_r(obsdat,iobs,OBSN,CFRAC(iclass,k))
                  iclass=iclass+1
                end if
              end do

              iclass=1
              ichan=1
              do iobs=OBS_M1C1,OBS_M7C6
                if(obs_columnActive_RH(obsdat,iobs)) then
                  call obs_headSet_r(obsdat,iobs,OBSN,RADMOY(iclass,ichan,k))
                  ichan=ichan+1
                  if (ichan>obs_getNchanAvhrr()) then
                    ichan=1
                    iclass=iclass+1
                  end if
                end if
              end do

              iclass=1
              ichan=1
              do iobs=OBS_S1C1,OBS_S7C6
                if(obs_columnActive_RH(obsdat,iobs)) then
                  call obs_headSet_r(obsdat,iobs,OBSN,radstd(iclass,ichan,k))
                  ichan=ichan+1
                  if (ichan>obs_getNchanAvhrr()) then
                    ichan=1
                    iclass=iclass+1
                  end if
                end if
              end do

            end if
            !============ IASI =====================================

          END IF

          if (allocated(TRINFO))  then
            IF ( NDATA > 0.or.NDATA_SF > 0 ) then
              call writeInfo(obsdat,familytype, TRINFO,LISTE_INFO,NELE_INFO  )
            END IF
          end if

        end do

        !---------UPPER AIR---------------------------
        if ( allocated(obsvalue) ) then
          deallocate ( obsvalue,VCOORD,VCORD,observ)
        end if
        if ( allocated(qcflag) ) then
          deallocate (qcflag,qcflags)
        end if
        if ( allocated(hiresTimeFlag) ) then
          deallocate (hiresTimeFlag)
        end if
        if ( allocated(hiresLatFlag) ) then
          deallocate (hiresLatFlag)
        end if
        if ( allocated(qi1val) ) then
          deallocate (qi1val)
        end if
        if ( allocated(qi2val) ) then
          deallocate (qi2val)
        end if
        if ( allocated(mtval) ) then
          deallocate (mtval)
        end if
        if ( allocated(lsval) ) then
          deallocate (lsval)
        end if
        if ( allocated(haval) ) then
          deallocate (haval)
        end if
        if ( allocated(gaval) ) then
          deallocate (gaval)
        end if
        if ( allocated(EMIS) ) then
          deallocate (EMIS,SURF_EMIS)
        end if
        if ( allocated(dataQcFlag2) ) then
          deallocate (dataQcFlag2)
        end if
        if ( allocated(dataQcFlagLEV) ) then
          deallocate (dataQcFlagLEV)
        end if
        if ( allocated(dataCloudFrac) ) then
          deallocate (dataCloudFrac)
        end if
        if ( allocated(dataCloudFracLEV) ) then
          deallocate (dataCloudFracLEV)
        end if
        if (allocated(azimuth)) then
          deallocate(azimuth)
        end if
        if ( allocated(BCOR) ) then
          deallocate (BCOR,BiasCorrection)
        end if
        if ( allocated(BCOR2) ) then
          deallocate (BCOR2,BiasCorrection2)
        end if        

        !---------SURFACE-----------------------------
        if ( allocated(obsvalue_sfc) ) then
          DEallocate(obsvalue_sfc,vcoord_sfc,OBSERV_SFC,BiasCorrection_sfc,BCOR_SFC)
        end if

        if ( allocated(qcflag_sfc) ) then
          DEallocate(  qcflag_sfc, qcflags_SFC)
        end if
        !--------SURFACE------------------------------
        
        if (  allocated(lat) ) then
          deallocate (lat,lon,date,hhmm,glbflag)
        end if
        if (  allocated(hlat) ) then
          deallocate (hlat,hlon,htime)
        end if
        if (  allocated(hlat_sfc) ) then
          deallocate (hlat_sfc,hlon_sfc,htime_sfc)
        end if
        if ( allocated(rinfo) ) then
          deallocate (rinfo,trinfo)
        end if
        if ( allocated(RADMOY) ) then
          deallocate (RADMOY,CFRAC,radstd)
        end if
        if ( allocated(phase) ) then
          deallocate (phase)
        end if


      end do REPORTS

    end if

    Deallocate(address)

    if ( flag_passage1 == 1 ) then
      write(*,*)
      write(*,*) ' descriptor block for grouped data present '
    end if
    if ( flag_passage2 == 1 ) then
      write(*,*)
      write(*,*) '- info block Present '
    end if
    if ( flag_passage3 == 0 ) then
      write(*,*)
      write(*,*) 'ERROR - observation block not seen ? Verify btyp'
    end if
    if ( flag_passage4 == 0 ) then
      write(*,*)
      write(*,*) 'ERROR - flag block not seen ? Verify btyp'
    end if

    call cleanup()

    numHeader = obs_numHeader(obsdat)
    write(*,*)' file   Nobs SUM = ',trim(brp_file),numHeader,SUM

  contains

    !--------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(4)

      errors(:) = 0
      if (allocated(address)) deallocate(address, stat=errors(1))
      call BURP_Free(File_in, iostat=errors(2))
      call BURP_Free(Rpt_in, iostat=errors(3))
      call BURP_Free(Block_in, iostat=errors(4))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_readBurp: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !--------- handle_error -----
    subroutine handle_error(icode,errormessage)
      implicit none
      character (len=*) :: errormessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errormessage))
      end if
    end subroutine handle_error

  end subroutine brpr_readBurp


  FUNCTION WRITE_BODY(obsdat,FAMTYP, ELEV,VERTCOORD,VCOORD_TYPE, &
                      obsvalue,qcflag,NELE,NVAL,LISTE_ELE,dataQcFlagLEV, &
                      dataCloudFracLEV, ROLAT,ROLON, SURF_EMIS_opt,BiasCorrection_opt)

    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat
    INTEGER ::  WRITE_BODY,VCOORD_TYPE
    REAL   , allocatable          ::  OBSVALUE(:,:)
    INTEGER, allocatable          ::  QCFLAG(:,:)
    REAL   , allocatable          ::  VERTCOORD(:)
    REAL                          ::  ROLAT(:), ROLON(:)
    REAL   , allocatable,optional ::  SURF_EMIS_opt(:)
    REAL   , allocatable,optional ::  BiasCorrection_opt(:,:)
    integer, intent(in)           ::  dataQcFlagLEV(:)
    integer, intent(in)           ::  dataCloudFracLEV(:)

    ! Locals:
    CHARACTER(len=2)  :: FAMTYP
    REAL              :: ELEVFACT,VCOORD,ZEMFACT
    INTEGER           :: NELE,NVAL,VCO,NONELEV
    integer           :: LISTE_ELE(:),NOBS,VARNO,IL,J,COUNT,NLV
    INTEGER           :: IFLAG,IFLAG2,BITSflagoff,BITSflagon,cloudFrac
    REAL(pre_obsReal) :: MISG,OBSV,ELEV,ELEV_R,REMIS,BCOR,rolat1,rolon1
    LOGICAL           :: L_EMISS,L_BCOR,L_dataQcFlag2, L_dataCloudFrac
    
    L_EMISS = present( SURF_EMIS_opt )
    L_BCOR  = present( BiasCorrection_opt )
    L_dataQcFlag2 = any(dataQcFlagLEV(:) /= MPC_missingValue_INT)
    L_dataCloudFrac = any(dataCloudFracLEV(:) /= MPC_missingValue_INT)

    NONELEV  =-1
    MISG=MPC_missingValue_R4
    ZEMFACT=0.01

    REMIS = MISG

    BITSflagoff=0
    DO J = 1, Bnbitsoff
      BITSflagoff = IBSET ( BITSflagoff, 13-BBITOFF(J) )
    END DO

    BITSflagon=0
    DO J = 1, Bnbitson
       BITSflagon = IBSET ( BITSflagon,  13-BBITON(J)  )
    END DO
    !write(*,*) ' write body BITSFLAGON= ',BITSFLAGON, ' BITSFLAGOFF= ',BITSFLAGOFF

    NOBS =obs_numHeader(obsdat) +1
    COUNT=obs_numBody (obsdat)

    NLV=0

! Special test for GB-GPS 
! Reports with missing ZTD cause problems with output of subsequent reports (missing OMP, OER, FLAG bits)
! Since July 2019, reports with missing ZTD are removed in the GB-GPS dbase files.
    if ( trim(FAMTYP) == trim('GP') )  then
      DO il = 1, NELE
         varno = LISTE_ELE(il)
         DO j = 1, NVAL
            obsv = obsvalue(il,j)
            if ( varno == 15031 .and.  obsv == MISG  ) then
               print * , 'write body: report rejected no ZTD'
               WRITE_BODY = NLV
               return
            endif
         END DO
      END DO
    endif
! Special test for GB-GPS

    if (  trim(FAMTYP) == trim('PR') .OR. trim(FAMTYP) == trim('SF') .OR. trim(FAMTYP) == trim('GP')  ) then
      ELEVFACT=1.
    else
      ELEVFACT=0.
    end if
    if (  trim(FAMTYP) == trim('TO') ) then
      !ELEV=0.
    END IF

    SELECT CASE(FAMTYP)
      CASE ( 'UA' , 'SW' , 'AI')
        VCO=2 ! PRESSURE COORD
      CASE ( 'SF' , 'SC', 'PR', 'RO', 'GP', 'AL' )
        VCO=1 ! HEIGHT COORD
      CASE ( 'TO'  )
        VCO=3 ! CHANNEL NUMBER
      CASE ( 'CH'  ) ! CONSTITUENTS (consider different possibilities)
        IF (VCOORD_TYPE == 7006 .OR. VCOORD_TYPE == 7007) THEN
           VCO=1
           IF (VCOORD_TYPE == 7006) ELEVFACT=1.0
        ELSE IF (VCOORD_TYPE == 7004 .OR. VCOORD_TYPE == 7204) THEN
           VCO=2
        ELSE IF (VCOORD_TYPE == 2150) THEN
           VCO=3
        ELSE            
           ! Vertical coordinate not provided or not recognized.
           IF (NVAL == 1) THEN
              VCO=5 ! Initializes as surface point measurement
              DO il = 1, NELE
                 if ( obsvalue(il,1) /= MPC_missingValue_R4) then
                    if (liste_ele(il) == 15198.OR.liste_ele(il) == 15009.OR. &
                        liste_ele(il) == 15020.OR.liste_ele(il) == 15021.OR. &
                        liste_ele(il) == 15024.OR.liste_ele(il) == 15200.OR. &
                        liste_ele(il) == 15001.OR.liste_ele(il) == 15005) then
                        
                        VCO=4  ! Assumes this is a total column measurement
                        exit
                    end if
                 end if
              END DO
           ELSE   
              call utl_abort('write_body: Invalid BURP vertical coordinate type.')
           END IF
        END IF
    END SELECT

    !-------------------SPECIAL CASES--------------
    DO il = 1, NELE
      varno=LISTE_ELE(il)
      DO j = 1, NVAL
        VCOORD=VERTCOORD(j)
        OBSV  =   obsvalue(il,j)
        if( L_EMISS )  then
          if( SURF_EMIS_opt(j) /= MISG)  then
            REMIS =  SURF_EMIS_opt(j)*ZEMFACT
          else
            REMIS = MISG
          end if
        end if
        if ( L_BCOR )  then
          BCOR  =  BiasCorrection_opt(il,j)
        end if
        if ( L_dataQcFlag2 )  then
          IFLAG2  =  dataQcFlagLEV(j)
        end if
        if ( L_dataCloudFrac )  then
          cloudFrac  =  dataCloudFracLEV(j)
        end if
        IFLAG = INT(qCflag(il,j))

        if(iand(iflag,BITSflagoff) /= 0) cycle

        if ( obsv /= MPC_missingValue_R4 .and. VCOORD /= MPC_missingValue_R4   ) then
          count = count  + 1
          NLV= NLV +1
          IFLAG = IBCLR(IFLAG,12)

          call obs_bodySet_r(obsdat,OBS_VAR,count,OBSV)
          call obs_bodySet_i(obsdat,OBS_VNM,count,VARNO)
          call obs_bodySet_i(obsdat,OBS_VCO,count,VCO)
          ELEV_R=VCOORD + ELEV*ELEVFACT
          call obs_bodySet_r(obsdat,OBS_PPP,count, ELEV_R)
          call obs_bodySet_i(obsdat,OBS_FLG,count,IFLAG)
          if ( FAMTYP == 'RO' ) then
            rolat1 = rolat(j)*MPC_RADIANS_PER_DEGREE_R8
            rolon1 = rolon(j)*MPC_RADIANS_PER_DEGREE_R8
            call obs_bodySet_r(obsdat,OBS_LATD,count,rolat1)
            call obs_bodySet_r(obsdat,OBS_LOND,count,rolon1)
          else
            call obs_bodySet_r(obsdat,OBS_LATD,count,obs_missingValue_R)
            call obs_bodySet_r(obsdat,OBS_LOND,count,obs_missingValue_R)
          end if
          if ( L_BCOR .and. obs_columnActive_RB(obsdat,OBS_BCOR) ) then
            call obs_bodySet_r(obsdat,OBS_BCOR,count,BCOR)
          end if

          if ( L_dataQcFlag2 ) then
             call obs_bodySet_i(obsdat,OBS_QCF2,count,IFLAG2)
          end if

          if ( L_dataCloudFrac ) then
             call obs_bodySet_i(obsdat,OBS_CLA,count,cloudFrac)
          end if

          if ( REMIS /= MPC_missingValue_R4 .and. FAMTYP == 'TO') THEN
            call obs_bodySet_r(obsdat,OBS_SEM,count,REMIS)
          else
            if ( FAMTYP == 'TO') then
                call obs_bodySet_r(obsdat,OBS_SEM,count,tvs_defaultEmissivity)
             else
                call obs_bodySet_r(obsdat,OBS_SEM,count,MISG)
             end if
          end if

          call obs_bodySet_i(obsdat,OBS_VCO,count,VCO)

          if (.not. filt_bufrCodeAssimilated(varno) .and. &
              .not. ovt_bufrCodeSkipped(varno)) then
            ! Add a row for the destination transform variable
            call obs_bodySet_i(obsdat,OBS_VNM,count+1,ovt_getDestinationBufrCode(varno))
            call obs_bodySet_i(obsdat,OBS_FLG,count+1,0)
            ELEV_R=VCOORD + ELEV*ELEVFACT
            call obs_bodySet_r(obsdat,OBS_PPP, count+1,ELEV_R)
            call obs_bodySet_i(obsdat,OBS_VCO, count+1,VCO)
            call obs_bodySet_r(obsdat,OBS_VAR, count+1,MISG)
            call obs_bodySet_r(obsdat,OBS_LATD,count+1,MISG)
            call obs_bodySet_r(obsdat,OBS_LOND,count+1,MISG)
            count = count + 1
            NLV = NLV + 1
            if (ovt_isWindObs(varno)) then
              ! Add an extra row for the other wind component
              call obs_bodySet_i(obsdat,OBS_VNM, count+1,ovt_getDestinationBufrCode(varno,extra_opt=.true.))
              call obs_bodySet_i(obsdat,OBS_FLG, count+1,0)
              call obs_bodySet_r(obsdat,OBS_PPP, count+1,ELEV_R)
              call obs_bodySet_i(obsdat,OBS_VCO, count+1,VCO)
              call obs_bodySet_r(obsdat,OBS_VAR, count+1,MISG)
              call obs_bodySet_r(obsdat,OBS_LATD,count+1,MISG)
              call obs_bodySet_r(obsdat,OBS_LOND,count+1,MISG)
              count = count + 1
              NLV = NLV + 1
            end if
          end if

        end if

      END DO

    END DO

    WRITE_BODY=NLV

  END FUNCTION  WRITE_BODY


  SUBROUTINE WRITE_HEADER(obsdat, STNID,LAT,LON,DATE,TIME,CODTYP,STATUS,ELEV,FILENUMB,PHASE_Opt)

    implicit none
    type (struct_obs), intent(inout) :: obsdat
    CHARACTER(LEN=9)       :: STNID

    integer     ::    DATE,TIME,CODTYP,STATUS
    integer     ::    FILENUMB
    INTEGER, optional :: phase_opt

    REAL(pre_obsReal) :: ELEV,LAT,LON

    integer     ::   NOBS

    NOBS=obs_numHeader(obsdat)  +1

    call obs_headSet_i(obsdat,OBS_ONM,nobs,nobs)
    call obs_headSet_r(obsdat,OBS_LAT,nobs,LAT)
    call obs_headSet_r(obsdat,OBS_LON,nobs,LON)
    call obs_headSet_i(obsdat,OBS_DAT,nobs,DATE)
    call obs_headSet_i(obsdat,OBS_ETM,nobs,TIME)
    call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
    call obs_headSet_i(obsdat,OBS_ST1,nobs,STATUS)
    call obs_headSet_r(obsdat,OBS_ALT,nobs,ELEV)
    !call obs_headSet_i(obsdat,OBS_IDF,nobs,FILENUMB)
    call obs_headSet_i(obsdat,OBS_OTP,nobs,FILENUMB)
    call obs_set_c(obsdat,'STID',nobs,STNID )
    if ( present(phase_opt) .and. &
         obs_columnActive_IH(obsdat,OBS_PHAS) ) then
      call obs_headSet_i(obsdat,OBS_PHAS,nobs,phase_opt)
    end if

  END SUBROUTINE  WRITE_HEADER

!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------

  SUBROUTINE WRITE_QI(obsdat, QI1value, QI2value, MTvalue, LSvalue, HAvalue, GAvalue)

    implicit none
    type (struct_obs), intent(inout) :: obsdat
    integer     ::  MTvalue, HAvalue, GAvalue, QI1value, QI2value, LSvalue
    integer     ::  NOBS

    NOBS = obs_numHeader(obsdat)

    call obs_headSet_i(obsdat,OBS_SWQ1,nobs,QI1value)
    call obs_headSet_i(obsdat,OBS_SWQ2,nobs,QI2value)
    call obs_headSet_i(obsdat,OBS_SWMT,nobs,MTvalue)
    call obs_headSet_i(obsdat,OBS_SWLS,nobs,LSvalue)
    call obs_headSet_i(obsdat,OBS_SWGA,nobs,GAvalue)
    call obs_headSet_i(obsdat,OBS_SWHA,nobs,HAvalue)

  END SUBROUTINE  WRITE_QI

!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------

  subroutine write_al(obsdat, azimuth)
    implicit none
    type (struct_obs), intent(inout) :: obsdat
    real(pre_obsReal) :: azimuth
    integer :: nobs

    nobs = obs_numHeader(obsdat)

    call obs_headSet_r(obsdat,OBS_AZA,nobs,azimuth)
  end subroutine write_al

  !--------------------------------------------------------------------------
  ! writeInfo
  !--------------------------------------------------------------------------
  subroutine writeInfo(obsdat, FAMTYP, RINFO, LISTE_INFO, NELE_INFO)
    ! :Purpose: Write values in obsSpaceData related to the info block

    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat
    real        ::   RINFO(NELE_INFO)
    CHARACTER*2 ::   FAMTYP
    integer     ::   NELE_INFO
    integer     ::   LISTE_INFO(NELE_INFO)

    ! Locals:
    REAL*4      ::   INFOV
    integer     ::   CODTYP
    integer     ::   IL,NOBS
    integer     ::   SENSOR,ORBIT,ID_SAT,INSTRUMENT,LAND_SEA,CONSTITUENT_TYPE
    integer     ::   TERRAIN_TYPE,QCFLAG1,QCFLAG2,QCFLAG3,RAINFLAG
    integer     ::   IGQISFLAGQUAL,IGQISQUALINDEXLOC,IRO_QCFLAG
    integer     ::   IFOV,ORIGIN_CENTRE,RAOBSTYPE, LAUNCHTIME
    real        ::   RIGQISFLAGQUAL,RIGQISQUALINDEXLOC,RCONSTITUENT,RQCFLAG1,RQCFLAG2,RQCFLAG3,RRAINFLAG
    real        ::   RTERRAIN_TYPE,RLAND_SEA,RID_SAT,RSENSOR,RINSTRUMENT,RRO_QCFLAG,RORIGIN_CENTRE
    real        ::   RORBIT, RIWV
    REAL(pre_obsReal) ::   RTANGENT_RADIUS,RGEOID,RSOLAR_AZIMUTH,RCLOUD_COVER,RSOLAR_ZENITH,RZENITH,RAZIMUTH
    real        ::   RFOV
    REAL(pre_obsReal) ::   cloudLiquidWaterObs, cloudLiquidWaterFG

    NOBS=obs_numHeader(obsdat)
    CODTYP=obs_headElem_i(obsdat,OBS_ITY,NOBS)
    !write(*,*)'  DEBUT WRITE_INFO NOBS CODTYP ----> ',NOBS,CODTYP,size(liste_info),size(RINFO),liste_info
    LAND_SEA   = 0
    INSTRUMENT = 0
    ID_SAT     = 0
    SENSOR     = 0
    ORIGIN_CENTRE = 0
    RAOBSTYPE  = MPC_missingValue_INT
    LAUNCHTIME = MPC_missingValue_INT
    ORBIT      = 0
    QCFLAG1    = 0
    QCFLAG2    = 0
    QCFLAG3    = 0
    RAINFLAG   = 0

    IRO_QCFLAG=MPC_missingValue_INT
    IGQISQUALINDEXLOC=0
    IGQISFLAGQUAL=0

    RTANGENT_RADIUS=real(MPC_missingValue_R8,pre_obsReal)
    RGEOID=real(MPC_missingValue_R8,pre_obsReal)
    TERRAIN_TYPE=-1
    RIWV= MPC_missingValue_R4
    RCLOUD_COVER = MPC_missingValue_R4
    CONSTITUENT_TYPE = MPC_missingValue_INT
    IFOV = MPC_missingValue_INT
    RIGQISQUALINDEXLOC = MPC_missingValue_R4
    RIGQISFLAGQUAL = MPC_missingValue_R4
    RRO_QCFLAG = MPC_missingValue_R4
    RSOLAR_AZIMUTH = real(MPC_missingValue_R8,pre_obsReal)
    RSOLAR_ZENITH = real(MPC_missingValue_R8,pre_obsReal)
    RZENITH = 90.
    RAZIMUTH = 0.
    cloudLiquidWaterObs = real(MPC_missingValue_R8,pre_obsReal)
    cloudLiquidWaterFG = real(MPC_missingValue_R8,pre_obsReal)

    do il=1,NELE_INFO
      INFOV=rinfo(il)
      SELECT CASE( liste_info(il) )
        CASE( 1007)
          RID_SAT=INFOV
          IF (RID_SAT == MPC_missingValue_R4 ) THEN 
            ID_SAT=0
          ELSE
            ID_SAT=NINT(RID_SAT)
          END IF
        CASE( 1033)
          RORIGIN_CENTRE=INFOV
          IF (RORIGIN_CENTRE == MPC_missingValue_R4 ) THEN 
            ORIGIN_CENTRE=0
          ELSE
            ORIGIN_CENTRE=NINT(RORIGIN_CENTRE)
          END IF
        CASE( 2048)
          RSENSOR = INFOV
          if (RSENSOR == MPC_missingValue_R4 ) THEN
            SENSOR = MPC_missingValue_INT
          ELSE
            SENSOR = NINT(RSENSOR)
          END IF
        CASE( 5040)
          RORBIT = INFOV
          if (RORBIT == MPC_missingValue_R4 ) THEN
            ORBIT = MPC_missingValue_INT
          ELSE
            ORBIT = NINT(RORBIT)
          END IF
        CASE( 33078)
          RQCFLAG1 = INFOV
          if (RQCFLAG1 == MPC_missingValue_R4 ) THEN
            QCFLAG1 = MPC_missingValue_INT
          ELSE
            QCFLAG1 = NINT(RQCFLAG1)
          END IF
        CASE( 33079)
          RQCFLAG2 = INFOV
          if (RQCFLAG2 == MPC_missingValue_R4 ) THEN
            QCFLAG2 = MPC_missingValue_INT
          ELSE
            QCFLAG2 = NINT(RQCFLAG2)
          END IF
        CASE( 020029)
          RRAINFLAG = INFOV
          if (RRAINFLAG == MPC_missingValue_R4 ) THEN
            RAINFLAG = MPC_missingValue_INT
          ELSE
            RAINFLAG = NINT(RRAINFLAG)
          END IF
        CASE( 33080)
          RQCFLAG3 = INFOV
          if (RQCFLAG3 == MPC_missingValue_R4 ) THEN
            QCFLAG3 = MPC_missingValue_INT
          ELSE
            QCFLAG3 = NINT(RQCFLAG3)
          END IF

        CASE( 2019)
          RINSTRUMENT = INFOV
          if (RINSTRUMENT == MPC_missingValue_R4 ) THEN
            INSTRUMENT = 0
          ELSE
            INSTRUMENT = NINT(RINSTRUMENT)
          END IF
        CASE( 5043)
          RFOV = INFOV
          if (RFOV == MPC_missingValue_R4 ) THEN
            IFOV = 0
          ELSE
            IFOV = NINT(RFOV)
          END IF
        CASE( 7024)
          RZENITH = INFOV
          if (RZENITH == MPC_missingValue_R4 ) THEN
            RZENITH = 90.
          END IF
        CASE( 7025)
          RSOLAR_ZENITH = INFOV
        CASE( 5021)
          RAZIMUTH=INFOV
          if (RAZIMUTH == MPC_missingValue_R4 ) THEN
            RAZIMUTH = 0.
          END IF
        CASE( 33060)
          RIGQISFLAGQUAL=INFOV
          if (RIGQISFLAGQUAL == MPC_missingValue_R4 )  then
            IGQISFLAGQUAL=0
          ELSE
            IGQISFLAGQUAL=NINT ( RIGQISFLAGQUAL )
          END IF
        CASE( 33062)
          RIGQISQUALINDEXLOC=INFOV
          if (RIGQISQUALINDEXLOC == MPC_missingValue_R4 )  then
            IGQISQUALINDEXLOC=0
          ELSE
            IGQISQUALINDEXLOC=NINT ( RIGQISQUALINDEXLOC )
          END IF
        CASE( 5022)
          RSOLAR_AZIMUTH=INFOV
        CASE( 8012)
          RLAND_SEA=INFOV
          if (RLAND_SEA == MPC_missingValue_R4 ) THEN
            LAND_SEA=99
          ELSE
            LAND_SEA=NINT ( RLAND_SEA )
          END IF
        CASE( 13095)
          RIWV=INFOV
          if (RIWV == MPC_missingValue_R4 ) THEN
            RIWV=0.
          END IF
        CASE( 13039)
          RTERRAIN_TYPE=INFOV
          if (RTERRAIN_TYPE == MPC_missingValue_R4 ) THEN
            TERRAIN_TYPE=-1
          ELSE
            TERRAIN_TYPE=NINT ( RTERRAIN_TYPE )
          END IF
        CASE( 20010)
          RCLOUD_COVER=INFOV
        CASE( 10035)
          RTANGENT_RADIUS=INFOV
        CASE( 10036)
          RGEOID=INFOV
        CASE( 33039)
          RRO_QCFLAG=INFOV
          if (RRO_QCFLAG == MPC_missingValue_R4 ) THEN
            IRO_QCFLAG=MPC_missingValue_INT
          ELSE
            IRO_QCFLAG=NINT ( RRO_QCFLAG )
          END IF
        CASE( 08046)          
          IF (trim(FAMTYP) == 'CH') THEN
             RCONSTITUENT=INFOV
             IF (RCONSTITUENT == MPC_missingValue_R4) THEN
                call utl_abort('writeInfo: Missing 08046 element for the CH family.')
             ELSE
                CONSTITUENT_TYPE=NINT(RCONSTITUENT)
             END IF
          END IF
        CASE(13209)
          cloudLiquidWaterObs = INFOV
        CASE(2011)
          raobsType = nint(infov)
        CASE(4197)
          launchTime = nint(infov)
      END SELECT
      if (liste_info(il) == clwFgElementId ) cloudLiquidWaterFG = INFOV
    end do

    !-------------------SPECIAL CASES--------------

    ! INSTRUMENT
    IF ( SENSOR == MPC_missingValue_INT) then
      IF ( INSTRUMENT == MPC_missingValue_INT) then
        INSTRUMENT=0
      END IF
    ELSE
      INSTRUMENT = obsu_cvt_obs_instrum(sensor)
    END IF

    ! AIRS
    IF ( INSTRUMENT == 420 ) ID_SAT = 784

    ! CrIS FSR
    if (codtyp == 202 .and. INSTRUMENT == 620) then
      INSTRUMENT = 2046 
    end if

    if (  trim(FAMTYP) == trim('GO') ) then
      call utl_abort('writeInfo: unknown familyType : ' // trim(FAMTYP))
    END IF
   

    !-------------------SPECIAL CASES--------------

    if ( obs_columnActive_IH(obsdat,OBS_TTYP)) call obs_headSet_i(obsdat,OBS_TTYP,nobs,TERRAIN_TYPE)
    if ( obs_columnActive_IH(obsdat,OBS_STYP)) call obs_headSet_i(obsdat,OBS_STYP,nobs,LAND_SEA)
    if ( obs_columnActive_IH(obsdat,OBS_ORBI)) call obs_headSet_i(obsdat,OBS_ORBI,nobs,ORBIT)
    if ( obs_columnActive_IH(obsdat,OBS_AQF1)) call obs_headSet_i(obsdat,OBS_AQF1,nobs,QCFLAG1)
    if ( obs_columnActive_IH(obsdat,OBS_AQF2)) call obs_headSet_i(obsdat,OBS_AQF2,nobs,QCFLAG2)
    if ( obs_columnActive_IH(obsdat,OBS_AQF3)) call obs_headSet_i(obsdat,OBS_AQF3,nobs,QCFLAG3)
    if ( obs_columnActive_IH(obsdat,OBS_RAIN)) call obs_headSet_i(obsdat,OBS_RAIN,nobs,RAINFLAG)
    if ( obs_columnActive_IH(obsdat,OBS_INS) ) call obs_headSet_i(obsdat,OBS_INS,nobs,INSTRUMENT  )
    if ( obs_columnActive_IH(obsdat,OBS_FOV) ) call obs_headSet_i(obsdat,OBS_FOV,nobs,IFOV )
    if ( obs_columnActive_IH(obsdat,OBS_SAT) ) call obs_headSet_i(obsdat,OBS_SAT,nobs,ID_SAT)
    if ( obs_columnActive_IH(obsdat,OBS_ORI) ) call obs_headSet_i(obsdat,OBS_ORI,nobs,ORIGIN_CENTRE)
    if ( obs_columnActive_IH(obsdat,OBS_TEC) ) call obs_headSet_i(obsdat,OBS_TEC,nobs,0)
    if ( obs_columnActive_IH(obsdat,OBS_GQF) ) call obs_headSet_i(obsdat,OBS_GQF,nobs,IGQISFLAGQUAL)
    if ( obs_columnActive_IH(obsdat,OBS_GQL) ) call obs_headSet_i(obsdat,OBS_GQL,nobs,IGQISQUALINDEXLOC)
    if ( obs_columnActive_IH(obsdat,OBS_ROQF) ) call obs_headSet_i(obsdat,OBS_ROQF,nobs,IRO_QCFLAG)
    if ( obs_columnActive_IH(obsdat,OBS_RTP) ) call obs_headSet_i(obsdat,OBS_RTP,nobs,raobsType)
    if ( obs_columnActive_IH(obsdat,OBS_LCH) ) call obs_headSet_i(obsdat,OBS_LCH,nobs,launchTime)
    if ( obs_columnActive_RH(obsdat,OBS_CLF) ) call obs_headSet_r(obsdat,OBS_CLF,nobs,RCLOUD_COVER )
    if ( obs_columnActive_RH(obsdat,OBS_SUN) ) call obs_headSet_r(obsdat,OBS_SUN,nobs,RSOLAR_ZENITH )
    if ( obs_columnActive_RH(obsdat,OBS_SAZ) ) call obs_headSet_r(obsdat,OBS_SAZ,nobs,RSOLAR_AZIMUTH )
    if ( obs_columnActive_RH(obsdat,OBS_SZA) ) call obs_headSet_r(obsdat,OBS_SZA,nobs,RZENITH )
    if ( obs_columnActive_RH(obsdat,OBS_AZA) ) call obs_headSet_r(obsdat,OBS_AZA,nobs,RAZIMUTH )
    if ( obs_columnActive_RH(obsdat,OBS_TRAD) ) call obs_headSet_r(obsdat,OBS_TRAD,nobs,RTANGENT_RADIUS)
    if ( obs_columnActive_RH(obsdat,OBS_IWV))  call obs_headSet_r(obsdat,OBS_IWV,nobs,RIWV)
    if ( obs_columnActive_RH(obsdat,OBS_GEOI) ) call obs_headSet_r(obsdat,OBS_GEOI,nobs,RGEOID)
    if (trim(FAMTYP) == trim('CH')) then
        if ( obs_columnActive_IH(obsdat,OBS_CHM) ) call obs_headSet_i(obsdat,OBS_CHM,nobs,CONSTITUENT_TYPE)
    else
        if ( obs_columnActive_IH(obsdat,OBS_CHM) ) call obs_headSet_i(obsdat,OBS_CHM,nobs,-1)
    end if
    if ( obs_columnActive_RH(obsdat,OBS_CLWO) ) call obs_headSet_r(obsdat,OBS_CLWO,nobs,cloudLiquidWaterObs)
    if ( obs_columnActive_RH(obsdat,OBS_CLWB) ) call obs_headSet_r(obsdat,OBS_CLWB,nobs,cloudLiquidWaterFG)

  END SUBROUTINE  writeInfo

  !--------------------------------------------------------------------------
  ! setInfoToMissing
  !--------------------------------------------------------------------------
  subroutine setInfoToMissing(obsdat)
    ! :Purpose: Set the obsSpaceData column related to the info block with
    !           missing values

    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nobs

    nobs = obs_numHeader(obsdat)

    if ( obs_columnActive_IH(obsdat,OBS_STYP)) call obs_headSet_i(obsdat,OBS_STYP,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_INS) ) call obs_headSet_i(obsdat,OBS_INS,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_FOV) ) call obs_headSet_i(obsdat,OBS_FOV,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_SAT) ) call obs_headSet_i(obsdat,OBS_SAT,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_ORI) ) call obs_headSet_i(obsdat,OBS_ORI,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_TEC) ) call obs_headSet_i(obsdat,OBS_TEC,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_GQF) ) call obs_headSet_i(obsdat,OBS_GQF,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_GQL) ) call obs_headSet_i(obsdat,OBS_GQL,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_ROQF) ) call obs_headSet_i(obsdat,OBS_ROQF,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_RTP) ) call obs_headSet_i(obsdat,OBS_RTP,nobs,mpc_missingValue_int)
    if ( obs_columnActive_IH(obsdat,OBS_LCH) ) call obs_headSet_i(obsdat,OBS_LCH,nobs,mpc_missingValue_int)

    if ( obs_columnActive_RH(obsdat,OBS_CLF) ) call obs_headSet_r(obsdat,OBS_CLF,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_SUN) ) call obs_headSet_r(obsdat,OBS_SUN,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_SAZ) ) call obs_headSet_r(obsdat,OBS_SAZ,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_SZA) ) call obs_headSet_r(obsdat,OBS_SZA,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_AZA) ) call obs_headSet_r(obsdat,OBS_AZA,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_TRAD) ) call obs_headSet_r(obsdat,OBS_TRAD,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_GEOI) ) call obs_headSet_r(obsdat,OBS_GEOI,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_CLWO) ) call obs_headSet_r(obsdat,OBS_CLWO,nobs,obs_missingValue_r)
    if ( obs_columnActive_RH(obsdat,OBS_CLWB) ) call obs_headSet_r(obsdat,OBS_CLWB,nobs,obs_missingValue_r)

  end subroutine  setInfoToMissing

  !--------------------------------------------------------------------------
  ! find_index
  !--------------------------------------------------------------------------
  integer  FUNCTION FIND_INDEX(LIST,ELEMENT)
    implicit none
    integer LIST(:)
    integer I,ELEMENT
    FIND_INDEX=-1
    do I=1,size (LIST)
      if (list(i) == element) THEN
        FIND_INDEX=i
        exit
      end if
    end do
    return
  END FUNCTION FIND_INDEX

  !--------------------------------------------------------------------------
  ! brpr_addCloudParametersandEmissivity
  !--------------------------------------------------------------------------
  subroutine brpr_addCloudParametersandEmissivity( obsSpaceData, fileIndex, burpFile )
    !
    ! :Purpose: Add to the input BURP file number fileIndex cloud parameters and emissivity.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData ! obsSpacedata structure
    integer, intent(in)              :: fileIndex    ! number of the burp file to update
    character (len=*), intent(in)    :: burpFile
    
    ! Locals:
    type(BURP_FILE)        :: inputFile
    type(BURP_RPT)         :: inputReport,copyReport
    type(BURP_BLOCK)       :: inputBlock
    character(len=9)       :: opt_missing
    integer                :: btyp10, btyp, bfam, error
    integer                :: btyp10des, btyp10inf, btyp10obs, btyp10flg, btyp10omp
    integer                :: nb_rpts, ref_rpt, ref_blk, count
    integer, allocatable   :: address(:), goodprof(:)
    real(8), allocatable   :: btobs(:,:)
    real(8)                :: emisfc
    integer                :: nbele,nvale,nte
    integer, allocatable   :: glbflag(:)
    integer                :: headerIndex, valIndex, tIndex, reportIndex, bodyIndex, headElem_i
    integer                :: ind008012,ind012163,ind055200,indEmis,indchan,ichn,ichnb
    integer                :: ind5021, ind7024, ind13039
    integer                :: ind14213, ind14214, ind14215, ind14216, ind14217, ind14218
    integer                :: ind14219, ind14220, ind14221, ind13214, ind59182
    integer                :: ind13209, indClwFG, ind13208, ind13095, ind25174, indtmp
    integer                :: idata2,idata3,idata,idatend
    integer                :: flag_passage1,flag_passage2,flag_passage3
    integer                :: flag_passage4,flag_passage5
    integer                :: idatyp
    real                   :: val_option_r4
    character(len=9)       :: station_id
    
    write(*,*) '----------------------------------------------------------'
    write(*,*) '------- Begin brpr_addCloudParametersandEmissivity -------'
    write(*,*) '----------------------------------------------------------'

    ! Initialisation

    flag_passage1 = 0
    flag_passage2 = 0
    flag_passage3 = 0
    flag_passage4 = 0
    flag_passage5 = 0
    
    opt_missing = 'MISSING'
    val_option_r4  = -7777.77

    call BURP_Set_Options(                  &                
         REAL_OPTNAME       = opt_missing,  &
         REAL_OPTNAME_VALUE = val_option_r4,&
         iostat             = error )
    call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Options")

    call BURP_Init(inputFile, iostat=error)
    call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_init inputFile")
    call BURP_Init(inputReport,copyReport)
    call BURP_Init(inputBlock)

    ! Opening file
    write(*,*) 'OPENED FILE = ', trim(burpFile)

    call BURP_New(inputFile, &
         FILENAME = burpFile, &
         MODE     = FILE_ACC_APPEND, &
         iostat   = error )
    call handle_error(error, "brpr_addCloudParametersandEmissivity: problem opening input file")
    ! Obtain input burp file number of reports

    call BURP_Get_Property(inputFile, NRPTS=nb_rpts)

    ! Scan input burp file to get all reports address

    allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0

    do
      ref_rpt = BURP_Find_Report(inputFile, &
           report      = inputReport,       &
           SEARCH_FROM = ref_rpt,           &
           iostat      = error)
      call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Find_Report")
      if (ref_rpt < 0) Exit

      call BURP_Get_Property(inputReport, STNID=station_id)
      if (station_id(1:2)==">>") cycle

      count = count + 1
      address(count) = ref_rpt
    end do

    write(*,*) 
    write(*,*) 'NUMBER OF REPORTS WITH OBSERVATIONS = ',count
    write(*,*) 
    
    if ( count > 0 ) then

      ! Create a new report
      
      call BURP_New(copyReport, ALLOC_SPACE=20000000, iostat=error)
      call handle_error(error, "brpr_addCloudParametersandEmissivity: problem allocating copyReport")

      ! Loop on reports

      REPORTS: do reportIndex = 1, count

        call BURP_Get_Report(inputFile,        &
             report    = inputReport,          &
             REF       = address(reportIndex), &
             iostat    = error)
        call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Get_Report")
        if (reportIndex == 1) then
          call BURP_Get_Property(inputReport, IDTYP=idatyp)
          write(*,*) "brpr_addCloudParametersandEmissivity idatyp ", idatyp
          idata2 = -1
          call obs_set_current_header_list(obsSpaceData, 'TO')
          HEADER: do
            headerIndex = obs_getHeaderIndex(obsSpaceData)
            if (headerIndex < 0) exit HEADER  
            if  ( obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == idatyp .and.  &
                  obs_headElem_i(obsSpaceData,OBS_OTP,headerIndex) == fileIndex) then
              idata2 = headerIndex
              exit HEADER
            end if
          end do HEADER
          if (idata2 == -1) then
            write(*,*) "datyp ",idatyp," not found in input file !"
            write(*,*) "Nothing to do here ! Exiting ..."
            call  cleanup()
            return
          end if
          idata3 = idata2
        end if

        ! First loop on blocks

        ! Find bad profiles not in CMA. This occurs if :
        !  - all observations are -1 and/or have a quality flag not zero

        ref_blk = 0

        BLOCKS1: do

          ref_blk = BURP_Find_Block(inputReport, &
               BLOCK       = inputBlock,         &
               SEARCH_FROM = ref_blk,            &
               iostat      = error)
          call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Find_Block #1")
          if (ref_blk < 0) EXIT BLOCKS1

          call BURP_Get_Property(inputBlock, &
               NELE   = nbele,               &
               NVAL   = nvale,               &
               NT     = nte,                 &
               BFAM   = bfam,                &
               BTYP   = btyp,                &
               iostat = error)
          call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Get_Property #2")

          ! observation block (btyp = 0100 100011X XXXX)
          ! 0100 1000110 0000 = 9312
          btyp10    = ishft(btyp,-5)
          btyp10obs = 291
           
          if ( btyp10 == btyp10obs .and. bfam == 0 ) then

            allocate(goodprof(nte), btobs(nvale,nte))

            goodprof(:) = 0
            btobs(:,:)  = 0.

            ind012163  = BURP_Find_Element(inputBlock, ELEMENT=012163)
            if (ind012163 == -1) call handle_error(ind012163, "brpr_addCloudParametersandEmissivity: cannot find element 12163 in inputBlock")

            do tIndex=1,nte
              do valIndex=1,nvale
                btobs(valIndex,tIndex) = BURP_Get_Rval(inputBlock, &
                     NELE_IND = ind012163,                         &
                     NVAL_IND = valIndex,                          &
                     NT_IND   = tIndex )
                if ( btobs(valIndex,tIndex) > 0. ) goodprof(tIndex) = 1
              end do
            end do
            
          end if

        end do BLOCKS1

        call BURP_copy_Header(TO=copyReport, FROM=inputReport)

        call BURP_Init_Report_Write(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Init_Report_Write")

        ! Second loop on blocks

        ! add new informations

        ref_blk = 0
        
        BLOCKS2: do

          if ( .not. allocated(goodprof) ) then
            write(*,*)
            write(*,*) 'Resume report is position # ',reportIndex
            EXIT BLOCKS2
          end if

          ref_blk = BURP_Find_Block(inputReport, &
               BLOCK       = inputBlock, &
               SEARCH_FROM = ref_blk, &
               iostat      = error)
          call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Find_Block #2")
          if (ref_blk < 0) EXIT BLOCKS2

          call BURP_Get_Property(inputBlock, &
               NELE   = nbele,               &
               NVAL   = nvale,               &
               NT     = nte,                 &
               BFAM   = bfam,                &
               BTYP   = btyp,                &
               iostat = error)
          call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Get_Property #3")
          
          ! descriptor block (btyp = 0010 100000X XXXX) 
          ! 0010 1000000 0000==5120 )
          !    if profile contains rejected observations (apart from blacklisted channels),
          !     set bit 6 in global flags.

          btyp10    = ishft(btyp,-5)
          btyp10des = 160

          if ( btyp10 == btyp10des ) then

            flag_passage1 = 1

            allocate(glbflag(nte))

            ind055200  = BURP_Find_Element(inputBlock, ELEMENT=055200)
            if (ind055200 == -1) call handle_error(ind055200, "brpr_addCloudParametersandEmissivity: cannot find element 55200 in inputBlock")
            do tIndex = 1, nte
              glbflag(tIndex) =  BURP_Get_Tblval(inputBlock, &
                   NELE_IND = ind055200,                     &
                   NVAL_IND = 1,                             &
                   NT_IND   = tIndex )
            end do

            do tIndex = 1, nte
              if (goodprof(tIndex)/= 1) glbflag(tIndex) = ibset(glbflag(tIndex),6)            
            end do

            do tIndex = 1, nte
              call BURP_Set_Tblval(inputBlock, &
                   NELE_IND = ind055200,       &
                   NVAL_IND = 1,               &
                   NT_IND   = tIndex,          &
                   TBLVAL   = glbflag(tIndex), &
                   iostat   = error)
              call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Tblval 55200")
            end do
              
            deallocate(glbflag)

          end if
          ! For Hyperspectral data
          !
          if ( tvs_isIdBurpHyperSpectral(idatyp) ) then

            write(*,*) 'brpr_addCloudParametersandEmissivity for IR data'
            ! info block (btyp = 0001 100000X XXXX) 
            ! 0001 100000X XXXX = 3072
            btyp10    = ishft(btyp,-5)
            btyp10inf = 96
            if ( btyp10 == btyp10inf ) then
              flag_passage2 = 1
              ind14213 = BURP_Find_Element(inputBlock, ELEMENT=014213)
              ind14214 = BURP_Find_Element(inputBlock, ELEMENT=014214)
              ind14215 = BURP_Find_Element(inputBlock, ELEMENT=014215)
              ind14216 = BURP_Find_Element(inputBlock, ELEMENT=014216)
              ind14217 = BURP_Find_Element(inputBlock, ELEMENT=014217)
              ind14218 = BURP_Find_Element(inputBlock, ELEMENT=014218)
              ind14219 = BURP_Find_Element(inputBlock, ELEMENT=014219)
              ind14220 = BURP_Find_Element(inputBlock, ELEMENT=014220)
              ind14221 = BURP_Find_Element(inputBlock, ELEMENT=014221)
              ind13214 = BURP_Find_Element(inputBlock, ELEMENT=013214)
              ind59182 = BURP_Find_Element(inputBlock, ELEMENT=59182)
              if (ind14213 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=11, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #1")
                ind14213 = nbele+ 1
                ind14214 = nbele+ 2
                ind14215 = nbele+ 3
                ind14216 = nbele+ 4
                ind14217 = nbele+ 5
                ind14218 = nbele+ 6
                ind14219 = nbele+ 7
                ind14220 = nbele+ 8
                ind14221 = nbele+ 9
                ind13214 = nbele+ 10
                ind59182 = nbele+ 11
                call BURP_Set_Element(inputBlock, NELE_IND=ind14213, ELEMENT=014213, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14213")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14214, ELEMENT=014214, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14214")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14215, ELEMENT=014215, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14215")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14216, ELEMENT=014216, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14216")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14217, ELEMENT=014217, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14217")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14218, ELEMENT=014218, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14218")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14219, ELEMENT=014219, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14219")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14220, ELEMENT=014220, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14220")
                call BURP_Set_Element(inputBlock, NELE_IND=ind14221, ELEMENT=014221, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 14221")
                call BURP_Set_Element(inputBlock, NELE_IND=ind13214, ELEMENT=013214, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13214")
                call BURP_Set_Element(inputBlock, NELE_IND=ind59182, ELEMENT=059182, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 59182")
              end if

              ind008012 = BURP_Find_Element(inputBlock, ELEMENT  = 008012)
              if (ind008012 == -1) call handle_error(ind008012, "brpr_addCloudParametersandEmissivity: cannot find element 8012 in inputBlock #1")

              do tIndex = 1, nte

                if ( goodprof(tIndex) == 1 ) then

                  if ( obs_headElem_i(obsSpaceData,OBS_OTP,idata2)  /= fileIndex) then
                    headElem_i = obs_headElem_i(obsSpaceData,OBS_OTP,idata2)
                    write(*,*) "File Inconsistency ", headElem_i, fileIndex
                    write(*,*) "Should not happen..."
                    call utl_abort('brpr_addCloudParametersandEmissivity')
                    end if

                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ETOP,idata2)),ind14213,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_VTOP,idata2)),ind14214,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ECF,idata2)),ind14215,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_VCF,idata2)),ind14216,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_HE,idata2)),ind14217,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ZTSR,idata2)),ind14218,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_NCO2,idata2),ind14219,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ZTM,idata2)),ind14220,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ZTGM,idata2)),ind14221,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ZLQM,idata2)),ind13214,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_ZPS,idata2)),ind59182,1,tIndex)
                  call Insert_into_burp_i(tvs_ChangedStypValue(obsSpaceData,idata2),ind008012,1,tIndex)
                  idata2 = idata2 + 1

                else

                  call Insert_into_burp_r4(-1.0,ind14213,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14214,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14215,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14216,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14217,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14218,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14219,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14220,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind14221,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind13214,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind59182,1,tIndex)
                  call Insert_into_burp_i(-1,ind008012,1,tIndex)

                end if


              end do
 
            end if

            ! observation block (btyp = 0100 100011X XXXX)
            ! 0100 1000110 0000 = 9312
            btyp10    = ishft(btyp,-5)
            btyp10obs = 291

            if ( btyp10 == btyp10obs .and. bfam == 0 ) then
              flag_passage3 = 1

              indEmis  = BURP_Find_Element(inputBlock, ELEMENT=055043)
              if (indEmis < 0) then
                indEmis=nbele+1
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block data")
                call BURP_Set_Element(inputBlock, NELE_IND=indEmis, ELEMENT=055043, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 55043 #1")
                indEmis=nbele+1
              end if
              indchan  = BURP_Find_Element(inputBlock, ELEMENT=005042)
              if (indchan == -1) call handle_error(indchan, "brpr_addCloudParametersandEmissivity: cannot find element 5042 in inputBlock")
              do tIndex = 1, nte
                do valIndex = 1, nvale
                  call Insert_into_burp_i(-1,indEmis,valIndex,tIndex)
                end do

                if ( goodprof(tIndex) == 1 ) then

                  if ( obs_headElem_i(obsSpaceData,OBS_OTP,idata3)  /= fileIndex) then
                    headElem_i = obs_headElem_i(obsSpaceData,OBS_OTP,idata3)
                    write(*,*) "File Inconsistency emissivity block", &
                         headElem_i, fileIndex, idata3
                    write(*,*) "Should not happen..."
                    call utl_abort('brpr_addCloudParametersandEmissivity')
                  end if
                  idata   = obs_headElem_i(obsSpaceData,OBS_RLN,idata3)
                  idatend = obs_headElem_i(obsSpaceData,OBS_NLV,idata3) + idata - 1
                  do bodyIndex = idata, idatend
                    emisfc = 100.d0 * obs_bodyElem_r(obsspacedata,OBS_SEM,bodyIndex)
                    ichn = NINT(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
                    ichn = MAX(0,MIN(ichn,tvs_maxChannelNumber+1))
                    bl: do valIndex=1, nvale
                      ichnb=BURP_Get_Tblval(inputBlock, &
                           NELE_IND = indchan,          &
                           NVAL_IND = valIndex,         &
                           NT_IND   = tIndex)
                      if (ichn==ichnb) then
                        call Insert_into_burp_r4(sngl(emisfc), indEmis, valIndex, tIndex)
                        exit bl
                      end if
                    end do bl

                  end do

                  idata3 = idata3 + 1

                end if

              end do

            end if

            ! flag block (btyp = 0111 100011X XXXX)
            ! 0111 1000110 0000 = 15456
            btyp10    = ishft(btyp,-5)
            btyp10flg = 483

            if ( btyp10 == btyp10flg ) then
              flag_passage4 = 1

              indEmis  = BURP_Find_Element(inputBlock, ELEMENT=255043)
              if (indEmis < 0) then
                indEmis=nbele+1
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block marqueur")
                call BURP_Set_Element(inputBlock, NELE_IND=indEmis, ELEMENT=255043, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 255043")
              end if

              do tIndex = 1, nte
                do valIndex = 1, nvale
                  call BURP_Set_Tblval(inputBlock, &
                       NELE_IND = indEmis,         &
                       NVAL_IND = valIndex,        &
                       NT_IND   = tIndex,          &
                       TBLVAL   = 0, &
                       iostat   = error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Tblval 255043")
                end do
              end do
            end if

            ! O-P block (btyp = 0100 100011X XXXX)
            ! 0100 1000110 0000 = 9312
            btyp10    = ishft(btyp,-5)
            btyp10omp = 291

            if ( btyp10 == btyp10omp .and. bfam == 14 ) then
              flag_passage5 = 1

              indEmis  = BURP_Find_Element(inputBlock, ELEMENT=055043)
              if (indEmis < 0) then
                indEmis=nbele+1
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block O-P")
                call BURP_Set_Element(inputBlock, NELE_IND=nbele+1, ELEMENT=055043, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 55043 #2")
              end if

              do tIndex = 1, nte
                do valIndex = 1, nvale
                  call Insert_into_burp_i(-1,indEmis,valIndex,tIndex)
                end do
              end do

            end if

          end if ! hyper Spectral

          if ( (tvs_isIdBurpInst(idatyp,'atms' )) .or. &
               (tvs_isIdBurpInst(idatyp,'amsua')) ) then 
            ! info block (btyp = 0001 100000X XXXX) 
            ! 0001 100000X XXXX = 3072
            btyp10    = ishft(btyp,-5)
            btyp10inf = 96
            if ( btyp10 == btyp10inf ) then
              flag_passage2 = 1
              indtmp = nbele
              ! CLW
              ind13209 = BURP_Find_Element(inputBlock, ELEMENT=013209)
              if (ind13209 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #2")
                ind13209 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13209, ELEMENT=013209, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13209 #1")
                indtmp = indtmp + 1
              end if
              ! clwFG
              if ( tvs_mwAllskyAssim .and. &
                   tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                indClwFG = BURP_Find_Element(inputBlock, ELEMENT=clwFgElementId)
                if (indClwFG < 0) then
                  call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #3")
                  indClwFG = indtmp + 1
                  call BURP_Set_Element(inputBlock, NELE_IND=indClwFG, ELEMENT=clwFgElementId, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element clwFG #1")
                  indtmp = indtmp + 1
                end if
              end if
              ! SCATERING INDEX
              ind13208 = BURP_Find_Element(inputBlock, ELEMENT=013208)
              if (ind13208 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #4")
                ind13208 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13208, ELEMENT=013208, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 013208")
                indtmp = indtmp + 1
              end if

              ind008012 = BURP_Find_Element(inputBlock, ELEMENT  = 008012)
              if (ind008012 == -1) call handle_error(ind008012, "brpr_addCloudParametersandEmissivity: cannot find element 8012 in inputBlock #2")

              do tIndex = 1, nte

                if ( goodprof(tIndex) == 1 ) then

                  if ( obs_headElem_i(obsSpaceData,OBS_OTP,idata2)  /= fileIndex) then
                    headElem_i = obs_headElem_i(obsSpaceData,OBS_OTP,idata2)
                    write(*,*) "File Inconsistency ", headElem_i, fileIndex
                    write(*,*) "Should not happen..."
                    call utl_abort('brpr_addCloudParametersandEmissivity')
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWO,idata2)),ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWB,idata2)),indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_SIO,idata2)),ind13208,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_STYP,idata2),ind008012,1,tIndex)
                  idata2 = idata2 + 1

                else
                  call Insert_into_burp_r4(-1.0,ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(-1.0,indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(-1.0,ind13208,1,tIndex)
                  call Insert_into_burp_i(-1,ind008012,1,tIndex)

                end if

              end do
 
            end if
          end if ! tvs_isIdBurpInst(idatyp,'atms')) .or. (tvs_isIdBurpInst(idatyp,'amsua')

          if (tvs_isIdBurpInst(idatyp,'ssmis' )) then
            ! info block (btyp = 0001 100000X XXXX) 
            ! 0001 100000X XXXX = 3072
            btyp10    = ishft(btyp,-5)
            btyp10inf = 96
            if ( btyp10 == btyp10inf ) then
              flag_passage2 = 1
              indtmp = nbele

              ! LAND SEA QUALIFIER ELE 8012
              ind008012 = BURP_Find_Element(inputBlock, ELEMENT  = 008012)
              if (ind008012 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #5")
                ind008012 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind008012, ELEMENT=008012, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 8012 ssmis")
                indtmp = indtmp + 1
              end if

              ! TERRAIN TYPE ELE 13039
              ind13039 = BURP_Find_Element(inputBlock, ELEMENT  = 13039)
              if (ind13039 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #6")
                ind13039 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13039, ELEMENT=13039, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13039")
                indtmp = indtmp + 1
              end if
              
              ! SAT ZENITH ANGLE ELE 7024
              ind7024 = BURP_Find_Element(inputBlock, ELEMENT  = 7024)
              if (ind7024 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #7")
                ind7024 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind7024, ELEMENT=7024, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 7024")
                indtmp = indtmp + 1
              end if
              
              ! SAT AZIMUTH ANGLE ELE 5021
              ind5021 = BURP_Find_Element(inputBlock, ELEMENT  = 5021)
              if (ind5021 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #8")
                ind5021 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind5021, ELEMENT=5021, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 5021")
                indtmp = indtmp + 1
              end if

              ! CLW
              ind13209 = BURP_Find_Element(inputBlock, ELEMENT=013209)
              if (ind13209 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #9")
                ind13209 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13209, ELEMENT=013209, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13209 #2")
                indtmp = indtmp + 1
              end if
              ! clwFG
              if ( tvs_mwAllskyAssim .and. &
                   tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                indClwFG = BURP_Find_Element(inputBlock, ELEMENT=clwFgElementId)
                if (indClwFG < 0) then
                  call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #10")
                  indClwFG = indtmp + 1
                  call BURP_Set_Element(inputBlock, NELE_IND=indClwFG, ELEMENT=clwFgElementId, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element clwFG #2")
                  indtmp = indtmp + 1
                end if
              end if
              ! SSMIS INTEGRATED WATER VAPOR 
              ind13095 = BURP_Find_Element(inputBlock, ELEMENT=013095)
              if (ind13095 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #11")
                ind13095 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13095, ELEMENT=013095, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13095")
                indtmp = indtmp + 1
              end if
              ! SCATERING INDEX
              ind13208 = BURP_Find_Element(inputBlock, ELEMENT=013208)
              if (ind13208 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #12")
                ind13208 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13208, ELEMENT=013208, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13208")
                indtmp = indtmp + 1
              end if

              
              do tIndex = 1, nte

                if ( goodprof(tIndex) == 1 ) then

                  if ( obs_headElem_i(obsSpaceData,OBS_OTP,idata2)  /= fileIndex) then
                    headElem_i = obs_headElem_i(obsSpaceData,OBS_OTP,idata2)
                    write(*,*) "File Inconsistency ", headElem_i, fileIndex
                    write(*,*) "Should not happen..."
                    call utl_abort('brpr_addCloudParametersandEmissivity')
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWO,idata2)),ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWB,idata2)),indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_SIO,idata2)),ind13208,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_IWV,idata2)),ind13095,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_STYP,idata2),ind008012,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_TTYP,idata2),ind13039,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_SZA,idata2)),ind7024,1,tIndex)
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_AZA,idata2)),ind5021,1,tIndex)
                  idata2 = idata2 + 1
                else
                  call Insert_into_burp_r4(-1.0,ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(-1.0,indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(-1.0,ind13208,1,tIndex)
                  !call Insert_into_burp_r4(-1.0,ind13095,1,tIndex)
                  call Insert_into_burp_i(-1,ind008012,1,tIndex)
                  call Insert_into_burp_i(-1,ind13039,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind7024,1,tIndex)
                  call Insert_into_burp_r4(-1.0,ind5021,1,tIndex)
                end if
              end do
 
            end if
          end if ! tvs_isIdBurpInst(idatyp,'ssmis')

          if (tvs_isIdBurpInst(idatyp,'mwhs2' )) then
            ! info block (btyp = 0001 100000X XXXX)
            ! 0001 100000X XXXX = 3072
            btyp10    = ishft(btyp,-5)
            btyp10inf = 96
            if ( btyp10 == btyp10inf ) then
              flag_passage2 = 1
              indtmp = nbele

              ! LAND SEA QUALIFIER ELE 8012
              ind008012 = BURP_Find_Element(inputBlock, ELEMENT  = 008012)
              if (ind008012 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #5")
                ind008012 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind008012, ELEMENT=008012, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 8012 ssmis")
                indtmp = indtmp + 1
              end if

              ! TERRAIN TYPE ELE 13039
              ind13039 = BURP_Find_Element(inputBlock, ELEMENT  = 13039)
              if (ind13039 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #6")
                ind13039 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13039, ELEMENT=13039, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13039")
                indtmp = indtmp + 1
              end if

              ! SPECIAL QC FLAG INTEGER
              ind25174 = BURP_Find_Element(inputBlock, ELEMENT=025174)
              if (ind25174 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #13")
                ind25174 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind25174, ELEMENT=025174, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 25174")
                indtmp = indtmp + 1
              end if

              ! CLW
              ind13209 = BURP_Find_Element(inputBlock, ELEMENT=013209)
              if (ind13209 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #9")
                ind13209 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13209, ELEMENT=013209, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13209 #2")
                indtmp = indtmp + 1
              end if

              ! clwFG
              if ( tvs_mwAllskyAssim .and. &
                   tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                indClwFG = BURP_Find_Element(inputBlock, ELEMENT=clwFgElementId)
                if (indClwFG < 0) then
                  call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #10")
                  indClwFG = indtmp + 1
                  call BURP_Set_Element(inputBlock, NELE_IND=indClwFG, ELEMENT=clwFgElementId, iostat=error)
                  call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element clwFG #2")
                  indtmp = indtmp + 1
                end if
              end if

              ! SCATERING INDEX
              ind13208 = BURP_Find_Element(inputBlock, ELEMENT=013208)
              if (ind13208 < 0) then
                call BURP_Resize_Block(inputBlock, ADD_NELE=1, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Resize_Block info #12")
                ind13208 = indtmp + 1
                call BURP_Set_Element(inputBlock, NELE_IND=ind13208, ELEMENT=013208, iostat=error)
                call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Set_Element 13208")
                indtmp = indtmp + 1
              end if

              do tIndex = 1, nte

                if ( goodprof(tIndex) == 1 ) then

                  if ( obs_headElem_i(obsSpaceData,OBS_OTP,idata2)  /= fileIndex) then
                    headElem_i = obs_headElem_i(obsSpaceData,OBS_OTP,idata2)
                    write(*,*) "File Inconsistency ", headElem_i, fileIndex
                    write(*,*) "Should not happen..."
                    call utl_abort('brpr_addCloudParametersandEmissivity')
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWO,idata2)),ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_CLWB,idata2)),indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(sngl(obs_headElem_r(obsSpaceData,OBS_SIO,idata2)),ind13208,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_STYP,idata2),ind008012,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_TTYP,idata2),ind13039,1,tIndex)
                  call Insert_into_burp_i(obs_headElem_i(obsSpaceData,OBS_INFG,idata2),ind25174,1,tIndex)
                  idata2 = idata2 + 1
                else
                  call Insert_into_burp_r4(-1.0,ind13209,1,tIndex)
                  if ( tvs_mwAllskyAssim .and. &
                       tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
                    call Insert_into_burp_r4(-1.0,indClwFG,1,tIndex)
                  end if
                  call Insert_into_burp_r4(-1.0,ind13208,1,tIndex)
                  call Insert_into_burp_i(-1,ind008012,1,tIndex)
                  call Insert_into_burp_i(-1,ind13039,1,tIndex)
                  call Insert_into_burp_i(-1,ind25174,1,tIndex)
                end if
              end do

            end if
          end if ! tvs_isIdBurpInst(idatyp,'mwhs2')

          ! Add block into new report

          if ( btyp == 5120 ) then
            call BURP_Write_Block(copyReport, inputBlock, &
                 ENCODE_BLOCK  = .true., &
                 iostat        = error)
          else
            call BURP_Write_Block(copyReport, inputBlock, &
                 ENCODE_BLOCK  = .true., &
                 CONVERT_BLOCK = .true., &
                 iostat        = error)
          end if
          call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Write_Block")
        end do BLOCKS2

        if ( allocated(goodprof) ) then
          deallocate (goodprof,btobs)
        end if
        ! Write new report into file        
        call BURP_Delete_Report(inputFile, inputReport, iostat=error)
        call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Delete_Report")
        call BURP_Write_Report(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_addCloudParametersandEmissivity: BURP_Write_Report")
      end do REPORTS
      if ( tvs_isIdBurpHyperSpectral(idatyp) ) then
        if ( flag_passage1 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - descriptor block not seen ? Verify btyp'
        end if
        if ( flag_passage2 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - info block not seen ? Verify btyp'
        end if
        if ( flag_passage3 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - observation block not seen ? Verify btyp'
        end if
        if ( flag_passage4 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - flag block not seen ? Verify btyp'
        end if
        if ( flag_passage5 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - O-P block not seen ? Verify btyp'
        end if
      else if ( (tvs_isIdBurpInst(idatyp,'atms' )) .or. &
                (tvs_isIdBurpInst(idatyp,'amsua')) .or. & 
                (tvs_isIdBurpInst(idatyp,'ssmis')) ) then 
        if ( flag_passage1 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - descriptor block not seen ? Verify btyp'
        end if
        if ( flag_passage2 == 0 ) then
          write(*,*)
          write(*,*) 'ERROR - info block not seen ? Verify btyp'
        end if
      end if

    end if !! End of 'if ( count > 0 )'

    call  cleanup()

  contains

    !--------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(5)

      errors(:) = 0
      if (allocated(address)) deallocate(address, stat = errors(1))
      call BURP_Free(InputFile, iostat=errors(2))
      call BURP_Free(InputReport, iostat=errors(3))
      call BURP_Free(InputReport, iostat=errors(4))
      call BURP_Free(InputBlock, iostat=errors(5))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_addCloudParametersandEmissivity: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !--------- handle_error -----
    subroutine handle_error(icode,errormessage)
      implicit none
      character (len=*) :: errormessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errormessage))
      end if
    end subroutine handle_error

    !--------- insert_into_burp_r4 -----
    subroutine Insert_into_burp_r4( r4val, pele, pval, pt )

      implicit none

      real(4), intent(in) :: r4val
      integer, intent(in) :: pele
      integer, intent(in) :: pval
      integer, intent(in) :: pt

      integer :: error
      
      if ( r4val >= 0. ) then
        call BURP_Set_Rval(inputBlock, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = r4val, &
             iostat   = error)
      else
        call BURP_Set_Rval(inputBlock, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = val_option_r4, &
             iostat   = error)
      end if
      if (error/=burp_noerr) then
        write(*,*) "Insert_into_burp_r4: r4val,pele,pval,pt", r4val, pele, pval, pt
        call handle_error(error, "brpr_addCloudParametersandEmissivity: Insert_into_burp_r4")
      end if

    end subroutine Insert_into_burp_r4

    !--------- insert_into_burp_i -----
    subroutine Insert_into_burp_i( ival, pele, pval, pt )
      !
      implicit none

      integer, intent(in) :: ival
      integer, intent(in) :: pele
      integer, intent(in) :: pval
      integer, intent(in) :: pt

      integer :: error

      if ( ival >= 0 ) then
        call BURP_Set_Rval(inputBlock, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = real(ival), &
             iostat   = error)
      else
        call BURP_Set_Rval(inputBlock, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = val_option_r4, &
             iostat   = error)
      end if
      
      if (error/=burp_noerr) then
        write(*,*) "Insert_into_burp_i: ival,pele,pval,pt", ival, pele, pval, pt
        call handle_error(error, "brpr_addCloudParametersandEmissivity: Insert_into_burp_i")
      end if
    end subroutine Insert_into_burp_i

  end subroutine brpr_addCloudParametersandEmissivity

  !--------------------------------------------------------------------------
  ! brpr_updateMissingObsFlags
  !--------------------------------------------------------------------------
  subroutine brpr_updateMissingObsFlags( burpFile )
    !
    ! :Purpose: Open burp file and set missing data flags to 2048.
    !
    implicit none

    ! Arguments:
    character (len=*), intent(in)    :: burpFile

    ! Locals:
    type(BURP_FILE)        :: inputFile
    type(BURP_RPT)         :: inputReport,copyReport
    type(BURP_BLOCK)       :: inputBlock
    integer                :: btyp10, btyp, bfam, error
    integer                :: btyp10obs, btyp10flg
    integer                :: nb_rpts, ref_rpt, ref_blk, count
    integer, allocatable   :: address(:)
    integer, allocatable   :: btobs(:,:)
    logical, allocatable   :: goodTB(:,:)
    integer                :: nbele,nvale,nte
    integer                :: valIndex, tIndex, reportIndex
    integer                :: ind012163,ind212163
    integer                :: flag_passage, flagval
    character(len=9)       :: station_id

    write(*,*) '----------------------------------------------------------'
    write(*,*) '-- Begin subroutine brpr_updateMissingObsFlags----'
    write(*,*) '----------------------------------------------------------'

    ! Initialisation
    ! Opening file
    write(*,*) 'OPENED FILE = ', trim(burpFile)

    call BURP_New(inputFile, &
         FILENAME = burpFile, &
         MODE     = FILE_ACC_APPEND, &
         iostat   = error )
    call handle_error(error, "brpr_updateMissingObsFlags: BURP_New inputFile")
    ! Obtain input burp file number of reports

    call BURP_Get_Property(inputFile, NRPTS=nb_rpts)

    ! Scan input burp file to get all reports address

    allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0

    do
      ref_rpt = BURP_Find_Report(inputFile, &
           report      = inputReport,       &
           SEARCH_FROM = ref_rpt,           &
           iostat      = error)
      call handle_error(error, "brpr_updateMissingObsFlags: BURP_Find_Report")
      if (ref_rpt < 0) Exit

      call BURP_Get_Property(inputReport, STNID=station_id)
      if (station_id(1:2)==">>") cycle

      count = count + 1
      address(count) = ref_rpt
    end do

    write(*,*)
    write(*,*) 'NUMBER OF REPORTS WITH OBSERVATIONS = ',count
    write(*,*)

    if ( count > 0 ) then

      ! Create a new report

      call BURP_New(copyReport, ALLOC_SPACE=20000000, iostat=error)
      call handle_error(error, "brpr_updateMissingObsFlags: Error creating new directory")

      ! Loop on reports

      REPORTS: do reportIndex = 1, count

        call BURP_Get_Report(inputFile,        &
             report    = inputReport,          &
             REF       = address(reportIndex), &
             iostat    = error)
        call handle_error(error, "brpr_updateMissingObsFlags: BURP_Get_Report")
        ! Find bad/missing TB. 

        ref_blk = 0

        BLOCKS1: do

          ref_blk = BURP_Find_Block(inputReport, &
               BLOCK       = inputBlock,         &
               SEARCH_FROM = ref_blk,            &
               convert = .false.,            &
               iostat      = error)
          call handle_error(error, "brpr_updateMissingObsFlags: BURP_Find_Block #1")
          if (ref_blk < 0) exit BLOCKS1

          call BURP_Get_Property(inputBlock, &
               NELE   = nbele,               &
               NVAL   = nvale,               &
               NT     = nte,                 &
               BFAM   = bfam,                &
               BTYP   = btyp,                &
               iostat = error)
          call handle_error(error, "brpr_updateMissingObsFlags: BURP_Get_Property #1")

          ! observation block (btyp = 0100 100011X XXXX)
          ! 0100 1000110 0000 = 9312
          btyp10    = ishft(btyp,-5)
          btyp10obs = 291

          if ( btyp10 == btyp10obs .and. bfam == 0 ) then

           ind012163  = BURP_Find_Element(inputBlock, ELEMENT=012163)
           if ( ind012163 < 0 ) exit BLOCKS1
           allocate(btobs( nvale,nte))
           allocate(goodTB(nvale,nte))
            goodTB(:,:) = .false.
            btobs(:,:)  = 0
            do tIndex=1,nte
              do valIndex=1,nvale
                btobs(valIndex,tIndex) = BURP_Get_Tblval(inputBlock, &
                     NELE_IND = ind012163,                         &
                     NVAL_IND = valIndex,                          &
                     NT_IND   = tIndex )
                if ( btobs(valIndex,tIndex) /= -1 ) goodTB(valIndex,tIndex) = .true.
              end do
            end do

          end if

        end do BLOCKS1

        call BURP_copy_Header(TO=copyReport, FROM=inputReport)

        call BURP_Init_Report_Write(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_updateMissingObsFlags: BURP_Init_Report_Write")        

        ! Second loop on blocks

        ! to set missingFlag when not goodTB

        ref_blk = 0

        BLOCKS2: do

          if ( .not. allocated(goodTB) ) then
            write(*,*)
            write(*,*) 'Resume report is position # ',reportIndex 
            exit BLOCKS2
          end if
          ref_blk = BURP_Find_Block(inputReport, &
               BLOCK       = inputBlock, &
               SEARCH_FROM = ref_blk, &
               convert = .false., &
               iostat      = error)
          call handle_error(error, "brpr_updateMissingObsFlags: BURP_Find_Block #2")
          if (ref_blk < 0) exit BLOCKS2

          call BURP_Get_Property(inputBlock, &
               NELE   = nbele,               &
               NVAL   = nvale,               &
               NT     = nte,                 &
               BFAM   = bfam,                &
               BTYP   = btyp,                &
               iostat = error)
          call handle_error(error, "brpr_updateMissingObsFlags: BURP_Get_Property #2")

          ! flag block (btyp = 0111 100011X XXXX)
          ! 0111 1000110 0000 = 15456
          btyp10    = ishft(btyp,-5)
          btyp10flg = 483

          if ( btyp10 == btyp10flg ) then
            flag_passage = 1

            ind212163  = BURP_Find_Element(inputBlock, ELEMENT=212163)
            if (ind212163 == -1) call handle_error(ind212163, "brpr_updateMissingObsFlags: cannot find element 212163 in inputBlock")

            do tIndex = 1, nte
              do valIndex = 1, nvale
                if ( .not. goodTB(valIndex, tIndex) ) then
                  flagval = BURP_Get_Tblval(inputBlock, &
                          NELE_IND = ind212163,         &
                          NVAL_IND = valIndex,        &
                          NT_IND   = tIndex,          &
                          iostat   = error)
                  call handle_error(error, "brpr_updateMissingObsFlags: BURP_Get_Tblval pour ind212163")
                  flagval = ibset(flagval,11)
                  flagval = ibset(flagval,7)
                  flagval = ibset(flagval,9)
                  call BURP_Set_Tblval(inputBlock, &
                       NELE_IND = ind212163,         &
                       NVAL_IND = valIndex,        &
                       NT_IND   = tIndex,          &
                       TBLVAL   = flagval, &
                       iostat   = error)
                  call handle_error(error, "brpr_updateMissingObsFlags: BURP_Set_Tblval pour ind212163")
                end if
              end do
            end do
          end if
          ! Add block into new report
          call BURP_Write_Block(copyReport, inputBlock, &
               ENCODE_BLOCK  = .false., &
               CONVERT_BLOCK = .false., &
               iostat        = error)
          if (error/=burp_noerr) then
            write(*,*)"Btyp= ",btyp
            call handle_error(error, "brpr_updateMissingObsFlags: BURP_Write_Block")
          end if
        end do BLOCKS2
        if (allocated(goodTB) ) then
          deallocate(goodTB)
          deallocate(btobs)
        end if
        ! Write new report into file        
        call BURP_Delete_Report(inputFile, inputReport, iostat=error)
        call handle_error(error, "brpr_updateMissingObsFlags: BURP_Delete_Report")
        call BURP_Write_Report(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_updateMissingObsFlags: BURP_Write_Report")
      end do REPORTS
    end if !! End of 'if ( count > 0 )'

    call  cleanup()

  contains

    !--------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(5)

      errors(:) = 0
      if (allocated(address)) deallocate(address, stat=errors(1))
      call BURP_Free(InputFile, iostat=errors(2))
      call BURP_Free(InputReport, iostat=errors(3))
      call BURP_Free(CopyReport, iostat=errors(4))
      call BURP_Free(InputBlock, iostat=errors(5))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_updateMissingObsFlags: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !--------- handle_error -----
    subroutine handle_error(icode,errormessage)
      implicit none
      character (len=*) :: errormessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errormessage))
      end if
    end subroutine handle_error
    
  end subroutine brpr_updateMissingObsFlags


  !-----------------------------------------------------------------------
  ! brpr_addElementsToBurp
  !-----------------------------------------------------------------------
  subroutine brpr_addElementsToBurp(inputFileName, familyType, beSilent_opt)
    !
    !:Purpose: to add element(s) for bias correction to data block of DERIALT BURP file
    !
    implicit none
    !Arguments:
    character(len=*), intent(in)  :: inputFileName
    character(len=*), intent(in)  :: familyType
    logical, optional             :: beSilent_opt
    !Locals:
    type(burp_file)             :: inputFile
    type(burp_rpt)              :: inputReport, copyReport
    type(burp_block)            :: inputBlock
    integer                     :: nb_rpts, ref_rpt, ref_blk, count
    integer, allocatable        :: address(:)
    integer                     :: nbele, nvale, nte
    integer                     :: valIndex, tIndex, reportIndex, btyp, idatyp, bfam, error
    integer                     :: ind, indele, nsize, iun_burpin
    integer                     :: nulnam
    character(len=9)            :: station_id
    character(len=7), parameter :: opt_missing='MISSING'
    character(len=codtyp_name_length) :: instrumName
    integer                     :: instrumID, previousInstrumID
    integer                     :: icodele(3)
    integer                     :: icodeleMrq(3) 
    integer                     :: btClearMrqElementID
    real, parameter             :: val_option = -9999.0
    integer, external           :: mrfmxl
    logical                     :: isDerialt, isInstrumUsingCLW
    logical                     :: beSilent
    integer                     :: nElem, elemIndex, burpIndex

    namelist /NAMADDTOBURP/ addBtClearToBurp, clwFgElementId, btClearElementId

    write(*,*) '-----------------------------------------------'
    write(*,*) '- begin brpr_addElementsToBurp -'
    write(*,*) '-----------------------------------------------'
 
    icodele(:)    = 0
    icodeleMrq(:) = 0
    
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    select case(familyType)
    case("TO")
      nElem = 1
      icodele(1) = 12233
    case("GP")
      nElem = 1
      icodele(1) = 15234
    case("UA")
      nElem = 2
      icodele(1:2) = (/12204,12243/)
    case default
      return
    end select

    icodeleMrq(1:nElem) =  200000 + icodele(1:nElem)

    ! Read the NAMADDTOBURP namelist (if it exists)
    
    addBtClearToBurp = .false.
    clwFgElementId = -1 
    btClearElementId = -1
    btClearMrqElementID = -200001
    
    if ( familyType == "TO" ) then
      if (utl_isNamelistPresent('NAMADDTOBURP','./flnml')) then
        ! read the namelist
        nulnam = 0
        error = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
        read(nulnam, nml=NAMADDTOBURP, iostat=error)
        if ( error /= 0 ) call utl_abort('brpr_addElementsToBurp: Error reading namelist')
        write(*,nml=NAMADDTOBURP)
        error = fclos(nulnam)
      else
        write(*,*)
        write(*,*) 'brpr_addElementsToBurp: Namelist block NAMADDTOBURP is missing in the namelist.'
        write(*,*) '                               The default value will be taken.'
        write(*,nml=NAMADDTOBURP)
      end if

      ! check clear-sky radiance element is in the namelist
      if ( addBtClearToBurp .and. btClearElementId < 0 ) then
        call utl_abort('brpr_addElementsToBurp: btClearElementId missing in the namelist')
      end if

      btClearMrqElementID = 200000 + btClearElementId
    end if

    ! initialisation
    ! --------------
    call burp_set_options(                 &
         real_optname       = opt_missing, &
         real_optname_value = val_option,  &
         iostat             = error )
    call handle_error(error, "brpr_addElementsToBurp: burp_set_options")

    call burp_init(inputFile,iostat=error)
    call handle_error(error, "brpr_addElementsToBurp: burp_init inputFile")
    call burp_init(inputReport,copyReport)
    call burp_init(inputBlock)

    ! opening file
    ! ------------
    write(*,*) 'opened file = ', trim( inputFileName )

    call burp_new(inputFile,         &
         filename = inputFileName,   &
         mode     = file_acc_append, &
         iostat   = error )
  
    if (error /= burp_noerr) then
      write(*,*) "cannot open BURP input file ", inputFileName
      call utl_abort('brpr_addElementsToBurp')
    end if

    ! obtain input burp file number of reports
    ! ----------------------------------------
    call burp_get_property(inputFile, nrpts=nb_rpts, io_unit= iun_burpin)

    nsize = mrfmxl(iun_burpin)
    if ( addBtClearToBurp ) then
      nsize = 4 * nsize
    else
      nsize = 3 * nsize
    end if

    write(*,*) "nsize= ", nsize
    write(*,*) 
    write(*,*) 'number of reports with observations in input file = ', nb_rpts - 1
    write(*,*) 

    ! scan input burp file to get all reports address
    ! -----------------------------------------------
    allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0
    isDerialt = .false.
    previousInstrumID = -1
    do
      ref_rpt = burp_find_report(inputFile, &
           report      = inputReport,       &
           search_from = ref_rpt,           &
           iostat      = error)
      if (ref_rpt < 0) exit
      count = count + 1
      address(count) = ref_rpt

      call burp_get_property(inputReport, stnid = station_id, idtyp = idatyp )
      if (station_id == ">>DERIALT") isDerialt = .true.
      if (station_id(1:2) == '>>') cycle
      
      if ( familyType == "TO" ) then
        if ( .not. beSilent ) then
          if ( count == 1 ) then
            write(*,*) 'brpr_addElementsToBurp: tvs_mwAllskyAssim =', tvs_mwAllskyAssim
            write(*,*) 'brpr_addElementsToBurp: clwFgElementId =', clwFgElementId 
          end if

          instrumName = codtyp_get_name(idatyp)
          instrumID = tvs_getInstrumentId(instrumName)
          isInstrumUsingCLW = tvs_isInstrumUsingCLW(tvs_getInstrumentId(instrumName))
          if (instrumID /= previousInstrumID) then
            write(*,*) 'brpr_addElementsToBurp: for report count =', count, &
              ', instrumentName=', instrumName, &
              ', instrumentId =', instrumID, &
              ', isInstrumUsingCLW =', isInstrumUsingCLW
            previousInstrumID = instrumID
          end if
        end if

        ! check clwFG element is in the namelist in all-sky mode.
        if ( tvs_mwAllskyAssim .and. clwFgElementId < 0 .and. &
          tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
          call utl_abort('brpr_addElementsToBurp: clwFgElementId missing in the namelist')
        end if
      end if
      
    end do

    if ( count > 0 .and. isDerialt) then
      write(*,*) "brpr_addElementsToBurp: modifying file..."
     


      ! create a new report
      ! ------------------     
      call burp_new(copyReport, alloc_space=nsize, iostat=error)
      call handle_error(error, "brpr_addElementsToBurp: burp_new problem allocating copyReport")

      ! loop on reports
      ! ---------------
      reports: do reportIndex = 1, count
    
        call burp_get_report(inputFile,        &
             report    = inputReport,          &
             ref       = address(reportIndex), &
             iostat    = error)      

        call burp_copy_header(to=copyReport,from=inputReport)
        
        call burp_get_property(inputReport, stnid = station_id, idtyp = idatyp )

        call burp_init_report_write(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_addElementsToBurp: burp_init_report_write")

        ! loop on blocks
        ! --------------------
        ref_blk = 0
      
        blocks: do

          ref_blk = burp_find_block(inputReport, &
               block       = inputBlock,         &
               search_from = ref_blk,            &
               iostat      = error)

          if (ref_blk < 0) exit blocks
         
          call burp_get_property(inputBlock, &
               nele   = nbele,               &
               nval   = nvale,               &
               nt     = nte,                 &
               bfam   = bfam,                &
               btyp   = btyp,                & 
               iostat = error)
          call handle_error(error, "brpr_addElementsToBurp: burp_get_property")

          if ( isObsBlock(familyType,btyp) .and. bfam == 0 ) then 

            do elemIndex = 1, nElem
              burpIndex = burp_find_element(inputBlock, element=icodele(elemIndex), iostat=error)
              if ( burpIndex <= 0 ) then
                nbele = nbele + 1
                call burp_resize_block(InputBlock, ADD_NELE = 1, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_resize_block #1")
                call burp_set_element(InputBlock, NELE_IND = nbele, ELEMENT = icodele(elemIndex), IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_set_element #1")
                do valIndex = 1,nvale
                  do tIndex = 1,nte
                    call burp_set_rval( inputBlock, &
                       nele_ind = nbele,            &
                       nval_ind = valIndex,         &
                       nt_ind   = tIndex,           &
                       rval   = val_option, iostat=error)
                    call handle_error(error, "brpr_addElementsToBurp: burp_set_rval #1")
                  end do
                end do
              end if
            end do
        
            ! Adding clear-sky radiance to data block for instrument in all-sky mode.
            if ( tvs_mwAllskyAssim .and. addBtClearToBurp .and. &
                 tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then
              
              indele = burp_find_element(inputBlock, element=btClearElementId)

              if ( indele <= 0 ) then
                nbele = nbele + 1
                call burp_resize_block(InputBlock, ADD_NELE = 1, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_resize")
                Call burp_set_element(InputBlock, NELE_IND = nbele, ELEMENT = btClearElementId, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_set_element #2")
                do valIndex = 1,nvale
                  do tIndex = 1,nte
                    call burp_set_Rval( inputBlock, &
                         nele_ind = nbele,            &
                         nval_ind = valIndex,         &
                         nt_ind   = tIndex,           &
                         Rval = MPC_missingValue_R4, iostat=error)
                    call handle_error(error, "brpr_addElementsToBurp: burp_set_rval #2")
                  end do
                end do
              end if
          
            end if

            call burp_write_block(copyReport, block  = inputBlock,  &
                 convert_block =.true., encode_block=.true., iostat=error)
            call handle_error(error, "brpr_addElementsToBurp: burp_write_block #1")

          else if ( isFlagBlock(familyType,btyp) .and. bfam == 0 ) then     !  MRQ block 
            do elemIndex = 1, nElem
              ind = burp_find_element(inputBlock, element=icodeleMrq(elemIndex), iostat=error)
              if ( ind <= 0 ) then
                nbele = nbele + 1
                call burp_resize_block(InputBlock, ADD_NELE = 1, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_resize_block #2")
                Call burp_set_element(InputBlock, NELE_IND = nbele, ELEMENT = icodeleMrq(elemIndex), IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_set_element #3")
                do valIndex = 1,nvale
                  do tIndex = 1, nte
                    call burp_set_tblval( inputBlock, &
                         nele_ind = nbele,            &
                         nval_ind = valIndex,         &
                         nt_ind   = tIndex,           &
                         tblval   = 0, iostat=error)
                        call handle_error(error, "brpr_addElementsToBurp: burp_set_tblval #1")
                  end do
                end do
              end if
            end do
        
            ! Adding clear-sky radiance to MRQ block for instrument in all-sky mode.
            if ( tvs_mwAllskyAssim .and. addBtClearToBurp .and. &
                 tvs_isInstrumUsingCLW(tvs_getInstrumentId(codtyp_get_name(idatyp))) ) then

              indele = burp_find_element(inputBlock, element=btClearMrqElementID)
              if ( indele <= 0 ) then
                nbele = nbele + 1
                call burp_resize_block(InputBlock, ADD_NELE = 1, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_resize_block #3")
                Call burp_set_element(InputBlock, NELE_IND = nbele, ELEMENT = btClearMrqElementID, IOSTAT = error)
                call handle_error(error, "brpr_addElementsToBurp: burp_set_element #4")
                do valIndex = 1,nvale
                  do tIndex = 1, nte
                    call burp_set_tblval( inputBlock, &
                         nele_ind = nbele,            &
                         nval_ind = valIndex,         &
                         nt_ind   = tIndex,           &
                         tblval   = 0, iostat=error)
                    call handle_error(error, "brpr_addElementsToBurp: burp_set_tblval #2")
                  end do
                end do
              end if
          
            end if

            call burp_write_block(copyReport, block  = inputBlock,  &
                 convert_block =.false., encode_block=.true.,iostat=error)
            call handle_error(error, "brpr_addElementsToBurp: burp_write_block #2")

          else !other blocks

            call burp_write_block(copyReport, block  = inputBlock,  &
                 convert_block = ( btyp /= 5120 .and. btyp /= 2080), iostat=error)
            call handle_error(error, "brpr_addElementsToBurp: burp_write_block #3")

          end if
  
        end do blocks
      
        ! write new report into file
        ! --------------------------
        call BURP_Delete_Report(inputFile, inputReport, iostat=error)
        call handle_error(error, "brpr_addElementsToBurp: burp_delete_report")
        call burp_write_report(inputFile,copyReport, iostat=error)
        call handle_error(error, "brpr_addElementsToBurp: burp_write_report")
    
      end do reports

    end if

    call  cleanup()

    write(*,*) '---------------------------------------------'
    write(*,*) '- end brpr_addElementsToBurp -'
    write(*,*) '---------------------------------------------'

  contains

    !-------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(5)

      errors(:) = 0
      if (allocated(address)) deallocate(address, stat=errors(1))
      call burp_free(inputFile, iostat=errors(2))
      call burp_free(inputReport, iostat=errors(3))
      call burp_free(copyReport, iostat=errors(4))
      call burp_free(inputBlock, iostat=errors(5))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_addElementsToBurp: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !-------- handle_error -----
    subroutine handle_error(icode, errorMessage)
      implicit none
      character (len=*) :: errorMessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errorMessage))
      end if
    end subroutine handle_error
    
  end subroutine brpr_addElementsToBurp

  !-----------------------------------------------------------------------
  ! brpr_burpClean
  !-----------------------------------------------------------------------
  subroutine brpr_burpClean(inputFileName, familyType)
    !
    !:Purpose: to remove observations that are flagged not to be assimilated
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: inputFileName
    character(len=*), intent(in) :: familyType

    ! Locals:
    type(burp_file)             :: inputFile
    type(burp_rpt)              :: inputReport, copyReport
    type(burp_block)            :: inputBlock, inputBlock2
    integer, allocatable        :: addresses(:), elementIdsRead(:), elementIdsBlock(:)
    integer, allocatable        :: flagValues(:,:,:)
    real(4), allocatable        :: obsValues(:,:,:)
    logical, allocatable        :: rejectObs(:,:)
    integer                     :: numReports, refBlock, intBurpValue
    integer                     :: numElem, numLevels, numObsProfiles, numElem2
    integer                     :: numLevels2, newNumLevels
    integer                     :: numObsProfiles2, newNumObsProfiles
    integer                     :: elemIndex, levelIndex, levelIndexGood
    integer                     :: obsProfIndex, obsProfIndexGood
    integer                     :: reportIndex, btyp, datyp, error
    integer                     :: nsize, iun_burpin, numReject, numRejectTotal
    character(len=7), parameter :: opt_missing='MISSING'
    real, parameter             :: missingValue = -9999.0
    integer, external           :: mrfbfl
    logical                     :: groupedData, foundFlags, foundObs, emptyReport
    logical                     :: resumeReport, cleanLevels, checkBlock
    character(len=2)            :: familyTypesToDo(7) = (/'AI','SW','TO','SC','GP','UA','SF'/)
    character(len=9)            :: stnid
    logical                     :: debug = .false.

    write(*,*)
    write(*,*) 'brpr_burpClean: starting'

    ! only apply for certain obs families for now
    if ( all( trim(familyType) /= familyTypesToDo(:) ) ) then
      write(*,*) 'brpr_burpClean: not applied to obs family = ', trim(familyType)
      return
    end if

    ! for some obs types we will remove levels/channels, not just complete profiles
    cleanLevels = (trim(familyType) == 'UA')

    ! obtain input burp file report addresses and number of reports
    call getBurpReportAddresses(inputFileName, addresses)
    numReports = size(addresses)

    ! initialisation
    call burp_set_options(                   &
         real_optname       = opt_missing,   &
         real_optname_value = missingValue,  &
         iostat             = error )
    call handle_error(error, "brpr_burpClean: burp_set_options")

    ! initialize burp objects
    call burp_init(inputFile,iostat=error)
    call handle_error(error, "brpr_burpClean: burp_init inputFile")
    call burp_init(inputReport)
    call burp_init(copyReport)
    call burp_init(inputBlock)
    call burp_init(inputBlock2)

    ! opening file
    write(*,*) 'brpr_burpClean: opening burp file = ', trim(inputFileName), ', ', trim(familyType)

    call burp_new(inputFile,         &
         filename = inputFileName,   &
         mode     = file_acc_append, &
         iostat   = error )
    call handle_error(error, "brpr_burpClean: problem opening BURP input file " // trim(inputFileName))

    call burp_get_property(inputFile, nrpts=numReports, io_unit= iun_burpin)
    nsize = mrfbfl(iun_burpin)
    if (debug) write(*,*) 'nsize= ', nsize

    ! determine if file contains grouped data
    groupedData = isGroupedData(inputFile, addresses)
    write(*,*) 'brpr_burpClean: Grouped data = ', groupedData

    ! get list of element ids used for checking flags
    call getElementIdsRead(familyType, elementIdsRead)

    ! ignore some elements when checking the flags
    do elemIndex = 1, size(elementIdsRead)
      if ( (elementIdsRead(elemIndex) == 10194) ) then
        elementIdsRead(elemIndex) = -1
      end if
    end do
    write(*,*) 'brpr_burpClean: Flags checked for these element IDs: ', elementIdsRead(:)

    numRejectTotal = 0 ! summed over entire file

    ! loop on reports
    reports: do reportIndex = 1, numReports

      numReject = 0

      call burp_get_report(inputFile,          &
           report    = inputReport,            &
           ref       = addresses(reportIndex), &
           iostat    = error) 
      call handle_error(error, "brpr_burpClean: burp_get_report")
      
      call burp_get_property(inputReport, stnid=stnid)
      resumeReport = (stnid(1:2) == ">>")

      ! create and initialize a new report
      if (reportIndex == 1) then
        call burp_new(copyReport, alloc_space=nsize, iostat=error)
        call handle_error(error, "brpr_burpClean: burp_new copyReport creation")
      end if
      call burp_copy_header(to=copyReport, from=inputReport)
      call burp_init_report_write(inputFile, copyReport, iostat=error)
      call handle_error(error, "brpr_burpClean: burp_init_report_write inputFile copyReport")

      ! get the flags block, check if obs missing, and list of element IDs
      foundFlags = .false.
      foundObs = .false.
      refBlock = 0      
      blocks: do

        refBlock = burp_find_block(inputReport, &
             block       = inputBlock,          &
             search_from = refBlock,            &
             iostat      = error)
        call handle_error(error, "brpr_burpClean: burp_find_block #1")
        if (refBlock < 0) exit blocks
         
        call burp_get_property(inputBlock, &
             nele   = numElem,             &
             nval   = numLevels,           &
             nt     = numObsProfiles,      &
             btyp   = btyp,                & 
             iostat = error)
        call handle_error(error, "brpr_burpClean: burp_get_property #1")

        ! only check "multi" blocks when cleanLevels is true
        if (cleanLevels) then
          checkBlock = btest(btyp,13)
        else
          checkBlock = .true.
        end if
          
        if (isFlagBlock(familyType, btyp) .and. checkBlock) then
          foundFlags = .true.
          if (debug) write(*,*) 'Found a block with flags: ', reportIndex, familyType, btyp
          if (debug) write(*,*) 'numObsProfiles, numLevels, numElem = ', numObsProfiles, numLevels, numElem 
          if (allocated(flagValues)) deallocate(flagValues)
          allocate(flagValues(numElem,numLevels,numObsProfiles))
          if (allocated(elementIdsBlock)) deallocate(elementIdsBlock)
          allocate(elementIdsBlock(numElem))
          do elemIndex = 1, numElem
            elementIdsBlock(elemIndex) = burp_get_element(inputBlock, index=elemIndex)
            do obsProfIndex = 1, numObsProfiles
              do levelIndex = 1, numLevels
                flagValues(elemIndex,levelIndex,obsProfIndex) = burp_get_tblval(inputBlock, &
                     nele_ind=elemIndex, nval_ind=levelIndex, nt_ind=obsProfIndex)
              end do
            end do
          end do
        end if

        if (isObsBlock(familyType, btyp) .and. checkBlock) then
          foundObs = .true.
          if (debug) write(*,*) 'Found a block with obs: ', reportIndex, familyType, btyp
          if (debug) write(*,*) 'numObsProfiles, numLevels, numElem = ', numObsProfiles, numLevels, numElem 
          if (allocated(obsValues)) deallocate(obsValues)
          allocate(obsValues(numElem,numLevels,numObsProfiles))
          do elemIndex = 1, numElem
            do obsProfIndex = 1, numObsProfiles
              do levelIndex = 1, numLevels
                obsValues(elemIndex,levelIndex,obsProfIndex) = burp_get_rval(inputBlock, &
                     nele_ind=elemIndex, nval_ind=levelIndex, nt_ind=obsProfIndex)
              end do
            end do
          end do
        end if

        if (foundFlags .and. foundObs) exit blocks

      end do blocks

      if (.not.foundFlags .or. .not.foundObs) then

        if (debug) write(*,*) 'No flag or obs block found for this report ', reportIndex

      else

        if (debug) write(*,*) 'numObsProfiles, numLevels, numElem = ', numObsProfiles, numLevels, numElem 
        ! determine which observations to keep
        if (allocated(rejectObs)) deallocate(rejectObs)
        allocate(rejectObs(numLevels,numObsProfiles))
        rejectObs(:,:) = .true.
        obsProfiles: do obsProfIndex = 1, numObsProfiles
          obsLevels: do levelIndex = 1, numLevels
            elements: do elemIndex = 1, numElem
              ! skip this element if it is not normally read
              if ( all(elementIdsBlock(elemIndex) /= (200000 + elementIdsRead(:))) ) cycle elements

              ! if at least one element in profile is 'good', then cannot reject
              if ( (.not.btest(flagValues(elemIndex,levelIndex,obsProfIndex),11)) .and.  &
                   (obsValues(elemIndex,levelIndex,obsProfIndex) /= missingValue) ) then
                if (debug) write(*,*) 'found a GOOD    observation: ', levelIndex, obsProfIndex,  &
                     elementIdsBlock(elemIndex), flagValues(elemIndex,levelIndex,obsProfIndex),   &
                     obsValues(elemIndex,levelIndex,obsProfIndex)
                rejectObs(levelIndex,obsProfIndex) = .false.
              else if (obsValues(elemIndex,levelIndex,obsProfIndex) == missingValue) then
                if (debug) write(*,*) 'found a MISSING observation: ', levelIndex, obsProfIndex,  &
                     elementIdsBlock(elemIndex), flagValues(elemIndex,levelIndex,obsProfIndex),  &
                     obsValues(elemIndex,levelIndex,obsProfIndex)
              else if (btest(flagValues(elemIndex,levelIndex,obsProfIndex),11)) then
                if (debug) write(*,*) 'found a BAD     observation: ', levelIndex, obsProfIndex,  &
                     elementIdsBlock(elemIndex), flagValues(elemIndex,levelIndex,obsProfIndex),  &
                     obsValues(elemIndex,levelIndex,obsProfIndex)
              end if

            end do elements
            if (cleanLevels) then
              ! count number of individual rejected levels
              if (rejectObs(levelIndex,obsProfIndex)) numReject = numReject + 1
            end if
          end do obsLevels
          if (debug) write(*,*) 'rejectObs = ',obsProfIndex,rejectObs(1,obsProfIndex)
          if (.not. cleanLevels) then
            ! count number of rejected complete profiles (all levels must be rejected)
            if (all(rejectObs(:,obsProfIndex))) numReject = numReject + 1
          end if
        end do obsProfiles

        numRejectTotal = numRejectTotal + numReject

      end if
        
      ! copy reduced blocks output report
      emptyReport = .false.
      refBlock = 0      
      blocks2: do

        refBlock = burp_find_block(inputReport, &
             block       = inputBlock2,         &
             search_from = refBlock,            &
             convert     = .false.,             &
             iostat      = error)
        call handle_error(error, "brpr_burpClean: burp_find_block #2")
        if (refBlock < 0) exit blocks2
         
        call burp_get_property(inputBlock2, &
             nele   = numElem2,             &
             nval   = numLevels2,           &
             nt     = numObsProfiles2,      &
             btyp   = btyp,                 &
             datyp  = datyp,                &
             iostat = error)
        call handle_error(error, "brpr_burpClean: burp_get_property #2")

        ! only modify "multi" blocks when cleanLevels is true
        if (cleanLevels) then
          checkBlock = btest(btyp,13)
        else
          checkBlock = .true.
        end if
          
        if (checkBlock) then

          if (debug) write(*,*) 'btyp, datyp = ', btyp, datyp
          if (cleanLevels) then
            newNumLevels = numLevels2 - numReject
            if (debug) write(*,*) 'ReportIndex = ', reportIndex
            if (debug) write(*,*) 'Reducing the number of levels from ', numLevels2, ' to ', newNumLevels

            if (newNumLevels >= 1) then

              ! shuffle the data
              do obsProfIndex = 1, numObsProfiles2
                levelIndexGood = 0
                do levelIndex = 1, numLevels2
                  if (rejectObs(levelIndex,obsProfIndex)) cycle
                  levelIndexGood = levelIndexGood + 1
                  do elemIndex = 1, numElem2
                    intBurpValue = burp_get_tblval(inputBlock2, nele_ind=elemIndex,  &
                         nval_ind=levelIndex, nt_ind=obsProfIndex, iostat=error)
                    call handle_error(error, "brpr_burpClean: burp_get_tblval")
                    call burp_set_tblval(inputBlock2, nele_ind=elemIndex,  &
                         nval_ind=levelIndexGood, nt_ind=obsProfIndex, tblval=intBurpValue, iostat=error)
                    call handle_error(error, "brpr_burpClean: burp_set_tblval") 
                    if (debug .and. levelIndex /= levelIndexGood) then
                      write(*,*) 'shuffling data: ', elemIndex, levelIndex, levelIndexGood, intBurpValue, reportIndex
                    end if
                  end do
                end do
              end do

              ! reduce the size of the block
              call burp_reduce_block(inputBlock2, new_nval=newNumLevels, iostat=error)
              call handle_error(error, "brpr_burpClean: burp_reduce_block")

            else ! newNumLevels < 1

              if (debug) write(*,*) 'All observation levels rejected for this report: ', reportIndex
              emptyReport = .true.

            end if

          else

            newNumObsProfiles   = numObsProfiles2 - numReject
            if (debug) write(*,*) 'ReportIndex = ', reportIndex
            if (debug) write(*,*) 'Reducing the number of observation profiles from ', numObsProfiles2, ' to ', newNumObsProfiles

            if (newNumObsProfiles >= 1) then

              ! shuffle the data
              obsProfIndexGood = 0
              do obsProfIndex = 1, numObsProfiles2
                if (all(rejectObs(:,obsProfIndex))) cycle ! all levels rejected
                obsProfIndexGood = obsProfIndexGood + 1
                do elemIndex = 1, numElem2
                  do levelIndex = 1, numLevels2
                    intBurpValue = burp_get_tblval(inputBlock2, nele_ind=elemIndex,  &
                         nval_ind=levelIndex, nt_ind=obsProfIndex, iostat=error)
                    call handle_error(error, "brpr_burpClean: burp_get_tblval #2")
                    call burp_set_tblval(inputBlock2, nele_ind=elemIndex,  &
                         nval_ind=levelIndex, nt_ind=obsProfIndexGood, tblval=intBurpValue, iostat=error)
                    call handle_error(error, "brpr_burpClean: burp_set_tblval #2") 
                    if (debug .and. obsProfIndex /= obsProfIndexGood) then
                      write(*,*) 'shuffling data: ', obsProfIndex, obsProfIndexGood, intBurpValue, reportIndex
                    end if
                  end do
                end do
              end do

              ! reduce the size of the block
              call burp_reduce_block(inputBlock2, new_nt=newNumObsProfiles, iostat=error)
              call handle_error(error, "brpr_burpClean: burp_reduce_block #2")

            else ! newNumObsProfiles < 1

              if (debug) write(*,*) 'All observation profiles rejected for this report: ', reportIndex
              emptyReport = .true.

            end if

          end if ! cleanLevels

        end if ! checkBlock

        if (.not. emptyReport) then
          call burp_write_block(copyReport, block=inputBlock2,  &
               convert_block=.false., iostat=error)
          call handle_error(error, "brpr_burpClean: burp_write_block")
        end if

      end do blocks2
      
      ! delete existing report and write new report into file
      call burp_delete_report(inputFile, inputReport, iostat=error)
      call handle_error(error, "brpr_burpClean: burp_delete_report")
      if (.not. emptyReport) then
        ! for grouped data modify "elev" to new number of obs profiles
        if (groupedData .and. .not.resumeReport) then
          call burp_set_property(copyReport,             &
                                 elev=newNumObsProfiles, &
                                 iostat=error)
          call handle_error(error, "brpr_burpClean: burp_set_property")
        end if
        call burp_write_report(inputFile, copyReport, iostat=error)
        call handle_error(error, "brpr_burpClean: burp_write_report")
      end if

    end do reports

    write(*,*) 'brpr_burpClean: finished - total number of obs profiles cleaned:', numRejectTotal
    write(*,*)

    call cleanup()

  contains

    !-------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(11)

      errors(:) = 0
      if (allocated(addresses))       deallocate(addresses, stat=errors(1))
      if (allocated(obsValues))       deallocate(obsValues, stat=errors(2))
      if (allocated(flagValues))      deallocate(flagValues, stat=errors(3))
      if (allocated(rejectObs))       deallocate(rejectObs, stat=errors(4))
      if (allocated(elementIdsBlock)) deallocate(elementIdsBlock, stat=errors(5))
      if (allocated(elementIdsRead))  deallocate(elementIdsRead, stat=errors(6))
      call burp_free(inputFile, iostat=errors(7))
      call burp_free(inputReport, iostat=errors(8))
      call burp_free(copyReport, iostat=errors(9))
      call burp_free(inputBlock, iostat=errors(10))
      call burp_free(inputBlock2, iostat=errors(11))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "brpr_burpClean: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !-------- handle_error -----
    subroutine handle_error(icode, errorMessage)
      implicit none
      character (len=*) :: errorMessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errorMessage))
      end if
    end subroutine handle_error
   
  end subroutine brpr_burpClean


  subroutine getBurpReportAddresses(fileName, addresses)
    !
    !:Purpose: Initial scan of file to get number of reports. Store address
    !          of each report in array addresses(numReports).
    !
    implicit none 

    ! Arguments:
    character(len=*), intent(in)        :: fileName
    integer, allocatable, intent(inout) :: addresses(:)

    ! Locals:
    type(burp_file) :: burpFile
    type(BURP_RPT)  :: report
    integer         :: numReports, refReport, error

    ! initialisation
    call burp_init(burpFile, iostat=error)
    call handle_error(error, "getBurpReportAddresses: burp_init inputfile")
    call burp_init(report)

    ! ouverture du fichier burp
    call burp_new(burpFile, filename=fileName, mode=file_acc_read, iostat=error)
    call handle_error(error, "getBurpReportAddresses: problem opening " // trim(fileName) )

    ! number of reports and maximum report size from BURP file
    call burp_get_property(burpFile, nrpts=numReports)
    if (numReports <= 1) then
      write(*,*) 'getBurpReportAddresses: BURP file ', trim(fileName)
      write(*,*) 'getBurpReportAddresses: WARNING: no observations in file'
    end if
    write(*,*)
    write(*,*) 'getBurpReportAddresses: Total number of reports = ', numReports
    write(*,*)

    if (allocated(addresses)) deallocate(addresses)
    allocate(addresses(numReports))
  
    addresses(:) = 0
    refReport = 0
    numReports = 0
    do
      refReport = burp_find_report(burpFile, report=report, search_from=refReport, iostat=error)
      call handle_error(error, "getBurpReportAddresses: burp_find_report error finding next burp report")
      if (refReport < 0) exit
      numReports = numReports+1
      addresses(numReports) = refReport
    end do

    call cleanup()

  contains

    !-------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(2)

      errors(:) = 0 
      call burp_free(report, iostat=errors(1))
      call burp_free(burpFile, iostat=errors(2))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "getBurpReportAddresses: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !-------- handle_error -----
    subroutine handle_error(icode, errorMessage)
      implicit none
      character (len=*) :: errorMessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errorMessage))
      end if
    end subroutine handle_error

  end subroutine getBurpReportAddresses


  function isGroupedData(burpFile,address) result(isGrouped)
    implicit none
    ! Arguments:
    type(burp_file) :: burpFile
    integer :: address(:)
    logical :: isGrouped

    ! Locals:
    type(burp_rpt)   :: report
    type(burp_block) :: block
    integer          :: reportIndex, refBlock, btyp, btyp10, error

    call burp_init(report)
    call burp_init(block)

    ! determine if this file contains grouped data
    isGrouped = .false.
    reports: do reportIndex = 1, size(address)
    
      call burp_get_report(burpFile,         &
           report    = report,          &
           ref       = address(reportIndex), &
           iostat    = error)      
      call handle_error(error, "isGroupedData: burp_get_report error getting next burp report")
      ! loop on blocks
      refBlock = 0      
      blocks: do

        refBlock = burp_find_block(report, &
             block       = block,         &
             search_from = refBlock,            &
             iostat      = error)
        call handle_error(error, "isGroupedData: burp_find_block error finding next burp block")
        if (refBlock < 0) exit blocks

        call burp_get_property(block, btyp=btyp)
        btyp10 = ishft(btyp,-5)
        if ( btyp10 == 160 ) then
          isGrouped = .true.
          exit reports
        end if

      end do blocks

    end do reports

    call cleanup()

  contains

    !-------- cleanup -----
    subroutine cleanup()
      implicit none
      integer :: errors(2)

      errors(:) = 0
      call burp_free(report, iostat=errors(1))
      call burp_free(block, iostat=errors(2))
      !Should we abort here ?
      if (any(errors /= 0)) write(*,*) "isGroupedData: error while deallocating memory: ", errors(:)
    end subroutine cleanup

    !-------- handle_error -----
    subroutine handle_error(icode, errorMessage)
      implicit none
      character (len=*) :: errorMessage
      integer :: icode
      if ( icode /= burp_noerr ) then
        write(*,*) 'error code', icode
        write(*,*) BURP_STR_ERROR()
        write(*,*) "history"
        call BURP_STR_ERROR_HISTORY()
        call cleanup()
        call utl_abort(trim(errorMessage))
      end if
    end subroutine handle_error

  end function isGroupedData
  

  function isFlagBlock(familyType, btyp) result(isFlag)
    implicit none

    ! Arguments:
    character(len=*) :: familyType
    integer :: btyp
    logical :: isFlag

    ! Locals:
    integer :: btyp10, btyp10flg, offset

    isFlag = .false.

    btyp10 = ishft(btyp,-5)
    btyp10flg = 483

    if (trim(familyType)=='UA') then
      offset = 288
      isFlag = (btyp10 == btyp10flg .or. btyp10 == btyp10flg-offset)
      if ( .not.isFlag ) then
        isFlag = (btyp10 == btyp10flg-2 .or. btyp10 == btyp10flg-offset-2)
      end if
    else
      select case(trim(familyType))
      case('AI','SW','SC')
        offset = 256
      case('RO','TO')
        offset = 0
      case('SF','GP')
        offset = 288
      end select
      isFlag = (btyp10 == btyp10flg-offset .or. btyp10 == btyp10flg-offset-2)
    end if
    
  end function isFlagBlock


  function isObsBlock(familyType, btyp) result(isObs)
    implicit none

    ! Arguments:
    character(len=*) :: familyType
    integer :: btyp
    logical :: isObs

    ! Locals:
    integer :: btyp10, btyp10obs, offset

    isObs = .false.

    btyp10 = ishft(btyp,-5)
    btyp10obs = 291

    if (trim(familyType)=='UA') then
      offset = 288
      isObs = (btyp10 == btyp10obs .or. btyp10 + offset == btyp10obs)
      if ( .not.isObs ) then
        isObs = (btyp10 == btyp10obs-2 .or. btyp10 == btyp10obs-offset-2)
      end if
    else
      select case(trim(familyType))
      case('AI','SW','SC')
        offset = 256
      case('RO','TO')
        offset = 0
      case('SF','GP')
        offset = 288
      end select
      isObs = (btyp10 == btyp10obs-offset .or. btyp10 == btyp10obs-offset-2)
    end if
    
  end function isObsBlock

  
  subroutine getElementIdsRead(familyType, elementIds)
    implicit none

    ! Arguments:
    character(len=*) :: familyType
    integer, allocatable :: elementIds(:)

    ! Locals:
    integer :: elementIndex, elementCount

    if (allocated(elementIds)) deallocate(elementIds)

    select case(trim(familyType))

    case('UA','AI','AL','SW','SC','PR','RO')
      call brpacma_nml('namburp_conv', beSilent_opt=.true.)
      allocate(elementIds(nelems))
      elementIds(:) = blistelements(1:nelems)

    case('SF')
      call brpacma_nml('namburp_sfc', beSilent_opt=.true.)
      allocate(elementIds(nelems_sfc))
      elementIds(:) = blistelements_sfc(1:nelems_sfc)

    case('GP')
      call brpacma_nml('namburp_sfc', beSilent_opt=.true.)
      ! do not include "formal error", since it was removed from obsSpaceData
      if (any(liste_ele_gps(1:nelems_gps) == bufr_nefe)) then
        allocate(elementIds(nelems_gps-1))
        elementCount = 0
        do elementIndex = 1, nelems_gps
          if (liste_ele_gps(elementIndex) /= bufr_nefe) then
            elementCount = elementCount + 1
            elementIds(elementCount) = liste_ele_gps(elementIndex)
          end if
        end do
      else
        allocate(elementIds(nelems_gps))
        elementIds(:) = liste_ele_gps(1:nelems_gps)
      end if
      
    case('TO')
      call brpacma_nml('namburp_tovs', beSilent_opt=.true.)
      allocate(elementIds(nelems))
      elementIds(:) = blistelements(1:nelems)

    case default
      call utl_abort('getElementIdsRead: unknown familyType: ' // trim(familyType))

    end select

  end subroutine getElementIdsRead

end module burpread_mod
