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
!! MODULE burpRead (prefix='brpr' category='6. Observation input/output')
!!
!! *Purpose*: To read and update BURP observation files. Data is stored in 
!!            obsSpaceData object.
!!
!--------------------------------------------------------------------------
module  burpread_mod

use codePrecision_mod
use bufr_mod
use burp_module
use ObsSpaceData_mod
use MathPhysConstants_mod
use earthconstants_mod
use utilities_mod
use obsUtil_mod

implicit none
save

private

! public procedures
public :: brpr_readBurp, brpr_updateBurp, brpr_getTypeResume


! MODULE CONSTANTS ...
!
!   These variables are set during object initialization (variational_init) and
!   are not changed thereafter.
! bits to verify in Quality Control Flag:


INTEGER*4              :: NELEMS,NELEMS_SFC,BLISTELEMENTS(20),BLISTELEMENTS_SFC(20)
INTEGER*4              :: BN_ITEMS
CHARACTER *3           :: BITEMLIST(20)
CHARACTER *7           :: TYPE_RESUME

INTEGER*4              :: BNBITSOFF,BNBITSON,BBITOFF(15),BBITON(15)
LOGICAL                :: ENFORCE_CLASSIC_SONDES,UA_HIGH_PRECISION_TT_ES,READ_QI_GA_MT_SW


CONTAINS

  character(len=7) function brpr_getTypeResume
    brpr_getTypeResume=TYPE_RESUME
  end function brpr_getTypeResume

  SUBROUTINE brpr_updateBurp(obsdat,familytype,brp_file,FILENUMB)

    !******************************************************************************** 
    !
    !**ID brpr_updateBurp -- UPDATE VARIABLES RELATIVE TO ASSIMILATION IN BURP FILES
    !
    !       AUTHOR:   P. KOCLAS (CMDA/SMC) Feb 2013
    !
    !       REVISION:
    !                 S. MACPHERSON (ARMA) Oct 2013
    !                  -- add 'GP' family (ground-based GPS) as type "SFC"
    !                 P. KOCLAS (CMDA/SMC)) Aug 2014
    !                  -- Fix for the Burp 24 bit header markers of the UA4D
    !                 P. KOCLAS (CMDA/SMC)) Oct 2014
    !                  -- CHANGE ADDSIZE FOR SW ( GROUPED SATWINDS)
    !                 Ping Du (CMDA), Nov/Dec 2014
    !                  -- Added CASE('CH') option with intitial settings.
    !                  -- Added use of FAMILYTYPE2 and the 'namburp_chm_sfc' and 
    !                     'namburp_chm' namelists for the CH family 
    !                 Y.J. Rochon (ARQI), Dec 2014
    !                  -- Added 08046 (constituent code table) in LISTE_INFO.
    !                     Corresponds to new ObsSpaceData integer-header element OBS_CHM.
    !                  -- Added 15008,15011,15012,15020,15024,15026-29,15199,
    !                     and 15230 under LISTE_ELE for CASE('CH'). This
    !                     denotes the/a possible set of elements.
    !                  -- Modified 'vcord_type' to an array of dimension 10.
    !                  -- Extended vcord_type from 1 to 8 elements for CH
    !                  -- Setting of ELEVFACT to 1.0 for 7006 coordinate.
    !                 Y.J. Rochon (ARQI), May 2015
    !                  -- Made use of the existing WINDS logical to determine 
    !                     when the wind-related fields are to be included in
    !                     the reports with uni-level data blocks.
    !                 P. KOCLAS (CMDA/SMC)) May 2015
    !                  -- CHANGES  for GROUPED scatterometres or surface obs
    !                 P. KOCLAS (CMDA/SMC)) Jan 2016
    !                  -- bug fix for GROUPED scatterometres or surface obs(element doubling)
    !                 Y.J. Rochon (ARQI), March 2016
    !                  -- Account for pre-existing FGE and OER blocks. See LBLOCK_OER_CP
    !                     ,NBLOCK_OER_CP, LBLOCK_FGE_CP and BLOKC_FGE_CP. The OER block
    !                     may contain additional data such as averaging kernel matrices 
    !                     and obs error correlation matrices at least for the CH family.
    !                     As well, the FGE values are not calculated at the analysis phase
    !                     so would not be available if requested as output.
    !                 M. Sitwell (ARQI), Aug 2016
    !                  -- Modified update of resume record so that it is written once instead
    !                     of twice.
    !
    !       OBJECT: UPDATE CMC BURP FILE 
    !
    !          WHEN SEARCHING FOR A SPECIFIC BLOCK BY ITS BTYP, VALUES OF
    !          BIT 0 TO 3 ARE IRRELEVANT WHILE BIT 4 IS 0 FOR GLOBAL AND 1
    !          FOR REGIONAL MODEL. HERE, WE SEARCH BLOCK BY THEIR FIRST
    !          10 BITS (BIT 5 TO 14).
    !
    !       ARGUMENTS:
    !               INPUT:
    !                  -OBSDAT    : instance of obsspace_data module object
    !                  -FAMILYTYPE: TYPE of family ('UA' ,'SF','AI,''SW','TO')...
    !                  -BRP_FILE  : FILENAME OF BURP FILE 
    !
    !
    !****************************************************************************

    IMPLICIT NONE
    INTEGER,  PARAMETER    :: NBLOC_LIST=6
    CHARACTER(LEN=128)     :: BRP_FILE
    CHARACTER *2           :: FAMILYTYPE
    INTEGER                :: FILENUMB,LNMX
    type (struct_obs), intent(inout) :: obsdat

    TYPE(BURP_FILE)        :: FILE_IN
    TYPE(BURP_RPT)         :: RPT_IN,CP_RPT
    TYPE(BURP_BLOCK)       :: BLOCK_IN,BLOCK_OMA,BLOCK_OMP,BLOCK_OER,BLOCK_FGE,BLOCK_FLG,BLOCK_FSO
    TYPE(BURP_BLOCK)       :: BLOCK_OMA_SFC,BLOCK_OMP_SFC,BLOCK_OER_SFC,BLOCK_FGE_SFC,BLOCK_FLG_SFC,BLOCK_FSO_SFC
    TYPE(BURP_BLOCK)       :: Block_FLG_CP,BLOCK_OBS_MUL_CP,BLOCK_MAR_MUL_CP,BLOCK_OBS_SFC_CP,BLOCK_MAR_SFC_CP

    CHARACTER(LEN=5)       :: FAMILYTYPE2
    CHARACTER(LEN=9)       :: OPT_MISSING
    INTEGER                :: BTYP,BFAM,BTYP10,BTYP10_uni,BTYP10FLG_uni,BTYP10obs_uni 
    INTEGER                :: BTYP10DES,BTYP10INF,BTYP10OBS,BTYP10FLG

    INTEGER                :: NB_RPTS,REF_RPT,REF_BLK,COUNT
    INTEGER, ALLOCATABLE   :: ADDRESS(:)
    REAL                   :: VCOORD

    INTEGER                :: NBELE,NVALE,NTE
    INTEGER                :: I,J,JJ,K,KK,KI,IL,Jo,ERROR,OBSN,KOBSN,ITEM,IER
    INTEGER                :: info_elepos,IND_ELE,IND_VCOORD,IND_QCFLAG
    INTEGER                :: IND_ELE_MAR,IND_ELEU,IND_ELEF,IND_ELE_stat,IND_ELE_tth,IND_ELE_esh
    INTEGER                :: IND_LAT,IND_LON,IND_TIME

    INTEGER                :: vcord_type(10),FLAG,SUM
    REAL                   :: RELEV,ELEVFACT
    REAL                   :: XLAT,XLON,XTIME
    INTEGER                :: status ,idtyp,lati,long,dx,dy,elev, &
                              drnd,date_h,hhmm_h,oars,runn,YMD_DATE,HM
    INTEGER                :: IND055200,IND5002, IND6002,IND4208,IND4197

    INTEGER                :: iele,NELE,NELE_SFC,NVAL,NT,NELE_INFO,LN
    INTEGER                :: bit_alt,btyp_offset,btyp_offset_uni
    INTEGER                :: BKNAT,BKTYP,BKSTP
    character(len = 5)     :: BURP_TYP
    CHARACTER(LEN=9)       :: STNID,STN_RESUME,STID
    LOGICAL                :: HIRES,HIPCS
    INTEGER                :: NDATA,NDATA_SF
    INTEGER                :: IFLAG,BITSflagoff

    INTEGER                :: OBS_START,SAVE_OBS,OBS_HIRES,ASSIM
    INTEGER                :: IL_INDEX,IRLN,INLV,LK,VNM
    REAL                   :: PPP,OBS,OMA,OMP,OER,FSO,FGE,OBSVA,CONVFACT
    INTEGER                :: FLG,TIME,ILEMU,ILEMV,ILEMD,VCOORD_POS

    INTEGER                :: BLOCK_LIST(NBLOC_LIST),bl

    INTEGER                :: new_bktyp,post_bit,STATUS_HIRES,BIT_STATUS,FILEN
    LOGICAL                :: REGRUP,WINDS,OMA_SFC_EXIST,OMA_ALT_EXIST

    INTEGER                :: LISTE_INFO(18),LISTE_ELE(20),LISTE_ELE_SFC(20),is_in_list
    INTEGER                :: ADDSIZE,SIZE_DATA_BK

    LOGICAL                :: LBLOCK_OER_CP, LBLOCK_FGE_CP
    TYPE(BURP_BLOCK)       :: BLOCK_OER_CP, BLOCK_FGE_CP
    logical                :: FSOFound

    DATA LISTE_INFO  &
       /1007,002019,007024,007025 ,005021, 005022, 008012, 013039,020010,2048,2022,33060, &
        33062,33039,10035,10036,08046,5043/

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
    ELEVFACT=0.
    BNBITSOFF=0
    BNBITSON=0
    ENFORCE_CLASSIC_SONDES=.false.
    UA_HIGH_PRECISION_TT_ES=.false.
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
        CALL BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'UA'
        LISTE_ELE(1:5) = (/12001,11001,11002,12192,10194/)
        NELE=5
        ENFORCE_CLASSIC_SONDES=.false.
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=10000
      CASE('AI')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:4) = (/12001,12192,11001,11002/)
        NELE=4
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=5000
      CASE('AL')
        BURP_TYP='uni'
        vcord_type(1)=7071

        LISTE_ELE(1:4) = (/10004,40030,12001,5021/)
        NELE=4
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=5000
      CASE('SW')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        vcord_type(2)=-1
!       ADDSIZE=5000
!       ' AMVS-> ADDSIZE 600000 oct 2014 pik'
        ADDSIZE=600000
      CASE('SF','GP')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELE_SFC=7
        LISTE_ELE_SFC(1:7) = (/12004,11011,11012,10051,10004,12203,15031/)

        CALL BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'SFC'
        WINDS=.TRUE.
        IF (trim(FAMILYTYPE) == 'GP') WINDS=.FALSE.
        ILEMU=11215
        ILEMV=11216
        ILEMD=11011
        ADDSIZE=5000
      CASE('SC')
        BURP_TYP='uni'
        LISTE_ELE_SFC(1:2) = (/11012,11011/)
        NELE=2
        CALL BRPACMA_NML('namburp_conv')
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

        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.TRUE.
        ADDSIZE=10000
      CASE('RO')
        BURP_TYP='multi'
        vcord_type(1)=7007
        vcord_type(2)=7040

        LISTE_ELE(1:2) = (/15036,15037/)
        NELE=2

        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        WINDS=.FALSE.
        !================GPS-RO CANNOT BE FILTERED=======
        BNBITSOFF=0
        BNBITSON=0
        !================GPS-RO CANNOT BE FILTERED=======
        NELE_INFO=16
      CASE('GO','MI','TO')
        BURP_TYP='multi'
        vcord_type(1)=5042
        vcord_type(2)=2150
        
        LISTE_ELE(1:1) = (/12163/)
        NELE=1

        CALL BRPACMA_NML('namburp_tovs')
        NELE=NELEMS

        NELE_INFO=16
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
        CALL BRPACMA_NML('namburp_chm_sfc')
        NELE_SFC=NELEMS_SFC
        WINDS=.FALSE.
        ADDSIZE=10000

        FAMILYTYPE2='CH'
        LISTE_ELE(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                                15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE=19
        CALL BRPACMA_NML('namburp_chm')
        NELE=NELEMS
    END SELECT
    LISTE_ELE    (1:NELE    )=BLISTELEMENTS(1:NELE)
    LISTE_ELE_SFC(1:NELE_SFC)=BLISTELEMENTS_SFC(1:NELE_SFC)
    if(NELE     > 0)write(*,*)  ' LISTE_ELE =',LISTE_ELE
    if(NELE_SFC > 0)write(*,*)  ' LISTE_ELE_SFC =',LISTE_ELE_SFC
    if(BNBITSOFF > 0 .or. BNBITSON > 0)write(*,*)  ' BNBITSON BNBITSOFF SIZE OF CP_RPT   =',BNBITSON,BNBITSOFF,LNMX*8


    TYPE_RESUME='POSTALT'
    BN_ITEMS=1
    BITEMLIST(1)='OMA'
    CALL BRPACMA_NML('namburp_update')
    WRITE(*,*)  ' BN_ITEMS   =',BN_ITEMS
    WRITE(*,'(a12,x)') ' ITEMS TO ADD IN BURP FILE REPORTS =', BITEMLIST(1:BN_ITEMS)
    WRITE(*,'(x,a9)' ) ' BTYP OF UPDATED BURP FILE=', TYPE_RESUME

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


    WRITE(*,*) '----------------------------------------------------'
    WRITE(*,*) '-----------     BEGIN brpr_updateBurp   ------------'
    WRITE(*,*) 'FAMILYTYPE   =',FAMILYTYPE
    WRITE(*,*) 'BURP_TYP btyp_offset    =',BURP_TYP, btyp_offset
    WRITE(*,*) 'BURP_TYP btyp_offset_uni=',BURP_TYP, btyp_offset_uni
    WRITE(*,*) '----------------------------------------------------'


    ! initialisation

    SUM=0
    opt_missing = 'MISSING'


    Call BURP_Set_Options( &
       & REAL_OPTNAME       = opt_missing, &
       & REAL_OPTNAME_VALUE = MPC_missingValue_R4, &
       & CHAR_OPTNAME       = 'MSGLVL', &
       & CHAR_OPTNAME_VALUE = 'FATAL', &
       & IOSTAT             = error )

    Call BURP_Init(File_in      ,IOSTAT=error)
    !Call BURP_Init(Rpt_in,IOSTAT=error)
    Call BURP_Init(Rpt_in,CP_RPT,IOSTAT=error)
    Call BURP_Init(Block_in     ,IOSTAT=error)

    Call BURP_Init(BLOCK_OMA    ,IOSTAT=error)
    Call BURP_Init(BLOCK_OMP    ,IOSTAT=error)
    Call BURP_Init(BLOCK_OER    ,IOSTAT=error)
    Call BURP_Init(BLOCK_FGE    ,IOSTAT=error)
    Call BURP_Init(BLOCK_FSO    ,IOSTAT=error)
    Call BURP_Init(BLOCK_OMA_SFC,IOSTAT=error)
    Call BURP_Init(BLOCK_OMP_SFC,IOSTAT=error)
    Call BURP_Init(BLOCK_OER_SFC,IOSTAT=error)
    Call BURP_Init(BLOCK_FGE_SFC,IOSTAT=error)
    Call BURP_Init(BLOCK_FSO_SFC,IOSTAT=error)
    Call BURP_Init(BLOCK_FLG_SFC,IOSTAT=error)

    Call BURP_Init(BLOCK_FLG    ,IOSTAT=error)
    Call BURP_Init(Block_FLG_CP ,IOSTAT=error)

    Call BURP_Init(BLOCK_OBS_MUL_CP ,IOSTAT=error)
    Call BURP_Init(BLOCK_MAR_MUL_CP ,IOSTAT=error)

    Call BURP_Init(BLOCK_OBS_SFC_CP ,IOSTAT=error)
    Call BURP_Init(BLOCK_MAR_SFC_CP ,IOSTAT=error)

    ! opening file
    ! ------------
    write(*,*) 'OPENING BURP FILE FOR UPDATE = ', trim(brp_file)

    Call BURP_New(File_in, FILENAME = brp_file, &
                   & MODE    = FILE_ACC_APPEND, &
                   & IOSTAT  = error )

    ! obtain input burp file number of reports
    ! ----------------------------------------
    Call BURP_Get_Property(File_in, NRPTS=nb_rpts)
    Call BURP_Init(Rpt_in,IOSTAT=error)

    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'IOSTAT    =',error
    WRITE(*,*) 'NUMBER OF REPORTS IN FILE  = ',nb_rpts
    WRITE(*,*) '-----------------------------------------'

    ! scan input burp file to get all reports address
    ! -----------------------------------------------

    Allocate(address(nb_rpts))
    address(:) = 0
    count      = 0
    ref_rpt    = 0
    bit_alt    = 0

    do
      ref_rpt = BURP_Find_Report(File_in, &
                 & REPORT      = Rpt_in, &
                 & SEARCH_FROM = ref_rpt, &
                 & IOSTAT      = error)
      Call burp_get_property(Rpt_in, STNID = stnid )
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
      Call BURP_New(Cp_rpt, ALLOC_SPACE=12*LNMX, IOSTAT=error)

      ! LOOP ON REPORTS
      REPORTS: do kk = 1, count
            

        Call BURP_Get_Report(File_in, &
                & REPORT    = Rpt_in, &
                & REF       = address(kk), &
                & IOSTAT    = error) 
        Call burp_get_property(Rpt_in, &
               STNID = stnid ,TEMPS =hhmm_h,FLGS = status ,IDTYP =idtyp,LATI = lati &
               ,LONG = long ,DX = dx ,DY = dy,ELEV=elev,DRND =drnd,DATE =date_h &
               ,OARS =oars,RUNN=runn ,IOSTAT=error)

        IF ( stnid(1:2) == ">>" ) THEN
          write(*,*)  ' RESUME RECORD POSITION IN BURP FILE =',stnid,kk
          Call BURP_Copy_Header(TO=Cp_rpt,FROM=Rpt_in)
          Call BURP_Init_Report_Write(File_in,Cp_Rpt, IOSTAT=error)
          Call BURP_Set_Property(Cp_Rpt,STNID =">>"//TYPE_RESUME)
          Call BURP_Delete_Report(File_in,Rpt_in, IOSTAT=error)
          Call BURP_Write_Report(File_in,Cp_rpt, IOSTAT=error)
          cycle REPORTS
        ELSE
          !write(*,*) ' UPDATING STN IN BURP FILE =',  TRIM(FAMILYTYPE),KK,stnid,lati,LONG,dx,DY,elev,idtyp
        END IF
        Call BURP_Copy_Header(TO=Cp_rpt,FROM=Rpt_in)
        Call BURP_Init_Report_Write(File_in,Cp_Rpt, IOSTAT=error)

        ! FIRST LOOP ON BLOCKS

        ref_blk = 0

        LBLOCK_OER_CP=.false.
        LBLOCK_FGE_CP=.false.
        HIRES=.FALSE.
        HIPCS=.FALSE.
        REGRUP=.false.
        NDATA_SF=-1
        !WRITE(*,*)'  record number =',kk,' obs_start =',obs_start
        BLOCK_LIST(1:6)=-1
        BLOCKS0: do
          ref_blk = BURP_Find_Block(Rpt_in, &
                     & BLOCK       = Block_in, &
                     & SEARCH_FROM = ref_blk, &
                     & IOSTAT      = error)

          if (ref_blk < 0) EXIT BLOCKS0

          Call BURP_Get_Property(Block_in, &
                                 & NELE   = nbele, &
                                 & NVAL   = nvale, &
                                 & NT     = nte, &
                                 & BFAM   = bfam, &
                                 & BTYP   = btyp, &
                                 & BKTYP  = bktyp, &
                                 & BKNAT  = BKNAT, &
                                 & BKSTP  = BKSTP, &
                                 & IOSTAT = error)
          if(trim(familytype) == 'AL')then

            ! Fudge the block type, because the data are simulated
            if(btyp == 1024) then
              btyp=1152
            else if(btyp == 7168)then
              btyp=7296
            end if
          end if

          btyp10    = ishft(btyp,-5)
          if ( btyp10 - BTYP10des == 0 ) then
            Block_FLG_CP=BLOCK_IN
            BLOCK_LIST(1)=BTYP
            REGRUP=.TRUE.
          elseif ( btyp10 - btyp10obs_uni == 0 .and. bkstp <= 4 ) then
            BLOCK_OBS_SFC_CP=BLOCK_IN
            BLOCK_LIST(2)=BTYP
            NDATA_SF=0
          elseif ( btyp10 - btyp10flg_uni == 0 .and. bkstp <= 4) then
            BLOCK_MAR_SFC_CP=BLOCK_IN
            BLOCK_LIST(3)=BTYP
          elseif ( btyp10 - btyp10obs == 0 .and. bfam == 0 ) then
            BLOCK_LIST(4)=BTYP
            BLOCK_OBS_MUL_CP=BLOCK_IN
          elseif ( btyp10 - btyp10flg == 0 ) then
            BLOCK_LIST(5)=BTYP
            BLOCK_MAR_MUL_CP=BLOCK_IN
          elseif ( (btyp10 - btyp10inf == 0) .or. (btyp10 - btyp10inf == 1) ) then
            BLOCK_LIST(6)=BTYP
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

          if (ref_blk < 0) cycle BLOCKS1

          Call BURP_Get_Property(Block_in, &
                           & NELE   = nbele, &
                           & NVAL   = nvale, &
                           & NT     = nte, &
                           & BFAM   = bfam, &
                           & BTYP   = btyp, &
                           & BKTYP  = bktyp, &
                           & BKNAT  = BKNAT, &
                           & IOSTAT = error)
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
                Call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BKTYP =new_bktyp, BKSTP=0)
                Call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BKTYP =new_bktyp, BKSTP=0)
              else
                Call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BKTYP =new_bktyp)
                Call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BKTYP =new_bktyp)
              end if

            end if

            il_index=0
             
            !Call BURP_Delete_BLOCK(Rpt_in,BLOCK=Block_in)

            OMA_SFC_EXIST=.true.

            if (.not.WINDS) then

               Call BURP_New(BLOCK_OMA_SFC,NELE =NBELE,NVAL=nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_OMP_SFC,NELE =NBELE,NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_OER_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_FGE_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                             ,IOSTAT = error) 
              Call BURP_New(BLOCK_FSO_SFC, NELE =NBELE, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2 &
                             ,IOSTAT = error)
            else

               IND_eleu  = BURP_Find_Element(Block_in, ELEMENT=11215, IOSTAT=error)
               IND_elef  = BURP_Find_Element(Block_in, ELEMENT=11011, IOSTAT=error)

               Call BURP_New(BLOCK_OMA_SFC,NELE =NBELE+2,NVAL=nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_OMP_SFC,NELE =NBELE+2,NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_OER_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_FGE_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                             ,IOSTAT = error) 
               Call BURP_New(BLOCK_FSO_SFC, NELE =NBELE+2, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2 &
                             ,IOSTAT = error)
               ILEMU = 11215
               ILEMV = 11216
               call BURP_Set_Element( BLOCK_OMA_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error) 
               call BURP_Set_Element( BLOCK_OMA_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error) 

               call BURP_Set_Element( BLOCK_OMP_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error) 
               call BURP_Set_Element( BLOCK_OMP_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error) 
   
               call BURP_Set_Element( BLOCK_OER_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error) 
               call BURP_Set_Element( BLOCK_OER_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error) 

               call BURP_Set_Element( BLOCK_FGE_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error) 
               call BURP_Set_Element( BLOCK_FGE_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error) 

               call BURP_Set_Element( BLOCK_FSO_SFC,NELE_IND = 1,ElEMENT=ILEMU,IOSTAT=error)
               call BURP_Set_Element( BLOCK_FSO_SFC,NELE_IND = 2,ElEMENT=ILEMV,IOSTAT=error)

               do k =1,nte
                 Call BURP_Set_Rval(Block_OMA_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_OMA_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_OMP_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_OMP_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_OER_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_OER_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_FGE_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_FGE_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 Call BURP_Set_Rval(Block_FSO_SFC,NELE_IND =1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 )
                 Call BURP_Set_Rval(Block_FSO_SFC,NELE_IND =2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 )
               end do
               il_index=2

               IF (IND_eleu < 0 .and. IND_elef > 0) THEN
                 Call BURP_RESIZE_BLOCK(BLOCK_OBS_SFC_CP,ADD_NELE = 2 ,IOSTAT=error) 
                 call BURP_Set_Element( BLOCK_OBS_SFC_CP,NELE_IND = nbele+1,ElEMENT=ILEMU,IOSTAT=error) 
                 call BURP_Set_Element( BLOCK_OBS_SFC_CP,NELE_IND = nbele+2,ElEMENT=ILEMV,IOSTAT=error) 
                 do k =1,nte
                   Call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =nbele+1,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                   Call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =nbele+2,NVAL_IND =1,NT_IND = k,RVAL = MPC_missingValue_R4 ) 
                 end do

                 Call BURP_RESIZE_BLOCK(BLOCK_MAR_SFC_CP,ADD_NELE = 2 ,IOSTAT=error) 
                 call BURP_Set_Element( BLOCK_MAR_SFC_CP,NELE_IND = nbele+1,ElEMENT=ILEMU+200000,IOSTAT=error) 
                 call BURP_Set_Element( BLOCK_MAR_SFC_CP,NELE_IND = nbele+2,ElEMENT=ILEMV+200000,IOSTAT=error) 
                 do k =1,nte
                   Call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =nbele+1,NVAL_IND =1,NT_IND = k,tblval = 0 ) 
                   Call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =nbele+2,NVAL_IND =1,NT_IND = k,tblval = 0 ) 
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

                 IND_ELE_MAR= BURP_Find_Element(Block_MAR_SFC_CP, ELEMENT=iele+200000, IOSTAT=error)
                 if (IND_ele_mar <= 0 ) cycle

                 IND_ele  = BURP_Find_Element(BLOCK_OBS_SFC_CP, ELEMENT=iele, IOSTAT=error)

                 if ( k == 1 ) then
                    if( OMA_SFC_EXIST .eqv.  .true.   ) then
                      if (iele /= ILEMU .and. iele /= ILEMV) then
                        il_index=il_index +1
                        call BURP_Set_Element (BLOCK_OMA_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                        Call BURP_Set_Element (BLOCK_OMP_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                        Call BURP_Set_Element (BLOCK_OER_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                        Call BURP_Set_Element (BLOCK_FGE_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                        Call BURP_Set_Element (BLOCK_FSO_SFC,NELE_IND= il_index,ElEMENT=iele,IOSTAT=error)
                      end if
                    end if
                 end if

                 is_in_list=-1
                 is_in_list=FIND_INDEX(LISTE_ELE_SFC,iele)
                 if (is_in_list < 0   .and. iele  /= ILEMU .and. iele /= ILEMV) cycle ELEMS_SFC
                    IND_ele_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele, IOSTAT=error)
                    Call BURP_Set_Rval(Block_OMA_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                    Call BURP_Set_Rval(Block_OMP_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                    Call BURP_Set_Rval(Block_OER_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4)
                    Call BURP_Set_Rval(Block_FGE_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                    Call BURP_Set_Rval(Block_FSO_SFC, NELE_IND =IND_ELE_stat ,NVAL_IND =1 , NT_IND  = k , RVAL = MPC_missingValue_R4) 

                    IFLAG =  BURP_Get_Tblval(Block_MAR_SFC_CP,NELE_IND = IND_ele_mar,NVAL_IND = 1, NT_IND = k)
                    OBSVA =  BURP_Get_Rval  (Block_OBS_SFC_CP,NELE_IND = IND_ele    ,NVAL_IND = 1, NT_IND = k)
                    if (OBSVA == MPC_missingValue_R4 .and. iele /= ILEMU  .and. iele /= ILEMV ) cycle
                    if(iand(IFLAG,BITSflagoff) /= 0) cycle

                    if (OBSN > obs_numHeader(obsdat) ) write(*,*)  ' debordement  surface OBS_START=',OBS_START
                    if (OBSN > obs_numHeader(obsdat))  cycle

                    IRLN=obs_headElem_i(obsdat,OBS_RLN,OBSN )
                    INLV=obs_headElem_i(obsdat,OBS_NLV,OBSN )

                    IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele, IOSTAT=error)
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
                        FLG=obs_bodyElem_i(obsdat,OBS_FLG ,LK)
                        KOBSN= KOBSN + 1
                        SUM=SUM +1
                        Call BURP_Set_Rval( Block_OER_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OER   ) 
                        Call BURP_Set_Rval( Block_FGE_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FGE  ) 
                        Call BURP_Set_Rval( Block_FSO_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FSO  )

                        IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA_SFC, ELEMENT=iele, IOSTAT=error)
                        Call BURP_Set_Rval( Block_OMA_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OMA)

                        IND_ELE_stat  = BURP_Find_Element(BLOCK_OMP_SFC, ELEMENT=iele, IOSTAT=error)
                        Call BURP_Set_Rval( Block_OMP_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = OMP)

                        IND_ELE_stat  = BURP_Find_Element(BLOCK_FSO_SFC, ELEMENT=iele, IOSTAT=error)
                        Call BURP_Set_Rval( Block_FSO_SFC, NELE_IND =IND_ele_stat ,NVAL_IND =1,NT_IND = k , RVAL = FSO)
   
                        Call BURP_Set_tblval(Block_MAR_SFC_CP,NELE_IND =IND_ele_mar,NVAL_IND =1,NT_IND = k ,TBLVAL= FLG) 
   
                        OBS=obs_bodyElem_r(obsdat,OBS_VAR,LK)
                        IND_ele  = BURP_Find_Element(BLOCK_OBS_SFC_CP, ELEMENT=iele, IOSTAT=error)
                        Call BURP_Set_Rval(Block_OBS_SFC_CP,NELE_IND =IND_ele,NVAL_IND =1,NT_IND = k,RVAL = OBS ) 

                        exit
                      end if

                    end do  OBSDATA
!-------------------------------------------
               end do  elems_sfc
!-------------------------------------------

            if ( REGRUP .and. KOBSN > 0   )  then
                STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START )
                STATUS=IBSET(STATUS,BIT_STATUS)
                ind055200 = BURP_Find_Element(Block_FLG_CP, ELEMENT=055200, IOSTAT=error)
                Call BURP_Set_tblval( Block_FLG_CP, NELE_IND =ind055200,NVAL_IND =1,NT_IND = k ,TBLVAL = STATUS ) 
                OBSN=OBSN +1
                OBS_START=OBS_START +1
            end if
!==============================================
            end do  REGRUP_SFC
!=============================================

            do item=1,BN_ITEMS

              if ( BITEMLIST(item) == 'OMA') then
                Call BURP_Reduce_Block(BLOCK_OMA_SFC, NEW_NELE =il_index )
                Call BURP_Write_Block( CP_RPT, BLOCK_OMA_SFC,&
                ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                cycle
              end if
              if ( BITEMLIST(item) == 'OMP') then
                Call BURP_Reduce_Block(BLOCK_OMP_SFC, NEW_NELE =il_index )
                Call BURP_Write_Block( CP_RPT, BLOCK_OMP_SFC,&
                ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                cycle
              end if
              if ( BITEMLIST(item) == 'OER') then
                if (.not.LBLOCK_OER_CP) then
                   Call BURP_Reduce_Block(BLOCK_OER_SFC, NEW_NELE =il_index )
                   Call BURP_Write_Block( CP_RPT, BLOCK_OER_SFC,&
                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                else
                   Call BURP_Set_Property(BLOCK_OER_CP ,BKTYP =new_bktyp)
                   Call BURP_Write_Block( CP_RPT, BLOCK_OER_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)                 
                end if
                cycle
              end if
              if ( BITEMLIST(item) == 'FGE') then
                if (.not.LBLOCK_FGE_CP) then
                   Call BURP_Reduce_Block(BLOCK_FGE_SFC, NEW_NELE =il_index )
                   Call BURP_Write_Block( CP_RPT, BLOCK_FGE_SFC,&
                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                else
                   Call BURP_Set_Property(BLOCK_FGE_CP ,BKTYP =new_bktyp)
                   Call BURP_Write_Block( CP_RPT, BLOCK_FGE_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error)                 
                end if
                cycle
              end if

              if ( BITEMLIST(item) == 'FSO') then
                   Call BURP_Reduce_Block(BLOCK_FSO_SFC, NEW_NELE =il_index )
                   Call BURP_Write_Block( CP_RPT, BLOCK_FSO_SFC,&
                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
                cycle
              end if

            end do

            Call BURP_Set_Property(BLOCK_OBS_SFC_CP ,BFAM =0)
            Call BURP_Set_Property(BLOCK_MAR_SFC_CP ,BFAM =0)

            Call BURP_Write_Block( CP_RPT, BLOCK_OBS_SFC_CP,&
                                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
            Call BURP_Write_Block( CP_RPT, BLOCK_MAR_SFC_CP,&
                                   ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 
            IF ( KOBSN > 0 .and. .not. REGRUP ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
              STATUS=IBSET(STATUS,BIT_STATUS)
              Call BURP_Set_Property(CP_RPT ,FLGS =STATUS)
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
              Call BURP_Set_Property(BLOCK_OBS_MUL_CP ,BKTYP =new_bktyp)
              Call BURP_Set_Property(BLOCK_MAR_MUL_CP ,BKTYP =new_bktyp)
            end if
            OBSN=OBS_START
            NVAL=NVALE ;  NT=NTE

            il_index=1
            Call BURP_New(BLOCK_OMA, NELE =1, NVAL =nvale,NT=NTE,bfam=12,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                  ,IOSTAT = error) 
            Call BURP_New(BLOCK_OMP, NELE =1, NVAL =nvale,NT=NTE,bfam=14,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=10 &
                  ,IOSTAT = error) 
            Call BURP_New(BLOCK_OER, NELE =1, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=14 &
                  ,IOSTAT = error) 
            Call BURP_New(BLOCK_FGE, NELE =1, NVAL =nvale,NT=NTE,bfam=10,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=15 &
                  ,IOSTAT = error) 
            Call BURP_New(BLOCK_FSO, NELE =1, NVAL =nvale,NT=NTE,bfam=1,BKNAT=BKNAT,BKTYP=new_bktyp,BKSTP=2  &
                  ,IOSTAT = error)
            
            VCOORD_POS=0
            k=0
            IND_VCOORD=-1
            do while (vcord_type(k+1) /= -1 .and. IND_VCOORD == -1)
               k=k+1 
               IND_VCOORD  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=vcord_type(k),IOSTAT=error)
            end do
            IF ( IND_VCOORD > 0 ) then
              IF (trim(FAMILYTYPE) == trim('CH')) THEN
                 ELEVFACT=0.0
                 IF (vcord_type(k) == 7006) ELEVFACT=1.0
              END IF
              call BURP_Set_Element(BLOCK_OMA,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error) 
              call BURP_Set_Element(BLOCK_OMP,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error) 
              call BURP_Set_Element(BLOCK_OER,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error) 
              call BURP_Set_Element(BLOCK_FGE,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error) 
              call BURP_Set_Element(BLOCK_FSO,NELE_IND= 1,ElEMENT=vcord_type(k),IOSTAT=error)
              VCOORD_POS=1
            ELSE IF (IND_VCOORD == -1) then
              !WRITE(*,*) '  PAS DE COORDONNEE VERTICALE famille ',trim(FAMILYTYPE)
              il_index=0
            end if
            VCOORD = -999.

            IND_eleu  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=ILEMU, IOSTAT=error)
            IND_elef  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=ILEMD, IOSTAT=error)

            OMA_ALT_EXIST=.false.
            if(WINDS .and. IND_eleu < 0 .and. IND_elef > 0) then
              Call BURP_RESIZE_BLOCK(BLOCK_OMA,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OMA,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OMA,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error) 
              OMA_ALT_EXIST=.true.

              Call BURP_RESIZE_BLOCK(BLOCK_OMP,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OMP,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OMP,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error) 

              Call BURP_RESIZE_BLOCK(BLOCK_OER,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OER,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OER,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error) 

              Call BURP_RESIZE_BLOCK(BLOCK_FGE,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_FGE,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_FGE,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error) 

              Call BURP_RESIZE_BLOCK(BLOCK_FSO,ADD_NELE = 2 ,IOSTAT=error)
              call BURP_Set_Element( BLOCK_FSO,NELE_IND = VCOORD_POS+1,ElEMENT=ILEMU,IOSTAT=error)
              call BURP_Set_Element( BLOCK_FSO,NELE_IND = VCOORD_POS+2,ElEMENT=ILEMV,IOSTAT=error)

              Call BURP_RESIZE_BLOCK(BLOCK_OBS_MUL_CP,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OBS_MUL_CP,NELE_IND = nbele+1,ElEMENT=ILEMU,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_OBS_MUL_CP,NELE_IND = nbele+2,ElEMENT=ILEMV,IOSTAT=error) 

              Call BURP_RESIZE_BLOCK(BLOCK_MAR_MUL_CP,ADD_NELE = 2 ,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_MAR_MUL_CP,NELE_IND = nbele+1,ElEMENT=ILEMU+200000,IOSTAT=error) 
              call BURP_Set_Element( BLOCK_MAR_MUL_CP,NELE_IND = nbele+2,ElEMENT=ILEMV+200000,IOSTAT=error) 

              do k=1,nte
                do jj=1,nvale
                  Call BURP_Set_Rval( Block_OMA,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_OMA,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_OMP,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_OMP,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_OER,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_OER,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_FGE,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_FGE,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  ) 
                  Call BURP_Set_Rval( Block_FSO,  NELE_IND =VCOORD_POS+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  )
                  Call BURP_Set_Rval( Block_FSO,  NELE_IND =VCOORD_POS+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4  )

                  Call BURP_Set_Rval(  BLOCK_OBS_MUL_CP, NELE_IND =nbele+1 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                  Call BURP_Set_Rval(  BLOCK_OBS_MUL_CP, NELE_IND =nbele+2 ,NVAL_IND =jj , NT_IND  = k , RVAL = MPC_missingValue_R4) 
                  Call BURP_Set_tblval(BLOCK_MAR_MUL_CP, NELE_IND =nbele+1 ,NVAL_IND =jj , NT_IND  = k , TBLVAL = 0  ) 
                  Call BURP_Set_tblval(BLOCK_MAR_MUL_CP, NELE_IND =nbele+2 ,NVAL_IND =jj , NT_IND  = k , TBLVAL = 0  ) 
                end do
              end do
              nbele=nbele+2

              il_index=il_index+2
            end if
            !Call BURP_Delete_BLOCK(Rpt_in,BLOCK=Block_in)

            ! LAT LON TIME IN DATA BLOCK
            IND_LAT   = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=5001, IOSTAT=error)
            IND_LON   = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=6001, IOSTAT=error)
            IND_TIME  = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=4015, IOSTAT=error)
            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES=.true.
            if(ENFORCE_CLASSIC_SONDES .eqv. .true.) hires=.false.

            if( TRIM(FAMILYTYPE2) == 'UA' .and. UA_HIGH_PRECISION_TT_ES .eqv. .true. ) HIPCS=.true.

            !print * , ' hires =true ? ndata_sf ',stnid,hires,NDATA_SF

            if ( HIRES .AND. NDATA_SF  > 0 ) OBS_START =OBS_START +1
            OBSN=OBS_START
            STATUS_HIRES=0

            regrup_LOOP: do k=1,nte
              KOBSN=0

              levels: do j=1,nvale

                !pikpik
                if(HIRES) KOBSN=0
                !pikpik

                elems: do IL = 1, NBELE  

                  iele=-1
                  iele=BURP_Get_Element(BLOCK_OBS_MUL_CP,INDEX =il,IOSTAT= error) 

                  IND_ELE_MAR= BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=iele+200000, IOSTAT=error)
                  if (IND_ele_mar <  0 ) cycle

                  IND_ele       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele, IOSTAT=error)
                  if (IND_ele == IND_LAT  .and. hires ) cycle
                  if (IND_ele == IND_LON  .and. hires ) cycle
                  if (IND_ele == IND_TIME .and. hires ) cycle
                  IND_ELE_STAT=-1
                  IND_ele_STAT  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele, IOSTAT=error)

                  if(j == 1 .and. il /= ind_vcoord .and. IND_ELE_STAT < 1 ) then
                  il_index=il_index +1
                  Call BURP_RESIZE_BLOCK(BLOCK_OMA,ADD_NELE = 1 ,IOSTAT=error) 
                  call BURP_Set_Element (BLOCK_OMA,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)

                  Call BURP_RESIZE_BLOCK(BLOCK_OMP,ADD_NELE = 1 ,IOSTAT=error) 
                  call BURP_Set_Element (BLOCK_OMP,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)

                  Call BURP_RESIZE_BLOCK(BLOCK_OER,ADD_NELE = 1 ,IOSTAT=error) 
                  call BURP_Set_Element (BLOCK_OER,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)

                  Call BURP_RESIZE_BLOCK(BLOCK_FGE,ADD_NELE = 1 ,IOSTAT=error) 
                  call BURP_Set_Element (BLOCK_FGE,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)

                  Call BURP_RESIZE_BLOCK(BLOCK_FSO,ADD_NELE = 1 ,IOSTAT=error)
                  call BURP_Set_Element (BLOCK_FSO,NELE_IND = il_index,ElEMENT=iele,IOSTAT=error)

                  do ki=1,nte
                    do jj=1,nvale
                      Call BURP_Set_Rval( Block_OMA, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                      Call BURP_Set_Rval( Block_OMP, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                      Call BURP_Set_Rval( Block_OER, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                      Call BURP_Set_Rval( Block_FGE, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 ) 
                      Call BURP_Set_Rval( Block_FSO, NELE_IND =il_index ,NVAL_IND =jj , NT_IND  = ki, RVAL = MPC_missingValue_R4 )
                    end do
                  end do
                end if
                VCOORD =  BURP_Get_Rval(BLOCK_OBS_MUL_CP, &
                                      &   NELE_IND = IND_VCOORD, &
                                      &   NVAL_IND = j, &
                                      &   NT_IND   = k)
                IF (il == IND_VCOORD) THEN
                  Call BURP_Set_Rval( Block_OMA, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                  Call BURP_Set_Rval( Block_OMP, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                  Call BURP_Set_Rval( Block_OER, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                  Call BURP_Set_Rval( Block_FGE, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  ) 
                  Call BURP_Set_Rval( Block_FSO, NELE_IND =1 ,NVAL_IND =j , NT_IND  = k , RVAL = VCOORD  )
                END IF

                IND_ele       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele, IOSTAT=error)
                !if ( kk < obs_headElem_i(obsdat,OBS_IDO,OBSN) ) print *, ' kk OBS_IDO cycle =',kk, obs_headElem_i(obsdat,OBS_IDO,OBSN)
                !if ( kk < obs_headElem_i(obsdat,OBS_IDO,OBSN) ) cycle elems

                is_in_list=-1
                is_in_list=FIND_INDEX(LISTE_ELE,iele)
                if (is_in_list < 0   .and. iele  /= ILEMU .and. iele /= ILEMV)cycle

                  IFLAG =  BURP_Get_Tblval(Block_MAR_MUL_CP,NELE_IND = IND_ELE_MAR,NVAL_IND = J, NT_IND = k)
                  OBSVA =  BURP_Get_Rval  (Block_OBS_MUL_CP,NELE_IND = IND_ele    ,NVAL_IND = J, NT_IND = k)
                  if(iand(iflag,BITSflagoff) /= 0) CYCLE ELEMS     

                  OBSVA =  BURP_Get_Rval  (Block_OBS_MUL_CP,NELE_IND = IND_ele    ,NVAL_IND = J, NT_IND = k)
                  !if( idtyp == 168 ) print * , ' bobossmi avant vcoord obsva    stnid =',  VCOORD,OBSVA,stnid

                  if (VCOORD == MPC_missingValue_R4 ) CYCLE ELEMS
                  if (OBSVA == MPC_missingValue_R4 .and. iele /= ILEMU  .and. iele /= ILEMV ) CYCLE ELEMS      

                  if (OBSN > obs_numHeader(obsdat)) write(*,*) ' debordement  altitude OBSN=',OBSN
                  if (OBSN > obs_numHeader(obsdat))  cycle
                  IRLN=obs_headElem_i(obsdat,OBS_RLN,OBSN)
                  INLV=obs_headElem_i(obsdat,OBS_NLV,OBSN)
                  TIME=obs_headElem_i(obsdat,OBS_ETM,OBSN)
                  STID=obs_elem_c(obsdat,'STID',OBSN)
                  if ( STID /= stnid ) cycle

                  IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele,  IOSTAT=error)

                  if(HIPCS) then
                    IND_ELE_tth   = BURP_Find_Element(BLOCK_OMA, ELEMENT=12101, IOSTAT=error)
                    IND_ELE_esh   = BURP_Find_Element(BLOCK_OMA, ELEMENT=12239, IOSTAT=error)
                  endif

                  IND_ELE       = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=iele, IOSTAT=error)
                  IND_eleu      = BURP_Find_Element(BLOCK_OBS_MUL_CP, ELEMENT=ILEMU, IOSTAT=error)
                  convfact=1.
                  if (iele == 10194) convfact=1./RG
                  do LK=IRLN,IRLN+INLV-1
                    ASSIM=obs_bodyElem_i(obsdat,OBS_ASS,LK)
                    VNM  =obs_bodyElem_i(obsdat,OBS_VNM,LK)
                    PPP  =obs_bodyElem_r(obsdat,OBS_PPP,LK) - (ELEV-400.)*ELEVFACT
                    if( abs( VCOORD -  PPP) < .01  .and.  VNM == iele  ) then 
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
                      FLG=obs_bodyElem_i(obsdat,OBS_FLG,LK)
                      KOBSN= KOBSN + 1
                      IND_ELE_stat  = BURP_Find_Element(BLOCK_OMA, ELEMENT=iele, IOSTAT=error)
                      if  ( OMA /= MPC_missingValue_R4 ) then
                        OMA=OMA*convfact
                      end if
                      Call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )

                      if(HIPCS) then
                        if(iele == 12001) Call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )
                        if(iele == 12192) Call BURP_Set_Rval(Block_OMA,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OMA )
                      endif

                      !if(trim(familytype) == 'TO' )print *,' bingo  stnid kk vnm ppp flg omp ',stnid,kk,vnm,ppp,flg,omp,oma
                      SUM=SUM +1
                      IND_ELE_stat  = BURP_Find_Element(BLOCK_OMP, ELEMENT=iele, IOSTAT=error)
                      if  ( OMP /= MPC_missingValue_R4 ) then
                        OMP=OMP*convfact
                      end if
                      Call BURP_Set_Rval(  Block_OMP,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OMP)

                      if(HIPCS) then
                        if(iele == 12001) Call BURP_Set_Rval(Block_OMP,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OMP )
                        if(iele == 12192) Call BURP_Set_Rval(Block_OMP,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OMP )
                      endif

                      Call BURP_Set_Rval(  Block_OER,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = OER  ) 

                      if(HIPCS) then
                        if(iele == 12001) Call BURP_Set_Rval(Block_OER,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = OER )
                        if(iele == 12192) Call BURP_Set_Rval(Block_OER,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = OER )
                      endif

                      Call BURP_Set_Rval(  Block_FGE,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = FGE ) 

                      if(HIPCS) then
                        if(iele == 12001) Call BURP_Set_Rval(Block_FGE,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = FGE )
                        if(iele == 12192) Call BURP_Set_Rval(Block_FGE,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = FGE )
                      endif

                      Call BURP_Set_Rval(  Block_FSO,  NELE_IND =IND_ele_stat ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )

                      if(HIPCS) then
                        if(iele == 12001) Call BURP_Set_Rval(Block_FSO,  NELE_IND =IND_ele_tth ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )
                        if(iele == 12192) Call BURP_Set_Rval(Block_FSO,  NELE_IND =IND_ele_esh ,NVAL_IND =j , NT_IND  = k , RVAL = FSO )
                      endif

                      IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=iele+200000, IOSTAT=error)

                      Call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )

                      if(HIPCS) then
                        if(iele == 12001) then
                          IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=212101, IOSTAT=error)
                          Call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )
                        endif
                        if(iele == 12192) then
                          IND_ele_mar  = BURP_Find_Element(Block_MAR_MUL_CP, ELEMENT=212239, IOSTAT=error)
                          Call BURP_Set_tblval(Block_MAR_MUL_CP,NELE_IND =IND_ELE_MAR ,NVAL_IND =j  , NT_IND  = k,TBLVAL = FLG )
                        endif
                      endif

                      IND_ele  = BURP_Find_Element(Block_OBS_MUL_CP, ELEMENT=iele, IOSTAT=error)
                      Call BURP_Set_Rval(Block_OBS_MUL_CP,NELE_IND =IND_ele,NVAL_IND =j,NT_IND = k,RVAL = OBS) 

                      EXIT

                    end if
                  end do 

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
                ind055200 = BURP_Find_Element(Block_FLG_CP, ELEMENT=055200, IOSTAT=error)
                Call BURP_Set_tblval( Block_FLG_CP, NELE_IND =ind055200,NVAL_IND =1,NT_IND = k ,TBLVAL = STATUS ) 
                OBSN=OBSN +1
                OBS_START=OBS_START +1
              end if

            end do regrup_LOOP

            IF (HIRES .and. KOBSN > 0 .and. .not. regrup ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
!pik 8-2014   STATUS=IBSET(STATUS,BIT_STATUS)
              STATUS_HIRES=IBSET(STATUS_HIRES,BIT_STATUS)
              Call BURP_Set_Property(CP_RPT ,FLGS =STATUS_HIRES)
            END IF
            IF (HIRES )OBS_START=OBSN
            IF (HIRES )SAVE_OBS=OBS_START
            IF (REGRUP)OBS_START=OBSN

            IF ( .not. HIRES .and. .not. regrup .and. KOBSN > 0 ) THEN
              STATUS=obs_headElem_i(obsdat,OBS_ST1,OBS_START)
              STATUS=IBSET(STATUS,BIT_STATUS)
              OBS_START=OBSN +1
              OBSN=OBSN +1
              Call BURP_Set_Property(CP_RPT ,FLGS =STATUS)
            END IF
            IF ( .not. HIRES .and. .not. regrup .and. KOBSN == 0 ) THEN
              write(*,*) ' KOBSN=0  stnid=',stnid
            END IF
            IF ( .not. HIRES .and.  regrup .and. KOBSN == 0 ) THEN
              write(*,*)' KOBSN=0 regrup stnid=',kk,stnid,obsn,obs_numHeader(obsdat)
            END IF

            IF (REGRUP) SAVE_OBS=OBS_START

            Call BURP_Reduce_Block(BLOCK_OMA, NEW_NELE =il_index )
            Call BURP_Reduce_Block(BLOCK_OMP, NEW_NELE =il_index )
            Call BURP_Reduce_Block(BLOCK_OER, NEW_NELE =il_index )
            Call BURP_Reduce_Block(BLOCK_FGE, NEW_NELE =il_index )
            Call BURP_Reduce_Block(BLOCK_FSO, NEW_NELE =il_index )

            Call BURP_Set_Property(BLOCK_OBS_MUL_CP ,BFAM =0)
            Call BURP_Set_Property(BLOCK_MAR_MUL_CP ,BFAM =0)
            Call BURP_Write_Block( CP_RPT, BLOCK_OBS_MUL_CP,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
            Call BURP_Write_Block( CP_RPT, BLOCK_MAR_MUL_CP,&
                      ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 

            do item=1,BN_ITEMS

              if ( BITEMLIST(item) == 'OMA') then
                Call BURP_Write_Block( CP_RPT, BLOCK_OMA,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
              end if
              if ( BITEMLIST(item) == 'OMP') then
                Call BURP_Write_Block( CP_RPT, BLOCK_OMP,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
              end if
              if ( BITEMLIST(item) == 'OER') then
                if (.not.LBLOCK_OER_CP) then
                   Call BURP_Write_Block( CP_RPT, BLOCK_OER,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                else
                   Call BURP_Set_Property(BLOCK_OER_CP ,BKTYP =new_bktyp)
                   Call BURP_Write_Block( CP_RPT, BLOCK_OER_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 
                  end if
              end if
              if ( BITEMLIST(item) == 'FGE') then
                if (.not.LBLOCK_FGE_CP) then
                   Call BURP_Write_Block( CP_RPT, BLOCK_FGE,&
                          ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error) 
                else
                   Call BURP_Set_Property(BLOCK_FGE_CP ,BKTYP =new_bktyp)
                   Call BURP_Write_Block( CP_RPT, BLOCK_FGE_CP,&
                        ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 
                  end if
              end if
              if ( BITEMLIST(item) == 'FSO') then
                 Call BURP_Write_Block( CP_RPT, BLOCK_FSO,&
                        ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .TRUE., IOSTAT= error)
              end if

            end do

            if (regrup ) OBS_START = OBSN
            if( .not. hires ) SAVE_OBS = OBS_START

          end if  ! bl == 4

          if ( bl == 6 ) then
            Call BURP_Write_Block( CP_RPT, BLOCK_in, ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 
          end if

          ! descriptor block (btyp = 0010 100000X XXXX) 

          BTYP10des = 160

          !if ( BTYP10 - BTYP10des == 0 ) then
          if (  bl == 1 ) then
            OBS_START=SAVE_OBS
          end if

          !==================== IASI  SPECIAL BLOCK==================
          if ( (BTYP == 9217 .or. BTYP == 15361) .and.  IDTYP == 186   ) then
            Call BURP_Write_Block( CP_RPT, BLOCK_in, ENCODE_BLOCK = .FALSE., CONVERT_BLOCK = .FALSE., IOSTAT= error) 
          end if
          !==================== IASI  SPECIAL BLOCK==================

        end do BLOCKS1

        if (  REGRUP )  then
          Call BURP_Write_Block( CP_RPT, Block_FLG_CP, ENCODE_BLOCK = .TRUE., CONVERT_BLOCK = .FALSE.)
        end if

        Call BURP_Delete_Report(File_in,Rpt_in, IOSTAT=error)
        Call BURP_Write_Report(File_in,CP_RPT,IOSTAT= error) 

      end do REPORTS

    end if

    Deallocate(address)

    Call BURP_Free(Rpt_in,CP_RPT,IOSTAT=error)
    Call BURP_Free(Block_in,     IOSTAT=error)
    Call BURP_Free(Block_OMA,    IOSTAT=error)
    Call BURP_Free(Block_OMP,    IOSTAT=error)
    Call BURP_Free(Block_OER,    IOSTAT=error)
    Call BURP_Free(Block_FGE,    IOSTAT=error)
    Call BURP_Free(Block_FSO,    IOSTAT=error)
    Call BURP_Free(Block_OMA_SFC,IOSTAT=error)
    Call BURP_Free(Block_OMP_SFC,IOSTAT=error)
    Call BURP_Free(Block_OER_SFC,IOSTAT=error)
    Call BURP_Free(Block_FGE_SFC,IOSTAT=error)
    Call BURP_Free(Block_FSO_SFC,IOSTAT=error)
    Call BURP_Free(Block_FLG_SFC,IOSTAT=error)
    Call BURP_Free(Block_FLG    ,IOSTAT=error)
    Call BURP_Free(Block_FLG_CP ,IOSTAT=error)
    Call BURP_Free(Block_MAR_MUL_CP ,IOSTAT=error)
    Call BURP_Free(Block_MAR_SFC_CP ,IOSTAT=error)
    Call BURP_Free(Block_OBS_MUL_CP ,IOSTAT=error)
    Call BURP_Free(Block_OBS_SFC_CP ,IOSTAT=error)
    Call BURP_Free(File_in,      IOSTAT=error)

    write(*,*) ' BURPFILE  UPDATED SUM = ',trim(brp_file),SUM

  END SUBROUTINE brpr_updateBurp


  SUBROUTINE BRPACMA_NML(NML_SECTION)

    !******************************************************************************** 
    !
    !       REVISION:
    !                 Ping Du (CMDA), Nov 2014 
    !                  -- Added NAMBURP_FILTER_CHM and NAMBURP_FILTER_CHM_SFC namelists
    !                  -- Added CASE('namburp_chm') and CASE('namburp_chm_sfc')
    !
    !
    !****************************************************************************

    IMPLICIT NONE

    INTEGER*4          :: NULNAM,IER,FNOM,FCLOS
    CHARACTER *256     :: NAMFILE
    CHARACTER(len = *) :: NML_SECTION

    NAMELIST /NAMBURP_FILTER_CONV/NELEMS,    BLISTELEMENTS,    BNBITSOFF,BBITOFF,BNBITSON,BBITON,ENFORCE_CLASSIC_SONDES,UA_HIGH_PRECISION_TT_ES,READ_QI_GA_MT_SW
    NAMELIST /NAMBURP_FILTER_SFC/ NELEMS_SFC,BLISTELEMENTS_SFC,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_FILTER_TOVS/NELEMS,BLISTELEMENTS,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_FILTER_CHM_SFC/NELEMS_SFC,BLISTELEMENTS_SFC,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_FILTER_CHM/NELEMS,BLISTELEMENTS,BNBITSOFF,BBITOFF,BNBITSON,BBITON
    NAMELIST /NAMBURP_UPDATE/BN_ITEMS, BITEMLIST,TYPE_RESUME


    NAMFILE=trim("flnml")
    nulnam=0
    IER=FNOM(NULNAM,NAMFILE,'R/O',0)
    WRITE(*,*) ' READ NML_SECTION =',trim(NML_SECTION)

    SELECT CASE(trim(NML_SECTION))
      CASE( 'namburp_sfc')
        READ(NULNAM,NML=NAMBURP_FILTER_SFC)
        write(*,nml=NAMBURP_FILTER_SFC)
      CASE( 'namburp_conv')
        READ(NULNAM,NML=NAMBURP_FILTER_CONV)
        write(*,nml=NAMBURP_FILTER_CONV)
      CASE( 'namburp_tovs')
        READ(NULNAM,NML=NAMBURP_FILTER_TOVS)
        write(*,nml=NAMBURP_FILTER_TOVS)
      CASE( 'namburp_chm_sfc')
        READ(NULNAM,NML=NAMBURP_FILTER_CHM_SFC)
        write(*,nml=NAMBURP_FILTER_CHM_SFC)
      CASE( 'namburp_chm')
        READ(NULNAM,NML=NAMBURP_FILTER_CHM)
        write(*,nml=NAMBURP_FILTER_CHM)
      CASE( 'namburp_update')
        READ(NULNAM,NML=NAMBURP_UPDATE)
        write(*,nml=NAMBURP_UPDATE)
    END SELECT

    ier=FCLOS(NULNAM)
  END SUBROUTINE BRPACMA_NML


  SUBROUTINE brpr_readBurp(obsdat,familytype,brp_file,FILENUMB)

    !***********************************************************************
    !
    !**ID brpr_readBurp -- SELECT VARIABLES RELATIVE TO AIRS IN BURP FILE
    !
    !       AUTHOR:   P. KOCLAS (CMDA/SMC) May 2011
    !
    !       REVISION:
    !                 S. MACPHERSON (ARMA) Oct 2013
    !                  -- add 'GP' family (ground-based GPS)
    !                 P. KOCLAS (CMDA) Oct 2014
    !                  -- Changed VCOORD Dimensions 
    !                  -- add VCORD array  to allow grouped SW or AI
    !                 Ping Du (CMDA) Nov/Dec 2014 
    !                  -- Added initial CASE('CH') section with addition of 'namburp_chm_sfc'
    !                     and 'namburp_chm' namelists
    !                  -- Added UNI_FAMILYTYPE.
    !                  -- Replaced "NDATA_SF= WRITE_BODY(*" to use UNI_FAMILYTYPE
    !                  -- Added use of FAMILYTYPE2 for the CH family.
    !                 Y.J. Rochon (ARQI) Dec 2014 - May 2016
    !                  -- Added 08046 (constituent code table) in LISTE_INFO.
    !                     Corresponds to new ObsSpaceData integer-header element OBS_CHM.
    !                  -- Added 15008,15011,15012,15020,15024,15026-29,15199,
    !                     and 15230 under LISTE_ELE for CASE('CH') - even if not used.
    !                  -- Modified 'vcord_type' to an array of dimension 10.
    !                  -- Extended 'vcord_type' from 1 to 7 elements for TR
    !                  -- Addition/use of 'vcoord_type'. Added 'vcoord_type' 
    !                     argument to WRITE_BODY.
    !                  -- Added flag_passage*=1 for btyp10*_uni cases.
    !                     Needed when only *_uni data being read to avoid
    !                     incorrect output error message.
    !                  -- Added consideration of NDATA_SF for calling WRITE_INFO
    !                     accompanied by initialization of NDATA_SF and NDATA to zero.
    !
    !       OBJECT: READ CMC BURP FILE 
    !
    !          WHEN SEARCHING FOR A SPECIFIC BLOCK BY ITS BTYP, VALUES OF
    !          BIT 0 TO 3 ARE IRRELEVANT WHILE BIT 4 IS 0 FOR GLOBAL AND 1
    !          FOR REGIONAL MODEL. HERE, WE SEARCH BLOCK BY THEIR FIRST
    !          10 BITS (BIT 5 TO 14).
    !
    !       ARGUMENTS:
    !               INPUT:
    !                      -BRP_FILE  : NAME OF BURP FILE
    !
    !
    !***********************************************************************

    IMPLICIT NONE
    CHARACTER(LEN=128)     :: BRP_FILE
    INTEGER                :: FILENUMB
    CHARACTER *2           :: FAMILYTYPE
    CHARACTER *2           :: UNI_FAMILYTYPE
    
    type (struct_obs), intent(inout) :: obsdat
    
    TYPE(BURP_FILE)        :: FILE_IN
    TYPE(BURP_RPT)         :: RPT_IN
    TYPE(BURP_BLOCK)       :: BLOCK_IN
    
    CHARACTER(LEN=5)       :: FAMILYTYPE2
    CHARACTER(LEN=9)       :: OPT_MISSING
    INTEGER                :: BTYP,BFAM,BKSTP,BTYP10,BTYP10_uni,BTYP10FLG_uni,BTYP10obs_uni 
    INTEGER                :: BTYP10DES,BTYP10INF,BTYP10OBS,BTYP10FLG

    INTEGER                :: NB_RPTS,REF_RPT,REF_BLK,COUNT
    INTEGER, ALLOCATABLE   :: ADDRESS(:)

    REAL   , ALLOCATABLE   :: OBSVALUE(:,:,:),OBSVALUE_SFC(:,:,:)
    REAL   , ALLOCATABLE   :: OBSERV  (:,:),    OBSERV_SFC(:,:)

    INTEGER, ALLOCATABLE   :: QI1VAL(:) ,QI2VAL(:), QIVAL(:), MTVAL(:), LSVAL(:), HAVAL(:), GAVAL(:)
    INTEGER, ALLOCATABLE   :: azimuth(:)
    INTEGER, ALLOCATABLE   :: QCFLAG  (:,:,:),  QCFLAG_SFC(:,:,:)
    INTEGER, ALLOCATABLE   :: QCFLAGS (:,:),   QCFLAGS_SFC(:,:)

    REAL   , ALLOCATABLE   :: VCOORD  (:,:),  VCOORD_SFC(:)
    REAL   , ALLOCATABLE   :: VCORD   (:)

    INTEGER, ALLOCATABLE   :: LAT(:),LON(:),HHMM(:),DATE(:),GLBFLAG(:)
    REAL   , ALLOCATABLE   :: HLAT(:,:),  HLON(:,:),  HTIME(:,:)
    REAL   , ALLOCATABLE   :: HLAT_SFC(:),HLON_SFC(:),HTIME_SFC(:)

    REAL   , ALLOCATABLE   :: RINFO(:,:)
    REAL   , ALLOCATABLE   :: TRINFO(:)

    REAL   , ALLOCATABLE   :: EMIS(:,:),SURF_EMIS(:)

    INTEGER, ALLOCATABLE   :: CFRAC(:,:)

    REAL(OBS_REAL), ALLOCATABLE :: RADMOY(:,:,:)
    REAL(OBS_REAL), ALLOCATABLE :: radstd(:,:,:)
    INTEGER                :: LISTE_INFO(18),LISTE_ELE(20),LISTE_ELE_SFC(20)

    INTEGER                :: NBELE,NVALE,NTE
    INTEGER                :: I,J,JJ,K,KK,KL,IL,ERROR,OBSN
    INTEGER                :: info_elepos,IND_ELE,IND_VCOORD,IND_QCFLAG,IND_SW
    INTEGER                :: IND055200,IND4208,ind4197,IND5002,IND6002,ind_al
    INTEGER                :: IND_LAT,IND_LON,IND_TIME,IND_EMIS
    INTEGER                :: FLAG_PASSAGE1,FLAG_PASSAGE2,FLAG_PASSAGE3,FLAG_PASSAGE4

    INTEGER                :: vcord_type(10),FLAG,SUM,vcoord_type
    REAL(OBS_REAL)         :: RELEV,XLAT,XLON,RELEV2
    REAL                   :: XTIME,SECONDS
    INTEGER                :: status ,idtyp,lati,long,dx,dy,elev, &
                               drnd,date_h,hhmm_h,oars,runn,YMD_DATE,HM,kstamp,kstamp2,HM_SFC,YMD_DATE_SFC

    INTEGER                :: iele,NELE,NELE_SFC,NVAL,NT,NELE_INFO,LN
    INTEGER                :: bit_alt,btyp_offset,btyp_offset_uni
    character(len = 5)     :: BURP_TYP
    CHARACTER(LEN=9)       :: STNID,STN_RESUME
    LOGICAL                :: HIRES,HIRES_SFC,HIPCS
    INTEGER                :: NDATA,NDATA_SF
    INTEGER                :: IER,date2,DATE3,time2,time_sonde,NEWDATE
    REAL                   :: RAD_MOY,RAD_STD
    INTEGER                :: iclass,NCHANAVHRR,NCLASSAVHRR,ichan,iobs,inorm
    INTEGER                :: infot
    
    DATA LISTE_INFO  &
       /1007,002019,007024,007025 ,005021, 005022, 008012, &
        013039,020010,2048,2022,33060,33062,33039,10035,10036,08046,5043/

    RELEV2=0.0
    FAMILYTYPE2= 'SCRAP'
    vcord_type(:)=-1
    vcord_type(1)=0
    NELE_INFO=1
    NELE_SFC=0
    NELE=0
    BNBITSOFF=0
    BNBITSON=0
    ENFORCE_CLASSIC_SONDES=.false.
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
        CALL BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'UA'
        LISTE_ELE(1:5) = (/12001,11001,11002,12192,10194/)
        NELE=5
        ENFORCE_CLASSIC_SONDES=.false.
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('AI')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:4) = (/12001,12192,11001,11002/)
        NELE=4
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('AL')
        BURP_TYP='uni'
        vcord_type(1)=7071

        LISTE_ELE(1) = 40030
        NELE=1
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('SW')
        BURP_TYP='uni'
        vcord_type(1)=7004

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('SF','GP')
        BURP_TYP='uni'
        vcord_type(1)=0
        NELE_SFC=7
        LISTE_ELE_SFC(1:7) = (/12004,11011,11012,10051,10004,12203,15031/)

        CALL BRPACMA_NML('namburp_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2= 'SFC'
      CASE('SC')
        vcord_type(1)=0
        BURP_TYP='uni'
        LISTE_ELE_SFC(1:2) = (/11012,11011/)
        NELE=2
        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS

        FAMILYTYPE2= 'UASFC2'
      CASE('PR')
        BURP_TYP='multi'
        vcord_type(1)=7006

        LISTE_ELE(1:2) = (/11001,11002/)
        NELE=2

        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
      CASE('RO')
        BURP_TYP='multi'
        vcord_type(1)=7007
        vcord_type(2)=7040

        LISTE_ELE(1:2) = (/15036,15037/)
        NELE=2

        CALL BRPACMA_NML('namburp_conv')
        NELE=NELEMS
        !================GPS-RO CANNOT BE FILTERED=======
        BNBITSOFF=0
        BNBITSON=0
        !================GPS-RO CANNOT BE FILTERED=======
        NELE_INFO=18
      CASE('GO','MI','TO')
        BURP_TYP='multi'
        vcord_type(1)=5042
        vcord_type(2)=2150

        LISTE_ELE(1:1) = (/12163/)
        NELE=1

        CALL BRPACMA_NML('namburp_tovs')
        NELE=NELEMS

        NELE_INFO=18
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
        CALL BRPACMA_NML('namburp_chm_sfc')
        NELE_SFC=NELEMS_SFC

        FAMILYTYPE2='CH'
        LISTE_ELE(1:19) = (/15008,15009,15010,15020,15021,15022,15023,15024,15026,15027,15028, &
                          15029,15198,15199,15200,15230,13001,13002,08090/)
        NELE=19
        CALL BRPACMA_NML('namburp_chm')
        NELE=NELEMS 
    END SELECT
    LISTE_ELE    (1:NELE    )=BLISTELEMENTS(1:NELE)
    LISTE_ELE_SFC(1:NELE_SFC)=BLISTELEMENTS_SFC(1:NELE_SFC)
    if(NELE     > 0)write(*,*)  ' LISTE_ELE =',LISTE_ELE
    if(NELE_SFC > 0)write(*,*)  ' LISTE_ELE_SFC =',LISTE_ELE_SFC(1:NELE_SFC)
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

    WRITE(*,*) '-----------------------------------------------'
    WRITE(*,*) '-----------  BEGIN brpr_readBurp   ------------'
    WRITE(*,*) 'FAMILYTYPE vcord_type   =',FAMILYTYPE,vcord_type
    WRITE(*,*) 'BURP_TYP btyp_offset    =',BURP_TYP, btyp_offset
    WRITE(*,*) '-----------------------------------------------'


    ! initialisation
    ! --------------
    SUM=0
    flag_passage1 = 0
    flag_passage2 = 0
    flag_passage3 = 0
    flag_passage4 = 0


    opt_missing = 'MISSING'


    Call BURP_Set_Options( &
       & REAL_OPTNAME       = opt_missing, &
       & REAL_OPTNAME_VALUE = MPC_missingValue_R4, &
       & CHAR_OPTNAME       = 'MSGLVL', &
       & CHAR_OPTNAME_VALUE = 'FATAL', &
       & IOSTAT             = error )

    Call BURP_Init(File_in      ,IOSTAT=error)
    Call BURP_Init(Rpt_in       ,IOSTAT=error)
    Call BURP_Init(Block_in     ,IOSTAT=error)


    ! opening file
    write(*,*) 'OPENING BURP FILE FOR READING = ', trim(brp_file)

    Call BURP_New(File_in, FILENAME = brp_file, &
                   & MODE    = FILE_ACC_READ, &
                   & IOSTAT  = error )


    ! obtain input burp file number of reports

    Call BURP_Get_Property(File_in, NRPTS=nb_rpts)

    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'IOSTAT    =',error
    WRITE(*,*) 'NUMBER OF REPORTS IN FILE  = ',nb_rpts
    WRITE(*,*) '-----------------------------------------'


    ! scan input burp file to get all reports address

    Allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0
    bit_alt = 0

    do
      ref_rpt = BURP_Find_Report(File_in, &
                    & REPORT      = Rpt_in,  &
                    & SEARCH_FROM = ref_rpt, &
                    & IOSTAT      = error)
      Call burp_get_property(Rpt_in, STNID = stnid )
      IF ( stnid(1:2) == ">>" ) then
        STN_RESUME=stnid
        TYPE_RESUME=STN_RESUME(3:9)
        SELECT CASE(stnid)
          CASE(">>BGCKALT", ">>POSTALT")
            bit_alt=1
          CASE(">>DERIALT")
            bit_alt=2
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

    !pik  write(*,'(a9,1x,a16,1x,i2)' )STN_RESUME,' bit_alt==== >  ',bit_alt
    write(*, *)STN_RESUME,' bit_alt==== >  ',bit_alt

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
            

        Call BURP_Get_Report(File_in, &
                      & REPORT    = Rpt_in, &
                      & REF       = address(kk), &
                      & IOSTAT    = error) 
        Call burp_get_property(Rpt_in, &
               STNID = stnid ,TEMPS =hhmm_h,FLGS = status ,IDTYP =idtyp,LATI = lati &
               ,LONG = long ,DX = dx ,DY = dy,ELEV=elev,DRND =drnd,DATE =date_h &
               ,OARS =oars,RUNN=runn ,IOSTAT=error)
        IF ( stnid(1:2) == ">>" ) cycle


        !  LOOP ON BLOCKS

        ref_blk = 0

        HIRES=.FALSE.
        HIPCS=.FALSE.
        HIRES_SFC=.FALSE.
        BLOCKS1: do

          ref_blk = BURP_Find_Block(Rpt_in, &
                    & BLOCK       = Block_in, &
                    & SEARCH_FROM = ref_blk, &
                    & IOSTAT      = error)

          if (ref_blk < 0) EXIT BLOCKS1

          Call BURP_Get_Property(Block_in, &
                      & NELE   = nbele, &
                      & NVAL   = nvale, &
                      & NT     = nte, &
                      & BFAM   = bfam, &
                      & BTYP   = btyp, &
                      & BKSTP  = BKSTP, &
                      & IOSTAT = error)


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
          if ( btyp10 - btyp10obs_uni == 0 .and. bkstp <= 4 ) then

            flag_passage3 = 1

            ALLOCATE(obsvalue_sfc(NELE_SFC,1,nte))
            ALLOCATE(  OBSERV_SFC(NELE_SFC,1) )

            ALLOCATE( vcoord_sfc(1))

            vcoord_type=0
            vcoord_SFC(:)       = 0
            obsvalue_sfc(:,:,:) = MPC_missingValue_R4
            IND_LAT   = BURP_Find_Element(Block_in, ELEMENT=5001, IOSTAT=error)
            IND_LON   = BURP_Find_Element(Block_in, ELEMENT=6001, IOSTAT=error)
            IND_TIME  = BURP_Find_Element(Block_in, ELEMENT=4015, IOSTAT=error)
            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES_SFC=.true.
            if (HIRES_SFC) ALLOCATE(HLAT_SFC(nte),HLON_SFC(nte),HTIME_SFC(nte) )
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

              iele=LISTE_ELE_SFC(IL)
              IND_ele  = BURP_Find_Element(Block_in, ELEMENT=iele, IOSTAT=error)
              if (IND_ele < 0 ) cycle

              do k=1,nte
                obsvalue_sfc(IL,1,k) =  BURP_Get_Rval(Block_in, &
                                         & NELE_IND = IND_ele, &
                                         & NVAL_IND = 1, &
                                         & NT_IND   = k)
              end do

            end do
          end if

          if ( btyp10 - btyp10flg_uni == 0 .and. bkstp <= 4 ) then

            flag_passage4 = 1

            ALLOCATE( qcflag_sfc (NELE_SFC,1,nte))
            ALLOCATE( qcflags_SFC(NELE_SFC,1) )
            QCFLAGS_SFC(:,:)=0

            do IL = 1, NELE_SFC
              iele=LISTE_ELE_SFC(IL) + 200000
              IND_QCFLAG  = BURP_Find_Element(Block_in, ELEMENT=iele, IOSTAT=error)
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

          if ( btyp10 - btyp10obs == 0 .and. bfam == 0 ) then

            flag_passage3 = 1

            NVAL=NVALE ;  NT=NTE
            ALLOCATE(obsvalue(NELE,nvale,nte),VCOORD(nvale,nte))
            ALLOCATE(  OBSERV(NELE,nvale)    )
            ALLOCATE(   VCORD(nvale)    )
                           
            obsvalue(:,:,:) = MPC_missingValue_R4
            OBSERV  (:,:)   = MPC_missingValue_R4
            VCOORD  (:,:)   = 0.
            VCORD   (:)     = 0.

            k=0
            IND_VCOORD=-1
            do while (vcord_type(k+1) /= -1 .and. IND_VCOORD == -1)
               k=k+1 
               IND_VCOORD  = BURP_Find_Element(Block_in, ELEMENT=vcord_type(k),IOSTAT=error)
            end do
            vcoord_type=0
            IF (IND_VCOORD > 0) vcoord_type=vcord_type(k)
            !if (IND_VCOORD == -1)write(*,*) 'PAS DE COORDONNEE VERTICALE STNID=',STNID,trim(FAMILYTYPE)

            ! LAT LON TIME IN DATA BLOCK
            IND_LAT   = BURP_Find_Element(Block_in, ELEMENT=5001, IOSTAT=error)
            IND_LON   = BURP_Find_Element(Block_in, ELEMENT=6001, IOSTAT=error)
            IND_TIME  = BURP_Find_Element(Block_in, ELEMENT=4015, IOSTAT=error)
            IND_EMIS  = BURP_Find_Element(Block_in, ELEMENT=55043,IOSTAT=error)
            if (IND_LAT > 0 .and. IND_LON > 0 .and. IND_TIME > 0 ) HIRES=.true.

            if( TRIM(FAMILYTYPE2) == 'UA' .and. UA_HIGH_PRECISION_TT_ES .eqv. .true. ) HIPCS=.true.

            if(HIRES) ALLOCATE(HLAT(nvale,nte),HLON(nvale,nte),HTIME(nvale,nte) )

            ALLOCATE(EMIS(nvale,nte))
            ALLOCATE(SURF_EMIS(nvale))
            EMIS(:,:)       = MPC_missingValue_R4
            OBSVALUE(:,:,:) = MPC_missingValue_R4

            do IL = 1, NELE

              iele=LISTE_ELE(IL)
              IND_ele  = BURP_Find_Element(Block_in, ELEMENT=iele, IOSTAT=error)
              if (IND_ele < 0 ) cycle

              if(HIPCS .and. iele == 12001) IND_ele = BURP_Find_Element(Block_in, ELEMENT=12101, IOSTAT=error)
              if(HIPCS .and. iele == 12192) IND_ele = BURP_Find_Element(Block_in, ELEMENT=12239, IOSTAT=error)

              do k=1,nte
                do j=1,nvale
                  obsvalue(IL,j,k) =  BURP_Get_Rval(Block_in, &
                                       & NELE_IND = IND_ele, &
                                       & NVAL_IND = j, &
                                       & NT_IND   = k)
                  IF (HIRES) THEN
                    HLAT(j,k) =BURP_Get_Rval(Block_in, &
                                           & NELE_IND = IND_LAT, &
                                           & NVAL_IND = j, &
                                           & NT_IND   = k)
                    HLON(j,k) =BURP_Get_Rval(Block_in, &
                                           & NELE_IND = IND_LON, &
                                           & NVAL_IND = j, &
                                           & NT_IND   = k)
                    HTIME(j,k)=BURP_Get_Rval(Block_in, &
                                           & NELE_IND = IND_TIME, &
                                           & NVAL_IND = j, &
                                           & NT_IND   = k)
                  END IF

                  IF (IND_EMIS > 0) THEN
                    EMIS(j,k) =BURP_Get_Rval(Block_in, &
                                           & NELE_IND = IND_EMIS, &
                                           & NVAL_IND = j, &
                                           & NT_IND   = k)
                  END IF

                  if (IND_VCOORD <= 0) cycle
                    VCOORD(j,k) =  BURP_Get_Rval(Block_in, &
                              &   NELE_IND = IND_VCOORD, &
                              &   NVAL_IND = j, &
                              &   NT_IND   = k)
                end do
              end do

            end do

!==================================================================================
!
! Lire les divers elements de la famille SW
!
            ALLOCATE(qi1val(nte))
            ALLOCATE(qi2val(nte))
            ALLOCATE(qival(nte))
            ALLOCATE(mtval(nte))
            ALLOCATE(lsval(nte))
            ALLOCATE(haval(nte))
            ALLOCATE(gaval(nte))
            QI1VAL(:) = 0
            QI2VAL(:) = 0
            QIVAL (:) = 0
            MTVAL (:) = 0
            LSVAL (:) = 0
            HAVAL (:) = 0
            GAVAL (:) = 0

            if (TRIM(FAMILYTYPE) == 'SW' .and. READ_QI_GA_MT_SW) then

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=33007, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                QI1VAL(k)= BURP_Get_Tblval(Block_in, &
                                           NELE_IND = IND_SW, &
                                           NVAL_IND = 1, &
                                           NT_IND   = k, &
                                           IOSTAT   = error)
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=33194, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                QI2VAL(k)= BURP_Get_Tblval(Block_in, &
                                           NELE_IND = IND_SW, &
                                           NVAL_IND = 1, &
                                           NT_IND   = k, &
                                           IOSTAT   = error)
              end do

              do k = 1, nte
                QIVAL(k) = QI2VAL(k)
                if(QIVAL(k) < 0)  QIVAL(k) = QI1VAL(k)
              enddo

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=2023, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                MTVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=8012, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                LSVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=13039, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                GAVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
              end do

              IND_SW  = BURP_Find_Element(Block_in, ELEMENT=2163, IOSTAT=error)
              if (IND_SW <= 0 ) cycle
              do k = 1, nte
                HAVAL(k)= BURP_Get_Tblval(Block_in, &
                                          NELE_IND = IND_SW, &
                                          NVAL_IND = 1, &
                                          NT_IND   = k, &
                                          IOSTAT   = error)
              end do

            !====================================================================
            !
            ! Lire les divers elements de la famille AL
            !
            else if (trim(familytype) == 'AL') then
              if (.not. allocated(azimuth)) then
                allocate(azimuth(nte))
                azimuth(:) = 0
              end if

              ! Read in the azimuth)
              ind_al=burp_find_element(block_in, element=BUFR_NEAZ, iostat=error)
              if (ind_al <= 0 ) cycle

              do k = 1, nte
                azimuth(k)= BURP_Get_Tblval(Block_in, &
                                            nele_ind = ind_al, &
                                            nval_ind = 1, &
                                            nt_ind   = k, &
                                            iostat   = error)
              end do
            end if ! AL
            !
            !====================================================================

          end if


          ! flag block (btyp = 0111 100011X XXXX)
          if ( btyp10 - btyp10flg == 0 ) then

            flag_passage4 = 1
            ALLOCATE(qcflag( NELE,nvale,nte))
            ALLOCATE(qcflags(NELE,nvale) )
            QCFLAG (:,:,:)  = 0
            QCFLAGS(:,:)    = 0

            do IL = 1, NELE

              iele=LISTE_ELE(IL)

              IND_QCFLAG  = BURP_Find_Element(Block_in, ELEMENT=200000+iele, IOSTAT=error)
              if (IND_QCFLAG <= 0 ) cycle
              do k = 1, nte
                do j = 1, nvale
                  QCFLAG(IL,j,k)= BURP_Get_Tblval(Block_in, &
                                  &   NELE_IND = IND_QCFLAG, &
                                  &   NVAL_IND = j, &
                                  &   NT_IND   = k, &
                                  &   IOSTAT   = error)
                  SUM = SUM +1
                end do
              end do

            end do

          end if

          ! info block (btyp = 0001 100000X XXXX) 
          BTYP10inf = 96

          if ( (btyp10 - btyp10inf == 0) .or. (btyp10 - btyp10inf == 1) ) then

            ALLOCATE( RINFO(NELE_INFO,nte))
            ALLOCATE(TRINFO(NELE_INFO))

            flag_passage2 = 1

            do kl=1,NELE_INFO
              info_elepos = BURP_Find_Element(Block_in, &
                                 & ELEMENT  = LISTE_INFO(kl), &
                                 & IOSTAT   = error)
              if (  info_elepos >= 0 )then
                
                do k =1 , nte
                  RINFO(kl,k)= BURP_Get_rval(Block_in, &
                              &   NELE_IND = info_elepos, &
                              &   NVAL_IND = 1, &
                              &   NT_IND   = k, &
                              &   IOSTAT   = error)
                  if  (RINFO(kl,k) == MPC_missingValue_R4)  THEN    
                    infot= BURP_Get_tblval(Block_in, &
                            &   NELE_IND = info_elepos, &
                            &   NVAL_IND = 1, &
                            &   NT_IND   = k, &
                            &   IOSTAT   = error)
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

          if ( BTYP10 - BTYP10des == 0 ) then

            flag_passage1 = 1

            ALLOCATE(GLBFLAG(nte))
            ALLOCATE(    lat(nte))
            ALLOCATE(    lon(nte))
            ALLOCATE(   date(nte))
            ALLOCATE(   hhmm(nte))

            ! DATE  004208  HHMM 004197  STATUS 055200  LAT 005002  LON 006002  DELAY 004195
            ind055200 = BURP_Find_Element(Block_in, ELEMENT=055200, IOSTAT=error)
            ind5002   = BURP_Find_Element(Block_in, ELEMENT=5002  , IOSTAT=error)
            ind6002   = BURP_Find_Element(Block_in, ELEMENT=6002  , IOSTAT=error)
            ind4208   = BURP_Find_Element(Block_in, ELEMENT=4208  , IOSTAT=error)
            ind4197   = BURP_Find_Element(Block_in, ELEMENT=4197  , IOSTAT=error)

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
            CFRAC(:,:)=-999

            IASIQUAL: DO k = 1, nte
              iclass=1

              NVALS :do j=1,nvale
                DO  il=1,nbele
                  iele=BURP_Get_Element(Block_in,INDEX =il,IOSTAT= error) 
                  SELECT CASE(iele)
                    CASE(25085)
                      CFRAC(iclass,k)= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                    CASE(5042)
                      !ICHAN= BURP_Get_TBLVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      ICHAN= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                    CASE(25142)
                      !INORM= BURP_Get_TBLVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      INORM= BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                    CASE(14047)
                      RAD_MOY=BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      RADMOY(iclass,ICHAN,k)=RAD_MOY * 10.d0**(-1.d0 * INORM ) * 100000.d0
                    CASE(14048)
                      RAD_STD=BURP_Get_RVAL(Block_in,NELE_IND = il, NVAL_IND = j, NT_IND = k, IOSTAT = error)
                      RADSTD(iclass,ICHAN,k)= RAD_STD * 10.d0**(-1.d0 * INORM ) * 100000.d0
                      IF (ICHAN==NCHANAVHRR) iclass=iclass+1
                      IF ( iclass==(NCLASSAVHRR+1) ) EXIT NVALS
                  END SELECT
                END DO
              end do NVALS

            END DO IASIQUAL

          end if

        end do BLOCKS1


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
            RELEV   =REAL(ELEV,OBS_REAL) - 400.
          END IF

          if (allocated(RINFO))      TRINFO(1:NELE_INFO)     =RINFO       (1:NELE_INFO,k)


          if(ENFORCE_CLASSIC_SONDES .eqv. .true.) hires=.false.

          IF  (HIRES ) THEN

            if (allocated(EMIS))     SURF_EMIS(1:NVAL)    =EMIS      (1:NVAL,k)

            NDATA   =0
            NDATA_SF=0
            IF ( allocated(obsvalue_sfc) ) THEN
              OBSERV_SFC(1:NELE_SFC,1:1)=obsvalue_sfc(1:NELE_SFC,1:1,k)
              QCFLAGS_sfc(1:NELE_SFC,1:1)=qcflag_sfc  (1:NELE_SFC,1:1,k)
              IF ( HIRES_SFC) THEN
                XLAT=HLAT_SFC(k);XLON=HLON_SFC(k);XTIME=HTIME_SFC(k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON
                ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                CALL INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)
                time2=time_sonde/10000
                YMD_DATE_SFC=date2
                HM_SFC=time2
              END IF

              NDATA_SF= WRITE_BODY(obsdat,UNI_FAMILYTYPE,RELEV,vcoord_sfc,vcoord_type,OBSERV_sfc,qcflags_sfc,NELE_SFC,1,LISTE_ELE_SFC)
              IF ( NDATA_SF > 0) THEN
                call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE_SFC,HM_SFC,idtyp,STATUS,RELEV,FILENUMB)
                OBSN=obs_numHeader(obsdat)
                call obs_setFamily(obsdat,trim(FAMILYTYPE),  OBSN )
                call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA_SF)
                IF (OBSN > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                END IF
              END IF

            END IF

            IF ( allocated(obsvalue) ) THEN

              ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
              do  JJ =1,nval
                OBSERV(1:NELE,1:1)  =obsvalue  (1:NELE,jj:jj,k)
                if (allocated(qcflag))     QCFLAGS(1:NELE,1:1)    =qcflag  (1:NELE,jj:jj,k)
                XLAT=HLAT(jj,k);XLON=HLON(jj,k);XTIME=HTIME(jj,k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON

                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                CALL INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)

                time2=time_sonde/10000

                VCORD(1)=VCOORD(jj,k)
                NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type,OBSERV,qcflags,NELE,1,LISTE_ELE,SURF_EMIS)

                IF (NDATA > 0) THEN
                  call WRITE_HEADER(obsdat,STNID,XLAT,XLON,date2,time2,idtyp,STATUS,RELEV,FILENUMB)
!==================================================================================
!
! Ajoute qivals dans les argument de WRITE_QI

                  if (TRIM(FAMILYTYPE) == 'SW' .and. READ_QI_GA_MT_SW) call WRITE_QI(obsdat,QIVAL(k),MTVAL(k),LSVAL(k),HAVAL(k),GAVAL(k))

                  if(trim(familytype) == 'AL')call write_al(obsdat, azimuth(k))


                  OBSN=obs_numHeader(obsdat)
                  call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                  call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA)
                  IF (OBSN > 1 ) THEN
                    LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                    call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                    !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                  ELSE
                    call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                    !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                  END IF
                END IF

              end do

            END IF

          ELSE

            if (allocated(EMIS))      SURF_EMIS(1:NVAL)        =EMIS     (1:NVAL,k)

            NDATA   =0
            NDATA_SF=0
            IF ( allocated(obsvalue_sfc) ) THEN
              IF ( HIRES_SFC) THEN
                XLAT=HLAT_SFC(k);XLON=HLON_SFC(k);XTIME=HTIME_SFC(k)
                IF ( XLON  < 0. ) XLON  = 360. + XLON
                ier= NEWDATE(kstamp2,YMD_DATE,HM*10000,3)
                XLAT=XLAT*MPC_RADIANS_PER_DEGREE_R8
                XLON=XLON*MPC_RADIANS_PER_DEGREE_R8

                CALL INCDATR(kstamp, kstamp2, XTIME/60.d0  )
                IER=newdate(kstamp,date2,time_sonde,-3)
                time2=time_sonde/10000
                YMD_DATE=date2
                HM=time2
              END IF

              OBSERV_SFC (1:NELE_SFC,1:1)=obsvalue_sfc(1:NELE_SFC,1:1,k)
              QCFLAGS_sfc(1:NELE_SFC,1:1)=qcflag_sfc  (1:NELE_SFC,1:1,k)

              NDATA_SF= WRITE_BODY(obsdat,UNI_FAMILYTYPE,RELEV,vcoord_sfc,vcoord_type,OBSERV_sfc,qcflags_sfc,NELE_SFC,1,LISTE_ELE_SFC)
              IF ( NDATA_SF > 0) THEN
                call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE,HM,idtyp,STATUS,RELEV,FILENUMB)
                OBSN=obs_numHeader(obsdat) 
                call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                call obs_headSet_i(obsdat,OBS_NLV ,OBSN,NDATA_SF)
                IF (OBSN  > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                  !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                  !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                END IF
              END IF
            END IF
          
            IF ( allocated(obsvalue) ) THEN
              OBSERV(1:NELE,1:NVAL)    =obsvalue(1:NELE,1:NVAL,k)
              QCFLAGS(1:NELE,1:NVAL)   =qcflag  (1:NELE,1:NVAL,k)
              VCORD(1:NVAL)            =VCOORD  (1:NVAL,k)
              NDATA= WRITE_BODY(obsdat,familytype,RELEV,VCORD,vcoord_type,OBSERV,qcflags,NELE,NVAL,LISTE_ELE,SURF_EMIS)

              IF (NDATA > 0) THEN

                IF (NDATA_SF == 0) THEN
                  call WRITE_HEADER(obsdat,STNID,XLAT,XLON,YMD_DATE,HM,idtyp,STATUS,RELEV,FILENUMB)
!==================================================================================
!
! Ajoute qivals dans les argument de WRITE_QI

                  if (TRIM(FAMILYTYPE) == 'SW') call WRITE_QI(obsdat,QIVAL(k),MTVAL(k),LSVAL(k),HAVAL(k),GAVAL(k))

                  if(trim(familytype) == 'AL')call write_al(obsdat, azimuth(k))

       	          OBSN=obs_numHeader(obsdat) 
                  call obs_setFamily(obsdat,trim(FAMILYTYPE), OBSN )
                END IF
                OBSN=obs_numHeader(obsdat) 
                call obs_headSet_i(obsdat,OBS_NLV,OBSN,NDATA+NDATA_SF)
                IF (OBSN   > 1 ) THEN
                  LN= obs_headElem_i(obsdat,OBS_RLN,OBSN-1) + obs_headElem_i(obsdat,OBS_NLV,OBSN-1)
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,LN)
                  !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                ELSE
                  call obs_headSet_i(obsdat,OBS_RLN,OBSN,1)
                  !call obs_headSet_i(obsdat,OBS_IDO,OBSN,kk)
                END IF

              END IF

            END IF

            !============ IASI =====================================
            if ( allocated(RADMOY) .and. NDATA > 0 ) then
              OBSN=obs_numHeader(obsdat)

              iclass=1
              do iobs=OBS_CF1,OBS_CF7
                if(obs_columnActive_IH(obsdat,iobs)) then
                  call obs_headSet_i(obsdat,iobs,OBSN,CFRAC(iclass,k))
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
              call WRITE_INFO(obsdat,familytype, TRINFO,LISTE_INFO,NELE_INFO  )
            END IF
          end if

        end do

        !---------UPPER AIR---------------------------
        if ( allocated(obsvalue) ) then
          DEALLOCATE ( obsvalue,VCOORD,VCORD,observ)
        end if
        if ( allocated(qcflag) ) then
          DEALLOCATE (qcflag,qcflags)
        end if
        if ( allocated(qi1val) ) then
          DEALLOCATE (qi1val)
        end if
        if ( allocated(qi2val) ) then
          DEALLOCATE (qi2val)
        end if
        if ( allocated(qival) ) then
          DEALLOCATE (qival)
        end if
        if ( allocated(mtval) ) then
          DEALLOCATE (mtval)
        end if
        if ( allocated(lsval) ) then
          DEALLOCATE (lsval)
        end if
        if ( allocated(haval) ) then
          DEALLOCATE (haval)
        end if
        if ( allocated(gaval) ) then
          DEALLOCATE (gaval)
        end if
        if ( allocated(EMIS) ) then
          DEALLOCATE (EMIS,SURF_EMIS)
        end if
        if (allocated(azimuth)) then
          deallocate(azimuth)
        end if

        !---------SURFACE-----------------------------
        if ( allocated(obsvalue_sfc) ) then
          DEALLOCATE(obsvalue_sfc,vcoord_sfc,OBSERV_SFC)
        end if

        if ( allocated(qcflag_sfc) ) then
          DEALLOCATE(  qcflag_sfc, qcflags_SFC)
        end if

        !--------SURFACE------------------------------
        if (  allocated(lat) ) then
          DEALLOCATE (lat,lon,date,hhmm,glbflag)
        end if
        if (  allocated(hlat) ) then
          DEALLOCATE (hlat,hlon,htime)
        end if
        if (  allocated(hlat_sfc) ) then
          DEALLOCATE (hlat_sfc,hlon_sfc,htime_sfc)
        end if
        if ( allocated(rinfo) ) then
          DEALLOCATE (rinfo,trinfo)
        end if
        if ( allocated(RADMOY) ) then
          DEALLOCATE (RADMOY,CFRAC,radstd)
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


    Call BURP_Free(File_in,      IOSTAT=error)
    Call BURP_Free(Rpt_in,       IOSTAT=error)
    Call BURP_Free(Block_in,     IOSTAT=error)

    write(*,*)' file   Nobs SUM = ',trim(brp_file),obs_numHeader(obsdat),SUM
  END SUBROUTINE brpr_readBurp




  FUNCTION WRITE_BODY(obsdat,FAMTYP, ELEV,VERTCOORD,VCOORD_TYPE, &
                      obsvalue,qcflag,NELE,NVAL,LISTE_ELE,SURF_EMIS_opt)

    !******************************************************************************** 
    !
    !       REVISION:
    !                 Ping Du (CMDA), Nov 2014
    !                  -- Added  CASE ('CH') with VCO=4 
    !                 Y.J. Rochon (ARQI) and M. Sitwell (ARQI), Dec 2014, Feb 2015, June 2015
    !                  -- Added VCOORD_TYPE as input argument and applied in function.
    !                  -- Added new VCO cases for the CH family.
    !                  -- Added abort statement if constituent vertical coordinate not recognized.
    !
    !****************************************************************************

    implicit none
    type (struct_obs), intent(inout) :: obsdat

    INTEGER ::  WRITE_BODY,VCOORD_TYPE
    REAL   , allocatable          ::   OBSVALUE(:,:)
    REAL   , allocatable,optional ::   SURF_EMIS_opt(:)
    INTEGER, allocatable          ::     QCFLAG(:,:)
    REAL   , allocatable          ::  VERTCOORD(:)

    CHARACTER*2 ::   FAMTYP
    REAL        ::   ELEVFACT,VCOORD,ZFACT,INFOV
    INTEGER     ::   NELE,NVAL
    integer     ::   LISTE_ELE(:)

    INTEGER     ::   ID_OBS,ID_DATA

    INTEGER     ::   NOBS
    INTEGER     ::   VARNO,IL,J,COUNT,NLV


    INTEGER     ::   IFLAG,LN,BITSflagoff,BITSflagon
    REAL(OBS_REAL) :: MISG,OBSV,ELEV,ELEV_R,REMIS,emmissivite
    INTEGER     ::   VCO
    INTEGER     ::   NONELEV
    REAL        ::   ZEMFACT
    LOGICAL     ::   L_EMISS
    
    if(present(SURF_EMIS_opt)) then
      L_EMISS=.true.
    else
      L_EMISS=.false.
    end if

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

    id_obs=NOBS


    if (  trim(FAMTYP) == trim('PR') .OR. trim(FAMTYP) == trim('SF') ) then
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
        if( L_EMISS .eqv. .true.)  then
          if( SURF_EMIS_opt(j) /= MISG)  then
            REMIS =  SURF_EMIS_opt(j)*ZEMFACT
          else
            REMIS = MISG
          end if
        end if
        IFLAG = INT(qCflag(il,j))
        if(iand(iflag,BITSflagoff) /= 0) cycle
        !burpmodule IFLAG = IBCLR(IFLAG,12)

        !if (VARNO == .10194)OBSV=OBSV*RG
        if ( obsv /= MPC_missingValue_R4 .and. VCOORD /= MPC_missingValue_R4   ) then
          count = count  + 1
          NLV= NLV +1
          ID_DATA=count
          IFLAG = IBCLR(IFLAG,12)
          call obs_bodySet_r(obsdat,OBS_VAR,count,OBSV)
          call obs_bodySet_i(obsdat,OBS_VNM,count,VARNO)
          call obs_bodySet_i(obsdat,OBS_VCO,count,VCO)
          ELEV_R=VCOORD + ELEV*ELEVFACT
          call obs_bodySet_r(obsdat,OBS_PPP,count, ELEV_R)
          call obs_bodySet_i(obsdat,OBS_FLG,count,IFLAG)
          if ( REMIS /= MPC_missingValue_R4 .and. FAMTYP == 'TO') THEN
            call obs_bodySet_r(obsdat,OBS_SEM,count,REMIS)
          else
            if ( FAMTYP == 'TO') then
                emmissivite=0.95
                call obs_bodySet_r(obsdat,OBS_SEM,count,emmissivite)
             else
                call obs_bodySet_r(obsdat,OBS_SEM,count,MISG)
             end if
          end if
          !call obs_set_i(obsdat,'OBS',count,NOBS)

          call obs_bodySet_i(obsdat,OBS_VCO,count,VCO)
          !call obs_bodySet_i(obsdat,OBS_IDD,count,ID_DATA)
          !call obs_bodySet_i(obsdat,OBS_IDD,count,0)

          if ( varno == 11001 .or. varno == 11011) then
            call obs_bodySet_r(obsdat,OBS_VAR,count,OBSV)
            if ( varno == 11001) then
              call obs_bodySet_i(obsdat,OBS_VNM,count+1,11003)
              call obs_bodySet_i(obsdat,OBS_FLG,count+1,0)
              !call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
              ELEV_R=VCOORD + ELEV*ELEVFACT
              call obs_bodySet_r(obsdat,OBS_PPP,count+1,ELEV_R)
              call obs_bodySet_i(obsdat,OBS_VCO,count+1,VCO)

              call obs_bodySet_i(obsdat,OBS_VNM,count+2,11004)
              !call obs_set_i(obsdat,'OBS',count+2,NOBS)
              call obs_bodySet_i(obsdat,OBS_FLG,count+2,0)
              !call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
              call obs_bodySet_r(obsdat,OBS_PPP,count+2,ELEV_R)
              call obs_bodySet_i(obsdat,OBS_VCO,count+2,VCO)
            else
              call obs_bodySet_i(obsdat,OBS_VNM,count+1,11215)
              call obs_bodySet_i(obsdat,OBS_FLG,count+1,0)
              !call obs_bodySet_i(obsdat,OBS_IDD,count+1,-1)
              ELEV_R=VCOORD + ELEV*ELEVFACT
              call obs_bodySet_r(obsdat,OBS_PPP,count+1,ELEV_R)
              call obs_bodySet_i(obsdat,OBS_VCO,count+1,VCO)

              call obs_bodySet_i(obsdat,OBS_VNM,count+2,11216)
              call obs_bodySet_i(obsdat,OBS_FLG,count+2,0)
              !call obs_bodySet_i(obsdat,OBS_IDD,count+2,-1)
              call obs_bodySet_r(obsdat,OBS_PPP,count+2,ELEV_R)
              call obs_bodySet_i(obsdat,OBS_VCO,count+2,VCO)
            end if

            call obs_bodySet_r(obsdat,OBS_VAR,count+1,MISG)
            call obs_bodySet_r(obsdat,OBS_VAR,count+2,MISG)

            count = count + 2
            NLV = NLV + 2
          end if

        end if

      END DO

    END DO


    WRITE_BODY=NLV

  END FUNCTION  WRITE_BODY


  subroutine GET_HEADER(obsdat, LAT,LON,DATE,TIME,CODTYP,STATUS,ELEV,NOBS,FILENUMB)

    implicit none
    type (struct_obs), intent(inout) :: obsdat

    INTEGER     ::   DATE,TIME,CODTYP,STATUS,FILENUMB
    REAL(OBS_REAL)  :: ELEV,LAT,LON

    INTEGER     ::   LN,NOBS

    NOBS=obs_numHeader(obsdat) 

    LAT     = obs_headElem_r(obsdat,OBS_LAT,nobs)
    LON     = obs_headElem_r(obsdat,OBS_LON,nobs)
    DATE    = obs_headElem_i(obsdat,OBS_DAT,nobs)
    TIME    = obs_headElem_i(obsdat,OBS_ETM,nobs)
    CODTYP  = obs_headElem_i(obsdat,OBS_ITY,nobs)
    STATUS  = obs_headElem_i(obsdat,OBS_ST1,nobs)
    !NOBS    = obs_headElem_i(obsdat,OBS_IDO,nobs)
    FILENUMB= obs_headElem_i(obsdat,OBS_OTP,nobs)
    ELEV    = obs_headElem_r(obsdat,OBS_ALT,nobs)

    RETURN

  END SUBROUTINE  GET_HEADER


  SUBROUTINE WRITE_HEADER(obsdat, STNID,LAT,LON,DATE,TIME,CODTYP,STATUS,ELEV,FILENUMB)

    implicit none
    type (struct_obs), intent(inout) :: obsdat
    CHARACTER(LEN=9)       :: STNID

    INTEGER     ::    DATE,TIME,CODTYP,STATUS
    INTEGER     ::    FILENUMB
    REAL(OBS_REAL) :: ELEV,LAT,LON

    INTEGER     ::   LN,NOBS

    NOBS=obs_numHeader(obsdat)  +1

    call obs_headSet_i(obsdat,OBS_ONM,nobs,nobs)
    call obs_headSet_r(obsdat,OBS_LAT,nobs,LAT)
    call obs_headSet_r(obsdat,OBS_LON,nobs,LON)
    call obs_headSet_i(obsdat,OBS_DAT,nobs,DATE)
    call obs_headSet_i(obsdat,OBS_ETM,nobs,TIME)
    call obs_headSet_i(obsdat,OBS_ITY,nobs,CODTYP)
    call obs_headSet_i(obsdat,OBS_ST1,nobs,STATUS)
    !call obs_headSet_i(obsdat,OBS_IDO,nobs,NOBS)
    call obs_headSet_r(obsdat,OBS_ALT,nobs,ELEV)
    !call obs_headSet_i(obsdat,OBS_IDF,nobs,FILENUMB)
    call obs_headSet_i(obsdat,OBS_OTP,nobs,FILENUMB)
    call obs_set_c(obsdat,'STID',nobs,STNID )

  END SUBROUTINE  WRITE_HEADER

!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------

  SUBROUTINE WRITE_QI(obsdat, QIvalue, MTvalue, LSvalue, HAvalue, GAvalue)

    implicit none
    type (struct_obs), intent(inout) :: obsdat
    INTEGER     ::  QIvalue, MTvalue, LSvalue, HAvalue, GAvalue
    INTEGER     ::  NOBS

    NOBS = obs_numHeader(obsdat)

    call obs_headSet_i(obsdat,OBS_SZA,nobs,QIvalue)
    call obs_headSet_i(obsdat,OBS_TEC,nobs,MTvalue)
    call obs_headSet_i(obsdat,OBS_AZA,nobs,LSvalue)
    call obs_headSet_i(obsdat,OBS_SUN,nobs,GAvalue)
    call obs_headSet_i(obsdat,OBS_CLF,nobs,HAvalue)

  END SUBROUTINE  WRITE_QI

!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------

  subroutine write_al(obsdat, azimuth)
    implicit none
    type (struct_obs), intent(inout) :: obsdat
    integer :: azimuth, nobs

    nobs = obs_numHeader(obsdat)

    call obs_headSet_i(obsdat,OBS_AZA,nobs,azimuth)
  end subroutine write_al

!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------

  subroutine WRITE_INFO(obsdat,FAMTYP, RINFO,LISTE_INFO,NELE_INFO  )

    !******************************************************************************** 
    !
    !       REVISION:
    !                 Y.J. Rochon (ARQI), Dec 2014
    !                  -- Added consideraton of 08046 for the CH family.
    !                     (see *CONSTITUENT*)
    !
    !****************************************************************************

    implicit none
    type (struct_obs), intent(inout) :: obsdat

    !REAL   , allocatable  ::      RINFO(:)
    REAL        ::      RINFO(NELE_INFO)

    CHARACTER*2 ::   FAMTYP
    REAL*4      ::   INFOV
    INTEGER     ::   NELE_INFO
    integer     ::   LISTE_INFO(NELE_INFO)

    INTEGER     ::   CODTYP

    INTEGER     ::   IL,J,NOBS
    INTEGER     ::   SENSOR,ID_SAT,INSTRUMENT,LAND_SEA,CONSTITUENT_TYPE
    INTEGER     ::   TERRAIN_TYPE,ZENITH,SOLAR_ZENITH,AZIMUTH,CLOUD_COVER,SOLAR_AZIMUTH
    INTEGER     ::   IGQISFLAGQUAL,IGQISQUALINDEXLOC,ITANGENT_RADIUS,IGEOID,IRO_QCFLAG
    INTEGER     ::   IFOV

    REAL        ::   RSOLAR_ZENITH,RCLOUD_COVER,RZENITH,RAZIMUTH,RSOLAR_AZIMUTH
    REAL        ::   RIGQISFLAGQUAL,RIGQISQUALINDEXLOC,RCONSTITUENT
    REAL        ::   RTERRAIN_TYPE,RLAND_SEA,RID_SAT,RSENSOR,RINSTRUMENT,RRO_QCFLAG
    REAL(OBS_REAL)        ::   RTANGENT_RADIUS,RGEOID
    REAL        ::   RFOV

    NOBS=obs_numHeader(obsdat)

    CODTYP=obs_headElem_i(obsdat,OBS_ITY,NOBS)
    !write(*,*)'  DEBUT WRITE_INFO NOBS CODTYP ----> ',NOBS,CODTYP,size(liste_info),size(RINFO),liste_info
    LAND_SEA  =0
    INSTRUMENT=0
    ID_SAT    =0
    ZENITH    =0
    SENSOR    =0
    AZIMUTH   =0

    SOLAR_AZIMUTH = 0
    SOLAR_ZENITH = 0
    RTANGENT_RADIUS=real(MPC_missingValue_R8,OBS_REAL)
    RGEOID=real(MPC_missingValue_R8,OBS_REAL)
    TERRAIN_TYPE=99
    !CLOUD_COVER = 0
    CONSTITUENT_TYPE = -99
    IFOV = -99

    !if ( allocated(rinfo) ) then

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
        CASE( 2048)
          RSENSOR = INFOV
          if (RSENSOR == MPC_missingValue_R4 ) THEN
            SENSOR = -99
          ELSE
            SENSOR = NINT(RSENSOR)
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
            ZENITH = 9000
          ELSE
            ZENITH = NINT ( (90.0 + RZENITH)*100 )
          END IF
        CASE( 7025)
          RSOLAR_ZENITH = INFOV
          if (RSOLAR_ZENITH == MPC_missingValue_R4 )  THEN 
            SOLAR_ZENITH = 0
            SOLAR_ZENITH = -99
          ELSE
            SOLAR_ZENITH = NINT ( (90.0 + RSOLAR_ZENITH)*100 )
          END IF
        CASE( 5021)
          RAZIMUTH=INFOV
          if (RAZIMUTH == MPC_missingValue_R4 ) THEN
            AZIMUTH=0
          ELSE
            AZIMUTH=NINT ( (RAZIMUTH)*100 )
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
          if (RSOLAR_AZIMUTH == MPC_missingValue_R4 )  then
            SOLAR_AZIMUTH=0
            SOLAR_AZIMUTH=-99
          ELSE
            SOLAR_AZIMUTH=NINT ( (RSOLAR_AZIMUTH)*100 )
          END IF
        CASE( 8012)
          RLAND_SEA=INFOV
          if (RLAND_SEA == MPC_missingValue_R4 ) THEN
            LAND_SEA=99
          ELSE
            LAND_SEA=NINT ( RLAND_SEA )
          END IF
        CASE( 13039)
          RTERRAIN_TYPE=INFOV
          if (RTERRAIN_TYPE == MPC_missingValue_R4 ) THEN
            TERRAIN_TYPE=99
          ELSE
            TERRAIN_TYPE=NINT ( RTERRAIN_TYPE )
          END IF
        CASE( 20010)
          RCLOUD_COVER=INFOV
          if (RCLOUD_COVER == MPC_missingValue_R4 ) THEN
            CLOUD_COVER=0
          ELSE
            CLOUD_COVER=NINT ( RCLOUD_COVER )
          END IF
        CASE( 10035)
          RTANGENT_RADIUS=INFOV
        CASE( 10036)
          RGEOID=INFOV
        CASE( 33039)
          RRO_QCFLAG=INFOV
          if (RRO_QCFLAG == MPC_missingValue_R4 ) THEN
            IRO_QCFLAG=-99
          ELSE
            IRO_QCFLAG=NINT ( RRO_QCFLAG )
          END IF
        CASE( 08046)          
          IF (trim(FAMTYP) == 'CH') THEN
             RCONSTITUENT=INFOV
             IF (RCONSTITUENT == MPC_missingValue_R4) THEN
                call utl_abort('WRITE_INFO: Missing 08046 element for the CH family.')
             ELSE
                CONSTITUENT_TYPE=NINT(RCONSTITUENT)
             END IF
          END IF
      END SELECT
    end do

    !-------------------SPECIAL CASES--------------

    ! INSTRUMENT
    IF ( SENSOR == -99)then
      IF ( INSTRUMENT == -99)then
        INSTRUMENT=0
      END IF
    ELSE
      INSTRUMENT = obsu_cvt_obs_instrum(sensor)
    END IF

    ! AIRS
    IF ( INSTRUMENT == 420 ) ID_SAT = 784

    if (  trim(FAMTYP) == trim('GO') ) then
      LAND_SEA=0
      ZENITH=0
    END IF
   

    !-------------------SPECIAL CASES--------------

    ! Is terrain type sea ice (iterrain=0)?, If so, set imask=2.
    IF ( TERRAIN_TYPE ==  0       ) THEN
      LAND_SEA = 2
    END IF

    if ( obs_columnActive_IH(obsdat,OBS_OFL) ) call obs_headSet_i(obsdat,OBS_OFL,nobs,LAND_SEA)
    if ( obs_columnActive_IH(obsdat,OBS_INS) ) call obs_headSet_i(obsdat,OBS_INS,nobs,INSTRUMENT  )
    if ( obs_columnActive_IH(obsdat,OBS_SZA) ) call obs_headSet_i(obsdat,OBS_SZA,nobs,ZENITH )
    if ( obs_columnActive_IH(obsdat,OBS_FOV) ) call obs_headSet_i(obsdat,OBS_FOV,nobs,IFOV )
    if ( obs_columnActive_IH(obsdat,OBS_CLF) ) call obs_headSet_i(obsdat,OBS_CLF,nobs,CLOUD_COVER )
    if ( obs_columnActive_IH(obsdat,OBS_AZA) ) call obs_headSet_i(obsdat,OBS_AZA,nobs,AZIMUTH )
    if ( obs_columnActive_IH(obsdat,OBS_SUN) ) call obs_headSet_i(obsdat,OBS_SUN,nobs,SOLAR_ZENITH )
    if ( obs_columnActive_IH(obsdat,OBS_SAZ) ) call obs_headSet_i(obsdat,OBS_SAZ,nobs,SOLAR_AZIMUTH )
    if ( obs_columnActive_IH(obsdat,OBS_SAT) ) call obs_headSet_i(obsdat,OBS_SAT,nobs,ID_SAT)
    if ( obs_columnActive_IH(obsdat,OBS_GQF) ) call obs_headSet_i(obsdat,OBS_GQF,nobs,IGQISFLAGQUAL)
    if ( obs_columnActive_IH(obsdat,OBS_GQL) ) call obs_headSet_i(obsdat,OBS_GQL,nobs,IGQISQUALINDEXLOC)
    !if( trim(FAMTYP) == trim('RO'))print *, 'geoid   QCFLAG TANGENT_RADIUS GEOID=',IRO_QCFLAG,RTANGENT_RADIUS,RGEOID
    if ( obs_columnActive_IH(obsdat,OBS_ROQF) ) call obs_headSet_i(obsdat,OBS_ROQF,nobs,IRO_QCFLAG)
    if ( obs_columnActive_RH(obsdat,OBS_TRAD) ) call obs_headSet_r(obsdat,OBS_TRAD,nobs,RTANGENT_RADIUS)
    if ( obs_columnActive_RH(obsdat,OBS_GEOI) ) call obs_headSet_r(obsdat,OBS_GEOI,nobs,RGEOID)
    if (trim(FAMTYP) == trim('CH')) then
        if ( obs_columnActive_IH(obsdat,OBS_CHM) ) call obs_headSet_i(obsdat,OBS_CHM,nobs,CONSTITUENT_TYPE)
    else
        if ( obs_columnActive_IH(obsdat,OBS_CHM) ) call obs_headSet_i(obsdat,OBS_CHM,nobs,-1)
    end if

  END SUBROUTINE  WRITE_INFO


  INTEGER  FUNCTION FIND_INDEX(LIST,ELEMENT)
    implicit none
    INTEGER LIST(:)
    INTEGER I,ELEMENT
    FIND_INDEX=-1
    do I=1,size (LIST)
      if (list(i) == element) THEN
        FIND_INDEX=i
        exit
      end if
    end do
    RETURN
  END FUNCTION FIND_INDEX

end module burpread_mod
