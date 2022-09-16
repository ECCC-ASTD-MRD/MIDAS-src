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

module burpFiles_mod
  ! MODULE burpFiles_mod (prefix='brpf' category='3. Observation input/output')
  !
  ! :Purpose: To store the filenames of the burp observation files and call
  !           subroutines in readBurp to read and update burp files.
  !
  use codePrecision_mod
  use mathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use burpread_mod
  use bufr_mod
  use utilities_mod
  use obsSubSpaceData_mod
  use burp_module
  use obsUtil_mod
  use obsVariableTransforms_mod

  implicit none
  save
  private

  ! public procedures
  public :: brpf_getDateStamp, brpf_readfile, brpf_updatefile
  public :: brpf_obsSub_read, brpf_obsSub_update

contains

  !--------------------------------------------------------------------------
  ! brpf_getDataStamp
  !--------------------------------------------------------------------------
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

  !--------------------------------------------------------------------------
  ! brpf_readFile
  !--------------------------------------------------------------------------
  subroutine brpf_readFile(obsdat,fileName,familyType,fileIndex)
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*) :: fileName
    character(len=*) :: familyType
    integer :: fileIndex

    ! locals
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headerIndexBegin, headerIndexEnd, headerIndex
    integer :: numBody, numHeader
    logical :: burp_chem
    real(pre_obsReal)  :: missingValue

    write(*,*) ' '
    write(*,*) 'brpf_readFile: Starting'
    write(*,*) ' '
    missingValue = real(MPC_missingValue_R8,pre_obsReal)
    
    bodyIndexBegin   = obs_numbody(obsdat) + 1
    headerIndexBegin = obs_numheader(obsdat) + 1
    call brpr_readBurp(obsdat,                         & ! INOUT
                       familyType, fileName, fileIndex)  ! IN
    bodyIndexEnd   = obs_numbody(obsdat)
    headerIndexEnd = obs_numheader(obsdat)
 
    burp_chem = trim(familyType) == 'CH'

    if ( trim(familyType) /= 'TO' .and. .not.burp_chem ) then

      call ovt_transformObsValues      (obsdat, headerIndexBegin, headerIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headerIndexBegin, headerIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headerIndexBegin, headerIndexEnd )

    end if  

    do headerIndex = headerIndexBegin, headerIndexEnd

      call obs_headSet_i(obsdat, OBS_OTP, headerIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headerIndex)

      ! For CH family, apply scaling from the element BUFR_SCALE_EXPONENT when present.
      if (burp_chem) call brpf_setScaleCH(obsdat, headerIndex, forward=.true.)

    end do

    ! initializations
    do bodyIndex = bodyIndexBegin, bodyIndexEnd

      if ( obs_columnActive_RB(obsdat, OBS_OMA) )  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMA0))  call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMP) )  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMP6))  call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OER) )  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_HPHT))  call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_HAHT))  call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_WORK))  call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_SIGI))  call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_SIGO))  call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_ZHA ))  call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, missingValue )

    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD data (element 15031)
    if ( trim(familyType) == 'GP') then
      write(*,*)' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError(obsdat, headerIndexBegin, headerIndexEnd )
    end if

    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) 'brpf_readFile: obs_numheader', numHeader
    write(*,*) 'brpf_readFile: obs_numbody  ', numBody

  end subroutine brpf_readFile

  !--------------------------------------------------------------------------
  ! brpf_updateFile
  !--------------------------------------------------------------------------
  subroutine brpf_updateFile(obsSpaceData,fileName,familyType,fileIndex)
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*) :: fileName
    character(len=*) :: familyType
    integer :: fileIndex

    ! locals
    integer :: headerIndex

    call utl_tmg_start(12,'----UpdateBurpFile')

    write(*,*)
    write(*,*) 'brpf_updateFile: Starting'

    ! CH family: Scaling of the obs related values to be stored in the BURP files
    if (familytype == 'CH') then  
      call obs_set_current_header_list(obsSpaceData,'CH')
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER
        call brpf_setScaleCH(obsSpaceData,headerIndex,forward=.false.)
      end do HEADER
    end if
    
    call brpr_updateBurp(obsSpaceData,familyType,fileName,fileIndex)

    write(*,*) 'brpf_updateFile: Done'
    write(*,*)

    call utl_tmg_stop(12)

  end subroutine brpf_updateFile

  !--------------------------------------------------------------------------
  ! brpf_setScaleCH
  !--------------------------------------------------------------------------
  subroutine brpf_setScaleCH(obsdat,headerIndex,forward)
    !
    ! :Purpose:  Apply or unapply scaling to CH observations  by multiplying
    !            (or dividing) with 10^{exponent} where the exponent is from
    !            element BUFR_SCALE_EXPONENT if provided.           
    !
    !
      implicit none

      ! Arguments      
      type (struct_obs), intent(inout):: obsdat ! struct_obs instance
      integer, intent(in) :: headerIndex ! header index in obsdat
      logical, intent(in) :: forward     ! applies scaling if .true., unapplies scaling if .false.

      ! Locals
      integer  :: bodyIndex,rln,nlv

      real(pre_obsReal) :: obsv
      real(pre_obsReal) :: vomp, voma, voer, vhpht, scale
      integer        :: nexp,iobs,iexp

      real(pre_obsReal), allocatable :: expnt(:)

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

      if (nexp*2 /= nlv) then
         ! Skip over obs assuming mantissa was filtered out in brpr_readBurp 
         ! (not inserted in obsSpaceData) due to quality flags.
         ! Set exponent quality flag to that of a 'Suspicious element' 
         
         do bodyIndex = RLN, NLV + RLN -1
            call obs_bodySet_r(obsdat,OBS_VAR,bodyIndex, 0.0D0 )
            call obs_bodySet_i(obsdat,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsdat,OBS_FLG,bodyIndex),02) )
            call obs_bodySet_i(obsdat,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsdat,OBS_FLG,bodyIndex),04) )
            call obs_bodySet_i(obsdat,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsdat,OBS_FLG,bodyIndex),09) )
         end do
              
         ! write(*,*) 'NLV =',nlv,' Nexp=',nexp    
         ! call utl_abort('brpf_setScaleCH: Inconsistent number of exponents')
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

  end subroutine brpf_setScaleCH

  !--------------------------------------------------------------------------
  ! brpf_obsSub_read
  !--------------------------------------------------------------------------
  function brpf_obsSub_read(filename,stnid,varno,nlev,ndim,block_type,bkstp_opt, &
                            match_nlev_opt,codtyp_opt,numColumns_opt) result(burp_out)
    !
    !:Purpose: To retrieve information from observation BURP file. Returns the
    !          data in a struct_oss_obsdata object. Can retrieve either 1D or 2D
    !          data from a report.
    !
    !:Comments:
    !
    !    - BUFR power 10 exponent element (i.e. data with BUFR number
    !      BUFR_SCALE_EXPONENT) will be applied only to 1D data if present.
    !    - As burp_out is for a specific input stnid, burp_out%code contains
    !      only the (lat/long and time coord.) with 22 characters.
    !    - Exponent BUFR data (i.e. data with BUFR number BUFR_SCALE_EXPONENT)
    !      will be applied only to 1D data.
    !
    !:Arguments:
    !   :filename:        observation family name
    !   :stnid:           station ID of observation
    !   :varno:           BUFR code (if <=0, search through all codes to obtain first
    !                     between 10000 and 16000)
    !   :nlev:            number of levels in the observation (number of rows) 
    !   :ndim:            number of dimensions for the retrieved data in each report
    !                     (e.g. ndim=1 for std, ndim=2 for averaging kernels)
    !   :numColumns_opt:  Number of columns (if different from nlev and for ndim=2)
    !   :block_type:      block type indicated by the two rightmost bits of bknat.
    !                     Valid values are 'DATA', 'INFO', '3-D', and 'MRQR'.
    !   :match_nlev_opt: =.true. (default) causes filtering out of report if
    !                    the report number of levels is different from the input
    !                    argument nlev
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: filename ! BURP file name
    character(len=9), intent(in)  :: stnid ! station ID of observation
    integer, intent(in)           :: varno
    integer, intent(in)           :: nlev  ! number of levels in the observation
    integer, intent(in)           :: ndim
    character(len=4), intent(in)  :: block_type
    integer, intent(in), optional :: numColumns_opt ! Number of columns (if different from nlev and for ndim=2)
    integer, intent(in), optional :: bkstp_opt ! bkstp number of requested block
    logical, intent(in), optional :: match_nlev_opt
    integer, intent(in), optional :: codtyp_opt(:) ! optional CODTYP list for search
    type(struct_oss_obsdata) :: burp_out   ! struct_oss_obsdata object
    
    ! Locals:
    character(len=9)  :: rep_stnid
    type(burp_file)   :: brp
    type(burp_rpt)    :: rep
    type(burp_block)  :: blk
    integer           :: error,ref_rpt,nrep,ref_blk,varno_ivar
    integer           :: ref_bkstp,nval,ivar,ilev,icount,icodtyp
    integer           :: date,time,ilat,ilon,iele,nele,icol,ivar_previous
    logical           :: match_nlev
    real(8)           :: exponent

    if (present(match_nlev_opt)) then
       match_nlev = match_nlev_opt
    else
       match_nlev = .true.
    end if

    ! initialize burp file, report, and block
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_READ, IOSTAT=error)
    if (error /= 0) call utl_abort('brpf_obsSub_read: Could not find/open BURP file: ' // trim(filename))

    write(*,*) "brpf_obsSub_read: Reading file " // trim(filename)
    write(*,*) "brpf_obsSub_read: Selecting STNID = ",stnid," BUFR = ",varno," block type = ",block_type
    write(*,*) "brpf_obsSub_read:           nlev = ",nlev," match_nlev = ",match_nlev
    if (present(bkstp_opt))  write(*,*) "brpf_obsSub_read:           bkstp  = ",bkstp_opt
    if (present(codtyp_opt)) write(*,*) "brpf_obsSub_read:           codtyp = ",codtyp_opt(:)

    ! get number of reports in file
    call BURP_Get_Property(brp, NRPTS=nrep)

    ! allocate memory
    if (ndim == 1) then
      call oss_obsdata_alloc(burp_out,nrep,dim1=nlev)
    else
      if (present(numColumns_opt)) then
        call oss_obsdata_alloc(burp_out,nrep,dim1=nlev,dim2_opt=numColumns_opt)
      else
        call oss_obsdata_alloc(burp_out,nrep,dim1=nlev,dim2_opt=nlev)
      end if
    end if
    
    icount = 0  ! counter of reports with same stnid, number of levels, and varno as input 
    ref_rpt = 0
    
    ! loop through reports    
    REPORTS: do

       ref_rpt = BURP_Find_Report(brp, REPORT=rep, SEARCH_FROM=ref_rpt, IOSTAT=error)

       if (ref_rpt<0) exit REPORTS
       
       call BURP_Get_Property(rep, STNID=rep_stnid, DATE=date, TEMPS=time, LATI=ilat, LONG=ilon, IDTYP=icodtyp) 

       if (present(codtyp_opt)) then
          if (.not.any(codtyp_opt(:) == icodtyp)) cycle REPORTS
       end if

       if (.not.utl_stnid_equal(stnid,rep_stnid)) cycle REPORTS

       ! loop through blocks
       ref_blk = 0
       BLOCKS: do
          
          ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
          if (ref_blk<0) exit BLOCKS
          
          call BURP_Get_Property(blk, NELE=nele, NVAL=nval, BKSTP=ref_bkstp, IOSTAT=error)

          if (.not.IS_Burp_Btyp(trim(block_type),BLOCK=blk)) cycle BLOCKS
          if (match_nlev .and. nval /= nlev) cycle BLOCKS
          if (present(bkstp_opt)) then
             if (bkstp_opt /= ref_bkstp) cycle BLOCKS
          end if

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
             if (varno_ivar < 10000.or.varno_ivar >= 16000) call utl_abort('brpf_obsSub_read: No valid element found for STNID ' // rep_stnid )
          end if

          ! required block found if code reaches this point, retrieve data and store in burp_out
          
          if (nval > nlev) call utl_abort('brpf_obsSub_read: number of levels in the report (' // trim(utl_str(nval)) // &
                                         ') exceeds the specified maximum number of levels (' // trim(utl_str(nlev)) // &
                                         ') for STNID ' // rep_stnid )

          icount=icount+1
          burp_out%code(icount) = oss_obsdata_get_header_code(ilon,ilat,date,time,rep_stnid)  ! this code is a unique identifier for this report

          if (ndim == 1) then
             ! retrieve 1D data

             do ilev=1,nval                   
               burp_out%data1d(ilev,icount) = BURP_Get_Rval(blk, NELE_IND=ivar, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)                
             end do
             
             if (ivar < nele) then             
               if (BURP_Get_Element(blk, INDEX=ivar+1, IOSTAT=error) == BUFR_SCALE_EXPONENT) then
                 ! Apply exponent
                 do ilev=1,nval                   
                   exponent = BURP_Get_Rval(blk, NELE_IND=ivar+1, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   burp_out%data1d(ilev,icount) = burp_out%data1d(ilev,icount) * 10**exponent
                 end do
               end if
             end if
                   
          else if (ndim == 2) then
             ! retrieve 2D data

             icol = 0
             ivar_previous=0
             do iele=1,nele
                ivar = BURP_Get_Element(blk, INDEX=iele, IOSTAT=error)
                if (ivar == varno) then
                   icol = icol+1
                   ivar_previous=ivar
                   do ilev=1,nval
                      burp_out%data2d(ilev,icol,icount) = BURP_Get_Rval(blk, NELE_IND=iele, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   end do                   
                else if (ivar_previous == varno .and. ivar == BUFR_SCALE_EXPONENT) then
                   ivar_previous=0
                   do ilev=1,nval
                      exponent = BURP_Get_Rval(blk, NELE_IND=iele, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                      burp_out%data2d(ilev,icol,icount) = burp_out%data2d(ilev,icol,icount) * 10**exponent
                   end do   
                else
                   ivar_previous=0                
                end if  
             end do
             if (present(numColumns_opt)) then
               if (icol /= numColumns_opt) call utl_abort('brpf_obsSub_read: number of columns (' // trim(utl_str(icol)) // &
                                         ') is not equal to the required number (' // trim(utl_str(numColumns_opt)) // &
                                         ') for STNID ' // rep_stnid )
             else
               if (icol > nlev ) call utl_abort('brpf_obsSub_read: number of columns (' // trim(utl_str(icol)) // &
                                         ') exceeds the maximum number (' // trim(utl_str(nlev)) // &
                                         ') for STNID ' // rep_stnid )
             end if
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

    write(*,*) "brpf_obsSub_read: Reading of file complete. Number of reports found: ",burp_out%nrep
    
    ! deallocate
    Call BURP_Free(brp,iostat=error)
    Call BURP_Free(rep,iostat=error)
    Call BURP_Free(blk,iostat=error)
    
  end function brpf_obsSub_read

  !--------------------------------------------------------------------------
  ! brpf_obsSub_update
  !--------------------------------------------------------------------------
  function brpf_obsSub_update(obsdata,filename,varno,block_type,bkstp_opt, &
                              multi_opt) result(nrep_modified)
    !
    !:Purpose: To add or modify data in BURP files from data stored in a
    !          struct_oss_obsdata object. Provided data can be either 1D or 2D.
    !
    !:Comments:
    !   - Currently assumes that all elements of varno(:) are distinct from
    !     each other.
    !   - Blocks with new data will overwrite any existing data of the same
    !     varno. Otherwise the new data will be appended to the block.
    !
    !:Arguments:
    !   :varno:        BUFR descriptors. Number of elements must be 
    !                  max(1,obsdata%dim2)
    !   :block_type:   block type indicated by the two rightmost bits of bknat.
    !                  Valid values are 'DATA', 'INFO', '3-D', and 'MRQR'.
    !   :multi_opt:    Indicates if intended report are for 'UNI' or 'MULTI'
    !                  level data (description is not accurate)
    !     
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata ! Input struct_oss_obsdata object for varno
    character(len=*), intent(in) :: filename ! BURP file name
    character(len=4), intent(in) :: block_type
    integer, intent(in) :: varno(:)
    integer, intent(in), optional :: bkstp_opt ! bkstp number of requested block
    character(len=*), intent(in), optional :: multi_opt
    !Output:
    integer :: nrep_modified ! Number of modified reports

    ! Locals:
    integer :: ncount
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
       write(*,*) 'brpf_obsSub_update: Skipped due to absence of data to update.'
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
    
    if (size(varno) < dim2) call utl_abort('brpf_obsSub_update: Number of BUFR elements not sufficient. ' // &
                                          trim(utl_str(size(varno))) // ' vs ' // trim(utl_str(dim2)))

    if (code_len < oss_obsdata_code_len()) call utl_abort('brpf_obsSub_update: Length of code string' &
                                          // ' needs to be increased to ' // trim(utl_str(oss_obsdata_code_len())))
     
    ! initialize burp file, report, and block system resources
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, R2=rep_new, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file in append mode (to replace or add data in a block)
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_APPEND, IOSTAT=error)
    if (error /= 0) call utl_abort('brpf_obsSub_update: Could not open BURP file: ' // trim(filename))

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

             blk_found = IS_Burp_Btyp(trim(block_type),BLOCK=blk) .and. dim1 == nval
             if (present(bkstp_opt)) blk_found = blk_found .and. bkstp_opt == ref_bkstp

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
          
          ! The call to BURP_Write_Block has ENCODE_BLOCK and CONVERT_BLOCK set
          ! to .true. in all cases, including when the block has not been
          ! modified, due to problems that can occur when writing blocks
          ! containing negative integers with datyp=4.
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
    
  end function brpf_obsSub_update

end module burpFiles_mod
