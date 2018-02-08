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
!! MODULE utilities_mod (prefix="utl")
!!
!! *Purpose*: A place to collect numerous simple utlity routines
!!
!--------------------------------------------------------------------------
module utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: utl_EZUVINT, utl_EZUVINT2, utl_EZGDEF, utl_CXGAIG, utl_FSTLIR, utl_FSTECR
  public :: utl_EZSINT, utl_EZSINT2, utl_findArrayIndex, utl_matSqrt
  public :: utl_writeStatus, utl_getfldprm, utl_abort, utl_checkAllocationStatus
  public :: utl_open_asciifile, utl_stnid_equal, utl_resize, utl_str
  public :: utl_get_stringId, utl_get_Id
  public :: utl_readFstField
  public :: utl_varNamePresentInFile

  ! module interfaces
  ! -----------------

  ! interface for resizing arrays
  interface utl_resize
     module procedure utl_resize_1d_real
     module procedure utl_resize_1d_int
     module procedure utl_resize_1d_str
     module procedure utl_resize_2d_real
     module procedure utl_resize_3d_real
  end interface utl_resize

  ! interface for conversion to a left-justified string (useful for calls to utl_abort)
  interface utl_str
     module procedure utl_int2str
     module procedure utl_float2str
  end interface utl_str

contains

  !--------------------------------------------------------------------------
  ! utl_EZUVINT2
  !--------------------------------------------------------------------------
  function utl_EZUVINT2(duuout, dvvout, duuin, dvvin, nio, nii) result(ierr)
    IMPLICIT NONE

    real(8) :: duuout(nio), dvvout(nio)
    real(4) :: duuin(nii) , dvvin(nii)
    integer :: iun, nio, nii

    integer :: ikey, ierr, ileni, ileno, jk1
    real, allocatable :: bufuuout4(:), bufvvout4(:)

    integer :: ezuvint

    allocate(bufuuout4(nio))
    allocate(bufvvout4(nio))

    ierr = ezuvint(bufuuout4, bufvvout4, duuin, dvvin)

    do jk1 = 1,nio
       duuout(jk1) = bufuuout4(jk1)
       dvvout(jk1) = bufvvout4(jk1)
    enddo

    deallocate(bufuuout4)
    deallocate(bufvvout4)

  end function utl_EZUVINT2

  !--------------------------------------------------------------------------
  ! utl_EZUVINT
  !--------------------------------------------------------------------------
  function utl_EZUVINT(duuout, dvvout, duuin, dvvin, nio, nii) result(ierr)
    IMPLICIT NONE

    real(8) :: duuout(nio), dvvout(nio)
    real(8) :: duuin(nii) , dvvin(nii)
    integer :: iun, nio, nii

    integer :: ikey, ierr, ileni, ileno, jk1
    real, allocatable :: bufuuout4(:), bufvvout4(:)
    real, allocatable :: bufuuin4(:), bufvvin4(:)

    integer :: ezuvint

    allocate(bufuuout4(nio))
    allocate(bufvvout4(nio))
    allocate(bufuuin4(nii))
    allocate(bufvvin4(nii))

    do jk1 = 1,nii
      bufuuin4(jk1) = duuin(jk1)
      bufvvin4(jk1) = dvvin(jk1)
    enddo

    ierr = ezuvint(bufuuout4, bufvvout4, bufuuin4, bufvvin4)

    do jk1 = 1,nio
       duuout(jk1) = bufuuout4(jk1)
       dvvout(jk1) = bufvvout4(jk1)
    enddo

    deallocate(bufuuin4)
    deallocate(bufvvin4)
    deallocate(bufuuout4)
    deallocate(bufvvout4)

  end function utl_EZUVINT

  !--------------------------------------------------------------------------
  ! utl_EZGDEF
  !--------------------------------------------------------------------------
  FUNCTION utl_EZGDEF(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, ax, ay) result(vezgdef)
    IMPLICIT NONE

    integer :: vezgdef

    integer :: ni, nj, ig1, ig2, ig3, ig4
    real(8) :: ax(*), ay(*)
    character(len=*) :: grtyp, grtypref

    integer :: ier1,ier2,jk,ilenx,ileny
    real, allocatable :: bufax4(:), bufay4(:)

    integer :: ezgdef

    if      (grtyp .eq. 'Y') then
       ilenx=max(1,ni*nj)
       ileny=ilenx
    else if (grtyp .eq. 'Z') then
       ilenx=max(1,ni)
       ileny=max(1,nj)
    else
       call utl_abort('VEZGDEF: Grid type not supported')
    endif

    allocate(bufax4(ilenx))
    allocate(bufay4(ileny))

    do jk = 1,ilenx
       bufax4(jk) = ax(jk)
    enddo
    do jk = 1,ileny
       bufay4(jk) = ay(jk)
    enddo

    ier2 = ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, &
         bufax4, bufay4)

    deallocate(bufax4)
    deallocate(bufay4)

    vezgdef=ier2

  end FUNCTION utl_EZGDEF

  !--------------------------------------------------------------------------
  ! utl_CXGAIG
  !--------------------------------------------------------------------------
  SUBROUTINE utl_CXGAIG(grtyp, ig1, ig2, ig3, ig4, xlat0, xlon0, dlat, dlon) 
    IMPLICIT NONE

    integer :: ig1, ig2, ig3, ig4   
    real(8) :: xlat0, xlon0, dlat, dlon 
    character(len=*) :: grtyp 

    real(4) :: xlat04, xlon04, dlat4, dlon4

    xlat04=xlat0
    xlon04=xlon0
    dlat4=dlat
    dlon4=dlon

    call cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat04, xlon04, dlat4, dlon4)

  end SUBROUTINE utl_CXGAIG

  !--------------------------------------------------------------------------
  ! utl_fstlir
  !--------------------------------------------------------------------------
  FUNCTION utl_FSTLIR(fld8, iun, ni, nj, nk, datev, etiket, &
       ip1, ip2, ip3, typvar, nomvar) result(vfstlir)
    IMPLICIT NONE

    integer :: vfstlir
    real(8) :: fld8(*)
    integer :: iun, ni, nj, nk, datev, ip1, ip2, ip3
    character(len=*) :: etiket
    character(len=*) :: nomvar
    character(len=*) :: typvar

    integer :: key1,key2, ierr, ilen, jk1, jk2, jk3, la
    real(4), allocatable :: buffer4(:)

    integer :: fstluk, fstinf

    !     Get field dimensions and allow memory for REAL copy of fld8.
    key1 = fstinf(iun, ni, nj, nk, datev, etiket, &
         ip1, ip2, ip3, typvar, nomvar)

    if(key1 >= 0) then
       ilen = ni*nj*nk
       allocate(buffer4(ilen))
       !     Read field
       key2 = fstluk(buffer4, key1, ni, nj, nk)
       if(key2 >= 0) then
          do jk3 = 1,nk
             do jk2 = 1,nj
                do jk1 = 1,ni
                   la=jk1+(jk2-1)*ni+(jk3-1)*ni*nj
                   fld8(la) = buffer4(la)
                enddo
             enddo
          enddo
       endif

       deallocate(buffer4)
    endif

    vfstlir=key1

  end FUNCTION utl_FSTLIR

  !--------------------------------------------------------------------------
  ! utl_fstecr
  !--------------------------------------------------------------------------
  FUNCTION utl_FSTECR(fld8, npak, iun, dateo, deet, &
       npas, ni, nj, nk, ip1, ip2, ip3, typvar, &
       nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, & 
       datyp, rewrit) result(vfstecr)
    IMPLICIT NONE

    integer :: vfstecr
    real(4) :: work
    integer, intent(in) :: ni,nj,nk
    real(8) :: fld8(ni,nj,nk)
    integer :: iun, datev, ip1, ip2, ip3, ig1, ig2, ig3, ig4
    integer :: npak, dateo, deet, npas, datyp
    logical :: rewrit  
    character(len=*) :: etiket 
    character(len=*) :: typvar
    character(len=*) :: grtyp 
    character(len=*) :: nomvar            

    integer :: ikey, ierr, jk1, jk2, jk3
    real(4), allocatable :: buffer4(:,:,:)

    integer :: fstecr

    allocate(buffer4(ni,nj,nk))

    do jk3 = 1,nk
       do jk2 = 1,nj
          do jk1 = 1,ni
             buffer4(jk1,jk2,jk3) = fld8(jk1,jk2,jk3)
          enddo
       enddo
    enddo

    ikey = fstecr(buffer4, work, npak, iun, dateo, deet, &
         npas, ni, nj, nk, ip1, ip2, ip3, typvar, nomvar, & 
         etiket, grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

    deallocate(buffer4)

    vfstecr=ikey

  end FUNCTION UTL_FSTECR

  !--------------------------------------------------------------------------
  ! utl_ezsint
  !--------------------------------------------------------------------------
  function utl_EZSINT(zout8, zin8, nio, njo, nko, nii, nji, nki) result(ierr)
    IMPLICIT NONE

    integer :: nii, nji, nki, nio, njo, nko     
    real(8) :: zout8(nio,njo,nko),zin8(nii,nji,nki)

    integer :: ierr, jk1, jk2, jk3
    real(4), allocatable :: bufferi4(:,:,:), buffero4(:,:,:)

    integer :: ezsint

    allocate(bufferi4(nii,nji,nki))
    allocate(buffero4(nio,njo,nko))

    do jk3 = 1,nki
       do jk2 = 1,nji
          do jk1 = 1,nii
             bufferi4(jk1,jk2,jk3) = zin8(jk1,jk2,jk3)
          enddo
       enddo
    enddo

    ierr = ezsint(buffero4,bufferi4)

    do jk3 = 1,nko
       do jk2 = 1,njo
          do jk1 = 1,nio
             zout8(jk1,jk2,jk3) = buffero4(jk1,jk2,jk3)
          enddo
       enddo
    enddo

    deallocate(bufferi4)
    deallocate(buffero4)

  end function utl_EZSINT

  !--------------------------------------------------------------------------
  ! utl_ezsint2
  !--------------------------------------------------------------------------
  function utl_EZSINT2(zout8, zin4, nio, njo, nko, nii, nji, nki)  result(ierr)
    IMPLICIT NONE

    integer :: nii, nji, nki, nio, njo, nko     
    real(8) :: zout8(nio,njo,nko)
    real(4) :: zin4(nii,nji,nki)

    integer :: ierr, jk1, jk2, jk3
    real(4), allocatable :: buffero4(:,:,:)

    integer :: ezsint

    allocate(buffero4(nio,njo,nko))

    ierr = ezsint(buffero4,zin4)

    do jk3 = 1,nko
       do jk2 = 1,njo
          do jk1 = 1,nio
             zout8(jk1,jk2,jk3) = buffero4(jk1,jk2,jk3)
          enddo
       enddo
    enddo

    deallocate(buffero4)

  end function utl_EZSINT2

  !--------------------------------------------------------------------------
  ! utl_findArrayIndex
  !--------------------------------------------------------------------------
  FUNCTION utl_findArrayIndex(KLIST, KLEN, KENTRY) result(isrcheq)
    implicit none
    !
    ! Find entry in list.
    !
    ! Arguments
    !     i   KLIST   : List.
    !     i   KLEN    : Dimension of list.
    !     i   KENTRY  : Entry.
    !     O   ISRCHEQ : Index of entry: (0, not found, >0, found)

    INTEGER :: ISRCHEQ
    INTEGER :: KENTRY, KLEN, JI
    INTEGER :: KLIST(KLEN)

    ISRCHEQ = 0
    DO JI=1,KLEN
       IF ( KLIST(JI) .EQ. KENTRY ) THEN
          ISRCHEQ = JI
          RETURN
       ENDIF
    ENDDO

  END FUNCTION UTL_FINDARRAYINDEX

  !--------------------------------------------------------------------------
  ! utl_matSqrt
  !--------------------------------------------------------------------------
  SUBROUTINE utl_MATSQRT(PA,KN,ZSIGN,printInformation_opt)
    !
    !**s/r MATSQRT     - Calculate square root of an error covariance
    !     .              matrix
    !
    ! Arguments
    !     .  PA(KN,KN)     :  on entry, the original matrix
    !     .                   on exit,  the sqrt     matrix
    !     .  KN            : order of the matrix
    IMPLICIT NONE
    !
    ! Arguments
    !
    INTEGER, intent(in) :: KN
    REAL(8), intent(inout) :: PA(KN,KN)
    REAL(8), intent(in) ::ZSIGN
    LOGICAL, intent(in), optional :: printInformation_opt
    !
    ! Local variables
    !
    LOGICAL :: printInformation
    INTEGER :: JI, J, INFO, IER, IWORK
    REAL(8) :: size_zwork
    REAL(8), allocatable :: ZWORK(:), ZRESULT(:,:), ZEIGENV2(:,:), ZEIGEN(:,:), ZEIGENV(:)

    if (present(printInformation_opt)) then
       printInformation = printInformation_opt
    else
       printInformation = .false.
    end if

    if (printInformation) then
       WRITE(*,*)' MATSQRT-Sqrt matrix of a symmetric matrix'
       WRITE(*,*)' zsign= ',zsign
       WRITE(*,*)'  -----------------------------------------------'
    end if
    !
    !     1. Computation of eigenvalues and eigenvectors
    !
    allocate(ZRESULT(KN,KN))
    allocate(ZEIGEN(KN,KN))
    allocate(ZEIGENV2(KN,KN))
    allocate(ZEIGENV(KN))

    DO JI=1,KN
       DO J=1,KN
          ZEIGEN(JI,J)=PA(JI,J)
       END DO
    END DO

    ! query the size of the 'zwork' vector by calling 'DSYEV' with 'iwork=-1'
    iwork=-1
    info = -1
    CALL DSYEV('V','U',KN, ZEIGEN,KN, ZEIGENV,size_ZWORK, IWORK, INFO )

    iwork=int(size_zwork)
    allocate(zwork(iwork))
    ! compute the eigenvalues
    CALL DSYEV('V','U',KN, ZEIGEN,KN, ZEIGENV,ZWORK, IWORK, INFO )
    deallocate(zwork)
    !
    if (printInformation) then
       WRITE(*,'(1x,"ORIGINAL EIGEN VALUES: ")')
       WRITE(*,'(1x,10f7.3)') (ZEIGENV(JI),JI=1,KN)
    end if
    !
    !     2.  Take SQRT of eigenvalues
    !
    DO JI=1,KN
       DO J=1,KN
          ZEIGENV2(JI,J)= 0.0d0
       END DO
    END DO
    DO JI=1,KN
       ZEIGENV2(JI,JI)= ZEIGENV(JI)**(0.5d0*ZSIGN)
    END DO
    !
    if (printInformation) then
       WRITE(*,'(1x,"SQRT OF ORIGINAL EIGEN VALUES: ")')
       WRITE(*,'(1x,10f7.3)') (ZEIGENV2(JI,JI),JI=1,KN)
    end if
    !
    IF (ZSIGN < 0. .and. printInformation) THEN
       Write(*,'(A,1x,e14.6)') "Condition number:", &
            maxval(ZEIGENV)/minval(ZEIGENV)
    ENDIF
    CALL DGEMM('N','N',KN,KN,KN,1.0d0,ZEIGEN,KN,ZEIGENV2,KN, &
         0.0D0 ,ZRESULT,KN)

    CALL DGEMM('N','T',KN,KN,KN,1.0D0,ZRESULT,KN,ZEIGEN,KN, &
         0.0d0,PA,KN)


    !
    !     4. Deallocate local arrays
    !
    deallocate(ZRESULT,ZEIGEN,ZEIGENV2,ZEIGENV)

    if (printInformation) then
       WRITE(*,*)'MATSQRT-----------Done--------------- '
       WRITE(*,*)' '
    end if

  END SUBROUTINE UTL_MATSQRT

  !--------------------------------------------------------------------------
  ! utl_writeStatus
  !--------------------------------------------------------------------------
  subroutine utl_writeStatus(cmsg)
    implicit none
    INTEGER :: iulstatus,fnom,fclos, ierr
    character(len=*) :: cmsg
    character(len=22):: clmsg
    
    clmsg='VAR3D_STATUS='//cmsg
    iulstatus = 0
    IERR =  FNOM(iulstatus,'VAR3D_STATUS.dot','SEQ+FMT',0)
    rewind (iulstatus)
    WRITE(iulstatus,'(a22)') clmsg
    ierr = fclos(iulstatus)
    
  end subroutine utl_writeStatus

  !--------------------------------------------------------------------------
  ! utl_getfldprm
  !--------------------------------------------------------------------------
  subroutine utl_getfldprm(kip1s,kip2,kip3,knlev,cdetiket,cdtypvar,kgid, &
                           cdvar,kstampv,knmaxlev,kinmpg,kip1style,kip1kind, &
                           ktrials,koutmpg)
    implicit none

    integer :: kstampv,knmaxlev,knlev,kgid
    integer :: kip1s(knmaxlev),kip1style,kip1kind,kip2,kip3
    integer :: ktrials, koutmpg  
    integer :: kinmpg(ktrials)
    character(len=*) :: cdtypvar
    character(len=*) :: cdvar
    character(len=*) :: cdetiket

    !
    !*    Purpose: Get 3D grid parameters for a specific trial field
    !              and check for consitancies between grid parameters
    !              of the levels.
    !
    !Arguments
    !
    ! Input:
    !     cdvar   : variable name to get the vertical levels from
    !     kstampv : valid date time stamp of the variable
    !     knmaxlev: maximum number of levels
    !     kinmpg  : file unit of trial field
    !     ktrials :  number of trial files.  
    !
    ! Output:
    !     kip1s(knmaxlev) : list of ip1s of variable cdvar
    !     kip2            : ip2 for variable cdvar
    !     kip3            : ip3 for variable cdvar
    !     knlev           : number of levels of variable cdvar
    !     cdetiket        : etiket of field cdvar
    !     cdtypvar        : typvar of field cdvar
    !     kgid            : handle of the field descriptor
    !     kip1style       : style in which ip1 is encoded (15 or 31 bits)
    !     kip1kind        : kind of vertical coord encoded in ip1
    !     koutmpg         : the unit which contains the selected records.  
    !

    integer :: fstinl,fstprm,ezqkdef,newdate
    integer :: ini,inj,ink,jlev,ier
    integer :: idateo, idateo2, idatyp, idatyp2, ideet, ideet2, idltf, &
         iextra1, iextra2, iextra3, iig12, iig22, &
         iig32, iig42, ilng, inbits,iig1,iig2,iig3,iig4, &
         inpas,inpas2, iswa, iubc, iip2, iip3
    !
    integer :: ipmode,idate2,idate3,idatefull
    integer :: k,ier1 
    real(4) :: zlev_r4
    character(len=12) :: cletiket
    character(len=4) :: clnomvar
    character(len=3) :: clnomvar_3
    character(len=2) :: cltypvar
    character(len=1) :: clgrtyp2,clgrtyp,clstring
    logical :: llflag
    integer :: ikeys(knmaxlev)
    !
    knlev = 0
    !
    do k=1,ktrials
       if(cdvar.eq.'U1') then
          clnomvar_3='UT1'
          ier = fstinl(kinmpg(k),INI,INJ, INK, kstampv, ' ', -1, -1, -1, &
               ' ',clnomvar_3,ikeys, knlev, knmaxlev)
       else if(cdvar.eq.'V1') then
          clnomvar_3='VT1'
          ier = fstinl(kinmpg(k),ini,inj, ink, kstampv, ' ', -1, -1, -1, &
               ' ',clnomvar_3,ikeys, knlev, knmaxlev)
       else
          ier = fstinl(kinmpg(k),INI,INJ, INK, kstampv, ' ', -1, -1, -1, &
               ' ',cdvar,IKEYS, KNLEV, knmaxlev)
       endif
       !
       if(knlev > 0 ) then
          ier1   = newdate(kstampv,idate2,idate3,-3)

          idatefull = idate2*100 + idate3/1000000
          idateo = -9999
          ideet = -9999
          inpas = -9999
          cdetiket = '-9999999'
          clgrtyp = '-'
          kip2 = -9999
          kip3 = -9999
          cdtypvar = '-'
          idatyp = -9999
          iig1 = -9999
          iig2 = -9999
          iig3 = -9999
          iig4 = -9999
          llflag = .true.
          koutmpg = kinmpg(k) 
          exit 
       endif
    enddo ! End of loop k   
    !
    if (knlev.gt.0) then
       do jlev = 1, knlev
          ier = fstprm(ikeys(jlev), idateo2, ideet2, inpas2, ini, inj, &
               ink,inbits,idatyp2, kip1s(jlev),iip2, iip3, &
               cltypvar,clnomvar,cletiket,clgrtyp2, iig12, iig22,iig32 &
               ,iig42,iswa,ilng,idltf,iubc,iextra1, iextra2, iextra3)
          llflag = (llflag.and.(idateo.eq.idateo2.or.idateo.eq.-9999))
          llflag = (llflag.and.(ideet.eq.ideet2.or.ideet.eq.-9999))
          llflag = (llflag.and.(inpas.eq.inpas2.or.inpas.eq.-9999))
          !          llflag = (llflag.and.(cdetiket.eq.cletiket.or.cdetiket.eq.
          !     &         '-9999999'))
          llflag = (llflag.and.(clgrtyp.eq.clgrtyp2.or.clgrtyp.eq.'-'))
          llflag = (llflag.and.(kip2.eq.iip2.or.kip2.eq.-9999))
          llflag = (llflag.and.(kip3.eq.iip3.or.kip3.eq.-9999))
          llflag = (llflag.and.(cdtypvar.eq.cltypvar.or.cdtypvar.eq.'-'))
          llflag = (llflag.and.(idatyp.eq.idatyp2.or.idatyp.eq.-9999))
          llflag = (llflag.and.(iig1.eq.iig12.or.iig1.eq.-9999))
          llflag = (llflag.and.(iig2.eq.iig22.or.iig2.eq.-9999))
          llflag = (llflag.and.(iig3.eq.iig32.or.iig3.eq.-9999))
          llflag = (llflag.and.(iig4.eq.iig42.or.iig4.eq.-9999))
          if (llflag) then
             idateo = idateo2
             ideet = ideet2
             inpas = inpas2
             cdetiket = cletiket
             clgrtyp = clgrtyp2
             kip2 = iip2
             kip3 = iip3
             cdtypvar = cltypvar
             idatyp = idatyp2
             iig1 = iig12
             iig2 = iig22
             iig3 = iig32
             iig4 = iig42
          else
             write(*,*) &
                  '****** Unit ', kinmpg &
                  ,' contains mixed dateo,deet,npas,etiket,grtyp,ip2,ip3' &
                  ,',typvar,datyp,ig1,ig2,ig3 and/or ig4 ' &
                  ,'for variable ',cdvar,' and datev, ',kstampv
             call utl_abort('GETFLDPRM2')
          endif
       enddo
       !
       kgid = ezqkdef(ini,inj,clgrtyp,iig1,iig2,iig3,iig4,koutmpg)
       !
       !-------Determine the style in which ip1 is encoded (15bits or 31 bits)
       !       A value <= 32767 (2**16 -1)  means that ip1 is compacted in 15 bits
       !       Determine the type of P which was encoded in IP1
       !
       if(kip1s(1) .le. 32767) then
          kip1style = 3
       else
          kip1style = 2
       endif
       !
       !-------Determine the type of P  (see doc. of convip)
       !
       ipmode = -1
       call CONVIP(kip1s(1),zlev_r4,KIP1KIND, &
            ipmode,clstring, .false. )
    else
       do k=1,ktrials
          ier = fstinl(kinmpg(k),ini,inj, ink, -1, ' ', -1, -1, -1, &
               ' ',cdvar,ikeys, knlev, knmaxlev)
       enddo
       write(*,*) 'Error - getfldprm2: no record found at time ' &
            ,idatefull,' for field ',cdvar,' but',knlev, &
            ' records found in unit ',kinmpg(k)
       call utl_abort('GETFLDPRM2')
    endif
    !
  end subroutine utl_getfldprm

  !--------------------------------------------------------------------------
  ! utl_abort
  !--------------------------------------------------------------------------
  subroutine utl_abort(message)
    
    implicit none
    character(len=*) :: message
    integer :: comm, ierr, rpn_comm_comm
    
    write(6,9000) message
9000 format(//,4X,"!!!---ABORT---!!!",/,8X,"MIDAS stopped in ",A)
    call flush(6)

    comm = rpn_comm_comm("WORLD")
    call mpi_abort( comm, 1, ierr )

  end subroutine utl_abort

  !--------------------------------------------------------------------------
  ! utl_open_asciifile
  !--------------------------------------------------------------------------
  subroutine utl_open_asciifile(filename,unit)
  ! 
  !  Purpose: Opens an ascii file for output 
  !
  !  Author: M. Sitwell, April 2016
  !
  !  Input
  !  
  !           filename      filename
  !           unit          unit number to use or 0 to let fnom set value
  !  Output
  !           unit          unit number associated with file
  !
  !-------------------------------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(out) :: unit
    logical :: file_exists
    integer :: ier
    integer :: fnom
    character(len=20) :: mode
    
    inquire(file=trim(filename), exist=file_exists)
    
    if (file_exists) then
       mode = 'FTN+APPEND+R/W'
    else
       mode = 'FTN+R/W'
    end if

    unit=0
    
    !ier = fnom(unit,trim(filename),trim(mode),0)
    ier = utl_open_file(unit,trim(filename),trim(mode))

    if (ier.ne.0) call utl_abort('utl_open_messagefile: Error associating unit number')

  end subroutine utl_open_asciifile

  !--------------------------------------------------------------------------
  ! utl_open_file
  !--------------------------------------------------------------------------
  function utl_open_file(unit,filename,mode) result(ier)
  ! 
  ! Purpose: This is a temporary subroutine to open a file with fnom that is needed due to
  !           a bug in fnom that does not allow an ascii file to be opened in 'APPEND' mode.  
  !
  !  Author: M. Sitwell, Aug 2016
  !
  !-------------------------------------------------------------------------------------------

    implicit none

    integer, intent(inout) :: unit
    character(len=*) :: filename,mode
    integer :: ier
    character(len=10) :: position,action
    integer :: fnom

    if (index(mode,'APPEND').gt.0) then
       position = 'APPEND'
    else
       position = 'ASIS'
    end if
    
    if (index(mode,'R/W').gt.0) then
       action = 'READWRITE'
    else
       action = 'READ'
    end if

    ier = fnom(unit,filename,mode,0)
    
    close(unit=unit)
    open(unit=unit, file=filename, position=position, action=action)

  end function utl_open_file
    
  !--------------------------------------------------------------------------
  ! utl_stnid_equal
  !--------------------------------------------------------------------------
  function utl_stnid_equal(id1,id2) result(same)
  !
  ! Author  : Y. Rochon  Nov 2014
  ! Revision: 
  !           M. Sitwell, Feb 2015
  !           - Code set as a function.
  !           Y. Rochon, July 2015
  !           - Accounted for ilen1>ilen2 with the additional characters being *
  !           M. Sitwell, Aug 2015
  !           - Made function symmetric so that * is treated as a wildcard
  !             for both arguments
  !
  ! Purpose: Compares STNID values allowing for * as wildcards and trailing blanks 
  !
  ! Input
  !          id1         reference stnid
  !          id2         stnid being verified
  !
  ! Output
  !
  !          same        logical indicating if id1 and id2 match
  !     
  !-----------------------------------------------------------------------------------------    

    implicit none

    logical :: same
    CHARACTER(len=*), intent(in) :: id1,id2
    integer :: ilen1,ilen2,ji

    same=.true.
    ilen1=len_trim(id1)
    ilen2=len_trim(id2)  
              
    do ji=1,min(ilen1,ilen2)
       if ( id1(ji:ji).ne.'*' .and. id2(ji:ji).ne.'*' .and. id2(ji:ji).ne.id1(ji:ji) ) then
          same = .false.
          exit
       end if
    end do
    
    if (same.and.ilen1.gt.ilen2) then
       do ji=ilen2+1,ilen1
          if (id1(ji:ji).ne.'*') then
              same=.false.
              exit
          end if
       end do
    else if (same.and.ilen2.gt.ilen1) then
       do ji=ilen1+1,ilen2
          if (id2(ji:ji).ne.'*') then
              same=.false.
              exit
          end if
       end do
    end if
        
  end function utl_stnid_equal
 
  !--------------------------------------------------------------------------
  ! utl_int2str
  !--------------------------------------------------------------------------
  character(len=20) function utl_int2str(i)
  !
  ! Author  : M. Sitwell Oct 2015
  !
  ! Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort. 
  !
  !-------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: i
    
    write(utl_int2str,*) i
    utl_int2str = adjustl(utl_int2str)
    
  end function utl_int2str
            
  !--------------------------------------------------------------------------
  ! utl_float2str
  !--------------------------------------------------------------------------
  character(len=20) function utl_float2str(x)
  !
  ! Author:  M. Sitwell April 2016 
  !
  ! Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort.
  !
  !-------------------------------------------------------------------------------------------

    implicit none

    real(8), intent(in) :: x

    write(utl_float2str,*) x
    utl_float2str = adjustl(utl_float2str)

  end function utl_float2str

  !--------------------------------------------------------------------------
  ! utl_resise_1d_real
  !--------------------------------------------------------------------------
  subroutine utl_resize_1d_real(arr,dim1)
  !
  ! Author  : M. Sitwell  April 2015
  ! Revision:  
  !
  ! Purpose: Resize 1D array
  !
  !---------------------------------------------------------------------------------------
    implicit none

    real(8), pointer, intent(inout) :: arr(:)
    integer, intent(in) :: dim1
    real(8), pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = 0.0D0
    
    deallocate(arr)

    arr => tmp
    
    nullify(tmp)

  end subroutine utl_resize_1d_real

  !--------------------------------------------------------------------------
  ! utl_resise_1d_int
  !--------------------------------------------------------------------------
  subroutine utl_resize_1d_int(arr,dim1)
  !
  ! Author  : M. Sitwell  April 2015
  !
  ! Purpose: Resize 1D array 
  !
  !---------------------------------------------------------------------------------------
     
    implicit none

    integer, pointer, intent(inout) :: arr(:)
    integer, intent(in) :: dim1
    integer, pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = 0
    
    deallocate(arr)

    arr => tmp
    
    nullify(tmp)

  end subroutine utl_resize_1d_int

  !--------------------------------------------------------------------------
  ! utl_resise_1d_str
  !--------------------------------------------------------------------------  
  subroutine utl_resize_1d_str(arr,dim1)
  !
  ! Author  : M. Sitwell  April 2016
  !
  ! Purpose: Resize 1D array
  !
  !---------------------------------------------------------------------------------------
     
    implicit none

    character(len=*), pointer, intent(inout) :: arr(:)
    integer, intent(in) :: dim1
    character(len=len(arr(1))), pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = ""
    
    deallocate(arr)
    arr => tmp
    nullify(tmp)

  end subroutine utl_resize_1d_str

  !--------------------------------------------------------------------------
  ! utl_resise_2d_real
  !--------------------------------------------------------------------------
  subroutine utl_resize_2d_real(arr,dim1,dim2)
  !
  ! Author  : M. Sitwell  April 2015
  !
  ! Purpose: Resize 2D array
  !
  ! Revision: 
  !           Y. Rochon Feb 2016
  !           - Added option to increase array sizes
  !
  !---------------------------------------------------------------------------------------

    implicit none

    real(8), pointer, intent(inout) :: arr(:,:)
    integer, intent(in) :: dim1,dim2
    real(8), pointer :: tmp(:,:)
    integer :: dim1_in,dim2_in,d1,d2

    dim1_in = size(arr,dim=1)
    dim2_in = size(arr,dim=2)
    d1 = min(dim1_in, dim1)
    d2 = min(dim2_in, dim2)

    allocate(tmp(dim1,dim2))
    tmp(1:d1,1:d2) = arr(1:d1,1:d2)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1,:) = 0.0D0
    if (dim2.gt.dim2_in) tmp(:,d2+1:dim2) = 0.0D0
      
    deallocate(arr)

    arr => tmp
    
    nullify(tmp)

  end subroutine utl_resize_2d_real

  !--------------------------------------------------------------------------
  ! utl_resise_3d_real
  !--------------------------------------------------------------------------
  subroutine utl_resize_3d_real(arr,dim1,dim2,dim3)
  !
  ! Author  : M. Sitwell  May 2015
  ! Revision: 
  !
  ! Purpose: Resize 3D array
  !
  !---------------------------------------------------------------------------------------

    implicit none

    real(8), pointer, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: dim1,dim2,dim3
    real(8), pointer :: tmp(:,:,:)
    integer :: dim1_in,dim2_in,dim3_in,d1,d2,d3

    dim1_in = size(arr,dim=1)
    dim2_in = size(arr,dim=2)
    dim3_in = size(arr,dim=3)
    d1 = min(dim1_in, dim1)
    d2 = min(dim2_in, dim2)
    d3 = min(dim3_in, dim3)

    allocate(tmp(dim1,dim2,dim3))
    tmp(1:d1,1:d2,1:d3) = arr(1:d1,1:d2,1:d3)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1,:,:) = 0.0D0
    if (dim2.gt.dim2_in) tmp(:,d2+1:dim2,:) = 0.0D0
    if (dim3.gt.dim3_in) tmp(:,:,d3+1:dim3) = 0.0D0
    
    deallocate(arr)

    arr => tmp
    
    nullify(tmp)

  end subroutine utl_resize_3d_real

  !--------------------------------------------------------------------------
  ! utl_get_stringId
  !--------------------------------------------------------------------------
  subroutine utl_get_stringId(cstringin,nobslev,CList,NListSize,Nmax,elemId)

  ! 
  !   Purpose: Get element ID from a list of accumulating character 
  !            strings (e.g. stnids). 
  !
  !            Called by filt_topoChm in filterobs_mod.ftn90
  !
  !   Author: Y.J. Rochon, ARQI/AQRD, Feb 2015
  !    
  !   Revisions:
  !            Y.J. Rochon, ARQI/AQRD, July 2015
  !            - Account for wildcards (use if utl_stnid_equal) when present.         
  !
  !   Input:
  !
  !       Nmax            Max allowed dimension.
  !       NListSize       Input number of identified IDs (must be >=0 and <=Nmax)
  !       CList           Input list of accumulated character strings
  !                       for uni and multi-level data.
  !       cstringin       Input character string
  !       nobslev         Number of elements in profile associated to cstringin.
  !
  !   Output:
  !
  !       NListSize       Updated number of identified IDs
  !       CList           Updated list of accumulated character strings
  !       elemId          Index of cstringin within CList_chm
  !        
  !-------------------------------------------------------------------------------------------
 
    implicit none

    integer, intent(in)    :: Nmax,nobslev
    integer, intent(inout) :: NListSize
    integer, intent(out)   :: elemId
    character(len=*), intent(in)     :: cstringin
    character(len=*),  intent(inout) :: CList(Nmax)

    integer :: i
    character(len=120) :: cstring
    
    elemId=0
    if (NListSize.gt.Nmax-1) then
       call utl_abort('utl_get_stringId: Dimension error, NListSize > Nmax-1.')     
    else if (NListSize.gt.0) then
       if (nobslev.eq.1) then 
          cstring=trim(cstringin)//'U'
          do i=1,NListSize
             if (trim(cstring).eq.trim(CList(i))) then
                 elemId=i
                 exit
             end if
          end do
       else 
          cstring=trim(cstringin)       
          do i=1,NListSize
             if (trim(cstring).eq.trim(CList(i))) then
                 elemId=i
                 exit
             end if
          end do
       end if
       
       if (elemId.eq.0) then
          do i=1,NListSize
             if (utl_stnid_equal(trim(CList(i)),trim(cstring))) then
                elemId=i
                exit
             end if
          end do
       end if
    end if

    if (elemID.eq.0) then
        NListSize=NListSize+1
        elemId=NListSize
        if (nobslev.eq.1) then
           CList(NListSize)=trim(cstringin)//'U'
        else
           CList(NListSize)=trim(cstringin)
        end if
    end if
    
  end subroutine utl_get_stringId

  !--------------------------------------------------------------------------
  ! utl_get_Id
  !--------------------------------------------------------------------------
  subroutine utl_get_Id(id,IdList,NListSize,Nmax,elemId)

  ! 
  !   Purpose: Get element ID from list of accumulating integer IDs.
  !
  !   Author: Y.J. Rochon, ARQI/AQRD, Feb 2015
  !
  !   Input:
  !
  !       Nmax         Max allowed dimension.
  !       NListSize    Input number of IDs (must be >=0 and <=Nmax)
  !       IdList       Input list of accumulated IDs.
  !       id           Input id for individual obs 
  !
  !   Output:
  !
  !       NListSize    Updated number of IDs 
  !       IdList       Updated list of accumulated IDs.
  !       elemId       Index of id within List
  !     
  !-------------------------------------------------------------------------------------------
 
    implicit none

    integer, intent(in)    :: Nmax,id
    integer, intent(inout) :: NListSize,IdList(Nmax)
    integer, intent(out)   :: elemId
  
    integer :: i
    
    elemId=0
    if (NListSize.gt.Nmax-1) then
       call utl_abort('utl_get_Id: Dimension error, NListSize > Nmax-1.')     
    else if (NListSize.gt.0) then
       do i=1,NListSize
          if (id.eq.IdList(i)) then
              elemId=i
              exit
          end if
       end do
    end if

    if (elemID.eq.0) then
        NListSize=NListSize+1
        elemId=NListSize
        IdList(NListSize)=id
    end if
    
    
  end subroutine utl_get_Id
  
  !--------------------------------------------------------------------------
  ! utl_readFstField
  !------------------;-------------------------------------------------------
  subroutine utl_readFstField(fname,varName,iip1,iip2,iip3,etiketi, &
                               ni,nj,nkeys,array,xlat_opt,xlong_opt,lvls_opt,kind_opt)
  !
  ! Author  : Y. Rochon, ARQI/AQRD, Nov 2015
  !
  ! Revision: 
  !
  ! Purpose:  Read specified field from standard RPN/fst file. Could be one
  !           to all levels depending on the input iip1,iip2,iip3 values
  !
  !           Currently assumes lat/long (or Gaussian) type grids.
  !           See hco_SetupFromFile for example toward future generalizations.
  !           Generalization would require having xlat and xlong being 2D.
  !
  ! IN
  !
  !     fname    input filename
  !     varName  search nomvar
  !     iip1     search ip1
  !     iip2     search ip2
  !     iip3     search ip3
  !     etiketi  search etiket
  !
  ! OUT
  !
  !     ni,nj     ni,nj values
  !     nkeys     number of records satisfying search criteria
  !     array     (ni,nj,nkeys) data arrray
  !     xlat_opt  1D latitude array (optional)
  !     xlong_opt 1D longitude array (optional)
  !     lvls_opt  1D vertical coordinate array (optional)
  !     kind_opt  vertical coordinate type according to convip (optional)
  !-----------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: iip1,iip2,iip3
    character(len=*), intent(in) :: varName,fname,etiketi
    integer, intent(out) :: ni, nj, nkeys
    integer, intent(out), optional :: kind_opt
    real(8), intent(out), allocatable :: array(:,:,:)
    real(8), intent(out), allocatable, optional :: lvls_opt(:), xlat_opt(:),xlong_opt(:)

    integer, external :: fnom,fclos,fstouv,fstfrm,fstinl,fstlir,fstluk,fstprm
    
    real(4) :: lvl_r4
    
    logical :: Exists,levels
    character(len=1) :: string
     
    integer, parameter :: iun=0
    integer :: i,ier, kindi

    integer, parameter :: maxkeys=1000
    integer :: keys(maxkeys),ini,inj,nk
    integer :: dateo, deet, npas, nbits, datyp
    integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    integer :: ig1, ig2, ig3, ig4  
    character*1 clgrtyp
    character*2 cltypvar
    character*4 nomvar
    character*12 cletiket
    real(4), allocatable :: buffer(:,:)

    real :: xlat1_4, xlon1_4, xlat2_4, xlon2_4
    
    ! Open file
    
    inquire(file=trim(fname),exist=Exists)
    if(.not.Exists) then
      write(*,*) 'File missing=',fname
      call utl_abort('utl_read_fst_field: did not find file.')
    else
      ier=fnom(iun,trim(fname),'RND+OLD+R/O',0)
      ier=fstouv(iun,'RND+OLD')
    end if

    ! Find reports in file for specified varName and iip*.

    ier = fstinl(iun,ni,nj,nk,-1,etiketi,iip1,iip2,iip3,'',varName,keys,nkeys,maxkeys) 

    if(ier.lt.0.or.nkeys.eq.0) then
      write(*,*) 'Search field missing ',varName, ' from file ',fname
      call utl_abort('utl_read_fst_field: did not find field.')
    else if (nk.gt.1) then
      write(*,*) 'Unexpected size nk ',nk,' for ',varName,' of file ',fname 
      call utl_abort('utl_read_fst_field')      
    end if

    if (present(xlat_opt).and.present(xlong_opt)) then

       !  Get lat and long if available.

       if (allocated(xlat_opt)) deallocate(xlat_opt,xlong_opt)
       allocate(xlat_opt(nj),xlong_opt(ni),buffer(ni*nj,1))
       xlat_opt(:)=-999.
       xlong_opt(:)=-999.

       ier = fstprm(keys(1),dateo, deet, npas, ni, nj, nk, nbits,    &         
                    datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, &
                    clgrtyp, ig1, ig2, ig3,                           &
                    ig4, swa, lng, dltf, ubc, extra1, extra2, extra3)  
    
       if (ni.gt.1) then
          ier=fstlir(buffer,iun,ni,inj,nk,-1,'',ig1,ig2,ig3,'','>>')
          if (ier.ge.0) xlong_opt(:)=buffer(1:ni,1) 
       end if
       if (nj.gt.1) then
          ier=fstlir(buffer,iun,ini,nj,nk,-1,'',ig1,ig2,ig3,'','^^')
          if (ier.ge.0) xlat_opt(:)=buffer(1:nj,1)     
       end if
       deallocate(buffer)

       if ( trim(clgrtyp) == 'Z' ) then

          ! Check for rotated grid

          ier = fstprm(ier,                                         & ! IN
                  dateo, deet, npas, ni, nj, nk, nbits,             & ! OUT
                  datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, & ! OUT
                  clgrtyp, ig1, ig2, ig3,                           & ! OUT
                  ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 ) ! OUT

          call cigaxg (clgrtyp,                                     & ! IN
                     xlat1_4, xlon1_4, xlat2_4, xlon2_4,            & ! OUT
                     ig1, ig2, ig3, ig4 ) ! IN

          if ( xlat1_4 /= xlat2_4 .or. xlon1_4 /= xlon2_4 ) &
             call utl_abort('utl_readFstField: Cannot currently handle rotated grid')
          
       else if (trim(clgrtyp) /= 'G') then

         call utl_abort('utl_readFstField: Cannot currently handle grid type ' // trim(clgrtyp) )

       end if
       
    end if 
    
    ! Get vertical coordinate
    
    if (present(lvls_opt)) then
       if (allocated(lvls_opt)) deallocate(lvls_opt)
       allocate(lvls_opt(nkeys))

       do i=1,nkeys
          ier = fstprm(keys(i),dateo, deet, npas, ni, nj, nk, nbits,    &         
                    datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, &
                    clgrtyp, ig1, ig2, ig3,                           &
                    ig4, swa, lng, dltf, ubc, extra1, extra2, extra3)  
          call convip(ip1,lvl_r4,kindi,-1,string,.false.)
          lvls_opt(i)=lvl_r4
       end do
    end if
     
    if (present(kind_opt)) then
        if (present(lvls_opt)) then         
            kind_opt=kindi
        else
            kind_opt=-1
        end if
    end if

    ! Get field
    
    allocate(array(ni,nj,nkeys),buffer(ni,nj))

    do i=1,nkeys
       ier=fstluk(buffer,keys(i),ni,nj,nk)
       array(:,:,i)=buffer(:,:)
    end do
    
    deallocate(buffer)
    
    ier=fstfrm(iun)  
    ier=fclos(iun)  

  end subroutine utl_readFstField

  !--------------------------------------------------------------------------
  ! utl_checkAllocationStatus(status, message, alloc_opt)
  !--------------------------------------------------------------------------
  subroutine utl_checkAllocationStatus(status, message, alloc_opt)
    
    implicit none
    character(len=*),intent(in) :: message
    integer, intent(in) :: status(:)
    logical, optional, intent(in) :: alloc_opt 

    logical :: flag

    if ( present(alloc_opt) ) then
       flag = alloc_opt
    else
       flag = .true.
    end if

    if (any(status /= 0)) then
       if (flag) then
          Write(6,'(A)') "Memory allocation failure !"
       else
          Write(6,'(A)') "Memory deallocation failure !"
       end if
       Write(6,'(64(i8,1x))') status
       call utl_abort(message)
    end if

 end subroutine utl_checkAllocationStatus

 !--------------------------------------------------------------------------
 ! utl_varNamePresentInFile
 !--------------------------------------------------------------------------
  function utl_varNamePresentInFile(fileName,varName) result(found)
    IMPLICIT NONE

    character(len=*) :: fileName
    character(len=*) :: varName
    logical :: found

    integer :: fnom, fstouv, fstfrm, fclos, fstinf
    integer :: ni, nj, nk, key, ierr
    integer :: unit = 0

    ierr = fnom(unit,fileName,'RND+OLD+R/O',0)
    ierr = fstouv(unit,'RND+OLD')
    
    key = fstinf(unit, ni, nj, nk, -1 ,' ', -1, -1, -1, ' ', trim(varName))
    
    if ( key > 0 )  then
      found = .true.
    else
      found = .false.
    end if
    
    ierr =  fstfrm(unit)
    ierr =  fclos (unit)

  end function utl_varNamePresentInFile

end module utilities_mod
