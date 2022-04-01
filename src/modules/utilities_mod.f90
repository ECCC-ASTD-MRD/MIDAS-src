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

module utilities_mod
  ! MODULE utilities_mod (prefix='utl' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: A place to collect numerous simple utility routines
  !
  use clib_interfaces_mod
  use randomNumber_mod

  implicit none
  save
  private

  ! public procedures
  public :: utl_ezuvint, utl_ezgdef, utl_cxgaig, utl_fstlir,  utl_fstlir_r4, utl_fstecr
  public :: utl_ezsint, utl_findArrayIndex, utl_matSqrt, utl_matInverse, utl_eigenDecomp
  public :: utl_pseudo_inverse
  public :: utl_writeStatus, utl_getfldprm, utl_abort, utl_checkAllocationStatus
  public :: utl_open_asciifile, utl_stnid_equal, utl_resize, utl_str
  public :: utl_get_stringId, utl_get_Id, utl_isNamelistPresent
  public :: utl_readFstField
  public :: utl_varNamePresentInFile
  public :: utl_reAllocate
  public :: utl_heapsort2d, utl_splitString, utl_stringArrayToIntegerArray, utl_parseColumns
  public :: utl_copyFile, utl_allReduce, utl_findloc, utl_findlocs
  public :: utl_randomOrderInt
  public :: utl_tmg_start

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

  interface utl_ezsint
    module procedure utl_ezsint_r4_2d
    module procedure utl_ezsint_r4_3d
    module procedure utl_ezsint_r4_2dTo1d
    module procedure utl_ezsint_r8_2d
    module procedure utl_ezsint_r8_3d
    module procedure utl_ezsint_r8_2dTo1d
  end interface utl_ezsint

  interface utl_ezuvint
    module procedure utl_ezuvint_r4_2d
    module procedure utl_ezuvint_r4_2dTo1d
    module procedure utl_ezuvint_r8_1d
    module procedure utl_ezuvint_r8_2d
  end interface utl_ezuvint

  interface utl_reAllocate
    module procedure utl_reAllocate_char_1d
    module procedure utl_reAllocate_char_2d
    module procedure utl_reAllocate_char_3d
    module procedure utl_reAllocate_log_1d
    module procedure utl_reAllocate_log_2d
    module procedure utl_reAllocate_log_3d
    module procedure utl_reAllocate_int_1d
    module procedure utl_reAllocate_int_2d
    module procedure utl_reAllocate_int_3d
    module procedure utl_reAllocate_r4_1d
    module procedure utl_reAllocate_r8_1d
    module procedure utl_reAllocate_r4_2d
    module procedure utl_reAllocate_r8_2d
    module procedure utl_reAllocate_r4_3d
    module procedure utl_reAllocate_r8_3d
    module procedure utl_reAllocate_r4_4d
    module procedure utl_reAllocate_r8_4d
    module procedure utl_reAllocate_r4_5d
    module procedure utl_reAllocate_r8_5d
  end interface utl_reAllocate

  interface utl_findloc
    module procedure utl_findloc_char
    module procedure utl_findloc_int
  end interface utl_findloc

  interface utl_findlocs
    module procedure utl_findlocs_char
  end interface utl_findlocs

contains


  subroutine utl_setezopt(interpDegree, extrapDegree_opt)
    implicit none

    ! arguments
    character(len=*) :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    character(len=12) :: extrapDegree
    integer           :: ierr, ezsetopt, ezsetval

    if ( trim(interpDegree) /= 'LINEAR' .and. &
         trim(interpDegree) /= 'CUBIC' .and. &
         trim(interpDegree) /= 'NEAREST' ) then
      write(*,*) 'utl_setezopt: interpDegree = ', trim(interpDegree)
      call utl_abort('utl_setezopt: invalid interpolation degree')
    end if

    if ( present(extrapDegree_opt) ) then
      extrapDegree = extrapDegree_opt
    else
      extrapDegree = 'VALUE'
    end if

    ierr = ezsetopt('INTERP_DEGREE', interpDegree)
    if ( trim(extrapDegree) == 'VALUE' ) then
      ierr = ezsetval('EXTRAP_VALUE', 0.0)
    end if
    ierr = ezsetopt('EXTRAP_DEGREE', extrapDegree)

  end subroutine utl_setezopt


  function utl_ezsint_r4_3d(zout4, zin4, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: zout4(:,:,:), zin4(:,:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    ierr = ezsint(zout4,zin4)

  end function utl_ezsint_r4_3d


  function utl_ezsint_r4_2d(zout4, zin4, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: zout4(:,:), zin4(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    ierr = ezsint(zout4,zin4)

  end function utl_ezsint_r4_2d


  function utl_ezsint_r4_2dTo1d(zout4, zin4, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: zout4(:), zin4(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    ierr = ezsint(zout4,zin4)

  end function utl_ezsint_r4_2dTo1d


  function utl_ezsint_r8_3d(zout8, zin8, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: zout8(:,:,:), zin8(:,:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: nii, nji, nki, nio, njo, nko     
    integer :: jk1, jk2, jk3
    real(4), allocatable :: bufferi4(:,:,:), buffero4(:,:,:)
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(zin8,1)
    nji = size(zin8,2)
    nki = size(zin8,3)

    nio = size(zout8,1)
    njo = size(zout8,2)
    nko = size(zout8,3)

    allocate(bufferi4(nii,nji,nki))
    allocate(buffero4(nio,njo,nko))

    do jk3 = 1,nki
      do jk2 = 1,nji
        do jk1 = 1,nii
          bufferi4(jk1,jk2,jk3) = zin8(jk1,jk2,jk3)
        end do
      end do
    end do

    ierr = ezsint(buffero4,bufferi4)

    do jk3 = 1,nko
      do jk2 = 1,njo
        do jk1 = 1,nio
          zout8(jk1,jk2,jk3) = buffero4(jk1,jk2,jk3)
        end do
      end do
    end do

    deallocate(bufferi4)
    deallocate(buffero4)

  end function utl_ezsint_r8_3d


  function utl_ezsint_r8_2d(zout8, zin8, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: zout8(:,:), zin8(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: nii, nji, nio, njo     
    integer :: jk1, jk2
    real(4), allocatable :: bufferi4(:,:), buffero4(:,:)
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(zin8,1)
    nji = size(zin8,2)

    nio = size(zout8,1)
    njo = size(zout8,2)

    allocate(bufferi4(nii,nji))
    allocate(buffero4(nio,njo))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufferi4(jk1,jk2) = zin8(jk1,jk2)
      end do
    end do

    ierr = ezsint(buffero4,bufferi4)

    do jk2 = 1,njo
      do jk1 = 1,nio
        zout8(jk1,jk2) = buffero4(jk1,jk2)
      end do
    end do

    deallocate(bufferi4)
    deallocate(buffero4)

  end function utl_ezsint_r8_2d


  function utl_ezsint_r8_2dTo1d(zout8, zin8, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: zout8(:), zin8(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: nii, nji, nio
    integer :: jk1, jk2
    real(4), allocatable :: bufferi4(:,:), buffero4(:)
    integer :: ezsint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(zin8,1)
    nji = size(zin8,2)

    nio = size(zout8,1)

    allocate(bufferi4(nii,nji))
    allocate(buffero4(nio))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufferi4(jk1,jk2) = zin8(jk1,jk2)
      end do
    end do

    ierr = ezsint(buffero4,bufferi4)

    do jk1 = 1,nio
      zout8(jk1) = buffero4(jk1)
    end do

    deallocate(bufferi4)
    deallocate(buffero4)

  end function utl_ezsint_r8_2dTo1d


  function utl_ezuvint_r4_2d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: uuout(:,:), vvout(:,:)
    real(4) :: uuin(:,:) , vvin(:,:)
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: ezuvint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    ierr = ezuvint(uuout, vvout, uuin, vvin)

  end function utl_ezuvint_r4_2d


  function utl_ezuvint_r4_2dTo1d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: uuout(:), vvout(:)
    real(4) :: uuin(:,:) , vvin(:,:)
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: ezuvint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    ierr = ezuvint(uuout, vvout, uuin, vvin)

  end function utl_ezuvint_r4_2dTo1d


  function utl_ezuvint_r8_1d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: uuout(:), vvout(:)
    real(8) :: uuin(:) , vvin(:)
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: nio, nii
    integer :: jk1
    real, allocatable :: bufuuout4(:), bufvvout4(:)
    real, allocatable :: bufuuin4(:), bufvvin4(:)
    integer :: ezuvint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(uuin)
    nio = size(uuout)

    allocate(bufuuout4(nio))
    allocate(bufvvout4(nio))
    allocate(bufuuin4(nii))
    allocate(bufvvin4(nii))

    do jk1 = 1,nii
      bufuuin4(jk1) = uuin(jk1)
      bufvvin4(jk1) = vvin(jk1)
    end do

    ierr = ezuvint(bufuuout4, bufvvout4, bufuuin4, bufvvin4)

    do jk1 = 1,nio
      uuout(jk1) = bufuuout4(jk1)
      vvout(jk1) = bufvvout4(jk1)
    end do

    deallocate(bufuuin4)
    deallocate(bufvvin4)
    deallocate(bufuuout4)
    deallocate(bufvvout4)

  end function utl_ezuvint_r8_1d


  function utl_ezuvint_r8_2d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: uuout(:,:), vvout(:,:)
    real(8) :: uuin(:,:) , vvin(:,:)
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: nio, njo, nii, nji
    integer :: jk1, jk2
    real, allocatable :: bufuuout4(:,:), bufvvout4(:,:)
    real, allocatable :: bufuuin4(:,:), bufvvin4(:,:)
    integer :: ezuvint

    call utl_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(uuin,1)
    nji = size(uuin,2)

    nio = size(uuout,1)
    njo = size(uuout,2)

    allocate(bufuuout4(nio,njo))
    allocate(bufvvout4(nio,njo))
    allocate(bufuuin4(nii,nji))
    allocate(bufvvin4(nii,nji))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufuuin4(jk1,jk2) = uuin(jk1,jk2)
        bufvvin4(jk1,jk2) = vvin(jk1,jk2)
      end do
    end do

    ierr = ezuvint(bufuuout4, bufvvout4, bufuuin4, bufvvin4)

    do jk2 = 1,njo
      do jk1 = 1,nio
        uuout(jk1,jk2) = bufuuout4(jk1,jk2)
        vvout(jk1,jk2) = bufvvout4(jk1,jk2)
      end do
    end do

    deallocate(bufuuin4)
    deallocate(bufvvin4)
    deallocate(bufuuout4)
    deallocate(bufvvout4)

  end function utl_ezuvint_r8_2d


  function utl_ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, ax, ay) result(vezgdef)
    implicit none

    integer :: vezgdef

    integer :: ni, nj, ig1, ig2, ig3, ig4
    real(8) :: ax(*), ay(*)
    character(len=*) :: grtyp, grtypref

    integer :: ier2,jk,ilenx,ileny
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
    end if

    allocate(bufax4(ilenx))
    allocate(bufay4(ileny))

    do jk = 1,ilenx
       bufax4(jk) = ax(jk)
    end do
    do jk = 1,ileny
       bufay4(jk) = ay(jk)
    end do

    ier2 = ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, &
         bufax4, bufay4)

    deallocate(bufax4)
    deallocate(bufay4)

    vezgdef=ier2

  end function utl_ezgdef


  subroutine utl_cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat0, xlon0, dlat, dlon) 
    implicit none

    integer :: ig1, ig2, ig3, ig4   
    real(8) :: xlat0, xlon0, dlat, dlon 
    character(len=*) :: grtyp 

    real(4) :: xlat04, xlon04, dlat4, dlon4

    xlat04=xlat0
    xlon04=xlon0
    dlat4=dlat
    dlon4=dlon

    call cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat04, xlon04, dlat4, dlon4)

  end subroutine utl_cxgaig


  function utl_fstlir(fld8, iun, ni, nj, nk, datev, etiket, &
       ip1, ip2, ip3, typvar, nomvar) result(vfstlir)
    implicit none

    integer :: vfstlir
    real(8) :: fld8(*)
    integer :: iun, ni, nj, nk, datev, ip1, ip2, ip3
    character(len=*) :: etiket
    character(len=*) :: nomvar
    character(len=*) :: typvar

    integer :: key1,key2, ilen, jk1, jk2, jk3, la
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
                end do
             end do
          end do
       end if

       deallocate(buffer4)
    end if

    vfstlir=key1

  end function utl_fstlir

  function utl_fstlir_r4(fld_r4, iun, ni, nj, nk, datev, etiket, &
       ip1, ip2, ip3, typvar, nomvar) result(vfstlir)
    implicit none

    integer :: vfstlir
    real(4) :: fld_r4(*)
    integer :: iun, ni, nj, nk, datev, ip1, ip2, ip3
    character(len=*) :: etiket
    character(len=*) :: nomvar
    character(len=*) :: typvar

    integer :: key1,key2, ilen, jk1, jk2, jk3, la
    real(4), allocatable :: buffer_r4(:)

    integer :: fstluk, fstinf

    !     Get field dimensions.
    key1 = fstinf(iun, ni, nj, nk, datev, etiket, &
         ip1, ip2, ip3, typvar, nomvar)

    if(key1 >= 0) then
       ilen = ni*nj*nk
       allocate(buffer_r4(ilen))
       !     Read field
       key2 = fstluk(buffer_r4, key1, ni, nj, nk)
       if(key2 >= 0) then
          do jk3 = 1,nk
             do jk2 = 1,nj
                do jk1 = 1,ni
                   la=jk1+(jk2-1)*ni+(jk3-1)*ni*nj
                   fld_r4(la) = buffer_r4(la)
                end do
             end do
          end do
       end if

       deallocate(buffer_r4)
    end if

    vfstlir=key1

  end function utl_fstlir_r4

  function utl_fstecr(fld8, npak, iun, dateo, deet, &
       npas, ni, nj, nk, ip1, ip2, ip3, typvar, &
       nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, & 
       datyp, rewrit) result(vfstecr)
    implicit none

    integer :: vfstecr
    real(4) :: work
    integer, intent(in) :: ni,nj,nk
    real(8) :: fld8(ni,nj,nk)
    integer :: iun, ip1, ip2, ip3, ig1, ig2, ig3, ig4
    integer :: npak, dateo, deet, npas, datyp
    logical :: rewrit  
    character(len=*) :: etiket 
    character(len=*) :: typvar
    character(len=*) :: grtyp 
    character(len=*) :: nomvar            

    integer :: ikey, jk1, jk2, jk3
    real(4), allocatable :: buffer4(:,:,:)

    integer :: fstecr

    allocate(buffer4(ni,nj,nk))

    do jk3 = 1,nk
      do jk2 = 1,nj
        do jk1 = 1,ni
          buffer4(jk1,jk2,jk3) = fld8(jk1,jk2,jk3)
        end do
      end do
    end do

    ikey = fstecr(buffer4, work, npak, iun, dateo, deet, &
         npas, ni, nj, nk, ip1, ip2, ip3, typvar, nomvar, & 
         etiket, grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

    deallocate(buffer4)

    vfstecr=ikey

  end function utl_fstecr


  function utl_findArrayIndex( klist, klen, kentry ) result(isrcheq)
    !
    ! :Purpose: Find entry in list.
    !
    ! :Arguments:
    !           :klist: Input list.
    !           :klen: Dimension of input list.
    !           :kentry: Entry.
    !           :isrcheq: Index of entry: (0, not found, >0, found)
    !
    implicit none

    ! Arguments:
    INTEGER :: ISRCHEQ
    INTEGER :: KENTRY
    integer :: KLEN
    INTEGER :: KLIST(KLEN)

    ! locals:
    integer :: JI
    

    ISRCHEQ = 0
    DO JI=1,KLEN
       IF ( KLIST(JI) .EQ. KENTRY ) THEN
          ISRCHEQ = JI
          RETURN
       END IF
    END DO

  end function utl_findArrayIndex


  subroutine utl_matsqrt(matrix, rank, exponentSign, printInformation_opt )
    ! 
    ! :Purpose: Calculate square root of an error covariance matrix
    !
    implicit none

    integer, intent(in)           :: rank
    real(8), intent(inout)        :: matrix(rank,rank)
    real(8), intent(in)           :: exponentSign
    logical, intent(in), optional :: printInformation_opt ! switch to print be more verbose

    real(8), allocatable :: eigenValues(:)
    real(8), allocatable :: work(:)
    real(8), allocatable :: eigenVectors(:,:)
    integer :: sizework, info, index, index1, index2 
    logical :: printInformation

    if (present(printInformation_opt)) then
       printInformation = printInformation_opt
    else
       printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)
      write(*,*) 'utl_matsqrt: Starting...'
    end if

    sizework = 64 * rank
    allocate(work(sizework))

    allocate(eigenValues (rank))
    allocate(eigenVectors(rank,rank))

    !- Calculate EigenVectors (V) and EigenValues (D) of B matrix
    eigenVectors(:,:) = matrix(:,:)

    call dsyev('V','U',rank,   & ! IN
               eigenVectors,   & ! INOUT
               rank,           & ! IN
               eigenValues,    & ! OUT
               work, sizework, & ! IN
               info )            ! OUT

    if ( info /= 0 ) then
      write(*,*)
      write(*,*) 'dsyev: ',info
      call utl_abort('utl_matsqrt: DSYEV failed!')
    end if

    if (printInformation) then
      write(*,*)
      write(*,'(1x,"Original EIGEN VALUES: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index),index=1,rank)
      if (exponentSign < 0.d0) then
        write(*,*)
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValues(:))/minval(eigenValues(:))
      end if
    end if

    !- Calculate Matrix^0.5 = V D^0.5 V^t
    where(eigenValues(:) < 0.d0)
      eigenValues = 0.d0
    end where

    eigenValues(:) = eigenValues(:)**(0.5d0*exponentSign)

    do index1 = 1, rank
      do index2 = 1, rank
        matrix(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                    * eigenVectors (index2,1:rank)   &
                                    * eigenValues(1:rank) )
      end do
    end do

    deallocate(eigenVectors)
    deallocate(eigenValues)
    deallocate(work)
    
    if (printInformation) then
      write(*,*)
      write(*,*) 'utl_matsqrt: Ending...'
    end if

  end subroutine utl_matsqrt

  !--------------------------------------------------------------------------
  ! utl_matInverse
  !--------------------------------------------------------------------------
  subroutine utl_matInverse(matrix, rank, inverseSqrt_opt, printInformation_opt)
    !
    ! :Purpose: Calculate the inverse of a covariance matrix 
    !           and, optionally, also the inverse square-root.
    !
    implicit none

    ! Arguments
    integer, intent(in)              :: rank                 ! order of the matrix
    real(8), intent(inout)           :: matrix(:,:)          ! on entry, the original matrix; on exit, the inverse
    real(8), intent(inout), optional :: inverseSqrt_opt(:,:) ! if present, the inverse sqrt matrix on exit 
    logical, intent(in), optional    :: printInformation_opt ! switch to print be more verbose

    ! Local variables
    integer :: index1, index2, info, sizework
    real(8) :: sizework_r8
    real(8), allocatable :: work(:), eigenVectors(:,:), eigenValues(:)
    logical :: printInformation

    if (present(printInformation_opt)) then
      printInformation = printInformation_opt
    else
      printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)' utl_matInvers: Inverse matrix of a symmetric matrix'
    end if

    !     1. Computation of eigenvalues and eigenvectors

    allocate(eigenVectors(rank,rank))
    allocate(eigenValues(rank))

    do index2=1,rank
      do index1=1,rank
        eigenVectors(index1,index2)=matrix(index1,index2)
      end do
    end do

    ! query the size of the 'work' vector by calling 'DSYEV' with 'sizework=-1'
    sizework = -1
    info = -1
    call dsyev('V','U',rank, eigenVectors, rank, eigenValues, sizework_r8, sizework, info)

    ! compute the eigenvalues
    sizework=int(sizework_r8)
    allocate(work(sizework))
    call dsyev('V','U',rank, eigenVectors,rank, eigenValues,work, sizework, info)
    deallocate(work)

    if (printInformation) then
      write(*,'(1x,"Original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index1),index1=1,rank)

      if(minval(eigenValues) > 1.0d-10) then
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValues)/minval(eigenValues)
      end if
    end if

    !     2.  Take inverse of eigenvalues

    do index1=1,rank
      if(eigenValues(index1) > 1.0d-10) then
        eigenValues(index1)= 1.0d0/eigenValues(index1)
      else
        write(*,*) 'utl_matInverse: WARNING eigenvalue is too small = ', index1, eigenValues(index1)
        eigenValues(index1) = 0.0d0
      end if
    end do

    if (printInformation) then
      write(*,'(1x,"Inverse of original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index1),index1=1,rank)
    end if

    !     3.  Compute the inverse matrix

    do index2 = 1, rank
      do index1 = 1, rank
        matrix(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                    * eigenVectors (index2,1:rank)   &
                                    * eigenValues(1:rank) )
      end do
    end do

    !     4.  If requested, computed the inverse square-root also

    if (present(inverseSqrt_opt)) then
      do index1=1,rank
        if(eigenValues(index1) > 1.0d-10) then
          eigenValues(index1)= sqrt(eigenValues(index1))
        else
          eigenValues(index1) = 0.0d0
        end if
      end do
      do index2 = 1, rank
        do index1 = 1, rank
          inverseSqrt_opt(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                               * eigenVectors (index2,1:rank)   &
                                               * eigenValues(1:rank) )
        end do
      end do
    end if

    !     5. Deallocate local arrays
    deallocate(eigenVectors,eigenValues)

    if (printInformation) then
      write(*,*) 'utl_matInverse: done'
      write(*,*) ' '
    end if

  end subroutine utl_matInverse

  !--------------------------------------------------------------------------
  ! utl_eigenDecomp
  !--------------------------------------------------------------------------
  subroutine utl_eigenDecomp(matrix, eigenValues, eigenVectors, tolerance, numReturned, printInformation_opt)
    !
    ! :Purpose: Calculate eigenValues/Vectors and return only those with eigenValues
    !           whose magnitude is greater than the specified tolerance.
    !
    implicit none

    ! Arguments
    real(8), intent(inout)        :: matrix(:,:)          ! on entry, the original matrix; on exit, the inverse
    real(8), intent(out)          :: eigenValues(:)       ! computed eigenValues
    real(8), intent(out)          :: eigenVectors(:,:)    ! computed eigenVectors
    real(8), intent(in)           :: tolerance            ! threshold for eigenValue magnitude to be returned
    integer, intent(out)          :: numReturned          ! number of eigenValues/Vectors returned
    logical, intent(in), optional :: printInformation_opt ! switch to print be more verbose

    ! Local variables
    integer :: rank, index1, index2, info, sizework
    real(8) :: sizework_r8
    real(8), allocatable :: work(:), eigenVectorsOrig(:,:), eigenValuesOrig(:)
    logical :: printInformation

    if (present(printInformation_opt)) then
      printInformation = printInformation_opt
    else
      printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)' utl_eigenDecomp: Eigen decomposition of a symmetric matrix'
    end if

    !     1. Computation of eigenvalues and eigenvectors

    rank = size(matrix,1)
    allocate(eigenVectorsOrig(rank,rank))
    allocate(eigenValuesOrig(rank))

    do index2 = 1, rank
      do index1 = 1, rank
        eigenVectorsOrig(index1,index2)=matrix(index1,index2)
      end do
    end do

    ! Query the size of the 'work' vector by calling 'DSYEV' with 'sizework=-1'
    sizework = -1
    info = -1
    call dsyev('V', 'U', rank, eigenVectorsOrig, rank, eigenValuesOrig,  &
               sizework_r8, sizework, info)

    ! Compute the eigenvalues/vectors
    sizework = int(sizework_r8)
    allocate(work(sizework))
    call dsyev('V', 'U', rank, eigenVectorsOrig, rank, eigenValuesOrig,  &
               work, sizework, info)
    deallocate(work)

    if (printInformation) then
      write(*,'(1x,"Original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValuesOrig(index1),index1=1,rank)

      if(minval(eigenValuesOrig) > tolerance) then
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValuesOrig)/minval(eigenValuesOrig)
      end if
    end if

    !     2.  Determine which eigen values/vectors to return

    numReturned = 0
    do index1 = rank, 1, -1
      if (eigenValuesOrig(index1) > tolerance) then
        numReturned = numReturned + 1
      else
        exit
      end if
    end do

    if (printInformation) then
      write(*,*) 'Number of eigen values returned =', numReturned, ' out of', rank
    end if

    !     3.  Copy eigenValues/Vectors into output arrays with reversed order

    do index1 = 1, numReturned
      ! And set negative values to zero
      eigenValues(index1) = max(0.0D0,eigenValuesOrig(rank-index1+1))
    end do
    do index1 = numReturned+1, rank
      eigenValues(index1) = 0.0D0
    end do

    do index2 = 1, numReturned
      do index1 = 1, rank
        eigenVectors(index1,index2) = eigenVectorsOrig(index1,rank-index2+1)
      end do
    end do
    do index2 = numReturned+1, rank
      do index1 = 1, rank
        eigenVectors(index1,index2) = 0.0D0
      end do
    end do

    !     4. Deallocate local arrays
    deallocate(eigenVectorsOrig,eigenValuesOrig)

    if (printInformation) then
      write(*,*) 'utl_eigenDecomp: done'
      write(*,*) ' '
    end if

  end subroutine utl_eigenDecomp

  !-----------------------------------------
  ! utl_pseudo_inverse
  !-----------------------------------------
  subroutine utl_pseudo_inverse(inputMatrix, pseudoInverse, threshold_opt)
    !
    ! :Purpose: to calculate the More-Penrose pseudo inverse of the matrix inputMatrix
    !
    implicit none
    !Arguments:
    real(8), intent(in)           :: inputMatrix(:,:)   ! Input Matrix
    real(8), intent(out)          :: pseudoInverse(:,:) ! its Moore Penrose Pseudo-Inverse
    real(8), optional, intent(in) :: threshold_opt 
    !Locals:
    real(8), allocatable :: copyMatrix(:,:), leftSingularVector(:,:), rightSingularVectorT(:,:)
    real(8), allocatable :: singularValues(:)
    integer :: info, lwork, lineIndex
    integer :: lineDim, columnDim, minDim
    real(8) :: thresh
    real(8), allocatable :: work(:)
    character(len=80) :: errorMessage

    lineDim = size(inputMatrix, dim=1)
    columnDim = size(inputMatrix, dim=2)
    minDim = min(lineDim, columnDim)

    allocate( copyMatrix(lineDim,columnDim), leftSingularVector(lineDim,lineDim) )
    allocate( rightSingularVectorT(columnDim,columnDim), singularValues(minDim) )

    copyMatrix(:,:) = inputMatrix(:,:) ! Work with a copy because copyMatrix will be overwriten
    lwork = max(10000, max(1, 3 * min(lineDim,columnDim) + max(lineDim,columnDim), 5 * minDim ))
    allocate(work(lwork))
    call dgesvd("A", "A", lineDim, columnDim, copyMatrix, lineDim, singularValues, &
         leftSingularVector, lineDim, rightSingularVectorT, columnDim, work, lwork, info ) 
   
    if (info /= 0) then
      write(errorMessage,*) "utl_pseudo_inverse: Problem in DGESVD ! ",info
      call utl_abort(errorMessage)
    end if

    deallocate(work)

    if (present(threshold_opt)) then
      thresh = threshold_opt
    else
      !according to wikipedia... as in matlab or numpy
      thresh = epsilon(thresh) * max(lineDim, columnDim) * maxval(singularValues)
    end if
    print *,"utl_pseudo_inverse: threshold= ",thresh

    pseudoInverse(:,:)=0.d0
    do lineIndex = 1, minDim
      If (singularValues(lineIndex) > thresh) then
        pseudoInverse(lineIndex,:) = ( 1.d0 / singularValues(lineIndex) ) * leftSingularVector(:,lineIndex)
      end if
    end do

    pseudoInverse = matmul( transpose(rightSingularVectorT), pseudoInverse)

    deallocate( singularValues, rightSingularVectorT )
    deallocate( leftSingularVector, copyMatrix )

  end subroutine utl_pseudo_inverse


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


  subroutine utl_getfldprm(kip1s,kip2,kip3,knlev,cdetiket,cdtypvar,kgid, &
                           cdvar,kstampv,knmaxlev,kinmpg,kip1style,kip1kind, &
                           ktrials,koutmpg)
    !
    ! :Purpose:  Get 3D grid parameters for a specific trial field
    !            and check for consitancies between grid parameters
    !            of the levels.
    !
    ! :Arguments:
    !
    !  :Input:
    !     :cdvar: variable name to get the vertical levels from
    !     :kstampv: valid date time stamp of the variable
    !     :knmaxlev: maximum number of levels
    !     :kinmpg: file unit of trial field
    !     :ktrials:  number of trial files.  
    !
    !  :Output:
    !     :kip1s: list of ip1s of variable cdvar
    !     :kip2: ip2 for variable cdvar
    !     :kip3: ip3 for variable cdvar
    !     :knlev: number of levels of variable cdvar
    !     :cdetiket: etiket of field cdvar
    !     :cdtypvar: typvar of field cdvar
    !     :kgid: handle of the field descriptor
    !     :kip1style: style in which ip1 is encoded (15 or 31 bits)
    !     :kip1kind: kind of vertical coord encoded in ip1
    !     :koutmpg: the unit which contains the selected records.  
    !
    implicit none

    integer :: kstampv,knmaxlev,knlev,kgid
    integer :: kip1s(knmaxlev),kip1style,kip1kind,kip2,kip3
    integer :: ktrials, koutmpg  
    integer :: kinmpg(ktrials)
    character(len=*) :: cdtypvar
    character(len=*) :: cdvar
    character(len=*) :: cdetiket

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
       end if
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
       end if
    end do ! End of loop k   
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
          end if
       end do
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
       end if
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
       end do
       write(*,*) 'Error - getfldprm2: no record found at time ' &
            ,idatefull,' for field ',cdvar,' but',knlev, &
            ' records found in unit ',kinmpg(k)
       call utl_abort('GETFLDPRM2')
    end if
    !
  end subroutine utl_getfldprm


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


  subroutine utl_open_asciifile(filename,unit)
    ! 
    ! :Purpose: Opens an ascii file for output 
    !
    ! :Arguments:
    ! :Input:
    !       :filename: filename
    !       :unit: unit number to use or 0 to let fnom set value
    ! :Output:
    !       :unit: unit number associated with file
    !
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(out) :: unit
    logical :: file_exists
    integer :: ier
    character(len=20) :: mode
    
    inquire(file=trim(filename), exist=file_exists)
    
    if (file_exists) then
       mode = 'FTN+APPEND+R/W'
    else
       mode = 'FTN+R/W'
    end if

    unit=0
    
    ier = utl_open_file(unit,trim(filename),trim(mode))

    if (ier.ne.0) call utl_abort('utl_open_messagefile: Error associating unit number')

  end subroutine utl_open_asciifile


  function utl_open_file(unit,filename,mode) result(ier)
    ! 
    ! :Purpose: This is a temporary subroutine to open a file with fnom that is needed due to
    !           a bug in fnom that does not allow an ascii file to be opened in 'APPEND' mode.  
    !
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
    

  function utl_stnid_equal(id1,id2) result(same)
    !
    ! :Purpose: Compares STNID values allowing for * as wildcards and trailing blanks 
    !
    ! :Arguments:
    !           :id1: reference stnid
    !           :id2: stnid being verified
    !           :same: logical indicating if id1 and id2 match
    !     
    implicit none

    logical :: same
    character(len=*), intent(in) :: id1
    character(len=*), intent(in) :: id2

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
 

  character(len=20) function utl_int2str(i)
    !
    ! :Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort. 
    !
    implicit none

    integer, intent(in) :: i
    
    write(utl_int2str,*) i
    utl_int2str = adjustl(utl_int2str)
    
  end function utl_int2str
            

  character(len=20) function utl_float2str(x)
    !
    ! :Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort.
    !
    implicit none

    real(8), intent(in) :: x

    write(utl_float2str,*) x
    utl_float2str = adjustl(utl_float2str)

  end function utl_float2str


  subroutine utl_resize_1d_real(arr,dim1)
    !
    ! :Purpose: Resize 1D array
    !
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


  subroutine utl_resize_1d_int(arr,dim1)
    !
    ! :Purpose: Resize 1D array 
    !
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

 
  subroutine utl_resize_1d_str(arr,dim1)
    !
    ! :Purpose: Resize 1D array
    !
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


  subroutine utl_resize_2d_real(arr,dim1,dim2)
    !
    ! :Purpose: Resize 2D array
    !
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


  subroutine utl_resize_3d_real(arr,dim1,dim2,dim3)
    !
    ! :Purpose: Resize 3D array
    !
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


  subroutine utl_get_stringId(cstringin,nobslev,CList,NListSize,Nmax,elemId)
    ! 
    ! :Purpose: Get element ID from a list of accumulating character strings (e.g. stnids). 
    !           Called by filt_topoChm in filterobs_mod.ftn90
    !
    ! :Arguments:
    !           :Nmax: Max allowed dimension.
    !           :NListSize: Input number of identified IDs (must be >=0 and <=Nmax)
    !           :CList: Input list of accumulated character strings for uni and multi-level data.
    !           :cstringin: Input character string
    !           :nobslev: Number of elements in profile associated to cstringin.
    !           :NListSize: Updated number of identified IDs
    !           :CList: Updated list of accumulated character strings
    !           :elemId: Index of cstringin within CList_chm
    !        
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


  subroutine utl_get_Id(id,IdList,NListSize,Nmax,elemId)
    ! 
    ! :Purpose: Get element ID from list of accumulating integer IDs.
    !
    ! :Arguments:
    !      :Input:
    !            :Nmax: Max allowed dimension.
    !            :NListSize: Input number of IDs (must be >=0 and <=Nmax)
    !            :IdList: Input list of accumulated IDs.
    !            :id: Input id for individual obs 
    !      :Output:
    !            :NListSize: Updated number of IDs 
    !            :IdList: Updated list of accumulated IDs.
    !            :elemId: Index of id within List
    !     
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
  

  subroutine utl_readFstField( fname, varName, iip1, iip2, iip3, etiketi, &
                               ni, nj, nkeys, array, xlat_opt, xlong_opt, lvls_opt, kind_opt )
    !
    ! :Purpose:  Read specified field from standard RPN/fst file. Could be one
    !           to all levels depending on the input iip1,iip2,iip3 values.
    !
    !           Currently assumes lat/long (or Gaussian) type grids.
    !           See hco_SetupFromFile for example toward future generalizations.
    !           Generalization would require having xlat and xlong being 2D.
    !
    ! :Arguments:
    !           :Input:
    !                 :fname: input filename
    !                 :varName:  search nomvar
    !                 :iip1: search ip1
    !                 :iip2: search ip2
    !                 :iip3: search ip3
    !                 :etiketi: search etiket
    !           :Output:
    !                 :ni: ni values
    !                 :nj: nj values
    !                 :nkeys: number of records satisfying search criteria
    !                 :array: data arrray
    !                 :xlat_opt: 1D latitude array (optional)
    !                 :xlong_opt: 1D longitude array (optional)
    !                 :lvls_opt: 1D vertical coordinate array (optional)
    !                 :kind_opt: vertical coordinate type according to convip (optional)
    !
    implicit none

    integer, intent(in) :: iip1,iip2,iip3
    character(len=*), intent(in) :: varName,fname,etiketi
    integer, intent(out) :: ni, nj, nkeys
    integer, intent(out), optional :: kind_opt
    real(8), intent(out), allocatable :: array(:,:,:)
    real(8), intent(out), allocatable, optional :: lvls_opt(:), xlat_opt(:),xlong_opt(:)

    integer, external :: fnom,fclos,fstouv,fstfrm,fstinl,fstlir,fstluk,fstprm

    real(4) :: lvl_r4

    logical :: Exists
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


  function utl_varNamePresentInFile(varName, fileName_opt, fileUnit_opt, typvar_opt) result(found)
    implicit none

    ! arguments:
    character(len=*), intent(in) :: varName
    character(len=*), optional, intent(in) :: fileName_opt
    integer, optional, intent(in) :: fileUnit_opt
    character(len=*), optional, intent(in) :: typvar_opt
    logical :: found

    ! locals:
    integer :: fnom, fstouv, fstfrm, fclos, fstinf
    integer :: ni, nj, nk, key, ierr
    integer :: unit
    character(len=128) :: fileName
    character(len=2)   :: typvar
    logical :: openFile

    if ( present(fileUnit_opt) ) then
      unit = fileUnit_opt
      openFile = .false.
    else
      unit = 0
      openFile = .true.
      if ( present(fileName_opt) ) then
        fileName = fileName_opt
      else
        call utl_abort('utl_varNamePresentInFile: please provide and file name or unit')
      end if
    end if

    if ( present(typvar_opt) ) then
      typvar = trim(typvar_opt)
    else
      typvar = ' '
    end if

    if (openFile) then
      ierr = fnom(unit,fileName,'RND+OLD+R/O',0)
      ierr = fstouv(unit,'RND+OLD')
    end if

    key = fstinf(unit, ni, nj, nk, -1 ,' ', -1, -1, -1, typvar, trim(varName))
    
    if ( key > 0 )  then
      found = .true.
    else
      found = .false.
    end if

    if (openFile) then
      ierr =  fstfrm(unit)
      ierr =  fclos (unit)
    end if

  end function utl_varNamePresentInFile


  subroutine utl_reAllocate_char_1d(array,dim1)
    implicit none
    character(len=128), allocatable :: array(:)
    integer :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))

  end subroutine utl_reAllocate_char_1d

  subroutine utl_reAllocate_char_2d(array,dim1,dim2)
    implicit none
    character(len=128), allocatable :: array(:,:)
    integer :: dim1, dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))

  end subroutine utl_reAllocate_char_2d

  subroutine utl_reAllocate_char_3d(array,dim1,dim2,dim3)
    implicit none
    character(len=128), allocatable :: array(:,:,:)
    integer :: dim1, dim2, dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))

  end subroutine utl_reAllocate_char_3d

  subroutine utl_reAllocate_log_1d(array,dim1)
    implicit none
    logical, allocatable :: array(:)
    integer :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = .true.

  end subroutine utl_reAllocate_log_1d

  subroutine utl_reAllocate_log_2d(array,dim1,dim2)
    implicit none
    logical, allocatable :: array(:,:)
    integer :: dim1, dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = .true.

  end subroutine utl_reAllocate_log_2d

  subroutine utl_reAllocate_log_3d(array,dim1,dim2,dim3)
    implicit none
    logical, allocatable :: array(:,:,:)
    integer :: dim1, dim2, dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = .true.

  end subroutine utl_reAllocate_log_3d

  subroutine utl_reAllocate_int_1d(array,dim1)
    implicit none
    integer, allocatable :: array(:)
    integer :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0d0

  end subroutine utl_reAllocate_int_1d

  subroutine utl_reAllocate_int_2d(array,dim1,dim2)
    implicit none
    integer, allocatable :: array(:,:)
    integer :: dim1, dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0d0

  end subroutine utl_reAllocate_int_2d

  subroutine utl_reAllocate_int_3d(array,dim1,dim2,dim3)
    implicit none
    integer, allocatable :: array(:,:,:)
    integer :: dim1, dim2, dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0d0

  end subroutine utl_reAllocate_int_3d

  subroutine utl_reAllocate_r4_1d(array,dim1)
    implicit none
    real(4), allocatable :: array(:)
    integer :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0.0d0

  end subroutine utl_reAllocate_r4_1d

  subroutine utl_reAllocate_r8_1d(array,dim1)
    implicit none
    real(8), allocatable :: array(:)
    integer :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0.0d0

  end subroutine utl_reAllocate_r8_1d

  subroutine utl_reAllocate_r4_2d(array,dim1,dim2)
    implicit none
    real(4), allocatable :: array(:,:)
    integer :: dim1, dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_2d

  subroutine utl_reAllocate_r8_2d(array,dim1,dim2)
    implicit none
    real(8), allocatable :: array(:,:)
    integer :: dim1, dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_2d

  subroutine utl_reAllocate_r4_3d(array,dim1,dim2,dim3)
    implicit none
    real(4), allocatable :: array(:,:,:)
    integer :: dim1, dim2, dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_3d


  subroutine utl_reAllocate_r8_3d(array,dim1,dim2,dim3)
    implicit none
    real(8), allocatable :: array(:,:,:)
    integer :: dim1, dim2, dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_3d


  subroutine utl_reAllocate_r4_4d(array,dim1,dim2,dim3,dim4)
    implicit none
    real(4), allocatable :: array(:,:,:,:)
    integer :: dim1, dim2, dim3, dim4

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4))
    array(:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_4d


  subroutine utl_reAllocate_r8_4d(array,dim1,dim2,dim3,dim4)
    implicit none
    real(8), allocatable :: array(:,:,:,:)
    integer :: dim1, dim2, dim3, dim4

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4))
    array(:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_4d


  subroutine utl_reAllocate_r4_5d(array,dim1,dim2,dim3,dim4,dim5)
    implicit none
    real(4), allocatable :: array(:,:,:,:,:)
    integer :: dim1, dim2, dim3, dim4, dim5

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4*dim5 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4,dim5))
    array(:,:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_5d


  subroutine utl_reAllocate_r8_5d(array,dim1,dim2,dim3,dim4,dim5)
    implicit none
    real(8), allocatable :: array(:,:,:,:,:)
    integer :: dim1, dim2, dim3, dim4, dim5

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4*dim5 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4,dim5))
    array(:,:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_5d


  subroutine utl_heapsort2d(array)
    !
    ! :Purpose: Sort a real 2D array in ascending order according
    !           to the first column
    ! 
    implicit none
    real(4), intent(inout) :: array(:,:)

    real(4) :: values(2) ! temporary value
    integer :: i,j,nsize
    integer :: ileft,iright

    nsize = size(array,1)
    ileft=nsize/2+1
    iright=nsize

    if (nsize == 1) return                  

    do 
      if(ileft > 1)then
        ileft=ileft-1
        values(:) = array(ileft,:)
      else
        values(:) = array(iright,:)
        array(iright,:) = array(1,:)
        iright = iright-1
        if (iright == 1) then
          array(1,:) = values(:)
          return
        end if
      end if
      i = ileft
      j = 2*ileft
      do while (j <= iright) 
        if (j < iright) then
          if (array(j,1) < array(j+1,1)) j=j+1
        endif
        if (values(1) < array(j,1)) then
          array(i,:) = array(j,:)
          i = j
          j = j+j
        else
          j = iright+1
        end if
      end do
      array(i,:) = values(:)
    end do

  end subroutine utl_heapsort2d


  subroutine utl_splitString(string,separator,stringArray)
    implicit none
    character(len=*), intent(in) :: string
    character(len=*), intent(in) :: separator
    character(len=256), allocatable :: stringArray(:)

    integer :: stringArraySize

    stringArraySize = count(transfer(string, 'a', len(string)) == separator) + 1

    allocate(stringArray(stringArraySize))

    read(string, *) stringArray(1:stringArraySize)

    write(*,*)  'utl_splitString: stringArraySize = ', stringArraySize
    write(*,*)  'utl_splitString: stringArray     = ', stringArray(:)
    
  end subroutine utl_splitString


  subroutine utl_stringArrayToIntegerArray(stringArray,integerArray)
    implicit none
    character(len=256) :: stringArray(:)
    integer, allocatable :: integerArray(:)

    integer :: arraySize, arrayIndex

    arraySize = size(stringArray)

    allocate(integerArray(arraySize))

    do arrayIndex = 1, arraySize
      read(stringArray(arrayIndex),'(i5)')  integerArray(arrayIndex)
    end do

    write(*,*)  'utl_stringArrayToIntegerArray: integerArray = ', integerArray(:)

  end subroutine utl_stringArrayToIntegerArray

  !--------------------------------------------------------------------------
  ! utl_isNamelistPresent
  !--------------------------------------------------------------------------
  function utl_isNamelistPresent(namelistSectionName, namelistFileName) result(found)
    !
    ! :Purpose: To find if a namelist name tag is present in a namelist file
    ! 
    implicit none
    logical :: found
    character(len=*), intent(in) :: namelistSectionName
    character(len=*), intent(in) :: namelistFileName

    integer :: unit, fnom, fclos, ierr
    character (len=1000) :: text
    character (len=100)  :: word, namelistSectionNameUpper
    logical :: namelistExist

    ! Check if namelistFileName is present
    inquire(file=namelistFileName,exist=namelistExist)
    if (.not. namelistExist) then
      call utl_abort('utl_isNamelistPresent: namelist file is missing : '// namelistFileName)
    end if

    ! Open the namelist file
    unit=0
    ierr=fnom(unit,namelistFileName,'FTN+SEQ+R/O',0)

    ! Search for namelistSectionName
    found = .false.
    namelistSectionNameUpper = namelistSectionName
    ierr = clib_toUpper(namelistSectionNameUpper)
    namelistLoop : do
      read (unit,"(a)",iostat=ierr) text ! read line into character variable
      if (ierr /= 0) exit
      if (trim(text) == "") cycle ! skip empty lines
      read (text,*) word ! read first word of line
      ierr = clib_toUpper(word)
      if (trim(word) == '&'//trim(namelistSectionNameUpper)) then ! case insensitive 
        ! found search string at beginning of line
        found = .true.
        exit
      end if
    end do namelistLoop

    ! Close the namelist file
    ierr=fclos(unit)

  end function utl_isNamelistPresent

  !-----------------------------------------------------------------
  ! utl_parseColumns
  !-----------------------------------------------------------------
  subroutine utl_parseColumns(line, numColumns, stringArray_opt)
    !
    ! :Purpose: To return column values in array of strings and
    !           the number of space-delimited columns in a string
    ! 
    implicit none
    ! Arguments
    character(len=*), intent(in) :: line
    integer, intent(out) :: numColumns
    character(len=*), intent(out), optional :: stringArray_opt(:)
    ! Locals
    integer :: linePosition, wordPosition, lineLength

    linePosition = 1
    lineLength = len_trim(line)
    numColumns = 0
    
    do while(linePosition <= lineLength)

      do while(line(linePosition:linePosition) == ' ') 
        linePosition = linePosition + 1
        if (lineLength < linePosition) return
      end do

      numColumns = numColumns + 1
      wordPosition = 0
      if (present(stringArray_opt)) then
        stringArray_opt(numColumns) = ''
      end if

      do
        if (linePosition > lineLength) return
        if (line(linePosition:linePosition) == ' ') exit
        if (present(stringArray_opt)) then
          wordPosition = wordPosition + 1
          stringArray_opt(numColumns)(wordPosition:wordPosition) = line(linePosition:linePosition)
        end if
        linePosition = linePosition + 1
      end do

    end do
    
  end subroutine utl_parseColumns

  !--------------------------------------------------------------------------
  ! utl_copyFile
  !--------------------------------------------------------------------------
  function utl_copyFile(filein, fileout) result(status)
    !
    !:Purpose: Copy the specified file to the new location and/or name
    !          This function is very general, but was initially written to
    !          copy files from the disk to the ram disk
    !
    !
    implicit none
    character(len=*) :: filein
    character(len=*) :: fileout
    integer :: status

    integer :: ierr, unitin, unitout
    integer(8) :: numChar
    character :: bufferB
    integer, parameter :: bufferSizeKB = 1024
    character :: bufferKB(bufferSizeKB)
    integer, parameter :: bufferSizeMB = 1024*1024
    character :: bufferMB(bufferSizeMB)

    write(*,*) 'utl_copyFile: copy from ', trim(filein), ' to ', trim(fileout)

    call utl_tmg_start(175,'low-level--utl_copyFile')

    unitin=10
    open(unit=unitin, file=trim(filein), status='OLD', form='UNFORMATTED', &
         action='READ', access='STREAM')

    unitout=11
    open(unit=unitout, file=trim(fileout), status='REPLACE', form='UNFORMATTED', &
         action='WRITE', access='STREAM')

    numChar = 0
    do 
      read(unitin,iostat=ierr) bufferMB
      if (ierr < 0) exit
      numChar = numChar + bufferSizeMB
      write(unitout) bufferMB
    end do

    do 
      read(unitin,iostat=ierr,pos=numChar+1) bufferKB
      if (ierr < 0) exit
      numChar = numChar + bufferSizeKB
      write(unitout) bufferKB
    end do

    do 
      read(unitin,iostat=ierr,pos=numChar+1) bufferB
      if (ierr < 0) exit
      numChar = numChar + 1
      write(unitout) bufferB
    end do

    write(*,*) 'utl_copyFile: copied ', numChar, ' bytes'

    close(unit=unitin)
    close(unit=unitout)

    if (numChar > 0) then
      status = 0
    else
      status = -1
      if (numChar == 0) then
        call utl_abort('utl_copyFile: ERROR, zero bytes copied')
      else
        ! Note: If 'numChar' becomes negative then it means it got bigger
        !       than the maximum integer the 'integer' type and so the
        !       variable 'numChar' wraps around and becomes negative.
        call utl_abort('utl_copyFile: ERROR, overflow detected since number of bytes copied is negative!')
      end if
    end if

    call tmg_stop(175)

  end function utl_copyFile

  !--------------------------------------------------------------------------
  ! utl_allReduce
  !--------------------------------------------------------------------------
  subroutine utl_allReduce(localGlobalValue)
    ! :Purpose: Perform mpi_allReduce to sum integer values over all
    !           mpi tasks and copy result back to same variable.

    implicit none

    ! Arguments:
    integer, intent(inout) :: localGlobalValue

    ! Locals:
    integer :: localValue, globalValue, ierr

    localValue = localGlobalValue
    call rpn_comm_allReduce(localValue, globalValue, 1, 'mpi_integer', &
                            'mpi_sum', 'grid', ierr)
    localGlobalValue = globalValue
    
  end subroutine utl_allReduce

  !--------------------------------------------------------------------------
  ! utl_findloc_char
  !--------------------------------------------------------------------------
  function utl_findloc_char(charArray, value) result(location)
    !
    ! :Purpose: A modified version of the fortran function `findloc`.
    !           If multiple matches are found in the array, a warning
    !           message is printed to the listing.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: charArray(:)
    character(len=*), intent(in) :: value
    integer                      :: location

    ! Locals:
    integer :: numFound, arrayIndex

    numFound = 0
    LOOP: do arrayIndex = 1, size(charArray)
      if (trim(charArray(arrayIndex)) == trim(value)) then
        numFound = numFound + 1
        ! return the first location found
        if (numFound == 1) location = arrayIndex
      end if
    end do LOOP

    ! give warning if more than 1 found
    if (numFound > 1) then
      write(*,*) 'utl_findloc_char: found multiple locations of ', trim(value)
      write(*,*) 'utl_findloc_char: number locations found =  ', numFound    
    end if

    ! return zero if not found
    if (numFound == 0) then
      location = 0
    end if

  end function utl_findloc_char

  !--------------------------------------------------------------------------
  ! utl_findloc_int
  !--------------------------------------------------------------------------
  function utl_findloc_int(intArray, value) result(location)
    !
    ! :Purpose: A modified version of the fortran function `findloc`.
    !           If multiple matches are found in the array, a warning
    !           message is printed to the listing.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: intArray(:)
    integer, intent(in) :: value
    integer             :: location

    ! Locals:
    integer :: numFound, arrayIndex

    numFound = 0
    LOOP: do arrayIndex = 1, size(intArray)
      if (intArray(arrayIndex) == value) then
        numFound = numFound + 1
        ! return the first location found
        if (numFound == 1) location = arrayIndex
      end if
    end do LOOP

    ! give warning if more than 1 found
    if (numFound > 1) then
      write(*,*) 'utl_findloc_int: found multiple locations of ', value
      write(*,*) 'utl_findloc_int: number locations found =  ', numFound    
    end if

    ! return zero if not found
    if (numFound == 0) then
      location = 0
    end if

  end function utl_findloc_int

  !--------------------------------------------------------------------------
  ! utl_findlocs_char
  !--------------------------------------------------------------------------
  function utl_findlocs_char(charArray, value) result(locations)
    !
    ! :Purpose: A modified version of the fortran function `findloc`.
    !           Returns an array of all matches found in the array.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: charArray(:)
    character(len=*), intent(in) :: value
    integer, allocatable         :: locations(:)

    ! Locals:
    integer :: numFound, arrayIndex

    if (allocated(locations)) deallocate(locations)

    ! count number of matches found
    numFound = 0
    do arrayIndex = 1, size(charArray)
      if (trim(charArray(arrayIndex)) == trim(value)) numFound = numFound + 1
    end do

    if (numFound > 0) then

      ! return all found locations
      allocate(locations(numFound))
      numFound = 0
      do arrayIndex = 1, size(charArray)
        if (trim(charArray(arrayIndex)) == trim(value)) then
          numFound = numFound + 1
          locations(numFound) = arrayIndex
        end if
      end do

    else

      ! return zero if not found
      allocate(locations(1))
      locations(1) = 0

    end if

  end function utl_findlocs_char

  !--------------------------------------------------------------------------
  ! utl_randomOrderInt
  !--------------------------------------------------------------------------
  subroutine utl_randomOrderInt(intArray,randomSeed)
    ! :Purpose: Randomly shuffle the order of the integer array elements.

    implicit none

    ! Arguments:
    integer, intent(inout) :: intArray(:)
    integer, intent(in)    :: randomSeed

    ! Locals:
    integer              :: arraySize, arrayIndex, arrayIndexMin
    integer, allocatable :: intArrayOut(:)
    real(8), allocatable :: realRandomArray(:)

    arraySize = size(intArray)
    allocate(realRandomArray(arraySize))
    allocate(intArrayOut(arraySize))

    call rng_setup(randomSeed)
    do arrayIndex = 1, arraySize
      realRandomArray(arrayIndex) = rng_uniform()
    end do

    do arrayIndex = 1, arraySize
      arrayIndexMin = minloc(realRandomArray,dim=1)
      realRandomArray(arrayIndexMin) = huge(1.0D0)
      intArrayOut(arrayIndex) = intArray(arrayIndexMin)
    end do

    intArray(:) = intArrayOut(:)

    deallocate(intArrayOut)
    deallocate(realRandomArray)

  end subroutine utl_randomOrderInt

  !--------------------------------------------------------------------------
  ! utl_tmg_start
  !--------------------------------------------------------------------------
  subroutine utl_tmg_start(blockIndex, blockLabel)
    ! :Purpose: Wrapper for rpnlib subroutine tmg_start

    implicit none

    ! Arguments:
    integer,          intent(in) :: blockIndex
    character(len=*), intent(in) :: blockLabel

    ! Locals:
    integer            :: labelLength
    integer, parameter :: labelPaddedLength = 40
    character(len=labelPaddedLength) :: blockLabelPadded

    blockLabelPadded = '........................................'
    labelLength = min(len_trim(blockLabel), labelPaddedLength)
    blockLabelPadded(1:labelLength) = blockLabel(1:labelLength)

    call tmg_start(blockIndex, blockLabelPadded)

  end subroutine utl_tmg_start
  
end module utilities_mod
