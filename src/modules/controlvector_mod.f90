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

module controlVector_mod
  ! MODULE controlVector_mod (prefix='cvm' category='6. High-level data objects')
  !
  ! :Purpose: The control vector and related information.  
  !
  use utilities_mod
  implicit none
  save
  private

  ! public variables
  public              :: cvm_nvadim, cvm_nvadim_mpiglobal
  ! public procedures
  public              :: cvm_setupSubVector, cvm_getSubVector, cvm_getSubVector_mpiglobal
  public              :: cvm_getSubVector_r4, cvm_getSubVector_mpiglobal_r4
  public              :: cvm_subVectorExists

  type struct_cvm
    private
    character(len=9) :: label                  = 'XXXXXXXXX'
    character(len=4) :: BmatrixType            = 'XXXX'
    integer          :: dimVector              = 0
    integer          :: subVectorBeg           = 1
    integer          :: subVectorEnd           = 0
    integer          :: dimVector_mpiglobal    = 0
    integer          :: subVectorBeg_mpiglobal = 1
    integer          :: subVectorEnd_mpiglobal = 0
  end type struct_cvm

  integer, parameter :: maxNumVectors = 50
  integer            :: numVectors = 0
  type(struct_cvm)   :: cvm_vector(maxNumVectors)

  integer             :: cvm_nvadim
  integer             :: cvm_nvadim_mpiglobal

contains

  subroutine cvm_setupSubVector(label, BmatrixType, dimVector)
    implicit none

    ! Arguments
    character(len=*) :: label
    character(len=*) :: BmatrixType
    integer :: dimVector

    ! Locals
    integer :: ierr, dimVector_mpiglobal

    if ( numVectors == maxNumVectors ) then
      call utl_abort('cvm_setupSubVector: number of allocated subvectors already at maximum allowed')
    end if

    call rpn_comm_allreduce(dimVector, dimVector_mpiglobal,  &
                            1, 'MPI_INTEGER', 'MPI_SUM', 'GRID', ierr)

    ! just return if subVector dimension is zero on all MPI tasks
    if ( dimVector_mpiglobal == 0 ) return

    numVectors = numVectors + 1

    if ( any(cvm_vector(:)%label == label) ) then
      write(*,*) 'cvm_setupSubVector: label = ', trim(label)
      call utl_abort('cvm_setupSubVector: this label is already present')
    end if

    cvm_vector(numVectors)%label = label
    cvm_vector(numVectors)%BmatrixType = BmatrixType
    cvm_vector(numVectors)%dimVector = dimVector
    cvm_vector(numVectors)%dimVector_mpiglobal = dimVector_mpiglobal

    if ( numVectors == 1 ) then
      cvm_vector(numVectors)%subVectorBeg = 1
      cvm_vector(numVectors)%subVectorEnd = cvm_vector(numVectors)%dimVector

      cvm_vector(numVectors)%subVectorBeg_mpiglobal = 1
      cvm_vector(numVectors)%subVectorEnd_mpiglobal = cvm_vector(numVectors)%dimVector_mpiglobal
    else
      cvm_vector(numVectors)%subVectorBeg = 1 + cvm_vector(numVectors-1)%subVectorEnd
      cvm_vector(numVectors)%subVectorEnd = cvm_vector(numVectors)%dimVector +   &
                                            cvm_vector(numVectors-1)%subVectorEnd

      cvm_vector(numVectors)%subVectorBeg_mpiglobal = 1 + cvm_vector(numVectors-1)%subVectorEnd_mpiglobal
      cvm_vector(numVectors)%subVectorEnd_mpiglobal = cvm_vector(numVectors)%dimVector_mpiglobal +  &
                                                      cvm_vector(numVectors-1)%subVectorEnd_mpiglobal
    end if

    cvm_nvadim = sum(cvm_vector(1:numVectors)%dimVector)
    cvm_nvadim_mpiglobal = sum(cvm_vector(1:numVectors)%dimVector_mpiglobal)

    write(*,*) 'cvm_setupSubVector: '
    write(*,*) '   added subVector with label                 = ', cvm_vector(numVectors)%label
    write(*,*) '   added subVector of type                    = ', cvm_vector(numVectors)%Bmatrixtype
    write(*,*) '   added subVector with dimension             = ', cvm_vector(numVectors)%dimVector
    write(*,*) '   added subVector with dimension (mpiglobal) = ', cvm_vector(numVectors)%dimVector_mpiglobal
    write(*,*) '   current total dimension             = ', cvm_nvadim
    write(*,*) '   current total dimension (mpiglobal) = ', cvm_nvadim_mpiglobal

  end subroutine cvm_setupSubVector


  function cvm_indexFromLabel(subVectorLabel) result(subVectorIndex)
    implicit none
    integer :: subVectorIndex

    ! Arguments
    character(len=*) :: subVectorLabel

    ! Locals
    logical :: found

    found = .false.
    index_loop: do subVectorIndex = 1, numVectors
      if ( trim(cvm_vector(subVectorIndex)%label) == trim(subVectorLabel) ) then
        found = .true.
        exit index_loop
      end if
    end do index_loop

    if ( .not.found ) then
      subVectorIndex = -1
    end if

  end function cvm_indexFromLabel


  function cvm_subVectorExists(subVectorLabel) result(exists)
    implicit none
    logical :: exists

    ! Arguments
    character(len=*) :: subVectorLabel

    ! Locals
    integer :: subVectorIndex

    subVectorIndex = cvm_indexFromLabel(subVectorLabel)

    if ( subVectorIndex < 0 ) then
      exists = .false.
      return
    end if

    if ( cvm_vector(subVectorIndex)%dimVector_mpiglobal > 0 ) then
      exists = .true.
    else
      exists = .false.
    end if

  end function cvm_subVectorExists


  function cvm_getSubVector(controlVector,subVectorLabel) result(subVector)
    implicit none
    real*8, pointer :: subVector(:)

    ! Arguments
    character(len=*) :: subVectorLabel
    real*8, target  :: controlVector(:)

    ! Locals
    integer         :: subVectorIndex, indexBeg, indexEnd

    subVectorIndex = cvm_indexFromLabel(subVectorLabel)

    if( subVectorIndex < 0 ) then
      call utl_abort('cvm_getSubVector: invalid subVector label')
    end if

    indexBeg = cvm_vector(subVectorIndex)%subVectorBeg
    indexEnd = cvm_vector(subVectorIndex)%subVectorEnd
    subVector => controlVector(indexBeg:indexEnd)

  end function cvm_getSubVector


  function cvm_getSubVector_r4(controlVector,subVectorLabel) result(subVector)
    implicit none
    real*4, pointer :: subVector(:)

    ! Arguments
    character(len=*) :: subVectorLabel
    real*4, target  :: controlVector(:)

    ! Locals
    integer         :: subVectorIndex, indexBeg, indexEnd

    subVectorIndex = cvm_indexFromLabel(subVectorLabel)

    if( subVectorIndex < 0 ) then
      call utl_abort('cvm_getSubVector_r4: invalid subVector label')
    end if

    indexBeg = cvm_vector(subVectorIndex)%subVectorBeg
    indexEnd = cvm_vector(subVectorIndex)%subVectorEnd
    subVector => controlVector(indexBeg:indexEnd)

  end function cvm_getSubVector_r4


  function cvm_getSubVector_mpiglobal(controlVector,subVectorLabel) result(subVector)
    implicit none
    real*8, pointer :: subVector(:)

    ! Arguments
    character(len=*) :: subVectorLabel
    real*8, target  :: controlVector(:)

    ! Locals
    integer         :: subVectorIndex, indexBeg, indexEnd

    subVectorIndex = cvm_indexFromLabel(subVectorLabel)

    if( subVectorIndex < 0 ) then
      call utl_abort('cvm_getSubVector_mpiglobal: invalid subVector label')
    end if

    indexBeg = cvm_vector(subVectorIndex)%subVectorBeg_mpiglobal
    indexEnd = cvm_vector(subVectorIndex)%subVectorEnd_mpiglobal
    subVector => controlVector(indexBeg:indexEnd)

  end function cvm_getSubVector_mpiglobal


  function cvm_getSubVector_mpiglobal_r4(controlVector,subVectorLabel) result(subVector)
    implicit none
    real*4, pointer :: subVector(:)

    ! Arguments
    character(len=*) :: subVectorLabel
    real*4, target  :: controlVector(:)

    ! Locals
    integer         :: subVectorIndex, indexBeg, indexEnd

    subVectorIndex = cvm_indexFromLabel(subVectorLabel)

    if( subVectorIndex < 0 ) then
      call utl_abort('cvm_getSubVector_mpiglobal_r4: invalid subVector label')
    end if

    indexBeg = cvm_vector(subVectorIndex)%subVectorBeg_mpiglobal
    indexEnd = cvm_vector(subVectorIndex)%subVectorEnd_mpiglobal
    subVector => controlVector(indexBeg:indexEnd)

  end function cvm_getSubVector_mpiglobal_r4

end module controlVector_mod
