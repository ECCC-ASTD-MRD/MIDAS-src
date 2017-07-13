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
!! MODULE controlVector (prefix="cvm")
!!
!! *Purpose*: The control vector and related information.  
!!
!! Revisions:
!!       Ping Du, June-Sept 2014
!!       - Introduced optional consideration of arrays/variables for
!!         constituents. See variables *Bchm*.
!!       M. Sitwell Aug 2015
!!       - Added cvm_numSubvector
!--------------------------------------------------------------------------
module controlVector_mod
  use utilities_mod
  implicit none
  save
  private

  ! public variables
  public              :: cvm_nvadim, cvm_nvadim_mpiglobal
  public              :: cvm_BHI, cvm_BEN, cvm_BCHM
  ! public procedures
  public              :: cvm_Setup, cvm_getSubVector, cvm_getSubVector_mpiglobal
  public              :: cvm_getSubVector_r4, cvm_getSubVector_mpiglobal_r4
  public              :: cvm_subVectorExists


  logical             :: initialized = .false.
  integer             :: cvm_dimBHI
  integer             :: cvm_dimBEN
  integer             :: cvm_nvadim
  integer             :: cvm_dimBchm
  integer             :: cvm_dimBHI_mpiglobal
  integer             :: cvm_dimBEN_mpiglobal
  integer             :: cvm_nvadim_mpiglobal
  integer             :: cvm_dimBchm_mpiglobal

  integer, parameter  :: cvm_numSubvector=3 ! total number of possible valid subvectors
  integer, parameter  :: cvm_BHI  = 1
  integer, parameter  :: cvm_BEN  = 2
  integer, parameter  :: cvm_BCHM = 3

contains

  subroutine cvm_setup(dimBhi_in,dimBen_in,dimBchm_in)
    implicit none

    integer           :: dimBHI_in,dimBEN_in,dimBCHM_in
    integer           :: ierr

    cvm_dimbhi = dimbhi_in
    cvm_dimben = dimben_in
    cvm_dimbchm = dimbchm_in 

    call rpn_comm_allreduce(cvm_dimbhi,cvm_dimbhi_mpiglobal,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(cvm_dimben,cvm_dimben_mpiglobal,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    call rpn_comm_allreduce(cvm_dimbchm,cvm_dimbchm_mpiglobal,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)

    cvm_nvadim = cvm_dimben + cvm_dimbhi + cvm_dimbchm 
    cvm_nvadim_mpiglobal = cvm_dimben_mpiglobal + cvm_dimbhi_mpiglobal + cvm_dimbchm_mpiglobal 

    write(*,*) 'cvm_setup: subvector dimensions            =',cvm_dimbhi,cvm_dimben, cvm_dimbchm,  &
                                                              cvm_nvadim
    write(*,*) 'cvm_setup: subvector dimensions (mpiglobal)=',cvm_dimbhi_mpiglobal,cvm_dimben_mpiglobal,       & 
                          cvm_dimbchm_mpiglobal, cvm_nvadim_mpiglobal

    initialized=.true.

  end subroutine cvm_setup


  function cvm_subVectorExists(subVectorIndex) RESULT(exists)
    implicit none

    integer :: subVectorIndex
    logical :: exists

    if( subVectorIndex < 0 .or. subVectorIndex > cvm_numSubvector ) then
      write(*,*) 'cvm_getSubVector: subVectorIndex = ',subVectorIndex
      call utl_abort('cvm_getSubVector: invalid value for subVectorIndex')
    endif

    exists = .false.

    if( subVectorIndex == cvm_BHI .and. cvm_dimbhi_mpiglobal > 0 ) then
      exists = .true.
    elseif( subVectorIndex == cvm_BEN .and. cvm_dimben_mpiglobal > 0 ) then
      exists = .true.
    elseif( subVectorIndex == cvm_BCHM .and. cvm_dimbchm_mpiglobal > 0 ) then
      exists = .true.
    endif

  end function cvm_subVectorExists


  function cvm_getSubVector(controlVector,subVectorIndex) RESULT(subVector)
    implicit none

    real*8, pointer :: subVector(:)
    real*8, target  :: controlVector(:)
    integer         :: subVectorIndex
    logical, save   :: firstCall(cvm_numSubvector)=.true.

    if( subVectorIndex < 0 .or. subVectorIndex > cvm_numSubvector ) then
       write(*,*) 'cvm_getSubVector: subVectorIndex = ',subVectorIndex
       call utl_abort('cvm_getSubVector: invalid value for subVectorIndex')
    end if

    nullify(subVector)

    if( subVectorIndex == cvm_BHI ) then
      subVector => controlVector(1:cvm_dimbhi)
    elseif( subVectorIndex == cvm_BEN ) then
      subVector => controlVector((cvm_dimbhi+1):(cvm_dimbhi+cvm_dimben))
    elseif( subVectorIndex == cvm_BCHM ) then
      subVector => controlVector((cvm_dimbhi+cvm_dimben+1):(cvm_dimbhi+cvm_dimben+cvm_dimbchm))
    else
      call utl_abort('cvm_getSubVector: unknown subVectorIndex!')
    endif

  end function cvm_getSubVector


  function cvm_getSubVector_r4(controlVector,subVectorIndex) RESULT(subVector)
    implicit none

    real*4, pointer :: subVector(:)
    real*4, target  :: controlVector(:)
    integer         :: subVectorIndex
    logical, save   :: firstCall(cvm_numSubvector)=.true.

    if( subVectorIndex < 0 .or. subVectorIndex > cvm_numSubvector ) then
       write(*,*) 'cvm_getSubVector_r4: subVectorIndex = ',subVectorIndex
       call utl_abort('cvm_getSubVector_r4: invalid value for subVectorIndex')
    end if

    nullify(subVector)

    if( subVectorIndex == cvm_BHI ) then
      subVector => controlVector(1:cvm_dimbhi)
    elseif( subVectorIndex == cvm_BEN ) then
      subVector => controlVector((cvm_dimbhi+1):(cvm_dimbhi+cvm_dimben))
    elseif( subVectorIndex == cvm_BCHM ) then
      subVector => controlVector((cvm_dimbhi+cvm_dimben+1):(cvm_dimbhi+cvm_dimben+cvm_dimbchm))
    else
      call utl_abort('cvm_getSubVector_r4: unknown subVectorIndex!')
    endif

  end function cvm_getSubVector_r4


  function cvm_getSubVector_mpiglobal(controlVector,subVectorIndex) RESULT(subVector)
    implicit none

    real*8, pointer :: subVector(:)
    real*8, target  :: controlVector(:)
    integer         :: subVectorIndex
    logical, save   :: firstCall(cvm_numSubvector)=.true.

    if( subVectorIndex < 0 .or. subVectorIndex > cvm_numSubvector ) then
       write(*,*) 'cvm_getSubVector_mpiglobal: subVectorIndex = ',subVectorIndex
       call utl_abort('cvm_getSubVector_mpiglobal: invalid value for subVectorIndex')
    end if

    nullify(subVector)

    if( subVectorIndex == cvm_BHI ) then
      subVector => controlVector(1:cvm_dimbhi_mpiglobal)
    elseif( subVectorIndex == cvm_BEN ) then
      subVector => controlVector((cvm_dimbhi_mpiglobal+1):(cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal))
    elseif( subVectorIndex == cvm_BCHM ) then
      subVector => controlVector((cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal+1):(cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal+cvm_dimbchm_mpiglobal))
    else
      call utl_abort('cvm_getSubVector_mpiglobal: unknown subVectorIndex!')
    endif

  end function cvm_getSubVector_mpiglobal


  function cvm_getSubVector_mpiglobal_r4(controlVector,subVectorIndex) RESULT(subVector)
    implicit none

    real*4, pointer :: subVector(:)
    real*4, target  :: controlVector(:)
    integer         :: subVectorIndex
    logical, save   :: firstCall(cvm_numSubvector)=.true.

    if( subVectorIndex < 0 .or. subVectorIndex > cvm_numSubvector ) then
       write(*,*) 'cvm_getSubVector_mpiglobal_r4: subVectorIndex = ',subVectorIndex
       call utl_abort('cvm_getSubVector_mpiglobal_r4: invalid value for subVectorIndex')
    end if

    nullify(subVector)

    if( subVectorIndex == cvm_BHI ) then
      subVector => controlVector(1:cvm_dimbhi_mpiglobal)
    elseif( subVectorIndex == cvm_BEN ) then
      subVector => controlVector((cvm_dimbhi_mpiglobal+1):(cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal))
    elseif( subVectorIndex == cvm_BCHM ) then
      subVector => controlVector((cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal+1):(cvm_dimbhi_mpiglobal+cvm_dimben_mpiglobal+cvm_dimbchm_mpiglobal))
    else
      call utl_abort('cvm_getSubVector_mpiglobal_r4: unknown subVectorIndex!')
    endif

  end function cvm_getSubVector_mpiglobal_r4

end module controlVector_mod
