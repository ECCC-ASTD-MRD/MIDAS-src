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

module ensembleStateVector_mod
  ! MODULE ensembleStateVector_mod (prefix='ens' category='2. High-level data objects')
  !
  ! :Purpose: Store and manipulate ensemble of state vectors and the ensemble
  !           mean.
  !
  use ramDisk_mod
  use mpi_mod
  use mpivar_mod
  use fileNames_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use mathPhysConstants_mod
  use utilities_mod
  use varNameList_mod
  implicit none
  save
  private

  ! public procedures
  public :: struct_ens, ens_allocate, ens_deallocate, ens_zero
  public :: ens_readEnsemble, ens_writeEnsemble, ens_copy, ens_copy4Dto3D, ens_add
  public :: ens_getOneLevMean_r8, ens_modifyVarName
  public :: ens_varExist, ens_getNumLev, ens_getNumMembers, ens_getNumSubEns
  public :: ens_computeMean, ens_removeMean, ens_recenter
  public :: ens_copyEnsMean, ens_copyToEnsMean, ens_copyMember, ens_insertMember
  public :: ens_computeStdDev, ens_copyEnsStdDev, ens_normalize
  public :: ens_getOneLev_r4, ens_getOneLev_r8
  public :: ens_getOffsetFromVarName, ens_getLevFromK, ens_getVarNameFromK 
  public :: ens_getNumK, ens_getKFromLevVarName, ens_getDataKind
  public :: ens_getVco, ens_getHco, ens_getLatLonBounds, ens_getNumStep
  public :: ens_varNamesList

  integer,external   :: get_max_rss

  type :: struct_oneLev_r4
    real(4), pointer :: onelevel(:,:,:,:) => null()
  end type struct_oneLev_r4

  type :: struct_oneLev_r8
    real(8), pointer :: onelevel(:,:,:,:) => null()
  end type struct_oneLev_r8

  type :: struct_ens
    private
    logical                       :: allocated = .false.
    integer                       :: numMembers
    integer                       :: dataKind = 4 ! default value
    type(struct_gsv)              :: statevector_work
    type(struct_oneLev_r8), allocatable :: allLev_ensMean_r8(:), allLev_ensStdDev_r8(:)
    type(struct_oneLev_r4), allocatable :: allLev_r4(:)
    type(struct_oneLev_r8), allocatable :: allLev_r8(:)
    logical                       :: meanIsComputed = .false.
    logical                       :: stdDevIsComputed = .false.
    logical                       :: meanIsRemoved = .false.
    integer, allocatable          :: subEnsIndexList(:), nEnsSubEns(:)
    integer                       :: numSubEns
    character(len=256)            :: enspathname
    character(len=12)             :: hInterpolateDegree='UNSPECIFIED' ! 'LINEAR' or 'CUBIC' or 'NEAREST'
  end type struct_ens

CONTAINS

  !--------------------------------------------------------------------------
  ! ens_allocate
  !--------------------------------------------------------------------------
  subroutine ens_allocate(ens, numMembers, numStep, hco_ens, vco_ens, &
                          dateStampList, varNames_opt, dataKind_opt, &
                          hInterpolateDegree_opt)
    !
    !:Purpose: Allocate an ensembleStateVector object
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens
    integer :: numMembers, numStep
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    integer :: dateStampList(:)
    character(len=*), optional :: varNames_opt(:)  ! allow specification of assigned variables
    integer, optional          :: dataKind_opt
    character(len=*), optional :: hInterpolateDegree_opt

    ! Locals:
    integer :: jk, lon1, lon2, lat1, lat2, k1, k2
    character(len=4), pointer :: varNames(:)

    if ( ens%allocated ) then
      write(*,*) 'ens_allocate: this object is already allocated, deallocating first.'
      call ens_deallocate( ens )
    end if

    if ( present(dataKind_opt) ) ens%dataKind = dataKind_opt

    if ( present(hInterpolateDegree_opt) ) then
      ! set the horizontal interpolation degree
      ens%hInterpolateDegree = trim(hInterpolateDegree_opt)
    else
      ens%hInterpolateDegree = 'LINEAR' ! default, for legacy purposes
    end if

    nullify(varNames)
    if (present(varNames_opt)) then
      allocate(varNames(size(varNames_opt)))
      varNames(:) = varNames_opt(:)
    else
      call gsv_varNamesList(varNames)
    end if

    call gsv_allocate( ens%statevector_work, &
                       numStep, hco_ens, vco_ens,  &
                       datestamplist_opt=dateStampList, mpi_local_opt=.true., &
                       varNames_opt=varNames, dataKind_opt=ens%dataKind, &
                       hInterpolateDegree_opt = hInterpolateDegree_opt)

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd

    if (ens%dataKind == 8) then
      allocate( ens%allLev_r8(k1:k2) )
      do jk = k1, k2
        allocate( ens%allLev_r8(jk)%onelevel(numMembers,numStep,lon1:lon2,lat1:lat2) )
      end do
    else if (ens%dataKind == 4) then
      allocate( ens%allLev_r4(k1:k2) )
      do jk = k1, k2
        allocate( ens%allLev_r4(jk)%onelevel(numMembers,numStep,lon1:lon2,lat1:lat2) )
      end do
    else
      call utl_abort('ens_allocate: unknown value of datakind')
    end if

    ens%allocated = .true.
    ens%numMembers = numMembers

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocate

  !--------------------------------------------------------------------------
  ! ens_allocateMean
  !--------------------------------------------------------------------------
  subroutine ens_allocateMean(ens)
    !
    !:Purpose: Allocate the ensemble mean arrays within an ensembleStateVector object
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens

    ! Locals:
    integer :: lon1, lon2, lat1, lat2, k1, k2, jk, numStep

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    allocate( ens%allLev_ensMean_r8(k1:k2) )
    do jk = k1, k2
      allocate( ens%allLev_ensMean_r8(jk)%onelevel(ens%numSubEns,numStep,lon1:lon2,lat1:lat2) )
      ens%allLev_ensMean_r8(jk)%onelevel(:,:,:,:) = 0.0d0
    end do

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocateMean

  !--------------------------------------------------------------------------
  ! ens_allocateStdDev
  !--------------------------------------------------------------------------
  subroutine ens_allocateStdDev(ens)
    !
    !:Purpose: Allocate the ensemble stddev arrays within an ensembleStateVector object
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens

    ! Locals:
    integer :: lon1, lon2, lat1, lat2, k1, k2, jk, numStep

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    allocate( ens%allLev_ensStdDev_r8(k1:k2) )
    do jk = k1, k2
      allocate( ens%allLev_ensStdDev_r8(jk)%onelevel(1,numStep,lon1:lon2,lat1:lat2) )
      ens%allLev_ensStdDev_r8(jk)%onelevel(:,:,:,:) = 0.0d0
    end do

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocateStdDev

  !--------------------------------------------------------------------------
  ! ens_deallocate
  !--------------------------------------------------------------------------
  subroutine ens_deallocate( ens )
    !
    !:Purpose: Deallocate an ensembleStateVector object
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens

    ! Locals:
    integer :: k1, k2, jk

    if ( .not. ens%allocated ) return

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd

    if (ens%dataKind == 8) then
      do jk = k1, k2
        deallocate( ens%allLev_r8(jk)%onelevel )
      end do
      deallocate( ens%allLev_r8 )
    else if (ens%dataKind == 4) then
      do jk = k1, k2
        deallocate( ens%allLev_r4(jk)%onelevel )
      end do
      deallocate( ens%allLev_r4 )
    end if

    if (ens%stdDevIsComputed) then
      do jk = k1, k2
        deallocate( ens%allLev_ensStdDev_r8(jk)%onelevel )
      end do
      deallocate( ens%allLev_ensStdDev_r8 )
    end if

    if (ens%meanIsComputed) then
      do jk = k1, k2
        deallocate( ens%allLev_ensMean_r8(jk)%onelevel )
      end do
      deallocate( ens%allLev_ensMean_r8 )
      deallocate( ens%subEnsIndexList )
      deallocate( ens%nEnsSubEns )
    end if

    ens%allocated = .false.

  end subroutine ens_deallocate

  !--------------------------------------------------------------------------
  ! ens_modifyVarName
  !--------------------------------------------------------------------------
  subroutine ens_modifyVarName(ens, oldVarName, newVarName) 
    !
    !:Purpose: Change an existing variable name within the ensemble.
    !          This is only used when the contents of a variable are
    !          transformed into another variable **in place**.
    !
    implicit none
    
    ! Arguments:
    type(struct_ens) :: ens
    character(len=*) :: oldVarName
    character(len=*) :: newVarName
    
    call gsv_modifyVarName(ens%statevector_work,oldVarName, newVarName)
    
  end subroutine ens_modifyVarName

  !--------------------------------------------------------------------------
  ! ens_copy
  !--------------------------------------------------------------------------
  subroutine ens_copy(ens_in,ens_out)
    !
    !:Purpose: Copy the contents of one ensembleStateVector object into another
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens_in, ens_out

    ! Locals:
    integer           :: lon1, lon2, lat1, lat2, k1, k2
    integer           :: jk, stepIndex, latIndex, lonIndex, memberIndex

    if (.not.ens_in%allocated) then
      call utl_abort('ens_copy: ens_in not yet allocated')
    end if
    if (.not.ens_out%allocated) then
      call utl_abort('ens_copy: ens_out not yet allocated')
    end if

    lon1 = ens_out%statevector_work%myLonBeg
    lon2 = ens_out%statevector_work%myLonEnd
    lat1 = ens_out%statevector_work%myLatBeg
    lat2 = ens_out%statevector_work%myLatEnd
    k1   = ens_out%statevector_work%mykBeg
    k2   = ens_out%statevector_work%mykEnd
 
    if ( ens_out%dataKind == 8 .and. ens_in%dataKind == 8 ) then

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens_out%statevector_work%numStep
              do memberIndex = 1, ens_out%numMembers
                ens_out%allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = &
                     ens_in %allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex)
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if ( ens_out%dataKind == 4 .and. ens_in%dataKind == 4 ) then

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens_out%statevector_work%numStep
              do memberIndex = 1, ens_out%numMembers
                ens_out%allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = &
                     ens_in %allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex)
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else
      call utl_abort('ens_copy: Data type must be the same for both ensembleStatevectors')
    end if

  end subroutine ens_copy

  !--------------------------------------------------------------------------
  ! ens_copy4Dto3D
  !--------------------------------------------------------------------------
  subroutine ens_copy4Dto3D(ens_in,ens_out)
    !
    !:Purpose: Copy contents of a 4D ensemble into a 3D ensemble object by
    !          extracting the middle time step.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens_in, ens_out

    ! Locals:
    integer           :: lon1, lon2, lat1, lat2, k1, k2
    integer           :: jk, latIndex, lonIndex, memberIndex
    integer           :: numStepIn, numStepOut, middleStepIndex

    if (.not.ens_in%allocated) then
      call utl_abort('ens_copy4Dto3D: ens_in not yet allocated')
    end if
    if (.not.ens_out%allocated) then
      call utl_abort('ens_copy4Dto3D: ens_out not yet allocated')
    end if

    lon1 = ens_out%statevector_work%myLonBeg
    lon2 = ens_out%statevector_work%myLonEnd
    lat1 = ens_out%statevector_work%myLatBeg
    lat2 = ens_out%statevector_work%myLatEnd
    k1   = ens_out%statevector_work%mykBeg
    k2   = ens_out%statevector_work%mykEnd
    numStepIn  =  ens_in%statevector_work%numStep
    numStepOut =  ens_out%statevector_work%numStep

    if (numStepOut /= 1) call utl_abort('ens_copy4Dto3D: output ensemble must have only 1 timestep')
    if (numStepIn == 1) then
      write(*,*) 'ens_copy4Dto3D: WARNING: input ensemble only has 1 timestep, will simply copy.'
    end if
    middleStepIndex = (numStepIn + 1) / 2

    if ( ens_out%dataKind == 8 .and. ens_in%dataKind == 8 ) then

      !$OMP PARALLEL DO PRIVATE (jk,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do memberIndex = 1, ens_out%numMembers
              ens_out%allLev_r8(jk)%onelevel(memberIndex,1,lonIndex,latIndex) = &
                   ens_in %allLev_r8(jk)%onelevel(memberIndex,middleStepIndex,lonIndex,latIndex)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if ( ens_out%dataKind == 4 .and. ens_in%dataKind == 4 ) then

      !$OMP PARALLEL DO PRIVATE (jk,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do memberIndex = 1, ens_out%numMembers
              ens_out%allLev_r4(jk)%onelevel(memberIndex,1,lonIndex,latIndex) = &
                   ens_in %allLev_r4(jk)%onelevel(memberIndex,middleStepIndex,lonIndex,latIndex)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else
      call utl_abort('ens_copy4Dto3D: Data type must be the same for both ensembleStatevectors')
    end if

  end subroutine ens_copy4Dto3D

  !--------------------------------------------------------------------------
  ! ens_add
  !--------------------------------------------------------------------------
  subroutine ens_add(ens_in, ens_inOut, scaleFactorIn_opt, scaleFactorInOut_opt)
    !
    !:Purpose: Add the contents of the ens_in ensemble to the ens_inOut ensemble.
    !
    implicit none

    ! arguments
    type(struct_ens)  :: ens_in, ens_inOut
    real(8), optional :: scaleFactorIn_opt
    real(8), optional :: scaleFactorInOut_opt

    ! locals
    integer           :: lon1, lon2, lat1, lat2, k1, k2
    integer           :: jk, stepIndex, latIndex, lonIndex, memberIndex
    real(4)           :: scaleFactorIn_r4, scaleFactorInOut_r4
    real(8)           :: scaleFactorIn, scaleFactorInOut

    if (.not.ens_in%allocated) then
      call utl_abort('ens_add: ens_in not yet allocated')
    end if
    if (.not.ens_inOut%allocated) then
      call utl_abort('ens_add: ens_inOut not yet allocated')
    end if

    lon1 = ens_inOut%statevector_work%myLonBeg
    lon2 = ens_inOut%statevector_work%myLonEnd
    lat1 = ens_inOut%statevector_work%myLatBeg
    lat2 = ens_inOut%statevector_work%myLatEnd
    k1   = ens_inOut%statevector_work%mykBeg
    k2   = ens_inOut%statevector_work%mykEnd
 
    if ( ens_inOut%dataKind == 8 .and. ens_in%dataKind == 8 ) then

      if (present(scaleFactorIn_opt)) then
        scaleFactorIn = scaleFactorIn_opt
      else
        scaleFactorIn = 1.0d0
      end if
      if (present(scaleFactorInOut_opt)) then
        scaleFactorInOut = scaleFactorInOut_opt
      else
        scaleFactorInOut = 1.0d0
      end if

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens_inOut%statevector_work%numStep
              do memberIndex = 1, ens_inOut%numMembers
                ens_inOut%allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = &
                     scaleFactorInOut * ens_inOut%allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) + &
                     scaleFactorIn    * ens_in%allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex)
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if ( ens_inOut%dataKind == 4 .and. ens_in%dataKind == 4 ) then

      if (present(scaleFactorIn_opt)) then
        scaleFactorIn_r4 = real(scaleFactorIn_opt,4)
      else
        scaleFactorIn_r4 = 1.0
      end if
      if (present(scaleFactorInOut_opt)) then
        scaleFactorInOut_r4 = real(scaleFactorInOut_opt,4)
      else
        scaleFactorInOut_r4 = 1.0
      end if

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens_inOut%statevector_work%numStep
              do memberIndex = 1, ens_inOut%numMembers
                ens_inOut%allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = &
                     scaleFactorInOut_r4 * ens_inOut%allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) + &
                     scaleFactorIn_r4    * ens_in%allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex)
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else
      call utl_abort('ens_add: Data type must be the same for both ensembleStatevectors')
    end if

  end subroutine ens_add

  !--------------------------------------------------------------------------
  ! ens_zero
  !--------------------------------------------------------------------------
  subroutine ens_zero(ens)
    !
    !:Purpose: Set the main contents of an ensembleStateVector object to zero
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens

    ! Locals:
    integer           :: lon1, lon2, lat1, lat2, k1, k2
    integer           :: jk, stepIndex, latIndex, lonIndex, memberIndex

    if (.not.ens%allocated) then
      call utl_abort('ens_zero: ens not yet allocated')
    end if

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1   = ens%statevector_work%mykBeg
    k2   = ens%statevector_work%mykEnd
 
    if ( ens%dataKind == 8 ) then

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens%statevector_work%numStep
              do memberIndex = 1, ens%numMembers
                ens%allLev_r8(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = 0.d0
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if ( ens%dataKind == 4 ) then

      !$OMP PARALLEL DO PRIVATE (jk,stepIndex,latIndex,lonIndex,memberIndex)    
      do jk = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do stepIndex = 1, ens%statevector_work%numStep
              do memberIndex = 1, ens%numMembers
                ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,lonIndex,latIndex) = 0.0
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end if

  end subroutine ens_zero

  !--------------------------------------------------------------------------
  ! ens_copyToStateWork (private routine)
  !--------------------------------------------------------------------------
  subroutine ens_copyToStateWork(ens, dataType, memberIndex_opt, &
                                 subEnsIndex_opt)
    !
    !:Purpose: Copy the selected contents of an ensembleStateVector into
    !          the 'work' stateVector within the object. The possible
    !          types of content are: 'member', 'mean', or 'stdDev'.
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens
    character(len=*) :: dataType
    integer, optional:: memberIndex_opt
    integer, optional:: subEnsIndex_opt

    ! Locals:
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    integer          :: k1, k2, jk, numStep, stepIndex
    integer          :: memberIndex, subEnsIndex

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    select case(trim(dataType))
    case ('member')
      if (present(memberIndex_opt)) then
        memberIndex = memberIndex_opt
      else
        call utl_abort('ens_copyToStateWork: memberIndex_opt must be provided with dataType=member')
      end if
    case ('mean','stdDev')
      if (present(subEnsIndex_opt)) then
        subEnsIndex = subEnsIndex_opt
      else
        subEnsIndex = 1
      end if
    case default
      write(*,*)
      write(*,*) 'Unsupported dataType: ', trim(dataType)
      write(*,*) '    please select either: member, mean or stdDev'
      call utl_abort('ens_copyToStateWork')
    end select

    if (ens%dataKind == 8) then
      ptr4d_r8 => gsv_getField_r8(ens%statevector_work)
      do stepIndex = 1, numStep
        do jk = k1, k2
          if (dataType == 'member') then
            ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_r8(jk)%onelevel(memberIndex,stepIndex,:,:)
          else if (dataType == 'mean') then
            ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:)
          else if (dataType == 'stdDev') then
            ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_ensStdDev_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:)
          end if
        end do
      end do
    else if (ens%dataKind == 4) then
      ptr4d_r4 => gsv_getField_r4(ens%statevector_work)
      do stepIndex = 1, numStep
        do jk = k1, k2
          if (dataType == 'member') then
            ptr4d_r4(:,:,jk,stepIndex) = ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:)
          else if (dataType == 'mean') then
            ptr4d_r4(:,:,jk,stepIndex) = real(ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:),4)
          else if (dataType == 'stdDev') then
            ptr4d_r4(:,:,jk,stepIndex) = real(ens%allLev_ensStdDev_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:),4)
          end if
        end do
      end do
    end if

  end subroutine ens_copyToStateWork

  !--------------------------------------------------------------------------
  ! ens_copyFromStateWork
  !--------------------------------------------------------------------------
  subroutine ens_copyFromStateWork(ens, dataType, memberIndex_opt, &
                                   subEnsIndex_opt)
    !
    !:Purpose: This is the inverse operation as ens_copyToStateWork.
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens
    character(len=*) :: dataType
    integer, optional:: memberIndex_opt
    integer, optional:: subEnsIndex_opt

    ! Locals:
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    integer          :: k1, k2, jk, numStep, stepIndex
    integer          :: memberIndex, subEnsIndex

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    select case(trim(dataType))
    case ('member')
      if (present(memberIndex_opt)) then
        memberIndex = memberIndex_opt
      else
        call utl_abort('ens_copyFromStateWork: memberIndex_opt must be provided with dataType=member')
      end if
    case ('mean','stdDev')
      if (present(subEnsIndex_opt)) then
        subEnsIndex = subEnsIndex_opt
      else
        subEnsIndex = 1
      end if
    case default
      write(*,*)
      write(*,*) 'Unsupported dataType: ', trim(dataType)
      write(*,*) '    please select either: member, mean or stdDev'
      call utl_abort('ens_copyToStateWork')
    end select

    if (ens%dataKind == 8) then
      ptr4d_r8 => gsv_getField_r8(ens%statevector_work)
      do stepIndex = 1, numStep
        do jk = k1, k2
          if (dataType == 'member') then
            ens%allLev_r8(jk)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r8(:,:,jk,stepIndex)
          else if (dataType == 'mean') then
            ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = ptr4d_r8(:,:,jk,stepIndex)
          else if (dataType == 'stdDev') then
            ens%allLev_ensStdDev_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = ptr4d_r8(:,:,jk,stepIndex)
          end if
        end do
      end do
    else if (ens%dataKind == 4) then
      ptr4d_r4 => gsv_getField_r4(ens%statevector_work)
      do stepIndex = 1, numStep
        do jk = k1, k2
          if (dataType == 'member') then
            ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r4(:,:,jk,stepIndex)
          else if (dataType == 'mean') then
            ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = real(ptr4d_r4(:,:,jk,stepIndex),8)
          else if (dataType == 'stdDev') then
            ens%allLev_ensStdDev_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = real(ptr4d_r4(:,:,jk,stepIndex),8)
          end if
        end do
      end do
    end if

  end subroutine ens_copyFromStateWork

  !--------------------------------------------------------------------------
  ! ens_getOneLev_r4
  !--------------------------------------------------------------------------
  function ens_getOneLev_r4(ens,kIndex) result(oneLevLevel)
    !
    !:Purpose: Return a 4D pointer to a single level of a real4 ensemble. The
    !          dimensions of the pointer are:
    !          (memberIndex, stepIndex, lonIndex, latIndex)
    !
    implicit none

    ! Arguments:
    real(4), pointer :: oneLevLevel(:,:,:,:)
    type(struct_ens) :: ens
    integer          :: kIndex

    ! Locals:
    integer          :: lon1, lat1

    lon1 = ens%statevector_work%myLonBeg
    lat1 = ens%statevector_work%myLatBeg

    oneLevLevel(1:,1:,lon1:,lat1:) => ens%allLev_r4(kIndex)%onelevel(:,:,:,:)

  end function ens_getOneLev_r4

  !--------------------------------------------------------------------------
  ! ens_getOneLev_r8
  !--------------------------------------------------------------------------
  function ens_getOneLev_r8(ens,kIndex) result(oneLevLevel)
    !
    !:Purpose: Return a 4D pointer to a single level of a real8 ensemble. The
    !          dimensions of the pointer are:
    !          (memberIndex, stepIndex, lonIndex, latIndex)
    !
    implicit none

    ! Arguments:
    real(8), pointer :: oneLevLevel(:,:,:,:)
    type(struct_ens) :: ens
    integer          :: kIndex

    ! Locals:
    integer          :: lon1, lat1

    lon1 = ens%statevector_work%myLonBeg
    lat1 = ens%statevector_work%myLatBeg

    oneLevLevel(1:,1:,lon1:,lat1:) => ens%allLev_r8(kIndex)%onelevel(:,:,:,:)

  end function ens_getOneLev_r8

  !--------------------------------------------------------------------------
  ! ens_getOneLevMean_r8
  !--------------------------------------------------------------------------
  function ens_getOneLevMean_r8(ens,subEnsIndex,kIndex) result(field)
    !
    !:Purpose: Return a 3D pointer to a single level the ensemble mean. The
    !          dimensions of the pointer are:
    !          (stepIndex, lonIndex, latIndex)
    !
    implicit none

    ! Arguments:
    real(8), pointer  :: field(:,:,:)
    type(struct_ens)  :: ens
    integer           :: subEnsIndex, kIndex

    ! Locals:
    integer           :: lon1, lat1

    lon1 = ens%statevector_work%myLonBeg
    lat1 = ens%statevector_work%myLatBeg

    field(1:, lon1:, lat1:) => ens%allLev_ensMean_r8(kIndex)%onelevel(subEnsIndex,:,:,:)

  end function ens_getOneLevMean_r8

  !--------------------------------------------------------------------------
  ! ens_copyEnsMean
  !--------------------------------------------------------------------------
  subroutine ens_copyEnsMean(ens, statevector, subEnsIndex_opt)
    !
    !:Purpose: Copy the ensemble mean into the supplied stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer, optional :: subEnsIndex_opt

    ! Locals:
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep, subEnsIndex
    character(len=4), pointer :: varNamesInEns(:)

    if( present(subEnsIndex_opt) ) then
      subEnsIndex = subEnsIndex_opt
    else
      subEnsIndex = 1
    end if

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. statevector%allocated) then
      nullify(varNamesInEns)
      call gsv_varNamesList(varNamesInEns,ens%statevector_work)
      call gsv_allocate(statevector, numStep,  &
                        ens%statevector_work%hco, ens%statevector_work%vco,  &
                        varNames_opt=varNamesInEns, datestamp_opt=tim_getDatestamp(), &
                        mpi_local_opt=.true., dataKind_opt=8 )
      deallocate(varNamesInEns)
    end if

    if (statevector%dataKind == 8) then
      ptr4d_r8 => gsv_getField_r8(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:)
        end do
      end do
    else
      ptr4d_r4 => gsv_getField_r4(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ptr4d_r4(:,:,jk,stepIndex) = real(ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:),4)
        end do
      end do
    end if

  end subroutine ens_copyEnsMean

  !--------------------------------------------------------------------------
  ! ens_copyToEnsMean
  !--------------------------------------------------------------------------
  subroutine ens_copyToEnsMean(ens, statevector, subEnsIndex_opt)
    !
    !:Purpose: Copy the supplied stateVector into the ensemble mean.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer, optional :: subEnsIndex_opt

    ! Locals:
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep, subEnsIndex
    character(len=4), pointer :: varNamesInEns(:)

    if( present(subEnsIndex_opt) ) then
      subEnsIndex = subEnsIndex_opt
    else
      subEnsIndex = 1
    end if

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. statevector%allocated) then
      call utl_abort('ens_copyToEnsMean: supplied stateVector must be allocated')
    end if

    if (statevector%dataKind == 8) then
      ptr4d_r8 => gsv_getField_r8(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = ptr4d_r8(:,:,jk,stepIndex)
        end do
      end do
    else
      ptr4d_r4 => gsv_getField_r4(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:) = real(ptr4d_r4(:,:,jk,stepIndex),8)
        end do
      end do
    end if

  end subroutine ens_copyToEnsMean

  !--------------------------------------------------------------------------
  ! ens_copyEnsStdDev
  !--------------------------------------------------------------------------
  subroutine ens_copyEnsStdDev(ens, statevector)
    !
    !:Purpose: Copy the ensemble StdDev into the supplied stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector

    ! Locals:
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep
    character(len=4), pointer :: varNamesInEns(:)

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. statevector%allocated) then
      nullify(varNamesInEns)
      call gsv_varNamesList(varNamesInEns,ens%statevector_work)
      call gsv_allocate(statevector, numStep,  &
                        ens%statevector_work%hco, ens%statevector_work%vco,  &
                        varNames_opt=varNamesInEns, datestamp_opt=tim_getDatestamp(), &
                        mpi_local_opt=.true., dataKind_opt=8 )
      deallocate(varNamesInEns)
    end if

    if (statevector%dataKind == 8) then
      ptr4d_r8 => gsv_getField_r8(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,:,:)
        end do
      end do
    else
      ptr4d_r4 => gsv_getField_r4(statevector)
      do stepIndex = 1, numStep
        do jk = k1, k2
          ptr4d_r4(:,:,jk,stepIndex) = real(ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,:,:),4)
        end do
      end do
    end if

  end subroutine ens_copyEnsStdDev

  !--------------------------------------------------------------------------
  ! ens_copyMember
  !--------------------------------------------------------------------------
  subroutine ens_copyMember(ens, statevector, memberIndex)
    !
    !:Purpose: Copy a selected ensemble member into the supplied stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer           :: memberIndex

    ! Locals:
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep, varIndex
    integer          :: gsvLevIndex, ensVarLevIndex, nLev
    character(len=4), pointer :: varNamesInEns(:)
    character(len=4), pointer :: varNamesInGsv(:)
    character(len=4) :: varName
    logical          :: sameVariables

    numStep = ens%statevector_work%numStep

    nullify(varNamesInEns)
    call gsv_varNamesList(varNamesInEns, ens%statevector_work)

    if (.not. statevector%allocated) then
      call gsv_allocate( statevector, numStep,  &
                         ens%statevector_work%hco, ens%statevector_work%vco,  &
                         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                         varNames_opt=varNamesInEns, dataKind_opt=8)
      varNamesInGsv => varNamesInEns
    else
      nullify(varNamesInGsv)
      call gsv_varNamesList(varNamesInGsv, statevector)
    end if

    sameVariables = .false.
    if (size(ens%statevector_work%varExistlist) == size(statevector%varExistlist)) then
      if (all(ens%statevector_work%varExistlist == statevector%varExistlist)) then
        sameVariables = .true.
      end if
    end if

    if (sameVariables) then

      k1 = ens%statevector_work%mykBeg
      k2 = ens%statevector_work%mykEnd

      if (ens%dataKind == 8) then
        ptr4d_r8 => gsv_getField_r8(statevector)
        do stepIndex = 1, numStep
          do jk = k1, k2
            ptr4d_r8(:,:,jk,stepIndex) = ens%allLev_r8(jk)%onelevel(memberIndex,stepIndex,:,:)
          end do
        end do
      else if (ens%dataKind == 4) then
        if (gsv_getDataKind(statevector) == 8) then
          ptr4d_r8 => gsv_getField_r8(statevector)
          do stepIndex = 1, numStep
            do jk = k1, k2
              ptr4d_r8(:,:,jk,stepIndex) = real(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:),8)
            end do
          end do
        else
          ptr4d_r4 => gsv_getField_r4(statevector)
          do stepIndex = 1, numStep
            do jk = k1, k2
              ptr4d_r4(:,:,jk,stepIndex) = ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:)
            end do
          end do
        end if
      end if

    else

      do varIndex = 1, size(varNamesInGsv)
        varName = varNamesInGsv(varIndex)
        if (.not. ens_varExist(ens,varName)) cycle
        nLev = gsv_getNumLev(statevector,vnl_varLevelFromVarname(varName))
        if (ens%dataKind == 8) then
          ptr4d_r8 => gsv_getField_r8(statevector,varName_opt=varName)
          do stepIndex = 1, numStep
            do gsvLevIndex = 1, nLev
              ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
              ptr4d_r8(:,:,gsvLevIndex,stepIndex) = ens%allLev_r8(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:)
            end do
          end do
        else if (ens%dataKind == 4) then
          if (gsv_getDataKind(statevector) == 8) then
            ptr4d_r8 => gsv_getField_r8(statevector,varName_opt=varName)
            do stepIndex = 1, numStep
              do gsvLevIndex = 1, nLev
                ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
                ptr4d_r8(:,:,gsvLevIndex,stepIndex) = real(ens%allLev_r4(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:),8)
              end do
            end do
          else
            ptr4d_r4 => gsv_getField_r4(statevector,varName_opt=varName)
            do stepIndex = 1, numStep
              do gsvLevIndex = 1, nLev
                ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
                ptr4d_r4(:,:,gsvLevIndex,stepIndex) = ens%allLev_r4(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:)
              end do
            end do
          end if
        end if
      end do

    end if

    if (associated(varNamesInGsv)) deallocate(varNamesInGsv)
    if (associated(varNamesInEns)) deallocate(varNamesInEns)

  end subroutine ens_copyMember

  !--------------------------------------------------------------------------
  ! ens_insertMember
  !--------------------------------------------------------------------------
  subroutine ens_insertMember(ens, statevector, memberIndex)
    !
    !:Purpose: Copy the supplied stateVector in the selected ensemble member.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer           :: memberIndex

    ! Locals:
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep, varIndex
    integer          :: gsvLevIndex, ensVarLevIndex, nLev
    character(len=4), pointer :: varNamesInEns(:)
    character(len=4), pointer :: varNamesInGsv(:)
    character(len=4) :: varName
    logical          :: sameVariables

    if (.not. ens%statevector_work%allocated) then
      call utl_abort('ens_insertMember: ens not allocated')
    end if

    numStep = ens%statevector_work%numStep

    nullify(varNamesInEns)
    call gsv_varNamesList(varNamesInEns, ens%statevector_work)
    nullify(varNamesInGsv)
    call gsv_varNamesList(varNamesInGsv, statevector)

    sameVariables = .false.
    if (size(ens%statevector_work%varExistlist) == size(statevector%varExistlist)) then
      if (all(ens%statevector_work%varExistlist == statevector%varExistlist)) then
        sameVariables = .true.
      end if
    end if

    if (sameVariables) then

      k1 = ens%statevector_work%mykBeg
      k2 = ens%statevector_work%mykEnd

      if (ens%dataKind == 8) then
        ptr4d_r8 => gsv_getField_r8(statevector)
        do stepIndex = 1, numStep
          do jk = k1, k2
            ens%allLev_r8(jk)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r8(:,:,jk,stepIndex)
          end do
        end do
      else if (ens%dataKind == 4) then
        if (gsv_getDataKind(statevector) == 8) then
          ptr4d_r8 => gsv_getField_r8(statevector)
          do stepIndex = 1, numStep
            do jk = k1, k2
              ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:) = real(ptr4d_r8(:,:,jk,stepIndex),4)
            end do
          end do
        else
          ptr4d_r4 => gsv_getField_r4(statevector)
          do stepIndex = 1, numStep
            do jk = k1, k2
              ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r4(:,:,jk,stepIndex)
            end do
          end do
        end if
      end if

    else

      do varIndex = 1, size(varNamesInGsv)
        varName = varNamesInGsv(varIndex)
        nLev = gsv_getNumLev(statevector,vnl_varLevelFromVarname(varName))
        if (ens%dataKind == 8) then
          ptr4d_r8 => gsv_getField_r8(statevector,varName_opt=varName)
          do stepIndex = 1, numStep
            do gsvLevIndex = 1, nLev
              ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
              ens%allLev_r8(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r8(:,:,gsvLevIndex,stepIndex)
            end do
          end do
        else if (ens%dataKind == 4) then
          if (gsv_getDataKind(statevector) == 8) then
            ptr4d_r8 => gsv_getField_r8(statevector,varName_opt=varName)
            do stepIndex = 1, numStep
              do gsvLevIndex = 1, nLev
                ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
                ens%allLev_r4(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:) = real(ptr4d_r8(:,:,gsvLevIndex,stepIndex),4)
              end do
            end do
          else
            ptr4d_r4 => gsv_getField_r4(statevector,varName_opt=varName)
            do stepIndex = 1, numStep
              do gsvLevIndex = 1, nLev
                ensVarLevIndex = gsvLevIndex + ens_getOffsetFromVarName(ens,varName)
                ens%allLev_r4(ensVarLevIndex)%onelevel(memberIndex,stepIndex,:,:) = ptr4d_r4(:,:,gsvLevIndex,stepIndex)
              end do
            end do
          end if
        end if
      end do ! varIndex

    end if

    if (associated(varNamesInGsv)) deallocate(varNamesInGsv)
    if (associated(varNamesInEns)) deallocate(varNamesInEns)

  end subroutine ens_insertMember

  !--------------------------------------------------------------------------
  ! ens_varExist
  !--------------------------------------------------------------------------
  function ens_varExist(ens,varName) result(varExist)
    !
    !:Purpose: Return true if the specified variable name exists in the ensemble.
    !
    implicit none
    logical                      :: varExist 

    ! Arguments:
    type(struct_ens)             :: ens
    character(len=*), intent(in) :: varName

    varExist = gsv_varExist(ens%statevector_work, varName)

  end function ens_varExist

  !--------------------------------------------------------------------------
  ! ens_varNamesList
  !--------------------------------------------------------------------------
  subroutine ens_varNamesList(varNames,ens)
    !
    !:Purpose: Return a list of the variable names that exist in the ensemble.
    !
    implicit none
    
    ! Arguments:
    type(struct_ens), optional :: ens
    character(len=4), pointer  :: varNames(:)

    if (associated(varNames)) then
      call utl_abort('ens_varNamesList: varNames must be NULL pointer on input')
    end if

    if (present(ens)) then
      call gsv_varNamesList(varNames, ens%statevector_work)
    else
      call gsv_varNamesList(varNames)
    end if

  end subroutine ens_varNamesList

  !--------------------------------------------------------------------------
  ! ens_getNumLev
  !--------------------------------------------------------------------------
  function ens_getNumLev(ens,varLevel) result(nlev)
    !
    !:Purpose: Return the number of vertical levels of the ensemble.
    !
    implicit none
    integer                       :: nlev

    ! Arguments:
    type(struct_ens), intent(in)  :: ens
    character(len=*), intent(in)  :: varLevel

    nlev = vco_getNumLev(ens%statevector_work%vco,varLevel)

  end function ens_getNumLev
  
  !--------------------------------------------------------------------------
  ! ens_getNumMembers
  !--------------------------------------------------------------------------
  function ens_getNumMembers(ens) result(numMembers)
    !
    !:Purpose: Return the number of members in the ensemble.
    !
    implicit none
    integer                       :: numMembers

    ! Arguments:
    type(struct_ens), intent(in)  :: ens

    numMembers = ens%numMembers

  end function ens_getNumMembers


  !--------------------------------------------------------------------------
  ! ens_getNumSubEns
  !--------------------------------------------------------------------------
  function ens_getNumSubEns(ens) result(numMembers)
    !
    !:Purpose: Return the number of sub-ensembles in the ensemble.
    !
    implicit none
    integer                       :: numMembers

    ! Arguments:
    type(struct_ens), intent(in)  :: ens

    numMembers = ens%numSubEns

  end function ens_getNumSubEns

  !--------------------------------------------------------------------------
  ! ens_getNumK
  !--------------------------------------------------------------------------
  function ens_getNumK(ens) result(numK)
    !
    !:Purpose: Return the number of kIndex (a.k.a. varLevs) values of the ensemble.
    !
    implicit none
    integer                       :: numK

    ! Arguments:
    type(struct_ens), intent(in)  :: ens

    numK = 1 + ens%statevector_work%mykEnd - ens%statevector_work%mykBeg

  end function ens_getNumK

  !--------------------------------------------------------------------------
  ! ens_getDataKind
  !--------------------------------------------------------------------------
  function ens_getDataKind(ens) result(dataKind)
    !
    !:Purpose: Return the floating point kind of the ensemble (4 or 8).
    !
    implicit none
    integer                       :: dataKind

    ! Arguments:
    type(struct_ens), intent(in)  :: ens

    dataKind = ens%dataKind

  end function ens_getDataKind

  !--------------------------------------------------------------------------
  ! ens_getOffsetFromVarName
  !--------------------------------------------------------------------------
  function ens_getOffsetFromVarName(ens,varName) result(offset)
    !
    !:Purpose: Return the offset of the kIndex for the specified variable name.
    !
    implicit none
    integer                      :: offset

    ! Arguments:
    type(struct_ens)             :: ens
    character(len=*), intent(in) :: varName

    if (.not. ens_varExist(ens,varName)) then
      call utl_abort('ens_getOffsetFromVarName: this varName is not present in ens: '//trim(varName))
    end if

    offset=gsv_getOffsetFromVarName(ens%statevector_work,varName)

  end function ens_getOffsetFromVarName

  !--------------------------------------------------------------------------
  ! ens_getLevFromK
  !--------------------------------------------------------------------------
  function ens_getLevFromK(ens,kIndex) result(levIndex)
    !
    !:Purpose: Return the level index from the kIndex value.
    !
    implicit none
    integer                      :: levIndex

    ! Arguments:
    type(struct_ens), intent(in) :: ens
    integer, intent(in)          :: kIndex

    levIndex = gsv_getLevFromK(ens%statevector_work,kIndex)

  end function ens_getLevFromK

  !--------------------------------------------------------------------------
  ! ens_getKFromLevVarName
  !--------------------------------------------------------------------------
  function ens_getKFromLevVarName(ens, levIndex, varName) result(kIndex)
    !
    !:Purpose: Return the kIndex value for the specified level index
    !          and variable name.
    !
    implicit none
    integer                      :: kIndex

    ! Arguments:
    type(struct_ens), intent(in) :: ens
    integer                      :: levIndex
    character(len=*)             :: varName

    kIndex = levIndex + gsv_getOffsetFromVarName(ens%statevector_work,trim(varName))

  end function ens_getKFromLevVarName

  !--------------------------------------------------------------------------
  ! ens_getVarNameFromK
  !--------------------------------------------------------------------------
  function ens_getVarNameFromK(ens,kIndex) result(varName)
    !
    !:Purpose: Return the variable name from the specified kIndex value.
    !
    implicit none
    character(len=4)             :: varName

    ! Arguments:
    type(struct_ens), intent(in) :: ens
    integer, intent(in)          :: kIndex

    varName = gsv_getVarNameFromK(ens%statevector_work,kIndex)

  end function ens_getVarNameFromK

  !--------------------------------------------------------------------------
  ! ens_getVco
  !--------------------------------------------------------------------------
  function ens_getVco(ens) result(vco_ptr)
    !
    !:Purpose: Return a pointer to the verticalCoord object associate with the ensemble.
    !
    implicit none
    type(struct_vco), pointer :: vco_ptr

    ! Arguments:
    type(struct_ens)          :: ens

    vco_ptr => ens%statevector_work%vco

  end function ens_getVco

  !--------------------------------------------------------------------------
  ! ens_getHco
  !--------------------------------------------------------------------------
  function ens_getHco(ens) result(hco_ptr)
    !
    !:Purpose: Return a pointer to the horizontalCoord object associate with the ensemble.
    !
    implicit none
    type(struct_hco), pointer :: hco_ptr

    ! Arguments:
    type(struct_ens)          :: ens

    hco_ptr => ens%statevector_work%hco

  end function ens_getHco

  !--------------------------------------------------------------------------
  ! ens_getLatLonBounds
  !--------------------------------------------------------------------------
  subroutine ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    !
    !:Purpose: Return the longitude and latitude index bounds for this mpi task.
    !
    implicit none

    ! Arguments:
    type(struct_ens)       :: ens
    integer, intent(out)   :: myLonBeg
    integer, intent(out)   :: myLonEnd
    integer, intent(out)   :: myLatBeg
    integer, intent(out)   :: myLatEnd

    myLonBeg = ens%statevector_work%myLonBeg
    myLonEnd = ens%statevector_work%myLonEnd
    myLatBeg = ens%statevector_work%myLatBeg
    myLatEnd = ens%statevector_work%myLatEnd

  end subroutine ens_getLatLonBounds

  !--------------------------------------------------------------------------
  ! ens_getNumStep
  !--------------------------------------------------------------------------
  function ens_getNumStep(ens) result(numStep)
    !
    !:Purpose: Return the number of time steps stored in the ensemble.
    !
    implicit none
    integer :: numStep

    ! Arguments:
    type(struct_ens) :: ens

    numStep = ens%statevector_work%numStep

  end function ens_getNumStep

  !--------------------------------------------------------------------------
  ! ens_computeMean
  !--------------------------------------------------------------------------
  subroutine ens_computeMean(ens, computeSubEnsMeans_opt, numSubEns_opt)
    !
    !:Purpose: Internally compute the ensemble mean.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    logical, optional :: computeSubEnsMeans_opt
    integer, optional :: numSubEns_opt

    ! Locals:
    logical           :: computeSubEnsMeans, lExists

    character(len=256), parameter :: subEnsIndexFileName = 'subEnsembleIndex.txt'

    integer           :: kulin, ierr, memberIndex, memberIndex2, stepIndex, subEnsIndex
    integer           :: k1, k2, jk, lon1, lon2, lat1, lat2, numStep, ji, jj
    integer           :: fnom, fclos

    if (present(computeSubEnsMeans_opt)) then
      computeSubEnsMeans = computeSubEnsMeans_opt
    else
      computeSubEnsMeans = .false.
    end if

    ! Read sub-ensemble index list from file, if it exists
    if (.not. allocated(ens%subEnsIndexList)) then
      allocate(ens%subEnsIndexList(ens%numMembers))
    end if
    if ( computeSubEnsMeans ) then
      write(*,*) 'ens_computeMean: checking in ensemble directory if file with sub-ensemble index list exists: ',subEnsIndexFileName
      inquire(file=trim(ens%enspathname) // trim(subEnsIndexFileName),exist=lExists)
      if ( lExists ) then
        kulin = 0
        ierr = fnom(kulin,trim(ens%enspathname) // trim(subEnsIndexFileName),'FMT+SEQ+R/O',0)
        do memberIndex = 1, ens%numMembers
          read(kulin,*) memberIndex2, ens%subEnsIndexList(memberIndex)
          write(*,*) 'read from sub-ensemble index list: ',memberIndex, memberIndex2, ens%subEnsIndexList(memberIndex)
        end do
        ierr   = fclos(kulin)
      else
        call utl_abort('ens_computeMean: could not find file with sub-ensemble index list')
      end if
    else
      ens%subEnsIndexList(:) = 1
    end if
    ens%numSubEns = maxval(ens%subEnsIndexList(:))
    if (.not. allocated(ens%nEnsSubEns)) then
      allocate(ens%nEnsSubEns(ens%numSubEns))
    end if
    ens%nEnsSubEns(:) = 0
    do memberIndex = 1, ens%numMembers
      ens%nEnsSubEns(ens%subEnsIndexList(memberIndex)) = ens%nEnsSubEns(ens%subEnsIndexList(memberIndex)) + 1
    end do
    write(*,*) 'ens_computeMean: number of sub-ensembles = ', ens%numSubEns
    write(*,*) 'ens_computeMean: number of members in each sub-ensemble = ', ens%nensSubEns(:)

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. allocated(ens%allLev_ensMean_r8)) then
      call ens_allocateMean(ens)
    else
      do jk = k1, k2
        ens%allLev_ensMean_r8(jk)%onelevel(:,:,:,:) = 0.0d0
      end do
    end if
    ens%meanIsComputed = .true.

    ! Compute ensemble mean(s)
    !$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex,subEnsIndex)
    do jk = k1, k2
      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, ens%statevector_work%numStep
            do memberIndex = 1, ens%numMembers
              ens%allLev_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj) = &
                   ens%allLev_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj) + &
                   dble(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj))
            end do
            do subEnsIndex = 1, ens%numSubEns
              ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,ji,jj) = &
                   ens%allLev_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,ji,jj) /  &
                   dble(ens%nEnsSubEns(subEnsIndex))
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! provide output argument value
    if ( present(numSubEns_opt) ) numSubEns_opt = ens%numSubEns

  end subroutine ens_computeMean

  !--------------------------------------------------------------------------
  ! ens_computeStdDev
  !--------------------------------------------------------------------------
  subroutine ens_computeStdDev(ens, containsScaledPerts_opt)
    !
    !:Purpose: Internally compute the ensemble stdDev.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    logical, intent(in), optional :: containsScaledPerts_opt

    ! Locals:
    integer           :: memberIndex, stepIndex, subEnsIndex
    integer           :: k1, k2, jk, lon1, lon2, lat1, lat2, numStep, ji, jj
    real(8), allocatable :: subEnsStdDev(:)
    logical           :: containsScaledPerts

    if ( present(containsScaledPerts_opt) ) then
      containsScaledPerts = containsScaledPerts_opt
    else
      containsScaledPerts = .false.
      if (.not.ens%meanIsRemoved .and. .not.ens%meanIsComputed) then
        if (mpi_myid == 0) write(*,*) 'ens_computeStdDev: compute Mean since it was not already done'
        call ens_computeMean( ens )
      end if
    end if

    ! Read sub-ensemble index list from file, if it exists
    ! The sub-ensembles should have been already read in routine 'ens_computeMean'
    write(*,*) 'ens_computeStdDev: number of sub-ensembles = ', ens%numSubEns
    write(*,*) 'ens_computeStdDev: number of members in each sub-ensemble = ', ens%nensSubEns(:)

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. allocated(ens%allLev_ensStdDev_r8)) then
      call ens_allocateStdDev(ens)
    else
      do jk = k1, k2
        ens%allLev_ensStdDev_r8(jk)%onelevel(:,:,:,:) = 0.0d0
      end do
    end if
    ens%StdDevIsComputed = .true.

    if (containsScaledPerts) then

      if (ens%numSubEns /= 1) then
        call utl_abort('ens_computeStdDev: sub-ensemble approach not compatible with scale perturbations')
      end if

      ! Compute the ensemble StdDev from previously scale ensemble perturbations 
      !  (i.e. pert = (fcst-mean)/(nEns-1) )

      !$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex)
      do jk = k1, k2
        do jj = lat1, lat2
          do ji = lon1, lon2
            do stepIndex = 1, ens%statevector_work%numStep
              
              ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) = 0.d0

              do memberIndex = 1, ens%numMembers
                ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) =      &
                     ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) + &
                     dble(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj))**2
              end do

              ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) = &
                   sqrt(ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj))

            end do
          end do
        end do
      end do
    !$OMP END PARALLEL DO

    else

      allocate(subEnsStdDev(ens%numSubEns))

      ! Compute global ensemble StdDev as the sqrt of the mean of each suchensemble variance
      !   var_subens(i) = sum( ( ens(j) - mean_subens(i) )**2, j=1..numEns_subens(i) ) / ( numEns_subens(i) - 1 )
      !   var_allensensemble = sum( numEns_subens(i) * var_subens(i), i=1..numSubEns)
      !   stddev = sqrt( var_allensensemble / numEnsTotal )

      !$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex,subEnsIndex,subEnsStdDev)
      do jk = k1, k2
        do jj = lat1, lat2
          do ji = lon1, lon2
            do stepIndex = 1, ens%statevector_work%numStep

              subEnsStdDev(:) = 0.0d0

              if (ens%meanIsRemoved) then
                do memberIndex = 1, ens%numMembers
                  subEnsStdDev(ens%subEnsIndexList(memberIndex)) =                      &
                       subEnsStdDev(ens%subEnsIndexList(memberIndex)) +                 &
                       (dble(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj)))**2
                end do
              else
                do memberIndex = 1, ens%numMembers
                  subEnsStdDev(ens%subEnsIndexList(memberIndex)) =                      &
                       subEnsStdDev(ens%subEnsIndexList(memberIndex)) +                 &
                       (dble(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj)) - &
                       ens%allLev_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj))**2
                end do
              end if

              do subEnsIndex = 1, ens%numSubEns
                ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) =      &
                     ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) + &
                     ens%nEnsSubEns(subEnsIndex)*subEnsStdDev(subEnsIndex)/(ens%nEnsSubEns(subEnsIndex)-1)
                     
              end do

              ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) =        &
                   sqrt( ens%allLev_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) / dble(ens%numMembers) )
              
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      deallocate(subEnsStdDev)

    end if

  end subroutine ens_computeStdDev

  !--------------------------------------------------------------------------
  ! ens_normalize
  !--------------------------------------------------------------------------
  subroutine ens_normalize(ens)
    !
    !:Purpose: Normalize the ensemble by the 3D ensemble stdDev.
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens

    ! Locals:
    integer :: lon1, lon2, lat1, lat2, k1, k2, numStep
    integer :: jk, jj, ji, stepIndex, memberIndex

    real(8) :: factor

    if (.not. ens%StdDevIsComputed) then
      if (mpi_myid == 0) write(*,*) 'ens_normalize: compute Std. Dev. since it was not already done'
      call ens_computeStdDev( ens )
    end if

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    !$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex,factor)
    do jk = k1, k2
      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, numStep
            do memberIndex = 1, ens%numMembers

              if (ens%allLev_ensStdDev_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj) > 0.0d0 ) then
                factor = 1.0d0/ens%allLev_ensStdDev_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj)
              else
                factor = 0.0d0
              endif

              ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj) =  &
                   real( real(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj),8) * factor, 4)
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine ens_normalize

  !--------------------------------------------------------------------------
  ! ens_removeMean
  !--------------------------------------------------------------------------
  subroutine ens_removeMean(ens)
    !
    !:Purpose: Subtract the ensemble mean from each member.
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens

    ! Locals:
    integer :: lon1, lon2, lat1, lat2, k1, k2, numStep
    integer :: jk, jj, ji, stepIndex, memberIndex

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    !$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex)
    do jk = k1, k2
      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, numStep
            do memberIndex = 1, ens%numMembers
              ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj) =  &
                   real( (real(ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj),8) -  &
                   ens%allLev_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj)), 4 )
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ens%meanIsRemoved = .true.

  end subroutine ens_removeMean

  !--------------------------------------------------------------------------
  ! ens_recenter
  !--------------------------------------------------------------------------
  subroutine ens_recenter(ens, recenteringMean, recenteringCoeff_opt,  &
                          recenteringCoeffArray_opt, &
                          alternativeEnsembleMean_opt, &
                          ensembleControlMember_opt, scaleFactor_opt, &
                          numMembersToRecenter_opt)
    !
    !:Purpose:
    !          To compute:
    !          ..math::
    !              x_recentered =
    !                      scaleFactor*x_original
    !                    + recenteringCoeff*(  x_recenteringMean
    !                                        - scaleFactor*x_ensembleMean
    !                                       )
    implicit none

    ! Arguments:
    type(struct_ens) :: ens
    type(struct_gsv) :: recenteringMean
    type(struct_gsv), optional :: alternativeEnsembleMean_opt, ensembleControlMember_opt
    real(8), optional :: recenteringCoeff_opt
    real(8), optional :: recenteringCoeffArray_opt(:)
    real(8), optional :: scaleFactor_opt(:)
    integer, optional :: numMembersToRecenter_opt

    ! Locals:
    integer,parameter    :: maxNumLevels=200
    real(8), pointer :: ptr4d_r8(:,:,:,:), alternativeEnsembleMean_r8(:,:,:,:), ptr4d_ensembleControlmember_r8(:,:,:,:)
    real(8) :: increment, scaleFactor(maxNumLevels), thisScaleFactor
    real(8) :: recenteringCoeffArray(ens%numMembers)
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer :: lon1, lon2, lat1, lat2, k1, k2, numStep, numMembersToRecenter
    integer :: jk, jj, ji, stepIndex, memberIndex, levIndex
    character(len=4) :: varLevel

    if ( present(recenteringCoeff_opt) ) then
      recenteringCoeffArray(:) = recenteringCoeff_opt
    else if ( present(recenteringCoeffArray_opt) ) then
      recenteringCoeffArray(:) = recenteringCoeffArray_opt(1:ens%numMembers)
    else
      call utl_abort('ens_recenter: Must specify recenteringCoeff_opt or recenteringCoeffArray_opt')
    end if

    if ( present(scaleFactor_opt) ) then
      ! scaleFactor cannot be used at the same time as a recenteringCoeff different from 1.0
      if ( any (abs(recenteringCoeffArray(:)  - 1.0D0) > 1.0D-5) ) then
        call utl_abort('ens_recenter: recenteringCoeff must be equal to 1.0 when using scaleFactor')
      end if
      scaleFactor = scaleFactor_opt
    else
      scaleFactor(:) = 1.0D0
    end if

    if (present(numMembersToRecenter_opt)) then
      numMembersToRecenter = numMembersToRecenter_opt
    else
      numMembersToRecenter = ens%numMembers
    end if

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    nullify(ptr4d_r4, ptr4d_r8)
    if (gsv_getDataKind(recenteringMean) == 8) then
      ptr4d_r8 => gsv_getField_r8(recenteringMean)
    else
      ptr4d_r4 => gsv_getField_r4(recenteringMean)
    end if
    if(present(alternativeEnsembleMean_opt)) then
      alternativeEnsembleMean_r8 => gsv_getField_r8(alternativeEnsembleMean_opt)
    else
      nullify(alternativeEnsembleMean_r8)
    end if

    if (present(ensembleControlMember_opt)) then
      ptr4d_ensembleControlmember_r8 => gsv_getField_r8(ensembleControlMember_opt)
    else
      nullify(ptr4d_ensembleControlmember_r8)
    end if

    !$OMP PARALLEL DO PRIVATE (jk,varLevel,levIndex,thisScaleFactor,jj,ji,stepIndex,memberIndex,increment)
    do jk = k1, k2

      ! define scaling factor as a function of vertical level and variable type
      varLevel = vnl_varLevelFromVarname(ens_getVarNameFromK(ens, jk))
      if ( trim(varLevel) == 'SF') then
        ! use lowest momentum level for surface variables
        levIndex = ens_getNumLev(ens, 'MM')
      else if ( (trim(varLevel) == 'MM') .and. (ens%statevector_work%vco%Vcode == 5002) ) then
        levIndex = ens_getLevFromK(ens, jk) + 1
      else
        levIndex = ens_getLevFromK(ens, jk)
      end if
      thisScaleFactor = scaleFactor(levIndex)

      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, numStep
            if(present(alternativeEnsembleMean_opt)) then
              if (associated(ptr4d_r8)) then
                increment = ptr4d_r8(ji,jj,jk,stepIndex) -  &
                            thisScaleFactor*alternativeEnsembleMean_r8(ji,jj,jk,stepIndex)
              else
                increment = real(ptr4d_r4(ji,jj,jk,stepIndex),8) -  &
                            thisScaleFactor*alternativeEnsembleMean_r8(ji,jj,jk,stepIndex)
              end if
            else
              if (associated(ptr4d_r8)) then
                increment = ptr4d_r8(ji,jj,jk,stepIndex) -  &
                            thisScaleFactor*ens%allLev_ensMean_r8(jk)%onelevel(1,stepIndex,ji,jj)
              else
                increment = real(ptr4d_r4(ji,jj,jk,stepIndex),8) -  &
                            thisScaleFactor*ens%allLev_ensMean_r8(jk)%onelevel(1,stepIndex,ji,jj)
              end if
            end if
            if (present(ensembleControlMember_opt)) then
              ptr4d_ensembleControlMember_r8(ji,jj,jk,stepIndex) =  &
                   thisScaleFactor*ptr4d_ensembleControlMember_r8(ji,jj,jk,stepIndex) +  &
                   recenteringCoeffArray(1)*increment
            else
              do memberIndex = 1, numMembersToRecenter
                ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj) =  &
                     real( real(thisScaleFactor*ens%allLev_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj),8) +  &
                     recenteringCoeffArray(memberIndex)*increment, 4)
              end do
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine ens_recenter

  !--------------------------------------------------------------------------
  ! ens_readEnsemble
  !--------------------------------------------------------------------------
  subroutine ens_readEnsemble(ens, ensPathName, biPeriodic, &
                              vco_file_opt, varNames_opt, checkModelTop_opt, &
                              containsFullField_opt)
    !
    !:Purpose: Read the ensemble from disk in parallel and do mpi communication
    !          so that all members for a given lat-lon tile are present on each
    !          mpi task.
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ens
    character(len=*) :: ensPathName
    logical          :: biPeriodic
    character(len=*), optional          :: varNames_opt(:)
    type(struct_vco), pointer, optional :: vco_file_opt
    logical, optional                   :: checkModelTop_opt
    logical, optional, intent(in)       :: containsFullField_opt

    ! Locals:
    type(struct_gsv) :: statevector_file_r4, statevector_hint_r4, statevector_member_r4
    type(struct_hco), pointer :: hco_file, hco_ens, hco_coregrid
    type(struct_vco), pointer :: vco_file, vco_ens
    real(4), allocatable :: gd_send_r4(:,:,:,:)
    real(4), allocatable :: gd_recv_r4(:,:,:,:)
    real(4), pointer     :: ptr3d_r4(:,:,:)
    integer,pointer :: dateStampList(:)
    integer :: batchnum, nsize, status, ierr
    integer :: yourid, youridx, youridy
    integer :: readFilePE(1000)
    integer :: memberIndexOffset, totalEnsembleSize
    integer :: length_envVariable
    integer :: lonPerPEmax, latPerPEmax, ni, nj, nk, numStep, numlevelstosend, numlevelstosend2
    integer :: memberIndex, memberIndex2, fileMemberIndex, stepIndex, jk, jk2, jk3
    character(len=256) :: ensFileName
    character(len=32)  :: envVariable
    character(len=2)   :: typvar
    character(len=12)  :: etiket
    character(len=4), pointer :: anlVar(:)
    logical :: verticalInterpNeeded, horizontalInterpNeeded, horizontalPaddingNeeded
    logical :: checkModelTop
    logical :: containsFullField
    character(len=4), pointer :: varNames(:)

    write(*,*) 'ens_readEnsemble: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. ens%allocated ) then
      call utl_abort('ens_readEnsemble: ensemble object not allocated!')
    end if

    !
    !- 1. Initial setup
    !
    lonPerPEmax = ens%statevector_work%lonPerPEmax
    latPerPEmax = ens%statevector_work%latPerPEmax
    ni          = ens%statevector_work%ni
    nj          = ens%statevector_work%nj
    nk          = ens%statevector_work%nk
    numStep     = ens%statevector_work%numStep
    dateStampList => ens%statevector_work%dateStampList

    ens%ensPathName = trim(ensPathName)

    ! Memory allocation
    numLevelsToSend = 10
    allocate(gd_send_r4(lonPerPEmax,latPerPEmax,numLevelsToSend,mpi_nprocs))
    allocate(gd_recv_r4(lonPerPEmax,latPerPEmax,numLevelsToSend,mpi_nprocs))
    gd_send_r4(:,:,:,:) = 0.0
    gd_recv_r4(:,:,:,:) = 0.0

    do memberIndex = 1, ens%numMembers
      readFilePE(memberIndex) = mod(memberIndex-1,mpi_nprocs)
    end do

    ! the default is to check whether output grid has a higher top than input grid during vertical interpolation
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if

    ! Retrieve environment variables related to doing an ensemble of perturbed analyses
    status = 0
    call get_environment_variable('envar_memberIndexOffset',envVariable,length_envVariable,status,.true.)
    if (status.gt.1) then
      write(*,*) 'ens_readEnsemble: Problem when getting the environment variable envar_memberIndexOffset'
      memberIndexOffset = 0
    else if (status == 1) then
      memberIndexOffset = 0
    else
      write(*,*) 'ens_readEnsemble: The environment variable envar_memberIndexOffset has been detected: ',envVariable
      read(envVariable,'(i8)') memberIndexOffset
      write(*,*) 'memberIndexOffset = ',memberIndexOffset
    end if

    status = 0
    call get_environment_variable('envar_totalEnsembleSize',envVariable,length_envVariable,status,.true.)
    if (status.gt.1) then
      write(*,*) 'ens_readEnsemble: Problem when getting the environment variable envar_totalEnsembleSize'
      totalEnsembleSize = ens%numMembers
    else if (status == 1) then
      totalEnsembleSize = ens%numMembers
    else
      write(*,*) 'ens_readEnsemble: The environment variable envar_totalEnsembleSize has been detected: ',envVariable
      read(envVariable,'(i8)') totalEnsembleSize
      write(*,*) 'totalEnsembleSize = ',totalEnsembleSize
    end if

    ! Set up hco and vco for ensemble files
    call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=1, copyToRamDisk_opt=.false.)

    nullify(anlVar)
    call gsv_varNamesList(anlVar)
    nullify(hco_file)
    call hco_SetupFromFile(hco_file, ensFileName, ' ', 'ENSFILEGRID', varName_opt=anlVar(1))
    if ( present(vco_file_opt) ) then
      ! use the input vertical grid provided
      vco_file => vco_file_opt
    else
      ! find the info from the ensemble files
      nullify(vco_file)
      if ( mpi_myid == 0 ) then
        call vco_SetupFromFile(vco_file, ensFileName)
      end if
      call vco_mpiBcast(vco_file)
    end if
    hco_ens  => gsv_getHco(ens%statevector_work)
    vco_ens  => gsv_getVco(ens%statevector_work)
    horizontalInterpNeeded = (.not. hco_equal(hco_ens, hco_file))
    verticalInterpNeeded   = (.not. vco_equal(vco_ens, vco_file))

    ! More efficient handling of common case where input is on Z grid, analysis in on G grid
    if ( hco_file%grtyp == 'Z' .and. hco_ens%grtyp == 'G' ) then
      if ( hco_file%ni == (hco_ens%ni+1) ) then
        write(*,*) 'ens_readEnsemble: no interpolation done for equivalent Gaussian grid stored as a Z grid'
        horizontalInterpNeeded = .false.
      end if
    end if
    ! In limited-area mode, avoid horizontal interpolation when the ensemble is on the coregrid
    horizontalPaddingNeeded = .false.
    if ( .not. hco_file%global ) then
      hco_coregrid => agd_getHco('CoreGrid')
      if ( hco_file%ni == hco_coregrid%ni .and. hco_file%nj == hco_coregrid%nj  ) then
        if (mpi_myid == 0) write(*,*) 'ens_readEnsemble: no interpolation needed for ensemble on the limited-area coregrid'
        horizontalInterpNeeded = .false.
        horizontalPaddingNeeded = .true.
      end if
    end if

    if (mpi_myid == 0) then
      write(*,*)
      write(*,*) 'ens_readEnsemble: dateStampList=',dateStampList(1:numStep)
      write(*,*)
      if (horizontalInterpNeeded )  write(*,*) 'ens_readEnsemble: HORIZONTAL interpolation is needed'
      if (verticalInterpNeeded   )  write(*,*) 'ens_readEnsemble: VERTICAL   interpolation is needed'
      if (horizontalPaddingNeeded)  write(*,*) 'ens_readEnsemble: HORIZONTAL padding       is needed'
    end if

    ! Input type
    if (present(containsFullField_opt)) then
      containsFullField = containsFullField_opt
    else
      containsFullField = .true.
    end if
    if (mpi_myid == 0) then
      write(*,*)
      write(*,*) 'ens_readEnsemble: containsFullField = ', containsFullField
    end if

    !
    !- 2.  Ensemble forecasts reading loop
    !

    nullify(varNames)
    if (present(varNames_opt)) then
      allocate(varNames(size(varNames_opt)))
      varNames(:) = varNames_opt(:)
    else
      call gsv_varNamesList(varNames)
    end if

    !- 2.1 Loop on time, ensemble member, variable, level
    do stepIndex = 1, numStep
      write(*,*) ' '
      write(*,*) 'ens_readEnsemble: starting to read time level ', stepIndex

      ! allocate the needed statevector objects
      call gsv_allocate(statevector_member_r4, 1, hco_ens, vco_ens,  &
                        datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                        varNames_opt = varNames, dataKind_opt = 4,  &
                        hInterpolateDegree_opt = ens%hInterpolateDegree)
      if (horizontalInterpNeeded .or. verticalInterpNeeded .or. horizontalPaddingNeeded) then
        call gsv_allocate(statevector_file_r4, 1, hco_file, vco_file,  &
                          datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                          varNames_opt = varNames, dataKind_opt = 4,  &
                          hInterpolateDegree_opt = ens%hInterpolateDegree)
      end if
      if (verticalInterpNeeded) then
        call gsv_allocate(statevector_hint_r4, 1, hco_ens, vco_file,  &
                          datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                          varNames_opt = varNames, dataKind_opt = 4, &
                          hInterpolateDegree_opt = ens%hInterpolateDegree)
      end if
      
      do memberIndex = 1, ens%numMembers

        if (mpi_myid == readFilePE(memberIndex)) then

          write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

          !  Read the file
          fileMemberIndex = 1+mod(memberIndex+memberIndexOffset-1, totalEnsembleSize)
          call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=fileMemberIndex)
          typvar = ' '
          etiket = ' '
          if (.not. horizontalInterpNeeded  .and. &
              .not. verticalInterpNeeded    .and. &
              .not. horizontalPaddingNeeded ) then
            call gsv_readFile(statevector_member_r4, ensFileName, etiket, typvar, &
                              containsFullField)
          else
            call gsv_readFile(statevector_file_r4, ensFileName, etiket, typvar, &
                              containsFullField)
          end if
          if (stepIndex == numStep) then
            ierr = ram_remove(ensFileName)
          end if

          ! do any required interpolation
          if (horizontalInterpNeeded .and. verticalInterpNeeded) then
            call gsv_hInterpolate_r4(statevector_file_r4, statevector_hint_r4)
            call gsv_vInterpolate_r4(statevector_hint_r4, statevector_member_r4,         &
                                     Ps_in_hPa_opt=.true.,checkModelTop_opt=checkModelTop)

          else if (horizontalInterpNeeded .and. .not. verticalInterpNeeded) then
            call gsv_hInterpolate_r4(statevector_file_r4, statevector_member_r4)

          else if (.not. horizontalInterpNeeded .and. verticalInterpNeeded) then
            if (horizontalPaddingNeeded) then
              call gsv_hPad(statevector_file_r4, statevector_hint_r4)
            else
              call gsv_copy(statevector_file_r4, statevector_hint_r4)
            end if
            call gsv_vInterpolate_r4(statevector_hint_r4, statevector_member_r4,         &
                                     Ps_in_hPa_opt=.true.,checkModelTop_opt=checkModelTop)

          else if (horizontalPaddingNeeded) then
            call gsv_hPad(statevector_file_r4, statevector_member_r4)
          end if

          ! unit conversion
          call gsv_fileUnitsToStateUnits( statevector_member_r4, containsFullField)

          !  Create bi-periodic forecasts when using scale-dependent localization in LAM mode
          if ( .not. hco_ens%global .and. biperiodic ) then
            ptr3d_r4 => gsv_getField3D_r4(statevector_member_r4)
            call agd_mach_r4(ptr3d_r4,    & ! INOUT
                             ni, nj, statevector_member_r4%nk)  ! IN
          end if

        end if ! locally read one member


        !  MPI communication: from 1 ensemble member per process to 1 lat-lon tile per process  
        if (readFilePE(memberIndex) == (mpi_nprocs-1) .or. memberIndex == ens%numMembers) then

          call tmg_start(13,'PRE_SUENS_COMM')
          batchnum = ceiling(dble(memberIndex)/dble(mpi_nprocs))

          do jk = 1, nk, numLevelsToSend
            jk2 = min(nk,jk+numLevelsToSend-1)
            numLevelsToSend2 = jk2 - jk + 1

            ptr3d_r4 => gsv_getField3D_r4(statevector_member_r4)
            !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
            do youridy = 0, (mpi_npey-1)
              do youridx = 0, (mpi_npex-1)
                yourid = youridx + youridy*mpi_npex
                gd_send_r4(1:ens%statevector_work%allLonPerPE(youridx+1),  &
                           1:ens%statevector_work%allLatPerPE(youridy+1), 1:numLevelsToSend2, yourid+1) =  &
                           ptr3d_r4(ens%statevector_work%allLonBeg(youridx+1):ens%statevector_work%allLonEnd(youridx+1),  &
                           ens%statevector_work%allLatBeg(youridy+1):ens%statevector_work%allLatEnd(youridy+1), jk:jk2)
              end do
            end do
            !$OMP END PARALLEL DO

            nsize = lonPerPEmax * latPerPEmax * numLevelsToSend2
            if (mpi_nprocs.gt.1) then
              call rpn_comm_alltoall(gd_send_r4(:,:,1:numLevelsToSend2,:),nsize,"mpi_real4",  &
                                     gd_recv_r4(:,:,1:numLevelsToSend2,:),nsize,"mpi_real4","GRID",ierr)
            else
              gd_recv_r4(:,:,1:numLevelsToSend2,1) = gd_send_r4(:,:,1:numLevelsToSend2,1)
            end if

            call tmg_start(110,'ENS_TO_ONELEV')
            !$OMP PARALLEL DO PRIVATE(jk3,memberIndex2,yourid)
            do jk3 = 1, numLevelsToSend2
              do memberIndex2 = 1+(batchnum-1)*mpi_nprocs, memberIndex
                yourid = readFilePE(memberIndex2)
                ens%allLev_r4(jk3+jk-1)%onelevel(memberIndex2,stepIndex, :, :) =  &
                     gd_recv_r4(1:ens%statevector_work%lonPerPE, 1:ens%statevector_work%latPerPE, jk3, yourid+1)
              end do
            end do
            !$OMP END PARALLEL DO
            call tmg_stop(110)

          end do ! jk
          call tmg_stop(13)

        end if ! MPI communication


      end do ! memberIndex

      ! deallocate the needed statevector objects
      call gsv_deallocate(statevector_member_r4)
      if (horizontalInterpNeeded .or. verticalInterpNeeded) call gsv_deallocate(statevector_file_r4)
      if (verticalInterpNeeded) call gsv_deallocate(statevector_hint_r4)

    end do ! time

    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)
    deallocate(datestamplist)
    call hco_deallocate(hco_file)
    if ( .not. present(vco_file_opt) ) then
      call vco_deallocate(vco_file)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'ens_readEnsemble: finished reading and communicating ensemble members...'

  end subroutine ens_readEnsemble

  !--------------------------------------------------------------------------
  ! ens_writeEnsemble
  !--------------------------------------------------------------------------
  subroutine ens_writeEnsemble(ens, ensPathName, ensFileNamePrefix, &
                               ctrlVarHumidity, etiket, typvar, &
                               etiketAppendMemberNumber_opt, varNames_opt, &
                               ip3_opt, containsFullField_opt, numBits_opt)
    !
    !:Purpose: Write the ensemble to disk by doing mpi transpose so that
    !          each mpi task can write a single member in parallel.
    !
    implicit none

    ! Arguments:
    type(struct_ens)  :: ens
    character(len=*)  :: ensPathName
    character(len=*)  :: ensFileNamePrefix
    character(len=*)  :: ctrlVarHumidity
    character(len=*)  :: etiket
    character(len=*)  :: typvar
    character(len=*), optional :: varNames_opt(:)  ! allow specification of variables
    integer, optional :: ip3_opt, numBits_opt
    logical, optional :: etiketAppendMemberNumber_opt
    logical, optional :: containsFullField_opt

    ! Locals:
    type(struct_gsv) :: statevector_member_r4
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    real(4), allocatable :: gd_send_r4(:,:,:,:)
    real(4), allocatable :: gd_recv_r4(:,:,:,:)
    real(4), pointer     :: ptr3d_r4(:,:,:)
    integer, allocatable :: dateStampList(:)
    integer :: batchnum, nsize, ierr
    integer :: yourid, youridx, youridy
    integer :: writeFilePE(1000)
    integer :: lonPerPE, lonPerPEmax, latPerPE, latPerPEmax, ni, nj, nk, numStep, numlevelstosend, numlevelstosend2
    integer :: memberIndex, memberIndex2, stepIndex, jk, jk2, jk3, ip3, ensFileExtLength, maximumBaseEtiketLength
    character(len=256) :: ensFileName
    character(len=12) :: etiketStr  ! this is the etiket that will be used to write files
    !! The two next declarations are sufficient until we reach 10^10 members
    character(len=10) :: memberIndexStr ! this is the member number in a character string
    character(len=10) :: ensFileExtLengthStr ! this is a string containing the same number as 'ensFileExtLength'
    character(len=4), pointer :: varNamesInEns(:)
    logical :: containsFullField

    write(*,*) 'ens_writeEnsemble: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. ens%allocated ) then
      call utl_abort('ens_writeEnsemble: ensemble object not allocated!')
    end if

    !- 1. Initial setup

    nullify(varNamesInEns)
    if (present(varNames_opt)) then
      allocate(varNamesInEns(size(varNames_opt)))
      varNamesInEns(:) = varNames_opt(:)
    else
      call gsv_varNamesList(varNamesInEns, ens%statevector_work)
    end if

    if (present(ip3_opt)) then
      ip3 = ip3_opt
    else
      ip3 = 0
    end if

    lonPerPE    = ens%statevector_work%lonPerPE
    latPerPE    = ens%statevector_work%latPerPE
    lonPerPEmax = ens%statevector_work%lonPerPEmax
    latPerPEmax = ens%statevector_work%latPerPEmax
    ni          = ens%statevector_work%ni
    nj          = ens%statevector_work%nj
    nk          = ens%statevector_work%nk
    numStep     = ens%statevector_work%numStep

    ens%ensPathName = trim(ensPathName)

    ! Memory allocation
    numLevelsToSend = 10
    allocate(gd_send_r4(lonPerPEmax,latPerPEmax,numLevelsToSend,mpi_nprocs))
    allocate(gd_recv_r4(lonPerPEmax,latPerPEmax,numLevelsToSend,mpi_nprocs))
    gd_send_r4(:,:,:,:) = 0.0
    gd_recv_r4(:,:,:,:) = 0.0

    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

    do memberIndex = 1, ens%numMembers
      writeFilePE(memberIndex) = mod(memberIndex-1,mpi_nprocs)
    end do

    hco_ens => gsv_getHco(ens%statevector_work)
    vco_ens => gsv_getVco(ens%statevector_work)

    if (mpi_myid == 0) then
      write(*,*)
      write(*,*) 'ens_writeEnsemble: dateStampList=',dateStampList(1:numStep)
      write(*,*)
    end if

    !
    !- 2.  Ensemble forecasts writing loop
    !

    !- 2.1 Loop on time, ensemble member, variable, level
    do stepIndex = 1, numStep
      write(*,*) ' '
      write(*,*) 'ens_writeEnsemble: starting to write time level ', stepIndex

      ! allocate the needed statevector objects
      call gsv_allocate(statevector_member_r4, 1, hco_ens, vco_ens,  &
                        datestamp_opt=dateStampList(stepIndex), mpi_local_opt=.false., &
                        varNames_opt=varNamesInEns, dataKind_opt=4)

      do memberIndex = 1, ens%numMembers

        !  MPI communication: from 1 lat-lon tile per process to 1 ensemble member per process
        if (writeFilePE(memberIndex) == 0) then

          call tmg_start(13,'PRE_SUENS_COMM')
          batchnum = ceiling(dble(memberIndex + mpi_nprocs - 1)/dble(mpi_nprocs))

          do jk = 1, nk, numLevelsToSend
            jk2 = min(nk,jk+numLevelsToSend-1)
            numLevelsToSend2 = jk2 - jk + 1

            if ( ens%dataKind == 8 ) then
              !$OMP PARALLEL DO PRIVATE(jk3,memberIndex2,yourid)
              do jk3 = 1, numLevelsToSend2
                do memberIndex2 = 1+(batchnum-1)*mpi_nprocs, min(ens%numMembers, batchnum*mpi_nprocs)
                  yourid = writeFilePE(memberIndex2)
                  gd_send_r4(1:lonPerPE,1:latPerPE,jk3,yourid+1) = real(ens%allLev_r8(jk3+jk-1)%onelevel(memberIndex2,stepIndex,:,:),4)
                end do
              end do
              !$OMP END PARALLEL DO
            else
              !$OMP PARALLEL DO PRIVATE(jk3,memberIndex2,yourid)
              do jk3 = 1, numLevelsToSend2
                do memberIndex2 = 1+(batchnum-1)*mpi_nprocs, min(ens%numMembers, batchnum*mpi_nprocs)
                  yourid = writeFilePE(memberIndex2)
                  gd_send_r4(1:lonPerPE,1:latPerPE,jk3,yourid+1) = ens%allLev_r4(jk3+jk-1)%onelevel(memberIndex2,stepIndex,:,:)
                end do
              end do
              !$OMP END PARALLEL DO
            end if

            nsize = lonPerPEmax * latPerPEmax * numLevelsToSend2
            if (mpi_nprocs > 1) then
              call rpn_comm_alltoall(gd_send_r4(:,:,1:numLevelsToSend2,:),nsize,"mpi_real4",  &
                                     gd_recv_r4(:,:,1:numLevelsToSend2,:),nsize,"mpi_real4","GRID",ierr)
            else
              gd_recv_r4(:,:,1:numLevelsToSend2,1) = gd_send_r4(:,:,1:numLevelsToSend2,1)
            end if

            ptr3d_r4 => gsv_getField3D_r4(statevector_member_r4)
            !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
            do youridy = 0, (mpi_npey-1)
              do youridx = 0, (mpi_npex-1)
                yourid = youridx + youridy*mpi_npex
                ptr3d_r4(ens%statevector_work%allLonBeg(youridx+1):ens%statevector_work%allLonEnd(youridx+1),  &
                         ens%statevector_work%allLatBeg(youridy+1):ens%statevector_work%allLatEnd(youridy+1), jk:jk2) = &
                         gd_recv_r4(1:ens%statevector_work%allLonPerPE(youridx+1),  &
                         1:ens%statevector_work%allLatPerPE(youridy+1), 1:numLevelsToSend2, yourid+1)

              end do
            end do
            !$OMP END PARALLEL DO

          end do ! jk
          call tmg_stop(13)

        end if ! MPI communication


        if (mpi_myid == writeFilePE(memberIndex)) then

          write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

          if ( typvar == 'A' .or. typvar == 'R' ) then
            if ( typvar == 'R' ) then
              call fln_ensAnlFileName( ensFileName, ensPathName, tim_getDateStamp(), memberIndex_opt=memberIndex,  &
                                       ensFileNamePrefix_opt=ensFileNamePrefix, ensFileNameSuffix_opt='inc' )
            else
              call fln_ensAnlFileName( ensFileName, ensPathName, tim_getDateStamp(), memberIndex_opt=memberIndex,  &
                                       ensFileNamePrefix_opt=ensFileNamePrefix )
            end if
            ensFileExtLength = 4
          else
            call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=memberIndex, ensFileNamePrefix_opt=ensFileNamePrefix, &
                                  shouldExist_opt=.false., ensembleFileExtLength_opt=ensFileExtLength )
          end if

          etiketStr = etiket
          !  Write the file
          if (present(etiketAppendMemberNumber_opt)) then
            if (etiketAppendMemberNumber_opt) then
              write(ensFileExtLengthStr,"(I1)") ensFileExtLength
              write(memberIndexStr,'(I0.' // trim(ensFileExtLengthStr) // ')') memberIndex
              !! 12 is the maximum length of an etiket for RPN fstd files
              maximumBaseEtiketLength = 12 - ensFileExtLength
              if ( len(trim(etiket)) >= maximumBaseEtiketLength ) then
                etiketStr = etiket(1:maximumBaseEtiketLength) // trim(memberIndexStr)
              else
                etiketStr = trim(etiket) // trim(memberIndexStr)
              end if
            end if
          end if

          ! Determine if ensemble is full fields (if yes, will be converted from K to C)
          if (present(containsFullField_opt)) then
            containsFullField = containsFullField_opt
          else
            containsFullField = (.not. ens%meanIsRemoved)
          end if

          call gsv_writeToFile( statevector_member_r4, ensFileName, etiketStr, ip3_opt = ip3, & 
                                typvar_opt = typvar, numBits_opt = numBits_opt,  &
                                containsFullField_opt = containsFullField )

        end if ! locally written one member

      end do ! memberIndex

      ! deallocate the needed statevector objects
      call gsv_deallocate(statevector_member_r4)

    end do ! time

    deallocate(varNamesInEns)
    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)
    deallocate(datestamplist)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'ens_writeEnsemble: finished communicating and writing ensemble members...'

  end subroutine ens_writeEnsemble

end module ensembleStateVector_mod
