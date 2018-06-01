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
!! MODULE ensembleStateVector (prefix="ens")
!!
!! *Purpose*: Store and manipulate ensemble of state vectors and 
!!            the ensemble mean.
!!
!--------------------------------------------------------------------------
MODULE ensembleStateVector_mod
  use ramDisk_mod
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
  public :: struct_ens, ens_allocate, ens_deallocate
  public :: ens_readEnsemble, ens_writeEnsemble
  public :: ens_copyToStateWork, ens_getRepackMean_r8
  public :: ens_varExist, ens_getNumLev
  public :: ens_computeMean, ens_removeMean, ens_copyEnsMean, ens_copyMember, ens_recenter, ens_recenterControlMember
  public :: ens_computeStdDev, ens_copyEnsStdDev
  public :: ens_getRepack_r4
  public :: ens_getOffsetFromVarName, ens_getLevFromK, ens_getVarNameFromK, ens_getNumK, ens_getKFromLevVarName

  integer,external   :: get_max_rss

  type :: struct_repack_r4
    real(4), pointer :: onelevel(:,:,:,:) => null()
  end type struct_repack_r4

  type :: struct_repack_r8
    real(8), pointer :: onelevel(:,:,:,:) => null()
  end type struct_repack_r8

  type :: struct_ens
    private
    logical                       :: allocated = .false.
    integer                       :: numMembers
    integer                       :: dataKind
    type(struct_gsv)              :: statevector_work
    type(struct_repack_r8), allocatable :: repack_ensMean_r8(:), repack_ensStdDev_r8(:)
    type(struct_repack_r4), allocatable :: repack_r4(:)
    logical                       :: meanIsComputed = .false.
    logical                       :: stdDevIsComputed = .false.
    integer, allocatable          :: subEnsIndexList(:), nEnsSubEns(:)
    integer                       :: numSubEns
    character(len=256)            :: enspathname
  end type struct_ens

CONTAINS

  subroutine ens_allocate(ens, numMembers, numStep, hco_ens, vco_ens, &
                          dateStampList, varNames_opt)
    implicit none

    ! arguments
    type(struct_ens) :: ens
    integer :: numMembers, numStep
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    integer :: dateStampList(:)
    character(len=*), optional :: varNames_opt(:)  ! allow specification of assigned variables
    ! locals
    integer :: memberIndex, ierr
    integer :: jk, lon1, lon2, lat1, lat2, k1, k2

    if ( ens%allocated ) then
      write(*,*) 'ens_allocate: this object is already allocated, deallocating first.'
      call ens_deallocate( ens )
    end if

    call gsv_allocate( ens%statevector_work, &
                       numStep, hco_ens, vco_ens,  &
                       datestamplist_opt=dateStampList, mpi_local_opt=.true., &
                       varNames_opt=VarNames_opt, dataKind_opt=4 )

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd

    allocate( ens%repack_r4(k1:k2) )
    do jk = k1, k2
      allocate( ens%repack_r4(jk)%onelevel(numMembers,numStep,lon1:lon2,lat1:lat2) )
    end do

    ens%allocated = .true.
    ens%numMembers = numMembers
    ens%dataKind = 4

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocate

  subroutine ens_allocateMean(ens)
    implicit none

    ! arguments
    type(struct_ens) :: ens

    ! locals
    integer :: subEnsIndex, ierr, lon1, lon2, lat1, lat2, k1, k2, jk, numStep

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    allocate( ens%repack_ensMean_r8(k1:k2) )
    do jk = k1, k2
      allocate( ens%repack_ensMean_r8(jk)%onelevel(ens%numSubEns,numStep,lon1:lon2,lat1:lat2) )
      ens%repack_ensMean_r8(jk)%onelevel(:,:,:,:) = 0.0d0
    end do

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocateMean


  subroutine ens_allocateStdDev(ens)
    implicit none

    ! arguments
    type(struct_ens) :: ens

    ! locals
    integer :: lon1, lon2, lat1, lat2, k1, k2, jk, numStep

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    allocate( ens%repack_ensStdDev_r8(k1:k2) )
    do jk = k1, k2
      allocate( ens%repack_ensStdDev_r8(jk)%onelevel(1,numStep,lon1:lon2,lat1:lat2) )
      ens%repack_ensStdDev_r8(jk)%onelevel(:,:,:,:) = 0.0d0
    end do

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine ens_allocateStdDev


  subroutine ens_deallocate( ens )
    implicit none

    ! arguments
    type(struct_ens) :: ens

    ! locals
    integer :: memberIndex, subEnsIndex, k1, k2, jk

    if ( .not. ens%allocated ) return

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd

    do jk = k1, k2
      deallocate( ens%repack_r4(jk)%onelevel )
    end do
    deallocate( ens%repack_r4 )

    if (ens%stdDevIsComputed) then
      do jk = k1, k2
        deallocate( ens%repack_ensStdDev_r8(jk)%onelevel )
      end do
      deallocate( ens%repack_ensStdDev_r8 )
    end if

    if (ens%meanIsComputed) then
      do jk = k1, k2
        deallocate( ens%repack_ensMean_r8(jk)%onelevel )
      end do
      deallocate( ens%repack_ensMean_r8 )
      deallocate( ens%subEnsIndexList )
      deallocate( ens%nEnsSubEns )
    end if

    ens%allocated = .false.

  end subroutine ens_deallocate


  subroutine ens_copyToStateWork(ens, memberIndex)
    implicit none

    ! arguments
    type(struct_ens) :: ens
    integer          :: memberIndex

    ! locals
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, numStep, stepIndex

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    ptr4d_r4 => gsv_getField_r4(ens%statevector_work)
    do stepIndex = 1, numStep
      do jk = k1, k2
        ptr4d_r4(:,:,jk,stepIndex) = ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,:,:) 
      end do
    end do

  end subroutine ens_copyToStateWork


  function ens_getRepack_r4(ens,kIndex) result(repackLevel)
    implicit none

    ! arguments
    type(struct_ens) :: ens
    integer          :: kIndex
    real(4),pointer  :: repackLevel(:,:,:,:)

    ! locals
    integer          :: lon1, lat1

    lon1 = ens%statevector_work%myLonBeg
    lat1 = ens%statevector_work%myLatBeg

    repackLevel(1:,1:,lon1:,lat1:) => ens%repack_r4(kIndex)%onelevel(:,:,:,:)

  end function ens_getRepack_r4


  function ens_getRepackMean_r8(ens,subEnsIndex,kIndex) result(field)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    integer           :: subEnsIndex, kIndex
    real(8),pointer   :: field(:,:,:)

    ! locals
    integer           :: lon1,lat1

    lon1 = ens%statevector_work%myLonBeg
    lat1 = ens%statevector_work%myLatBeg

    field(1:, lon1:, lat1:) => ens%repack_ensMean_r8(kIndex)%onelevel(subEnsIndex,:,:,:)

  end function ens_getRepackMean_r8


  subroutine ens_copyEnsMean(ens, statevector, subEnsIndex_opt)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer, optional :: subEnsIndex_opt

    ! locals
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep, subEnsIndex

    if( present(subEnsIndex_opt) ) then
      subEnsIndex = subEnsIndex_opt
    else
      subEnsIndex = 1
    end if

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. statevector%allocated) then
      call gsv_allocate(statevector, numStep,  &
                        ens%statevector_work%hco, ens%statevector_work%vco,  &
                        datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., dataKind_opt=8 )
    end if

    ptr4d_r8 => gsv_getField_r8(statevector)
    do stepIndex = 1, numStep
      do jk = k1, k2
        ptr4d_r8(:,:,jk,stepIndex) = ens%repack_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,:,:)
      end do
    end do

  end subroutine ens_copyEnsMean

  subroutine ens_copyEnsStdDev(ens, statevector)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector

    ! locals
    real(8), pointer :: ptr4d_r8(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    if (.not. statevector%allocated) then
      call gsv_allocate(statevector, numStep,  &
                        ens%statevector_work%hco, ens%statevector_work%vco,  &
                        datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., dataKind_opt=8 )
    end if

    ptr4d_r8 => gsv_getField_r8(statevector)
    do stepIndex = 1, numStep
      do jk = k1, k2
        ptr4d_r8(:,:,jk,stepIndex) = ens%repack_ensStdDev_r8(jk)%onelevel(1,stepIndex,:,:)
      end do
    end do

  end subroutine ens_copyEnsStdDev


  subroutine ens_copyMember(ens, statevector, memberIndex)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    type(struct_gsv)  :: statevector
    integer           :: memberIndex

    ! locals
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    integer          :: k1, k2, jk, stepIndex, numStep

    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    call gsv_allocate( statevector, numStep,  &
                       ens%statevector_work%hco, ens%statevector_work%vco,  &
                       datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., dataKind_opt=4 )

    ptr4d_r4 => gsv_getField_r4(statevector)
    do stepIndex = 1, numStep
      do jk = k1, k2
        ptr4d_r4(:,:,jk,stepIndex) = ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,:,:)
      end do
    end do

  end subroutine ens_copyMember


  function ens_varExist(ens,varName) result(varExist)
    implicit none

    ! arguments
    type(struct_ens)             :: ens
    character(len=*), intent(in) :: varName
    logical                      :: varExist 

    varExist = gsv_varExist(ens%statevector_work, varName)

  end function ens_varExist


  function ens_getNumLev(ens,varLevel) result(nlev)
    implicit none

    ! arguments
    type(struct_ens), intent(in)  :: ens
    character(len=*), intent(in)  :: varLevel
    integer                       :: nlev

    nlev = vco_getNumLev(ens%statevector_work%vco,varLevel)

  end function ens_getNumLev


  function ens_getNumK(ens) result(numK)
    implicit none

    ! arguments
    type(struct_ens), intent(in)  :: ens
    integer                       :: numK

    numK = 1 + ens%statevector_work%mykEnd - ens%statevector_work%mykBeg

  end function ens_getNumK


  function ens_getOffsetFromVarName(ens,varName) result(offset)
    implicit none
    type(struct_ens)             :: ens
    character(len=*), intent(in) :: varName
    integer                      :: offset

    offset=gsv_getOffsetFromVarName(ens%statevector_work,varName)

  end function ens_getOffsetFromVarName


  function ens_getLevFromK(ens,kIndex) result(levIndex)
    implicit none

    ! arguments
    type(struct_ens), intent(in) :: ens
    integer, intent(in)          :: kIndex
    integer                      :: levIndex

    levIndex = gsv_getLevFromK(ens%statevector_work,kIndex)

  end function ens_getLevFromK


  function ens_getKFromLevVarName(ens, levIndex, varName) result(kIndex)
    implicit none

    ! arguments
    type(struct_ens), intent(in) :: ens
    integer                      :: levIndex
    character(len=*)             :: varName
    integer                      :: kIndex

    kIndex = levIndex + gsv_getOffsetFromVarName(ens%statevector_work,trim(varName))

  end function ens_getKFromLevVarName


  function ens_getVarNameFromK(ens,kIndex) result(varName)
    implicit none

    ! arguments
    type(struct_ens), intent(in) :: ens
    integer, intent(in)          :: kIndex
    character(len=4)             :: varName

    varName = gsv_getVarNameFromK(ens%statevector_work,kIndex)

  end function ens_getVarNameFromK


  subroutine ens_computeMean(ens, computeSubEnsMeans_opt, numSubEns_opt)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    logical, optional :: computeSubEnsMeans_opt
    integer, optional :: numSubEns_opt

    ! locals
    logical           :: computeSubEnsMeans, lExists
    character(len=256), parameter :: subEnsIndexFileName = 'subEnsembleIndex.txt'
    integer           :: kulin, ierr, memberIndex, memberIndex2, stepIndex, subEnsIndex
    integer           :: k1, k2, jk, lon1, lon2, lat1, lat2, numStep, ji, jj
    integer           :: fnom, fclos
    real(4), pointer  :: ptr_repack_r4(:,:,:,:)
    real(8), pointer  :: ptr4d_r8(:,:,:,:)

    if (present(computeSubEnsMeans_opt)) then
      computeSubEnsMeans = computeSubEnsMeans_opt
    else
      computeSubEnsMeans = .false.
    end if

    ! Read sub-ensemble index list from file, if it exists
    allocate(ens%subEnsIndexList(ens%numMembers))
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
    allocate(ens%nEnsSubEns(ens%numSubEns))
    ens%nEnsSubEns(:) = 0
    do memberIndex = 1, ens%numMembers
      ens%nEnsSubEns(ens%subEnsIndexList(memberIndex)) = ens%nEnsSubEns(ens%subEnsIndexList(memberIndex)) + 1
    end do
    write(*,*) 'ens_computeMean: number of sub-ensembles = ', ens%numSubEns
    write(*,*) 'ens_computeMean: number of members in each sub-ensemble = ', ens%nensSubEns(:)

    call ens_allocateMean(ens)
    ens%meanIsComputed = .true.

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep
    ! Compute ensemble mean(s)
!$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex,subEnsIndex)
    do jk = k1, k2
      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, ens%statevector_work%numStep
            do memberIndex = 1, ens%numMembers
              ens%repack_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj) = &
                   ens%repack_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj) + &
                   dble(ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj))
            end do
            do subEnsIndex = 1, ens%numSubEns
              ens%repack_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,ji,jj) = &
                   ens%repack_ensMean_r8(jk)%onelevel(subEnsIndex,stepIndex,ji,jj) /  &
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


  subroutine ens_computeStdDev(ens)
    implicit none

    ! arguments
    type(struct_ens)  :: ens

    ! locals
    integer           :: kulin, ierr, memberIndex, memberIndex2, stepIndex, subEnsIndex
    integer           :: k1, k2, jk, lon1, lon2, lat1, lat2, numStep, ji, jj
    real(8), allocatable  :: subEnsStdDev(:)

    if (.not.ens%meanIsComputed) then
      if (mpi_myid == 0) write(*,*) 'ens_computeStdDev: compute Mean since it was not already done'
      call ens_computeMean( ens )
    end if

    ! Read sub-ensemble index list from file, if it exists
    ! The sub-ensembles should have been already read in routine 'ens_computeMean'
    write(*,*) 'ens_computeStdDev: number of sub-ensembles = ', ens%numSubEns
    write(*,*) 'ens_computeStdDev: number of members in each sub-ensemble = ', ens%nensSubEns(:)

    call ens_allocateStdDev(ens)
    ens%StdDevIsComputed = .true.

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

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
            do memberIndex = 1, ens%numMembers
              subEnsStdDev(ens%subEnsIndexList(memberIndex)) = subEnsStdDev(ens%subEnsIndexList(memberIndex)) + &
                   (dble(ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj))-ens%repack_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj))**2
            end do
            do subEnsIndex = 1, ens%numSubEns
              ens%repack_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) = ens%repack_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) + &
                   ens%nEnsSubEns(subEnsIndex)*subEnsStdDev(subEnsIndex)/(ens%nEnsSubEns(subEnsIndex)-1)
            end do
            ens%repack_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) = sqrt( ens%repack_ensStdDev_r8(jk)%onelevel(1,stepIndex,ji,jj) / dble(ens%numMembers) )
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

    deallocate(subEnsStdDev)

  end subroutine ens_computeStdDev


  subroutine ens_removeMean(ens)
    implicit none

    ! arguments
    type(struct_ens) :: ens

    ! locals
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
              ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj) =  &
                   real( (real(ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj),8) -  &
                   ens%repack_ensMean_r8(jk)%onelevel(ens%subEnsIndexList(memberIndex),stepIndex,ji,jj)), 4 )
            end do
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine ens_removeMean


  subroutine ens_recenter(ens,recenteringMean,recenteringCoeff,ensembleCenter_opt,ensembleControlMember_opt)
    implicit none

    !! We want to compute:
    !!    x_recentered = x_original + recenteringCoeff*(x_recenteringMean - x_ensembleMean)

    ! arguments
    type(struct_ens) :: ens
    type(struct_gsv) :: recenteringMean
    type(struct_gsv), optional :: ensembleCenter_opt, ensembleControlMember_opt
    real(8)          :: recenteringCoeff

    ! locals
    real(8), pointer :: ptr4d_r8(:,:,:,:), ensembleCenter_r8(:,:,:,:), ptr4d_ensembleControlmember_r8(:,:,:,:)
    real(8) :: increment
    integer :: lon1, lon2, lat1, lat2, k1, k2, numStep
    integer :: jk, jj, ji, stepIndex, memberIndex

    lon1 = ens%statevector_work%myLonBeg
    lon2 = ens%statevector_work%myLonEnd
    lat1 = ens%statevector_work%myLatBeg
    lat2 = ens%statevector_work%myLatEnd
    k1 = ens%statevector_work%mykBeg
    k2 = ens%statevector_work%mykEnd
    numStep = ens%statevector_work%numStep

    ptr4d_r8 => gsv_getField_r8(recenteringMean)
    if(present(ensembleCenter_opt)) then
      ensembleCenter_r8 => gsv_getField_r8(ensembleCenter_opt)
    else
      nullify(ensembleCenter_r8)
    end if

    if (present(ensembleControlMember_opt)) then
      ptr4d_ensembleControlmember_r8 => gsv_getField_r8(ensembleControlMember_opt)
    else
      nullify(ptr4d_ensembleControlmember_r8)
    end if

!$OMP PARALLEL DO PRIVATE (jk,jj,ji,stepIndex,memberIndex,increment)
    do jk = k1, k2
      do jj = lat1, lat2
        do ji = lon1, lon2
          do stepIndex = 1, numStep
            if(present(ensembleCenter_opt)) then
              increment = ptr4d_r8(ji,jj,jk,stepIndex) - ensembleCenter_r8(ji,jj,jk,stepIndex)
            else
              increment = ptr4d_r8(ji,jj,jk,stepIndex) - ens%repack_ensMean_r8(jk)%onelevel(1,stepIndex,ji,jj)
            end if
            if (present(ensembleControlMember_opt)) then
              ptr4d_ensembleControlMember_r8(ji,jj,jk,stepIndex) = ptr4d_ensembleControlMember_r8(ji,jj,jk,stepIndex) + recenteringCoeff*increment
            else
              do memberIndex = 1, ens%numMembers
                ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj) =  &
                     real( real(ens%repack_r4(jk)%onelevel(memberIndex,stepIndex,ji,jj),8) + recenteringCoeff*increment, 4)
              end do
            end if
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine ens_recenter


  subroutine ens_recenterControlMember(ens,hco_ens,vco_ens,ensFileName,ensPathName,recenteringMean,recenteringCoeff,HUcontainsLQ,ensembleCenter_opt)
    implicit none

    !! We want to compute:
    !!    x_recentered = x_original + recenteringCoeff*(x_recenteringMean - x_ensembleMean)

    ! arguments
    type(struct_ens) :: ens
    type(struct_vco), pointer :: vco_ens
    type(struct_hco), pointer :: hco_ens
    character(len=*) :: ensFileName, ensPathName
    type(struct_gsv) :: recenteringMean
    real(8)          :: recenteringCoeff
    logical          :: HUcontainsLQ
    type(struct_gsv), optional :: ensembleCenter_opt

    ! locals
    type(struct_gsv) :: statevector_ensembleControlMember
    integer          :: stepIndex, numStep

    numStep = ens%statevector_work%numStep

    call fln_ensFileName( ensFileName, ensPathName, memberIndex = 0)

    call gsv_allocate(statevector_ensembleControlMember, numStep, hco_ens, vco_ens, &
         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.)

    do stepIndex = 1, numStep
      if(mpi_myid == 0) write(*,*) ''
      if(mpi_myid == 0) write(*,*) 'ens_shiftEnsembleControlMember: reading ensemble control member for time step: ',stepIndex
      call gsv_readFromFile(statevector_ensembleControlMember, trim(ensFileName), ' ', ' ',  &
           stepIndex_opt=stepIndex, unitConversion_opt=.true., HUcontainsLQ_opt=HUcontainsLQ )
    end do

    call ens_recenter(ens,recenteringMean,recenteringCoeff,ensembleCenter_opt = ensembleCenter_opt, &
         ensembleControlMember_opt = statevector_ensembleControlMember)

    call gsv_deallocate(statevector_ensembleControlMember)

  end subroutine ens_recenterControlMember


  subroutine ens_readEnsemble(ens, ensPathName, biPeriodic, ctrlVarHumidity, &
                              vco_file_opt, varNames_opt)
    implicit none

    ! arguments
    type(struct_ens) :: ens
    character(len=*) :: ensPathName
    logical          :: biPeriodic
    character(len=*) :: ctrlVarHumidity
    character(len=*), optional :: varNames_opt(:)
    type(struct_vco), pointer, optional :: vco_file_opt

    ! locals
    type(struct_gsv) :: statevector_file_r4, statevector_hint_r4, statevector_member_r4
    type(struct_hco), pointer :: hco_file, hco_ens, hco_coregrid
    type(struct_vco), pointer :: vco_file, vco_ens
    real(4), allocatable :: gd_send_r4(:,:,:,:)
    real(4), allocatable :: gd_recv_r4(:,:,:,:)
    real(4), pointer     :: ptr3d_r4(:,:,:), ptr4d_r4(:,:,:,:)
    real(8)              :: multFactor
    integer,pointer :: dateStampList(:)
    integer :: batchnum, nsize, status, ierr
    integer :: yourid, youridx, youridy
    integer :: readFilePE(1000)
    integer :: memberIndexOffset, totalEnsembleSize
    integer :: length_envVariable
    integer :: lonPerPEmax, latPerPEmax, ni, nj, nk, numStep, numlevelstosend, numlevelstosend2
    integer :: memberIndex, memberIndex2, fileMemberIndex, stepIndex, jvar, jk, jk2, jk3
    character(len=256) :: ensFileName
    character(len=32)  :: envVariable
    character(len=2)   :: typvar
    character(len=12)  :: etiket
    character(len=4)   :: varName
    character(len=4)   :: ensVarNamesWanted_dummy(vnl_numvarmax)
    character(len=4), allocatable :: ensVarNamesWanted(:)
    logical :: verticalInterpNeeded, horizontalInterpNeeded, horizontalPaddingNeeded

    write(*,*) 'ens_readEnsemble: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. ens%allocated ) then
      call utl_abort('ens_readEnsemble: ensemble object not allocated!')
    end if

    !
    !- 1. Initial setup
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
    call fln_ensFileName(ensFileName, ensPathName, 1)
    nullify(hco_file)
    call hco_SetupFromFile(hco_file, ensFileName, ' ', 'ENSFILEGRID')
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

    !
    !- 2.  Ensemble forecasts reading loop
    !

    !- 2.1 Loop on time, ensemble member, variable, level
    do stepIndex = 1, numStep
      write(*,*) ' '
      write(*,*) 'ens_readEnsemble: starting to read time level ', stepIndex

      ! allocate the needed statevector objects
      call gsv_allocate(statevector_member_r4, 1, hco_ens, vco_ens,  &
                        datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                          varNames_opt = varNames_opt, dataKind_opt = 4)
      if (horizontalInterpNeeded .or. verticalInterpNeeded .or. horizontalPaddingNeeded) then
        call gsv_allocate(statevector_file_r4, 1, hco_file, vco_file,  &
                            datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                            varNames_opt = varNames_opt, dataKind_opt = 4)
      end if
      if (verticalInterpNeeded) then
        call gsv_allocate(statevector_hint_r4, 1, hco_ens, vco_file,  &
                          datestamp_opt = dateStampList(stepIndex), mpi_local_opt = .false., &
                          varNames_opt = varNames_opt, dataKind_opt = 4)
      end if
      
      do memberIndex = 1, ens%numMembers

        if (mpi_myid == readFilePE(memberIndex)) then

          write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

          !  Read the file
          fileMemberIndex = 1+mod(memberIndex+memberIndexOffset-1, totalEnsembleSize)
          call fln_ensFileName(ensFileName, ensPathName, fileMemberIndex)
          typvar = ' '
          etiket = ' '
          if (.not. horizontalInterpNeeded  .and. &
              .not. verticalInterpNeeded    .and. &
              .not. horizontalPaddingNeeded ) then
            call gsv_readFile(statevector_member_r4, ensFileName, etiket, typvar)
          else
            call gsv_readFile(statevector_file_r4, ensFileName, etiket, typvar)
          end if
          if (stepIndex == numStep) then
            ierr = ram_remove(ensFileName)
          end if

          ! do any required interpolation
          if (horizontalInterpNeeded .and. verticalInterpNeeded) then
            call gsv_hInterpolate_r4(statevector_file_r4, statevector_hint_r4)
            call gsv_vInterpolate_r4(statevector_hint_r4, statevector_member_r4, Ps_in_hPa_opt=.true.)

          else if (horizontalInterpNeeded .and. .not. verticalInterpNeeded) then
            call gsv_hInterpolate_r4(statevector_file_r4, statevector_member_r4)

          else if (.not. horizontalInterpNeeded .and. verticalInterpNeeded) then
            if (horizontalPaddingNeeded) then
              call gsv_hPad(statevector_file_r4, statevector_hint_r4)
            else
              call gsv_copy(statevector_file_r4, statevector_hint_r4)
            end if
            call gsv_vInterpolate_r4(statevector_hint_r4, statevector_member_r4, Ps_in_hPa_opt=.true.)

          else if (horizontalPaddingNeeded) then
            call gsv_hPad(statevector_file_r4, statevector_member_r4)
          end if

          ! unit conversion
          VAR_LOOP: do jvar = 1, vnl_numvarmax
            varName = vnl_varNameList(jvar)
            if ( .not. gsv_varExist(statevector_member_r4, varName) ) cycle VAR_LOOP
            ptr3d_r4 => gsv_getField3D_r4(statevector_member_r4, varName)

            if ( trim(varName) == 'UU' .or. trim(varName) == 'VV') then
              multFactor = MPC_M_PER_S_PER_KNOT_R8 ! knots -> m/s
            else if ( trim(varName) == 'P0' ) then
              multFactor = MPC_PA_PER_MBAR_R8 ! hPa -> Pa
            else
              multFactor = 1.0D0 ! no conversion
            end if
            ptr3d_r4(:,:,:) = real( multFactor * ptr3d_r4(:,:,:), 4 )

            if      (trim(varName) == 'HU' .and. ctrlVarHumidity == 'LQ') then
              ptr3d_r4(:,:,:) = sngl(log(max(real(ptr3d_r4(:,:,:),8),MPC_MINIMUM_HU_R8)))
            else if (trim(varName) == 'HU' .and. ctrlVarHumidity == 'HU') then
              ptr3d_r4(:,:,:) = sngl(    max(real(ptr3d_r4(:,:,:),8),MPC_MINIMUM_HU_R8))
            end if
          end do VAR_LOOP

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

            call tmg_start(110,'ENS_TO_REPACK')
!$OMP PARALLEL DO PRIVATE(jk3,memberIndex2,yourid)
            do jk3 = 1, numLevelsToSend2
              do memberIndex2 = 1+(batchnum-1)*mpi_nprocs, memberIndex
                yourid = readFilePE(memberIndex2)
                ens%repack_r4(jk3+jk-1)%onelevel(memberIndex2,stepIndex, :, :) =  &
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


  subroutine ens_writeEnsemble(ens, ensPathName, ensFileNamePrefix, ctrlVarHumidity, etiket, &
                               typvar, etiketAppendMemberNumber_opt, varNames_opt, ip3_opt, numBits_opt)
    implicit none

    ! arguments
    type(struct_ens)  :: ens
    character(len=*)  :: ensPathName
    character(len=*)  :: ensFileNamePrefix
    character(len=*)  :: ctrlVarHumidity
    character(len=*)  :: etiket
    character(len=*)  :: typvar
    character(len=*), optional :: varNames_opt(:)  ! allow specification of variables
    integer, optional :: ip3_opt, numBits_opt
    logical, optional :: etiketAppendMemberNumber_opt

    ! locals
    type(struct_gsv) :: statevector_member_r4
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    real(4), allocatable :: gd_send_r4(:,:,:,:)
    real(4), allocatable :: gd_recv_r4(:,:,:,:)
    real(4), pointer     :: ptr3d_r4(:,:,:), ptr4d_r4(:,:,:,:)
    real(8)              :: multFactor
    integer, allocatable :: dateStampList(:)
    integer :: batchnum, nsize, status, ierr
    integer :: yourid, youridx, youridy
    integer :: writeFilePE(1000)
    integer :: lonPerPE, lonPerPEmax, latPerPE, latPerPEmax, ni, nj, nk, numStep, numlevelstosend, numlevelstosend2
    integer :: memberIndex, memberIndex2, stepIndex, jvar, jk, jk2, jk3, ip3, ensFileExtLength, maximumBaseEtiketLength
    character(len=256) :: ensFileName
    character(len=12) :: etiketStr  ! this is the etiket that will be used to write files
    character(len=6) :: memberIndexStrFormat  !  will contain the character string '(I0.4)' to have 4 characters in the member extension
    character(len=10), allocatable :: memberIndexStr ! this is the member number in a character string

    write(*,*) 'ens_writeEnsemble: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. ens%allocated ) then
      call utl_abort('ens_writeEnsemble: ensemble object not allocated!')
    end if

    !- 1. Initial setup

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
                        varNames_opt=varNames_opt, dataKind_opt=4)

      do memberIndex = 1, ens%numMembers


        !  MPI communication: from 1 lat-lon tile per process to 1 ensemble member per process
        if (writeFilePE(memberIndex) == 0) then

          call tmg_start(13,'PRE_SUENS_COMM')
          batchnum = ceiling(dble(memberIndex + mpi_nprocs - 1)/dble(mpi_nprocs))

          do jk = 1, nk, numLevelsToSend
            jk2 = min(nk,jk+numLevelsToSend-1)
            numLevelsToSend2 = jk2 - jk + 1

!$OMP PARALLEL DO PRIVATE(jk3,memberIndex2,yourid)
            do jk3 = 1, numLevelsToSend2
              do memberIndex2 = 1+(batchnum-1)*mpi_nprocs, min(ens%numMembers, batchnum*mpi_nprocs)
                yourid = writeFilePE(memberIndex2)
                gd_send_r4(1:lonPerPE,1:latPerPE,jk3,yourid+1) = ens%repack_r4(jk3+jk-1)%onelevel(memberIndex2,stepIndex,:,:)
              end do
            end do
!$OMP END PARALLEL DO

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

          call fln_ensFileName( ensFileName, ensPathName, memberIndex, ensFileNamePrefix_opt = ensFileNamePrefix, &
               shouldExist_opt = .false., ensembleFileExtLength_opt = ensFileExtLength )

          etiketStr = etiket
          !  Write the file
          if (present(etiketAppendMemberNumber_opt)) then
            if (etiketAppendMemberNumber_opt) then
              allocate(character(len=ensFileExtLength) :: memberIndexStr)
              write(memberIndexStr,"(I1)") ensFileExtLength
              memberIndexStrFormat = '(I0.' // trim(memberIndexStr) // ')'
              write(memberIndexStr,memberIndexStrFormat) memberIndex
              !! 12 is the maximum length of an etiket for RPN fstd files
              maximumBaseEtiketLength = 12 - ensFileExtLength
              if ( len(trim(etiket)) >= maximumBaseEtiketLength ) then
                etiketStr = etiket(1:maximumBaseEtiketLength) // trim(memberIndexStr)
              else
                etiketStr = trim(etiket) // trim(memberIndexStr)
              end if
              deallocate(memberIndexStr)
            end if
          end if

          call gsv_writeToFile( statevector_member_r4, ensFileName, etiketStr, ip3_opt = ip3, typvar_opt = typvar , numBits_opt = numBits_opt)
        end if ! locally written one member


      end do ! memberIndex

      ! deallocate the needed statevector objects
      call gsv_deallocate(statevector_member_r4)

    end do ! time

    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)
    deallocate(datestamplist)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'ens_writeEnsemble: finished communicating and writing ensemble members...'

  end subroutine ens_writeEnsemble

end module ensembleStateVector_mod
