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

module gridBinning_mod
  ! MODULE gridBinning_mod (prefix='gbi' category='3. High-level transformations')
  !
  ! :Purpose: To compute categorical mean and standard deviation for gridded data
  !           contained in a gridStateVector or in an ensemble of
  !           gridStateVectors (e.g. the respective mean over land and sea)
  !
  use mpi_mod
  use mpivar_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use utilities_mod
  use horizontalCoord_mod
  use timeCoord_mod
  implicit none
  save
  private

  ! public procedures
  public :: struct_gbi
  public :: gbi_setup, gbi_deallocate
  public :: gbi_stdDev, gbi_mean

  type :: struct_gbi
    character(len=24) :: binningStrategy
    type(struct_gsv)  :: statevector_bin2d
    integer           :: numBins2d
  end type struct_gbi

  integer, external  :: get_max_rss

  ! Control parameter for the level of listing output
  logical, parameter :: verbose = .true.

  ! module interfaces
  interface gbi_stdDev
    module procedure gbi_stdDev_ens
  end interface gbi_stdDev

  interface gbi_mean
    module procedure gbi_mean_gsv
  end interface gbi_mean

contains

  !--------------------------------------------------------------------------
  ! gbi_setup
  !--------------------------------------------------------------------------
  subroutine gbi_setup(gbi, binningStrategy, statevector_template, hco_coregrid, &
                       mpi_distribution_opt, writeBinsToFile_opt)
    implicit none

    type(struct_gbi) :: gbi
    character(len=*),  intent(in) :: binningStrategy
    type(struct_gsv) :: statevector_template
    type(struct_hco), pointer :: hco_coregrid
    character(len=*), optional, intent(in) :: mpi_distribution_opt
    logical, optional, intent(in) :: writeBinsToFile_opt

    type(struct_gsv) :: statevector_landSeaTopo

    integer :: myLonBeg, myLonEnd
    integer :: myLatBeg, myLatEnd
    integer :: latIndex, lonIndex

    real(8) :: binCategory

    logical :: allocHeightSfc
    logical :: writeBinsToFile
    logical :: skip
    logical :: mpi_local

    real(4), pointer :: bin2d(:,:,:)
    real(4), pointer :: data2d(:,:,:)

    character(len=24) :: mpi_distribution

    !
    !- 1.  Allocate a 2D statevector that will contains the bin category for each 
    !      grid point
    !
    if (.not. gsv_allocated(statevector_template)) then
      call utl_abort('gbi_setup: the input template statevector was not allocated')
    end if

    mpi_local = statevector_template%mpi_local
    if (present(mpi_distribution_opt)) then
      mpi_distribution = mpi_distribution_opt
      if ( mpi_distribution ==  'None') then
        mpi_local = .false.
      end if
    else
      mpi_distribution = 'Tiles'
    end if

    if (trim(binningStrategy) == 'landSeaTopo') then
      allocHeightSfc = .true.
    else
      allocHeightSfc = .false.
    end if

    call gsv_allocate(gbi%statevector_bin2d, 1, statevector_template%hco,    &
                      statevector_template%vco, varNames_opt=(/'BIN'/),      &
                      mpi_distribution_opt=mpi_distribution, dataKind_opt=4, &
                      dateStamp_opt=statevector_template%dateStamp3d,        &
                      mpi_local_opt=mpi_local, allocHeightSfc_opt=allocHeightSfc)

    !
    !- 2.  Set the bins
    !
    call gsv_getField(gbi%statevector_bin2d,bin2d)

    myLonBeg = gbi%statevector_bin2d%myLonBeg
    myLonEnd = gbi%statevector_bin2d%myLonEnd
    myLatBeg = gbi%statevector_bin2d%myLatBeg
    myLatEnd = gbi%statevector_bin2d%myLatEnd

    select case(trim(binningStrategy))
    case ('GridPoint')
      if (mpi_myid == 0) then
        write(*,*) 'gbi_setup : BIN_TYPE = No horizontal averaging '
        write(*,*) '            > One bin per horizontal grid point'
      end if

      gbi%numBins2d = statevector_template%hco%ni * statevector_template%hco%nj
      binCategory = 0
      do latIndex = 1, statevector_template%hco%nj
        do lonIndex = 1, statevector_template%hco%ni
          binCategory = binCategory + 1
          if (lonIndex >= myLonBeg .and. lonIndex <= myLatEnd .and. &
              latIndex >= myLatBeg .and. latIndex <= myLatEnd ) then
            bin2d(lonIndex,latIndex,1) = real(binCategory,4)
          end if
        end do
      end do

    case ('YrowBand')
      if (mpi_myid == 0) then
        write(*,*)
        write(*,*) 'gbi_setup : BIN_TYPE = One bin per Y row'
      end if

      gbi%numBins2d = statevector_template%hco%nj
      binCategory = 0
      do latIndex = 1, statevector_template%hco%nj
        binCategory = binCategory + 1
        if (latIndex >= myLatBeg .and. latIndex <= myLatEnd ) then
          bin2d(:,latIndex,1) = real(binCategory,4)
        end if
      end do

    case ('HorizontalMean','FullDomain')
      if (mpi_myid == 0) then
        write(*,*)
        write(*,*) 'gbi_setup : BIN_TYPE = Average over all horizontal points'
      end if

      gbi%numBins2d = 1
      binCategory   = 1
      bin2d(:,:,1)  = real(binCategory,4)

    case ('landSeaTopo')
      if (mpi_myid == 0) then
        write(*,*)
        write(*,*) 'gbi_setup : BIN_TYPE = land/sea binning + storage of the topography field'
      end if

      call gsv_allocate(statevector_landSeaTopo, 1, statevector_template%hco,     &
                        statevector_template%vco, varNames_opt=(/'MG'/),          &
                        mpi_distribution_opt=mpi_distribution, dataKind_opt=4,    &
                        dateStamp_opt=statevector_template%dateStamp3d,           &
                        mpi_local_opt=mpi_local, hInterpolateDegree_opt='LINEAR', &
                        allocHeightSfc_opt=.true. )

      call gsv_readFromFile(statevector_landSeaTopo, './dataForGridBinning.fst', ' ', ' ', &
                            readHeightSfc_opt=.true.)

      call gsv_getField(statevector_landSeaTopo,data2d)

      gbi%numBins2d = 2 ! (land or sea)
      do latIndex= myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd

          skip = .false.
          if (.not. statevector_template%hco%global) then
            if (lonIndex > hco_coregrid%ni .or. latIndex > hco_coregrid%nj) then
              ! We are in the extension zone
              skip = .true.
            end if
          end if

          if (skip) then
            bin2d(lonIndex,latIndex,1) = -1.0
          else
            if (data2D(lonIndex,latIndex,1) >= 0.05) then
              bin2d(lonIndex,latIndex,1) = 1.0 ! land point
            else
              bin2d(lonIndex,latIndex,1) = 2.0 ! sea  point
            end if
          end if

        end do
      end do

      ! Store the topography in gbi%statevector_bin2d
      call gsv_copyHeightSfc(statevector_landSeaTopo, & ! IN
                         gbi%statevector_bin2d )    ! OUT

      call gsv_deallocate(statevector_landSeaTopo)

    case default
      call utl_abort('gbi_setup : Invalid Binning Strategy : '//trim(BinningStrategy))
    end select

    !
    !- 3.  Write the bins to a file (if desired)
    !
    if ( present(writeBinsToFile_opt) ) then
      writeBinsToFile = writeBinsToFile_opt
    else
      writeBinsToFile = .false.
    end if
    
    if (writeBinsToFile .and. (mpi_local .or. mpi_myid == 0) ) then
      call gsv_writeToFile(gbi%statevector_bin2d, './gridBinning.fst', & ! IN
                           'BINNING')                                    ! IN
    end if

  end subroutine gbi_setup
  
  !--------------------------------------------------------------------------
  ! gbi_deallocate
  !--------------------------------------------------------------------------
  subroutine gbi_deallocate(gbi)
    implicit none

    type(struct_gbi) :: gbi

    call gsv_deallocate(gbi%statevector_bin2d)

  end subroutine gbi_deallocate

  !--------------------------------------------------------------------------
  ! gbi_mean_gsv
  !--------------------------------------------------------------------------
  subroutine gbi_mean_gsv(gbi, statevector_in, statevector_out)
    implicit none

    type(struct_gbi) :: gbi
    type(struct_gsv) :: statevector_in
    type(struct_gsv) :: statevector_out

    integer :: myBinCount(gbi%numBins2d)
    integer ::   binCount(gbi%numBins2d)
    real(8) :: myBinSum  (gbi%numBins2d)
    real(8) ::   binMean (gbi%numBins2d)

    real(8), pointer :: field4d(:,:,:,:)
    real(8), pointer :: mean4d (:,:,:,:)

    real(4), pointer :: bin2d(:,:,:)

    integer :: myLonBeg, myLonEnd
    integer :: myLatBeg, myLatEnd
    integer :: latIndex, lonIndex, varLevIndex, stepIndex, binIndex
    integer :: nVarLev, nStep, ier

    if (mpi_myid == 0 .and. verbose) then
      write(*,*)
      write(*,*) 'gbi_mean_gsv: Starting...'
    end if

    nVarLev   = statevector_in%nk
    nStep     = statevector_in%numStep
    myLonBeg  = statevector_in%myLonBeg
    myLonEnd  = statevector_in%myLonEnd
    myLatBeg  = statevector_in%myLatBeg
    myLatEnd  = statevector_in%myLatEnd

    call gsv_getField(statevector_in,field4d)
    call gsv_getField(statevector_out,mean4d)

    call gsv_getField(gbi%statevector_bin2d,bin2d)

    do stepIndex = 1, nStep
      do varLevIndex = 1, nVarLev

        !- Sum the values per bin
        myBinCount(:) = 0
        myBinSum  (:) = 0.d0
        do latIndex= myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            binIndex = int(bin2d(lonIndex,latIndex,1))
            if (binIndex /= -1) then
              myBinCount(binIndex) = myBinCount(binIndex) + 1
              myBinSum(binIndex)   = myBinSum(binIndex)   + field4d(lonIndex,latIndex,varLevIndex,stepIndex)
            end if
          end do
        end do

        !- Compute the mean per bin
        call rpn_comm_allreduce(myBinCount,binCount,gbi%numBins2d,"MPI_INTEGER"         ,"MPI_SUM","GRID",ier)
        do binIndex = 1, gbi%numBins2d
          binMean(binIndex) = myBinSum(binIndex)
          call mpi_allreduce_sumreal8scalar(binMean(binIndex),"GRID")
        end do

        do binIndex = 1, gbi%numBins2d
          if (binCount(binIndex) /= 0 ) then
            binMean(binIndex) = binMean(binIndex) / real(binCount(binIndex),8)
          else
            binMean(binIndex) = 0.d0
          end if
        end do

        !- Distribute the bin mean values at each grid point
        do latIndex= myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            binIndex = int(bin2d(lonIndex,latIndex,1))
            if (binIndex /= -1) then
              mean4d(lonIndex,latIndex,varLevIndex,stepIndex) = binMean(binIndex)
            else
              mean4d(lonIndex,latIndex,varLevIndex,stepIndex) = 0.d0
            end if
          end do
        end do

      end do
    end do

    if (mpi_myid == 0 .and. verbose) then
      write(*,*)
      write(*,*) 'gbi_mean_gsv: Done!'
    end if

  end subroutine gbi_mean_gsv

  !--------------------------------------------------------------------------
  ! gbi_stdDev_ens
  !--------------------------------------------------------------------------
  subroutine gbi_stdDev_ens(gbi, ens, statevector)
    implicit none

    type(struct_gbi) :: gbi
    type(struct_ens) :: ens
    type(struct_gsv) :: statevector

    integer :: myBinCount(gbi%numBins2d)
    integer :: binCount  (gbi%numBins2d)
    real(8) :: myBinSum  (gbi%numBins2d)
    real(8) :: binStdDev (gbi%numBins2d)

    real(4), pointer :: ptr4d_r4 (:,:,:,:)
    real(8), pointer :: stdDev_r8(:,:,:,:)

    real(4), pointer :: bin2d(:,:,:)

    integer :: myLonBeg, myLonEnd
    integer :: myLatBeg, myLatEnd
    integer :: latIndex, lonIndex, varLevIndex, stepIndex, binIndex, memberIndex
    integer :: nVarLev, nStep, nEns, ier

    character(len=4), pointer :: varNamesList(:)

    if (mpi_myid == 0 .and. verbose) then
      write(*,*)
      write(*,*) 'gbi_stdDev_ens: Starting...'
    end if

    if (.not. gsv_allocated(statevector)) then
      nullify(varNamesList)
      call ens_varNamesList(varNamesList,ens) 
      call gsv_allocate(statevector, ens_getNumStep(ens),               &
           ens_getHco(ens), ens_getVco(ens), varNames_opt=varNamesList, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,      &
           dataKind_opt=8 )
    end if

    call ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    nVarLev = ens_getNumK(ens)
    nStep   = ens_getNumStep(ens)
    nEns    = ens_getNumMembers(ens)

    call gsv_getField(statevector,stdDev_r8)
    call gsv_getField(gbi%statevector_bin2d,bin2d)

    do varLevIndex = 1, nVarLev

      ptr4d_r4 => ens_getOneLev_r4(ens,varLevIndex)

      do stepIndex = 1, nStep

        !- Sum the square values per bin
        myBinCount(:) = 0
        myBinSum  (:) = 0.d0
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            binIndex = int(bin2d(lonIndex,latIndex,1))
            if (binIndex /= -1) then
              do memberIndex = 1, nEns
                myBinCount(binIndex) = myBinCount(binIndex) + 1
                myBinSum(binIndex)   = myBinSum(binIndex)   + real(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex),8)**2
              end do
            end if
          end do
        end do

        !- Compute the stdDev per bin
        call rpn_comm_allreduce(myBinCount,binCount ,gbi%numBins2d,"MPI_INTEGER"         ,"MPI_SUM","GRID",ier)
        call rpn_comm_allreduce(myBinSum  ,binStdDev,gbi%numBins2d,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ier)

        do binIndex = 1, gbi%numBins2d
          if (binCount(binIndex) /= 0 ) then
            binStdDev(binIndex) = sqrt(binStdDev(binIndex) / real(binCount(binIndex)-1,8))
          else
            binStdDev(binIndex) = 0.d0
          end if
        end do

        !- Distribute the bin stdDev values at each grid point
        do latIndex= myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            binIndex = int(bin2d(lonIndex,latIndex,1))
            if (binIndex /= -1) then
              stdDev_r8(lonIndex,latIndex,varLevIndex,stepIndex) = binStdDev(binIndex)
            else
              stdDev_r8(lonIndex,latIndex,varLevIndex,stepIndex) = 0.d0
            end if
          end do
        end do

      end do
    end do

    if (mpi_myid == 0 .and. verbose) then
      write(*,*)
      write(*,*) 'gbi_stdDev_ens: Done!'
    end if
    
  end subroutine gbi_stdDev_ens

end module gridBinning_mod
