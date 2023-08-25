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

module var1DIdealize_mod
    ! MODULE var1DIdealize_mod (prefix='var1D' category='4. Data Object transformations')
    !
    ! :Purpose: contains all 1Dvar-related methods.
    !
    use columnData_mod
    use columnVariableTransforms_mod
    use controlVector_mod
    use gridStatevector_mod
    use horizontalCoord_mod
    use midasMpi_mod 
    use obsSpaceData_mod
    use timeCoord_mod
    use utilities_mod
    use verticalCoord_mod
    use codeprecision_mod
    use mathphysconstants_mod
    use randomNumber_mod
    use bMatrix1Dvar_mod
    use innovation_mod
    use var1D_mod
    use gridStateVectorFileIO_mod
    use interpolation_mod
    use varNameList_mod
    use increment_mod
    use rMatrix_mod
    use humidityLimits_mod
    use obsoperators_mod
    use obsErrors_mod
  
    implicit none
    save
    private
  
    ! public procedures
    public :: var1DIdealize_simulateBackgroundState

    ! namelist variables
    integer  :: writeNumBits
    logical  :: writeHiresIncrement
    logical  :: imposeRttovHuLimits, useAnalIncMask
    character(len=12) :: etiket_anlm, etiket_rehm, etiket_rebm
    character(len=12) :: hInterpolationDegree
    logical :: applyLiebmann
    logical :: SSTSpread  
    integer :: SSTSpreadMaxBoxSize
    character(len=10) :: SSTSubgrid

  contains

  !--------------------------------------------------------------------------
  ! var1DIdealize_simulateBackgroundState
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_simulateBackgroundState(columnTruthOnTrlLev, columnSimTrlOnTrlLev, &
                                                   obsSpaceData, vco_anl, seed)

    !
    !:Purpose: Simulate the background state by adding a perturbation from the reference state (Truth)
    !
    implicit none

    ! arguments
    type(struct_columnData), target, intent(inout) :: columnTruthOnTrlLev
    type(struct_columnData), target, intent(out)   :: columnSimTrlOnTrlLev
    type(struct_obs),                intent(in)    :: obsSpaceData
    type(struct_vco), pointer,       intent(in)    :: vco_anl
    integer,                         intent(in)    :: seed
    
    ! locals:
    type(struct_columnData), target :: columnPertOnAnLev
    type(struct_columnData), target :: columnTruthOnAnlLev
    type(struct_columnData), target :: columnPertOnTrlLev
    real(8), allocatable            :: controlVector(:)
    integer                         :: cvIndex
    type(struct_gsv)                :: stateVectorPertOnAnLevTruth
    type(struct_gsv)                :: stateVectorPertOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevSim
    type(struct_gsv)                :: stateVectorTrlOnAnlLevTruth
    character(len=50)               :: prefixFileName
    logical                         :: containsFullField

    allocate(controlVector(cvm_nvadim))
    ! Generate perturbation sampling following gaussian distribution with zero mean and one std
    call rng_setup(abs(seed))
    do cvIndex = 1, cvm_nvadim
      controlVector(cvIndex) = rng_gaussian()
    end do

    ! Compute (B^1/2)*Pert (column)
    call col_setVco(columnPertOnAnLev, vco_anl)
    call col_allocate(columnPertOnAnLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData)

    call col_setVco(columnPertOnTrlLev, col_getVco(columnTruthOnTrlLev))
    call col_allocate(columnPertOnTrlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    
    ! Interpolate (B^1/2)*Pert from analysis to trial level
    call var1DIdealize_vInterpPertAnLev2TrlLev(columnPertOnAnLev, columnPertOnTrlLev, columnTruthOnTrlLev)
    
    call col_setVco(columnSimTrlOnTrlLev, col_getVco(columnTruthOnTrlLev))
    call col_allocate(columnSimTrlOnTrlLev, col_getNumCol(columnTruthOnTrlLev), &
                      setToZero_opt=.true.)
    call col_copy(columnTruthOnTrlLev, columnSimTrlOnTrlLev)

    ! Add the truth and (B^1/2)*Pert columns
    call col_add(columnPertOnTrlLev, columnSimTrlOnTrlLev)

    ! Compute the pressure levels
    call cvt_transform(columnSimTrlOnTrlLev, 'ZandP_nl')

    ! Restrict the simulated humidity background within physically reasonable values.
    call qlim_rttovLimit(columnSimTrlOnTrlLev)

    ! Interpolate the truth from trial to analysis increment levels
    call col_setVco(columnTruthOnAnlLev, vco_anl)
    call col_allocate(columnTruthOnAnlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    call inn_setupColumnsOnAnlIncLev(columnTruthOnTrlLev, columnTruthOnAnlLev)

    ! Write trial into standard files
    prefixFileName = 'SimTrialOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevSim, obsSpaceData, columnSimTrlOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnTrlLevSim, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevTruth, obsSpaceData, columnTruthOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnTrlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnAnlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnAnlLevTruth, obsSpaceData, columnTruthOnAnlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnAnlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'PertOnTrlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnTrlLevTruth, obsSpaceData, columnPertOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorPertOnTrlLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    prefixFileName = 'PertOnAnlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnAnLevTruth, obsSpaceData, columnPertOnAnLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorPertOnAnLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    if (mmpi_myId ==0) then
      call gsv_deallocate(stateVectorPertOnTrlLevTruth)
      call gsv_deallocate(stateVectorPertOnAnLevTruth)
      call gsv_deallocate(stateVectorTrlOnAnlLevTruth)
    end if

    call col_deallocate(columnPertOnTrlLev)
    call col_deallocate(columnPertOnAnLev)
    call col_deallocate(columnTruthOnAnlLev)
    deallocate(controlVector)

  end subroutine var1DIdealize_simulateBackgroundState

  !--------------------------------------------------------------------------
  ! var1DIdealize_vInterpPertAnLev2TrlLev
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_vInterpPertAnLev2TrlLev(columnAnlLev, columnTrlLev, columnPresRef)
    !
    ! :Purpose: Vertically Interpolate the generated perturbation from analysis to trial level. 
    !           An reference column on trial level is required to compute the pressure level, which 
    !           is used for the vertical interpolation
    !
    implicit none

    ! arguments
    type(struct_columnData), intent(in)     :: columnAnlLev  ! Column data in analysis level
    type(struct_columnData), intent(inout)  :: columnTrlLev  ! Column data in trial level
    type(struct_columnData), intent(in)     :: columnPresRef ! Column data where sfc pressure variables 
                                                             ! will be used for vertical interpolation

    ! locals:
    integer                    :: numColumns
    integer                    :: columnIndex, varIndex
    real(8), allocatable       :: pSfcRef(:,:)
    real(8), pointer           :: columnAnlLev_ptr(:), columnTrlLev_ptr(:)

    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Starting'

    ! Check the column size
    if (.not. (col_getNumCol(columnAnlLev) == col_getNumCol(columnTrlLev) .and.    &
        col_getNumCol(columnAnlLev) == col_getNumCol(columnPresRef))) then
      write(*,*) 'Column size columnAnlLev, columnTrlLev and columnPresRef', col_getNumCol(columnAnlLev), &
                  col_getNumCol(columnTrlLev), col_getNumCol(columnPresRef)
      call utl_abort('var1DIdealize_vInterpPertAnLev2TrlLev: The columnAnlLev, columnTrlLev and columnPresRef &
                                 do not have equal number of columns')
    end if

    numColumns = col_getNumCol(columnAnlLev)
    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Column size', numColumns

    ! Extract the surface pressure from the columnPresRef
    allocate(pSfcRef(1,numColumns))
    do columnIndex = 1, numColumns
      pSfcRef(1, columnIndex) = col_getElem(columnPresRef, 1, columnIndex, 'P0')
    end do

    ! Vertical Interpolation
    do varIndex = 1, vnl_numvarmax3D
      ! Check if varName is an analysis variable
      if (.not. varneed(vnl_varNameList3D(varIndex))) cycle
      if (.not. col_varExist(columnAnlLev, vnl_varNameList3D(varIndex)) ) cycle
      call int_vInterp_col(columnAnlLev, columnTrlLev, vnl_varNameList3D(varIndex), sfcPressureRef_opt=pSfcRef)
    end do
  
    ! Copy 2D surface variables
    do varIndex = 1, vnl_numvarmax2D
      if (.not. varneed(vnl_varNameList2D(varIndex))) cycle
      if (.not. col_varExist(columnAnlLev, vnl_varNameList2D(varIndex))) cycle
      if (col_getNumCol(columnAnlLev) > 0) then       
        do columnIndex = 1, col_getNumCol(columnAnlLev)
          columnTrlLev_ptr  => col_getColumn(columnTrlLev , columnIndex, vnl_varNameList2D(varIndex))
          columnAnlLev_ptr => col_getColumn(columnAnlLev, columnIndex, vnl_varNameList2D(varIndex))
          columnTrlLev_ptr(:) = columnAnlLev_ptr(:)
        end do
      end if
    end do

    deallocate(pSfcRef)
    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Finished'
    contains

    logical function varneed(varName)
      implicit none
      ! Arguements: 
      character(len=*) :: varName ! Variable Name

      ! Locals:
      integer          :: varIndex2

      varneed=.false.
      do varIndex2=1,VNL_NUMVARMAX
        if (trim(varName) == trim(bmat1D_includeAnlVar(varIndex2))) then
          varneed=.true.
       end if
      end do

    end function varneed
  end subroutine var1DIdealize_vInterpPertAnLev2TrlLev

  !--------------------------------------------------------------------------
  ! var1DIdealize_writeSimTrial
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_writeSimTrial(statevectorSim, prefixFileName, etiket, containsFullField)
    !
    ! :Purpose: Write the simulate background state from statevector strucure (1Dvar case) 
    !           into output standard file
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)     :: statevectorSim     ! Statevector to be written in file
    character(len=*), intent(in)     :: prefixFileName     ! Prefix of the filename
    character(len=*), intent(in)     :: etiket             ! Etiket of the filename
    logical,          intent(in)     :: containsFullField  ! Logical for full field values or Perturbation/Increments

    ! Locals:
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=100)   :: fileName
   
    if(mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: STARTING'

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevectorSim)) then
        dateStamp = gsv_getDateStamp(statevectorSim,stepIndex)
        if (mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: writing for time step: ',stepIndex, dateStamp

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if (nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if

        fileName = './'//trim(prefixFileName)//'_' // trim(coffset) // 'm'
        call gio_writeToFile( statevectorSim, fileName, trim(etiket), scaleFactor_opt = 1.0d0, &
                              ip3_opt = 0, stepIndex_opt = stepIndex, containsFullField_opt=containsFullField )
      end if
    end do

    if(mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: Finished'
  end subroutine var1DIdealize_writeSimTrial
end module var1DIdealize_mod

  

   
