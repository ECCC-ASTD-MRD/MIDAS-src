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
!! MODULE increment (prefix="inc")
!!
!! *Purpose*: To add a 4D increment to a given 4D background/reference state
!!            and to output the results
!!
!--------------------------------------------------------------------------
MODULE increment_mod
  use mpivar_mod
  use timeCoord_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use humidityLimits_mod
  use utilities_mod
  use variableTransforms_mod
  use BMatrix_mod
  use chem_postproc_mod
  implicit none
  save
  private

  ! public procedures
  public :: inc_computeAndWriteAnalysis, inc_getIncrement, inc_writeIncrement

  ! namelist variables
  integer  :: writeNumBits
  logical  :: writeHiresIncrement
  logical  :: imposeRttovHuLimits
  character(len=12) :: etiket_anlm, etiket_rehm, etiket_rebm
  character(len=12) :: hInterpolationDegree

  NAMELIST /NAMINC/ writeHiresIncrement, etiket_rehm, etiket_anlm, &
       etiket_rebm, writeNumBits, imposeRttovHuLimits, hInterpolationDegree

CONTAINS

  !--------------------------------------------------------------------------
  ! inc_computeAndWriteAnalysis
  !--------------------------------------------------------------------------
  subroutine inc_computeAndWriteAnalysis(statevector_incLowRes_opt)
    implicit none

    type(struct_gsv), intent(in), optional :: statevector_incLowRes_opt

    type(struct_gsv) :: statevector_incHighRes
    type(struct_gsv) :: statevector_trial, statevector_analysis
    type(struct_gsv) :: statevector_PsfcLowRes, statevector_Psfc
    type(struct_gsv) :: statevector_PsfcLowRes_varsLevs, statevector_Psfc_varsLevs

    type(struct_vco), pointer :: vco_trl => null()
    type(struct_vco), pointer :: vco_inc => null()
    type(struct_hco), pointer :: hco_trl => null()

    integer              :: fclos, fnom, fstopc, ierr
    integer              :: stepIndex, numStep, middleStep
    integer              :: nulnam, dateStamp
    integer, allocatable :: dateStampList(:)
    integer              :: get_max_rss

    character(len=256)  :: trialFileName, incFileName, anlFileName
    character(len=4)    :: coffset

    real(8)             :: deltaHours
    real(8), pointer    :: PsfcTrial(:,:,:,:), PsfcAnalysis(:,:,:,:)
    real(8), pointer    :: PsfcIncrement(:,:,:,:)
    real(8), pointer    :: PsfcIncLowResFrom3Dgsv(:,:,:,:), PsfcIncLowRes(:,:,:,:)
    real(8), pointer    :: GZsfc_increment(:,:), GZsfc_trial(:,:)

    logical  :: allocGZsfc, writeGZsfc, useIncLevelsOnly

    !
    !- Set/Read values for the namelist NAMINC
    !

    !- Setting default values
    writeHiresIncrement = .true.
    imposeRttovHuLimits = .true.
    etiket_rehm = 'INCREMENT'
    etiket_rebm = 'INCREMENT'
    etiket_anlm = 'ANALYSIS'
    writeNumBits = 16
    hInterpolationDegree = 'LINEAR'

    !- Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=naminc, iostat=ierr)
    if ( ierr /= 0) call utl_abort('inc_computeAndWriteAnalysis: Error reading namelist')
    if ( mpi_myid == 0 ) write(*,nml=naminc)
    ierr = fclos(nulnam)

    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    ! Setup timeCoord module (date read from trial file)
    numStep = tim_nstepobsinc
    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

    !- Initialize the trial state grid
    if (mpi_myid == 0) write(*,*) ''
    if (mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: Set hco parameters for trial grid'
    trialFileName = './trlm_01'
    call hco_setupFromFile( hco_trl, trim(trialFileName), ' ')
    if ( mpi_myid == 0 ) then
      call vco_setupFromFile( vco_trl, trim(trialFileName) )
    end if
    call vco_mpiBcast(vco_trl)

    !- Do we need to read all the vertical levels from the trial fields?
    if (present(statevector_incLowRes_opt)) then
      vco_inc => statevector_incLowRes_opt%vco
    else
      call difdatr(datestamplist(1),tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      incFileName = './rebm_' // trim(coffset) // 'm'

      if ( mpi_myid == 0 ) then
        call vco_setupFromFile(vco_inc, incFileName)
      end if
      call vco_mpiBcast(vco_inc)
    end if

    useIncLevelsOnly = vco_subsetOrNot(vco_inc, vco_trl)
    if ( useIncLevelsOnly ) then
      ! Read only the increment levels
      write(*,*)
      write(*,*) 'inc_computeAndWriteAnalysis: only the increment levels will be read in the trials'
      call  vco_deallocate(vco_trl)
      vco_trl => vco_inc
    else
      ! Read them all
      write(*,*)
      write(*,*) 'inc_computeAndWriteAnalysis: all the vertical levels will be read in the trials'
      if (.not. present(statevector_incLowRes_opt)) then
        call vco_deallocate(vco_inc)
      end if
    end if

    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    ! Read trial files
    !
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading background state for all time steps'
    call gsv_allocate(statevector_trial, tim_nstepobsinc, hco_trl, vco_trl,   &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocGZsfc_opt=.true., hInterpolateDegree_opt=hInterpolationDegree)
    call gsv_readTrials(statevector_trial)

    !
    !- Get the increment of Psfc
    !
    if (gsv_varExist(varName='P0')) then
      call gsv_allocate(statevector_Psfc, numStep, hco_trl, vco_trl, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        varNames_opt=(/'P0'/), allocGZsfc_opt=.true., &
                        hInterpolateDegree_opt=hInterpolationDegree)

      if (present(statevector_incLowRes_opt)) then
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: horiz interpolation of the Psfc increment'

        ! Extract Psfc inc at low resolution
        call gsv_allocate(statevector_PsfcLowRes, numStep, statevector_incLowRes_opt%hco, &
                          statevector_incLowRes_opt%vco,  &
                          dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                          varNames_opt=(/'P0'/))
        PsfcIncLowRes          => gsv_getField_r8(statevector_PsfcLowRes,'P0')
        PsfcIncLowResFrom3Dgsv => gsv_getField_r8(statevector_incLowRes_opt,'P0')
        PsfcIncLowRes(:,:,1,:) = PsfcIncLowResFrom3Dgsv(:,:,1,:)

        ! Interpolate
        call gsv_interpolate(statevector_PsfcLowRes,statevector_Psfc)
        call gsv_deallocate(statevector_PsfcLowRes)

      else
        !- Read from file
        do stepIndex = 1, numStep
          dateStamp = datestamplist(stepIndex)
          if(mpi_myid == 0) write(*,*) ''
          if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading Psfc increment for time step: ',stepIndex, dateStamp

          call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
          if(nint(deltaHours*60.0d0).lt.0) then
            write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
          else
            write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
          endif
          incFileName = './rebm_' // trim(coffset) // 'm'

          call gsv_readFromFile( statevector_Psfc, trim(incFileName), ' ', ' ', stepIndex,  &
                                 containsFullField_opt=.false. )
        end do

      end if

      !
      !- Compute analysis Psfc to use for interpolation of increment
      !
      PsfcTrial     => gsv_getField_r8(statevector_trial,'P0')
      PsfcIncrement => gsv_getField_r8(statevector_Psfc ,'P0')
      PsfcAnalysis  => gsv_getField_r8(statevector_Psfc ,'P0')

      PsfcAnalysis(:,:,1,:) = PsfcTrial(:,:,1,:) + PsfcIncrement(:,:,1,:)

      !
      !- Copy the surface GZ from trial into statevector_Psfc
      !
      GZsfc_increment => gsv_getGZsfc(statevector_Psfc)
      GZsfc_trial     => gsv_getGZsfc(statevector_trial)
      GZsfc_increment(:,:) = GZsfc_trial(:,:)
      allocGZsfc = .true.
    else
      allocGZsfc = .false.
    end if
    writeGZsfc = allocGZsfc

    !
    !- Compute the analysis
    !
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: compute the analysis'

    call gsv_allocate(statevector_incHighRes, numStep, hco_trl, vco_trl, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocGZsfc_opt=allocGZsfc, hInterpolateDegree_opt=hInterpolationDegree)
    call gsv_allocate(statevector_analysis, numStep, hco_trl, vco_trl, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocGZsfc_opt=allocGZsfc, hInterpolateDegree_opt=hInterpolationDegree)
    call gsv_copy(statevector_trial, statevector_analysis)

    if (present(statevector_incLowRes_opt)) then
      !- Interpolate and add the input increments
      if (gsv_varExist(varName='P0')) then
        call gsv_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis,&
                                   PsfcReference_opt=PsfcAnalysis(:,:,1,:))
      else
        call gsv_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis)
      end if
    else
      !- Read the increments from files
      do stepIndex = 1, numStep
        dateStamp = datestamplist(stepIndex)
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading increment for time step: ',stepIndex, dateStamp

        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        endif
        incFileName = './rebm_' // trim(coffset) // 'm'

        write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
        if (gsv_varExist(varName='P0')) then
          call gsv_readFromFile( statevector_incHighRes, trim(incFileName), ' ', ' ', stepIndex,  &
                                 PsfcReference_opt=PsfcAnalysis(:,:,1,stepIndex),  &
                                 containsFullField_opt=.false. )
        else
          call gsv_readFromFile( statevector_incHighRes, trim(incFileName), ' ', ' ', stepIndex, &
                                 containsFullField_opt=.false. )
        end if
      end do

      call gsv_add(statevector_incHighRes, statevector_analysis)

    end if

    !
    !- Impose limits on humidity analysis and recompute increment
    !
    write(*,*) 'inc_computeAndWriteAnalysis: calling qlim_gsvSaturationLimit'
    call qlim_gsvSaturationLimit(statevector_analysis)
    if( imposeRttovHuLimits ) call qlim_gsvRttovLimit(statevector_analysis)
    call gsv_copy(statevector_analysis, statevector_incHighRes)
    call gsv_add(statevector_trial, statevector_incHighRes, -1.0d0)

    !
    !- Write interpolated increment files
    !
    if( writeHiresIncrement ) then
      do stepIndex = 1, numStep
        dateStamp = datestamplist(stepIndex)
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: writing interpolated increment for time step: ',stepIndex, dateStamp

        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        endif
        incFileName = './rehm_' // trim(coffset) // 'm'

        write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
        call gsv_writeToFile( statevector_incHighRes, trim(incFileName), etiket_rehm,  &
                              stepIndex_opt=stepIndex, typvar_opt='R', numBits_opt=writeNumBits, &
                              containsFullField_opt=.false. )

        if (gsv_varExist(varName='P0')) then
          ! Also write analysis value of Psfc and surface GZ to increment file
          call gsv_writeToFile( statevector_Psfc, trim(incFileName), etiket_rehm,  &
                                stepIndex_opt=stepIndex, typvar_opt='A', writeGZsfc_opt=.true., &
                                numBits_opt=writeNumBits, containsFullField_opt=.false. )
        end if
      end do
    end if

    !
    !- Write analysis state to file
    !
    do stepIndex = 1, numStep
      dateStamp = datestamplist(stepIndex)
      if(mpi_myid == 0) write(*,*) ''
      if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: writing analysis for time step: ',stepIndex, dateStamp

      call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      anlFileName = './anlm_' // trim(coffset) // 'm'

      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
      call gsv_writeToFile( statevector_analysis, trim(anlFileName), etiket_anlm,  &
                            stepIndex_opt=stepIndex, typvar_opt='A', writeGZsfc_opt=writeGZsfc, &
                            numBits_opt=writeNumBits, containsFullField_opt=.true. )
    end do

    call gsv_deallocate(statevector_analysis)
    if (gsv_varExist(varName='P0')) then
      call gsv_deallocate(statevector_Psfc)
    end if
    call gsv_deallocate(statevector_incHighRes)
    call gsv_deallocate(statevector_trial)

  end subroutine inc_computeAndWriteAnalysis

  subroutine inc_getIncrement(incr_cv,statevector_incr,nvadim_mpilocal)

    implicit none

    ! arguments
    real(8) :: incr_cv(:)
    type(struct_gsv) :: statevector_incr
    integer :: nvadim_mpilocal

    ! compute increment from control vector (multiply by B^1/2)
    call bmat_sqrtB(incr_cv, nvadim_mpilocal, statevector_incr)

    ! Compute new diagnotics based on NAMSTATE
    if ( gsv_varExist(statevector_incr,'QR') .and. gsv_varExist(statevector_incr,'DD') ) then
       write(*,*)
       write(*,*) 'User is asking for Vort-Div analysis increment'
       call vtr_transform( statevector_incr, & ! INOUT
                           'UVtoVortDiv' )     ! IN
       if ( gsv_varExist(statevector_incr,'PP') .and. gsv_varExist(statevector_incr,'CC') ) then
          write(*,*)
          write(*,*) 'User is asking for Psi-Chi analysis increment'
          call vtr_transform( statevector_incr, & ! INOUT
                              'VortDivToPsiChi')  ! IN
       end if
    end if

    ! Adjust and or transform chemical consituent concentration increments as needed.
    ! This includes ensuring non-negative analysis values on the analysis/increment grid.
    if (gsv_varKindExist('CH')) call chm_transform_final_increments(statevector_incr)

  end subroutine inc_getIncrement

  subroutine inc_writeIncrement(statevector_incr)

    implicit none

    ! arguments
    type(struct_gsv)     :: statevector_incr
    ! locals
    integer              :: stepIndex, dateStamp, nulnam
    integer              :: fclos, fnom, ierr
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=2)     :: flnum
    character(len=30)    :: fileName

    if(mpi_myid == 0) write(*,*) 'inc_writeIncrement: STARTING'

    etiket_rebm = 'INCREMENT'

    !- Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=naminc, iostat=ierr)
    if ( ierr /= 0) call utl_abort('inc_computeAndWriteAnalysis: Error reading namelist')
    if ( mpi_myid == 0 ) write(*,nml=naminc)
    ierr = fclos(nulnam)

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc

      dateStamp = gsv_getDateStamp(statevector_incr,stepIndex)
      if(mpi_myid == 0) write(*,*) 'inc_writeIncrement: writing increment for time step: ',stepIndex, dateStamp

      ! write the increment file for this time step
      call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      fileName = './rebm_' // trim(coffset) // 'm'
      call gsv_writeToFile( statevector_incr, fileName, etiket_rebm, 1.0d0, 0,  &
                            stepIndex, containsFullField_opt=.false. )

    enddo

  end subroutine inc_writeIncrement

END MODULE increment_mod