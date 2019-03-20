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
!! MODULE enkf (prefix="enkf" category='1. High-level functionality')
!!
!! *Purpose*: Implementation of the EnKF in MIDAS.
!!
!--------------------------------------------------------------------------
MODULE enkf_mod
  use mpi_mod
  use gridStateVector_mod
  use mathPhysConstants_mod
  use utilities_mod
  use fileNames_mod
  use varNameList_mod
  use tt2phi_mod
  use obsSpaceData_mod
  use columnData_mod
  implicit none
  save
  private

  ! public procedures
  public :: enkf_computeColumnsMean, enkf_computeColumnsPerturbations
  public :: enkf_extractObsRealBodyColumn, enkf_extractObsIntBodyColumn
  public :: enkf_gatherHX

  integer, external :: get_max_rss

contains

  subroutine enkf_computeColumnsMean(column_mean, columns)
    implicit none

    ! arguments
    type(struct_columnData) :: column_mean, columns(:)

    ! locals
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8) :: multFactor
    real(8), pointer :: column_ptr(:)

    call tmg_start(146,'ENKF_COLSMEAN')

    nEns = size(columns)
    multFactor = 1.0d0 / real(nEns,8)
    write(*,*) 'enkf_computeColumnsMean: nEns =', nEns

    call col_zero(column_mean)

    do memberIndex = 1, nEns

        column_mean%all(:,:) = column_mean%all(:,:) +  &
                               multFactor * columns(memberIndex)%all(:,:)

        column_mean%gz_T(:,:) = column_mean%gz_T(:,:) +  &
                                multFactor * columns(memberIndex)%gz_T(:,:)
        column_mean%gz_M(:,:) = column_mean%gz_M(:,:) +  &
                                multFactor * columns(memberIndex)%gz_M(:,:)
        column_mean%gz_sfc(:,:) = column_mean%gz_sfc(:,:) +  &
                                multFactor * columns(memberIndex)%gz_sfc(:,:)

        column_mean%pressure_T(:,:) = column_mean%pressure_T(:,:) +  &
                                      multFactor * columns(memberIndex)%pressure_T(:,:)
        column_mean%pressure_M(:,:) = column_mean%pressure_M(:,:) +  &
                                      multFactor * columns(memberIndex)%pressure_M(:,:)

        column_mean%dP_dPsfc_T(:,:) = column_mean%dP_dPsfc_T(:,:) +  &
                                      multFactor * columns(memberIndex)%dP_dPsfc_T(:,:)
        column_mean%dP_dPsfc_M(:,:) = column_mean%dP_dPsfc_M(:,:) +  &
                                      multFactor * columns(memberIndex)%dP_dPsfc_M(:,:)

    end do

    !if (col_varExist('P0')) then
    !  call col_calcPressure(column_mean)
    !end if

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of column_mean:'
      column_ptr => col_getColumn(column_mean,1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(column_mean,levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(column_mean,levIndex,1,'MM')
      end do
    end if

    call tmg_stop(146)

  end subroutine enkf_computeColumnsMean


  subroutine enkf_computeColumnsPerturbations(columns, column_mean)
    implicit none

    ! arguments
    type(struct_columnData) :: columns(:), column_mean

    ! locals
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8), pointer :: column_ptr(:)

    call tmg_start(147,'ENKF_COLSPERTS')

    nEns = size(columns)
    write(*,*) 'enkf_computeColumnsPerturbations: nEns =', nEns

    !
    ! Remove ensemble mean from all variables, except: gz_sfc, dP_dPsfc_T/M, oltv
    !
    do memberIndex = 1, nEns

        columns(memberIndex)%all(:,:) = columns(memberIndex)%all(:,:) -  &
                                        column_mean%all(:,:)

        columns(memberIndex)%gz_T(:,:) = columns(memberIndex)%gz_T(:,:) -  &
                                         column_mean%gz_T(:,:)
        columns(memberIndex)%gz_M(:,:) = columns(memberIndex)%gz_M(:,:) -  &
                                         column_mean%gz_M(:,:)

        columns(memberIndex)%pressure_T(:,:) = columns(memberIndex)%pressure_T(:,:) -  &
                                      column_mean%pressure_T(:,:)
        columns(memberIndex)%pressure_M(:,:) = columns(memberIndex)%pressure_M(:,:) -  &
                                      column_mean%pressure_M(:,:)

    end do

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of columns(1):'
      column_ptr => col_getColumn(columns(1),1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(columns(1),levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(columns(1),levIndex,1,'MM')
      end do
    end if

    call tmg_stop(147)

  end subroutine enkf_computeColumnsPerturbations


  subroutine enkf_extractObsRealBodyColumn(outputVector, obsSpaceData, obsColumnIndex)
    implicit none

    ! arguments
    real(8)          :: outputVector(:)
    type(struct_obs) :: obsSpaceData
    integer          :: obsColumnIndex

    ! locals
    integer :: bodyIndex

    call tmg_start(148,'ENKF_EXTRACTBODY')

    do bodyIndex = 1, obs_numBody(obsSpaceData)
      outputVector(bodyIndex) = obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex)
    end do

    call tmg_stop(148)

  end subroutine enkf_extractObsRealBodyColumn


  subroutine enkf_extractObsIntBodyColumn(outputVector, obsSpaceData, obsColumnIndex)
    implicit none

    ! arguments
    integer          :: outputVector(:)
    type(struct_obs) :: obsSpaceData
    integer          :: obsColumnIndex

    ! locals
    integer :: bodyIndex

    call tmg_start(148,'ENKF_EXTRACTBODY')

    do bodyIndex = 1, obs_numBody(obsSpaceData)
      outputVector(bodyIndex) = obs_bodyElem_i(obsSpaceData,obsColumnIndex,bodyIndex)
    end do

    call tmg_stop(148)

  end subroutine enkf_extractObsIntBodyColumn


  subroutine enkf_gatherHX(HXens,HXensT_mpiglobal)
    implicit none

    ! arguments
    real(8) :: HXens(:,:)
    real(8),pointer :: HXensT_mpiglobal(:,:)

    ! locals
    integer :: ierr, nEns, numBody, procIndex, memberIndex, numBody_mpiglobal
    integer :: allNumBody(mpi_nprocs), displs(mpi_nprocs)

    numBody = size(HXens,1)
    nEns     = size(HXens,2)

    write(*,*) 'enkf_gatherHX: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call rpn_comm_gather( numBody, 1, 'mpi_integer', allNumBody, 1, 'mpi_integer', &
                          0, 'GRID', ierr )
    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allNumBody(procIndex-1)
      end do
    else
      displs(:) = 0
    end if

    numBody_mpiglobal = sum(allNumBody(:))
    if( mpi_myid == 0 ) then
      allocate(HXensT_mpiglobal(nEns,numBody_mpiglobal))
    else
      allocate(HXensT_mpiglobal(nEns,1))
    end if

    do memberIndex = 1, nEns
      call rpn_comm_gatherv( HXens(:,memberIndex), numBody, 'mpi_double_precision', &
                             HXensT_mpiglobal(memberIndex,:), allNumBody, displs, &
                             'mpi_double_precision', 0, 'GRID', ierr )
    end do

    write(*,*) 'enkf_gatherHX: finished'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_gatherHX


end module enkf_mod
