!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!! MODULE stateToColumn (prefix="s2c")
!!
!! *Purpose*: Tangent-linear and adjoint of bilinear interpolation to 
!!            between a gridStateVector object and a columnData object.
!!            This transformation has been optimized to improve MPI load
!!            balancing when applying the observation operators.
!!            (Replaces oda_L and oda_LT)
!!
!--------------------------------------------------------------------------
module stateToColumn_mod
  use MathPhysConstants_mod
  use mpivar_mod 
  use gridstatevector_mod
  use obsSpaceData_mod
  use columnData_mod
  use timeCoord_mod
  use tt2phi_mod
  use utilities_mod
  
  implicit none
  save
  private
  
  ! public routines
  public :: s2c_tl, s2c_ad, s2c_column_hbilin

  ! private module variables
  logical :: initialized = .false.
  integer :: numHeaderTile, numHeaderColumn
  integer, allocatable :: PE_sender_mpiglobal(:), PE_receiver_mpiglobal(:)
  integer, allocatable :: myPESourceColumn(:)
  real(8), allocatable :: myYposTile(:), myXposTile(:)
  real(8), allocatable :: myLatTile(:), myLonTile(:)
  real(8), allocatable :: myTimeInterpWeight(:,:)
  real(8), pointer     :: fieldsWithHalo(:,:,:,:)


CONTAINS 

  !---------------------------------------------------------
  ! S2C_SETUP
  !---------------------------------------------------------
  subroutine s2c_setup(statevector,column,obsSpaceData)
    implicit none
    ! Purpose: gather information needed concerning observations that live on
    !          other mpi tasks, but are associated with the horizontal tile
    !          on the local mpi task
    type(struct_gsv)        :: statevector
    type(struct_columnData) :: column
    type(struct_obs)        :: obsSpaceData

    integer :: procIndex, procIndex2, headerIndex, headerIndex2, stepIndex, ierr, nsize, status
    integer :: numToSend, PE_sender, PE_receiver
    integer :: numToSend_mpiglobal(mpi_nprocs), numToSend_tmp(mpi_nprocs)
    real(8), allocatable :: myYposToSend(:)
    real(8), allocatable :: myXposToSend(:)
    real(8), allocatable :: myLatToSend(:)
    real(8), allocatable :: myLonToSend(:)
    real(8), allocatable :: myTimeInterpWeightToSend(:,:)
    real(8) :: Lat, Lon, ypos, xpos, LatRot, LonRot

    if(mpi_myid == 0) write(*,*) 's2c_setup: gathering information for load balancing'

    ! Allocate some module arrays
    allocate(PE_sender_mpiglobal(mpi_nprocs))
    allocate(PE_receiver_mpiglobal(mpi_nprocs))

    ! Number of headers on local columns
    numHeaderColumn = col_getNumCol(column)

    write(*,*) ' ' 
    write(*,*) 's2c_setup: Column-related values:'
    write(*,*) 's2c_setup: numHeaderColumn              = ', numHeaderColumn

    ! determine for each receiver TILE which COLUMN is the sender
    PE_sender = -1
    do headerIndex = 1, numHeaderColumn
      if(obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex) /=  &
         obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)) then
        PE_sender = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
      end if
    end do
    call rpn_comm_allgather(PE_sender, 1, 'mpi_integer',  &
                            PE_sender_mpiglobal, 1, 'mpi_integer', 'GRID', ierr)

    ! determine for each sender TILE which COLUMN is the receiver
    do procIndex = 0, mpi_nprocs-1
      PE_receiver_mpiglobal(procIndex+1) = -1
      do procIndex2 = 0, mpi_nprocs-1
        if(PE_sender_mpiglobal(procIndex2+1) == procIndex) then
          PE_receiver_mpiglobal(procIndex+1) = procIndex2
        end if
      end do
    end do
    PE_receiver = PE_receiver_mpiglobal(mpi_myid+1)

    numToSend_tmp(:) = 0
    do headerIndex = 1, numHeaderColumn
      if(obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex) /=  &
         obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)) then
        numToSend_tmp(PE_sender+1) = numToSend_tmp(PE_sender+1) + 1
      end if
    end do

    call rpn_comm_allreduce(numToSend_tmp, numToSend_mpiglobal, mpi_nprocs, &
                            'MPI_INTEGER', 'MPI_SUM', 'GRID', ierr)
    numToSend = numToSend_mpiglobal(mpi_myid+1)

    ! Number of headers on local TILE
    if(PE_sender == -1) then
      ! I am a sender from TILE to COLUMN
      numHeaderTile = numHeaderColumn + numToSend_mpiglobal(mpi_myid+1)
    else
      ! I am a receiver from TILE to COLUMN
      numHeaderTile = numHeaderColumn - numToSend_mpiglobal(PE_sender+1)
    end if

    if(mpi_myid == 0) write(*,*) 's2c_setup: numToSend_mpiglobal   = ',numToSend_mpiglobal(:)
    if(mpi_myid == 0) write(*,*) 's2c_setup: PE_sender_mpiglobal   = ',PE_sender_mpiglobal(:)
    if(mpi_myid == 0) write(*,*) 's2c_setup: PE_receiver_mpiglobal = ',PE_receiver_mpiglobal(:)

    write(*,*) ' ' 
    write(*,*) 's2c_setup: numHeaderTile              = ', numHeaderTile

    ! Only do the following if this mpi task is either sender or receiver
    if(numHeaderTile /= numHeaderColumn) then
      allocate(myPEsourceColumn(numHeaderColumn))
      do headerIndex = 1, numHeaderColumn
        myPEsourceColumn(headerIndex) = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
      end do

      if(numHeaderTile > 0) then
        allocate(myYPosTile(numHeaderTile))
        allocate(myXPosTile(numHeaderTile))
        allocate(myLatTile(numHeaderTile))
        allocate(myLonTile(numHeaderTile))
        allocate(myTimeInterpWeight(numHeaderTile,statevector%numStep))
      end if

      if(PE_sender == -1) then
        ! I am a sender from TILE to COLUMN:
      
        ! copy lat-lon that are already on correct mpi task
        headerIndex2 = 0
        do headerIndex = 1, numHeaderColumn
          if(obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex) ==  &
             obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)) then
            headerIndex2 = headerIndex2 + 1
            call col_getLatLon( column, headerIndex,                 & ! IN
                                Lat, Lon, ypos, xpos, LatRot, LonRot ) ! OUT
            myYPosTile(headerIndex2) = ypos
            myXPosTile(headerIndex2) = xpos
            myLatTile(headerIndex2)  = lat
            myLonTile(headerIndex2)  = lon

            if(btest(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex),5) ) then
              myTimeInterpWeight(headerIndex2,:) = 0.0d0
            else
              do stepIndex = 1, statevector%numStep
                myTimeInterpWeight(headerIndex2,stepIndex) = &
                  tim_getTimeInterpWeight(headerIndex,stepIndex)
              end do
            end if
          end if
        end do

        write(*,*) 's2c_setup: headerIndex2, numToSend_mpiglobal(mpi_myid+1)=',headerIndex2, numToSend_mpiglobal(mpi_myid+1)

        ! receive lat-lon, etc from receiver
        nsize = numHeaderTile-numHeaderColumn
        if(nsize > 0) then
          call rpn_comm_recv(myYPosTile(numHeaderColumn+1:numHeaderTile),nsize, &
                            'mpi_double_precision',PE_receiver,PE_receiver*2000+mpi_myid,  &
                            'GRID',status,ierr)
          call rpn_comm_recv(myXPosTile(numHeaderColumn+1:numHeaderTile),nsize, &
                            'mpi_double_precision',PE_receiver,PE_receiver*2000+mpi_myid,  &
                            'GRID',status,ierr)
          call rpn_comm_recv(myLatTile(numHeaderColumn+1:numHeaderTile),nsize, &
                            'mpi_double_precision',PE_receiver,PE_receiver*2000+mpi_myid,  &
                            'GRID',status,ierr)
          call rpn_comm_recv(myLonTile(numHeaderColumn+1:numHeaderTile),nsize, &
                            'mpi_double_precision',PE_receiver,PE_receiver*2000+mpi_myid,  &
                            'GRID',status,ierr)
          call rpn_comm_recv(myTimeInterpWeight(numHeaderColumn+1:numHeaderTile,:),nsize*statevector%numStep, &
                            'mpi_double_precision',PE_receiver,PE_receiver*2000+mpi_myid,  &
                            'GRID',status,ierr)
        end if
      else
        ! I am a receiver from TILE to COLUMN:

        ! copy lat-lon that are already on correct mpi task
        headerIndex2 = 0
        do headerIndex = 1, numHeaderColumn
          if(obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex) ==  &
             obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)) then
            headerIndex2 = headerIndex2 + 1

            call col_getLatLon( column, headerIndex,                 & ! IN
                                Lat, Lon, ypos, xpos, LatRot, LonRot ) ! OUT
            myYPosTile(headerIndex2) = ypos
            myXPosTile(headerIndex2) = xpos
            myLatTile(headerIndex2)  = lat
            myLonTile(headerIndex2)  = lon

            if(btest(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex),5) ) then
              myTimeInterpWeight(headerIndex2,:) = 0.0d0
            else
              do stepIndex = 1, statevector%numStep
                myTimeInterpWeight(headerIndex2,stepIndex) = &
                  tim_getTimeInterpWeight(headerIndex,stepIndex)
              end do
            end if
          end if
        end do

        ! send lat-lon to sender
        nsize = numHeaderColumn - numHeaderTile
        if(nsize > 0) then
          allocate(myYposToSend(nsize))
          allocate(myXposToSend(nsize))
          allocate(myLatToSend(nsize))
          allocate(myLonToSend(nsize))
          allocate(myTimeInterpWeightToSend(nsize,statevector%numStep))
          headerIndex2 = 0
          do headerIndex = 1, numHeaderColumn
            if(obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex) /=  &
               obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)) then
              headerIndex2 = headerIndex2 + 1

              call col_getLatLon( column, headerIndex,                 & ! IN
                                  Lat, Lon, ypos, xpos, LatRot, LonRot ) ! OUT
              myYPosToSend(headerIndex2) = ypos
              myXPosToSend(headerIndex2) = xpos
              myLatToSend(headerIndex2)  = lat
              myLonToSend(headerIndex2)  = lon

              if(btest(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex),5) ) then
                myTimeInterpWeightToSend(headerIndex2,:) = 0.0d0
              else
                do stepIndex = 1, statevector%numStep
                  myTimeInterpWeightToSend(headerIndex2,stepIndex) = &
                    tim_getTimeInterpWeight(headerIndex,stepIndex)
                end do
              end if

            end if
          end do
          write(*,*) 's2c_setup: headerIndex2, nsize=',headerIndex2, nsize
          call rpn_comm_send(myYPosToSend,nsize, &
                             'mpi_double_precision',PE_sender,mpi_myid*2000+PE_sender,  &
                             'GRID',ierr)
          call rpn_comm_send(myXPosToSend,nsize, &
                             'mpi_double_precision',PE_sender,mpi_myid*2000+PE_sender,  &
                             'GRID',ierr)
          call rpn_comm_send(myLatToSend,nsize, &
                             'mpi_double_precision',PE_sender,mpi_myid*2000+PE_sender,  &
                             'GRID',ierr)
          call rpn_comm_send(myLonToSend,nsize, &
                             'mpi_double_precision',PE_sender,mpi_myid*2000+PE_sender,  &
                             'GRID',ierr)
          call rpn_comm_send(myTimeInterpWeightToSend,nsize*statevector%numStep, &
                             'mpi_double_precision',PE_sender,mpi_myid*2000+PE_sender,  &
                             'GRID',ierr)

          deallocate(myYposToSend)
          deallocate(myXposToSend)
          deallocate(myLatToSend)
          deallocate(myLonToSend)
          deallocate(myTimeInterpWeightToSend)
        end if

      end if

    else  ! here numHeaderTile  ==  numHeaderColumn, so no sending or receiving

      if (numHeaderTile > 0) then

        allocate(myYPosTile(numHeaderTile))
        allocate(myXPosTile(numHeaderTile))
        allocate(myLatTile(numHeaderTile))
        allocate(myLonTile(numHeaderTile))
        allocate(myTimeInterpWeight(numHeaderTile,statevector%numStep))

        ! copy lat-lon that are already on correct mpi task
        do headerIndex = 1, numHeaderTile
          call col_getLatLon( column, headerIndex,                 & ! IN
                              Lat, Lon, ypos, xpos, LatRot, LonRot ) ! OUT
          myYPosTile(headerIndex) = ypos
          myXPosTile(headerIndex) = xpos
          myLatTile(headerIndex)  = lat
          myLonTile(headerIndex)  = lon

          if(btest(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex),5) ) then
            myTimeInterpWeight(headerIndex,:) = 0.0d0
          else
            do stepIndex = 1, statevector%numStep
              myTimeInterpWeight(headerIndex,stepIndex) = &
                tim_getTimeInterpWeight(headerIndex,stepIndex)
            end do
          end if
        end do
   
      end if ! numHeaderTile > 0

    end if ! numHeaderTile  /=  numHeaderColumn

    ! allocate the 4D array used to store the statevector with a Halo
    allocate(fieldsWithHalo(statevector%nk,  &
                            statevector%myLonBeg:(statevector%myLonEnd+1),  &
                            statevector%myLatBeg:(statevector%myLatEnd+1),statevector%numStep))

    initialized = .true.

    if(mpi_myid == 0) write(*,*) 's2c_setup: END'

  end subroutine s2c_setup

  !---------------------------------------------------------
  ! Tangent linear operator (replaces oda_L)
  !---------------------------------------------------------
  subroutine s2c_tl(statevector,column,columng,obsSpaceData)
    implicit none
    !
    ! Purpose: Horizontal interpolation to transform a 
    !          statevector object into a column object
    !
    ! Author:
    ! Revisions:
    !           Yves Rochon and Mike Sitwell, Oct 2016
    !           - Addition of gsv_varExist test for calls to ltt2phi*
    !           
    type(struct_columnData) :: column, columng
    type(struct_obs)        :: obsSpaceData
    type(struct_gsv)        :: statevector

    integer :: ilev, ierr, myLonEndP1, myLatEndP1

    call tmg_start(101,'INTERP_BARR_TL')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(101)

    if(.not.initialized) call s2c_setup(statevector,column,obsSpaceData)

    myLonEndP1 = statevector%myLonEnd + 1
    myLatEndP1 = min(statevector%myLatEnd + 1, statevector%nj)

    !
    !- 1.  Interpolation to obs location
    !
    if(mpi_myid == 0) write(*,*) 's2c_tl - Horizontal interpolation StateVector --> ColumnData'

    !- 1.1 Communicate extra latitude needed for interpolation  
    call tmg_start(39,'INTERP_COMM')
    call commLatLon(statevector)
    call tmg_stop(39)

    call gd2mvo

    !
    !- 2.  Variable conversions
    !

    !- 2.1 Mass fields (TT,PS,HU) to hydrostatic geopotential
    if (col_getNumLev(columng,'MM') > 1 .and. gsv_varExist(statevector,'TT') .and. gsv_varExist(statevector,'HU') &
        .and. gsv_varExist(statevector,'P0') ) then
       call tmg_start(36,'INTERP_TT2PHI_TL')
       call tt2phi_tl(column,columng)
       call tmg_stop(36)
    end if

    !- 2.2 Rotated wind to Meteorological wind
    if ( gsv_varExist(statevector,'UU') .and. gsv_varExist(statevector,'VV') .and.  &
         statevector%hco%Rotated ) then
      write(*,*) 'uvrot2uv Active'
      call uvrot2uv('UU', 'VV', col_getNumLev(column,'MM')) ! IN
    end if
  
  CONTAINS

    !--------------------------------------------------------------------------
    ! GD2MVO
    !--------------------------------------------------------------------------
    subroutine gd2mvo
      !
      ! s/r GD2MVO  - Horizontal bilinear interpolation of the model variables
      !               in grid-point space to observation locations.
      !
      !     numLev   : number of levels
      !
      implicit none
    
      integer :: numLev
      integer :: jlev, headerIndex, stepIndex, ilon, ilat
      real(8) :: dldx, dldy, dlw1, dlw2, dlw3, dlw4
      real(8) :: xpos, ypos, lat, lon, latrot, lonrot
      real(8), pointer :: column_ptr(:,:), column_tile(:,:)
      integer :: get_max_rss

      if (numHeaderColumn > 0) then
        column_ptr => col_getAllColumns(column)
        column_ptr(:,:) = 0.0d0
        numLev = size(column_ptr,1)
      else
        numLev = 1
        write(*,*) 'gd2mvo: numHeaderTile not positive, set numLev to 1'
      end if

      if (numHeaderTile > 0) then
        allocate(column_tile(numLev,numHeaderTile))
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      else
        allocate(column_tile(numLev,1))
        write(*,*) 'gd2mvo: numHeaderTile not positive: ',numHeaderTile
      end if
      column_tile(:,:) = 0.0d0

      ! Note: We assume here that all the obs between the poles and the last grid points
      !       (i.e. outside the grid) have been moved within the grid by sugomobs

      do stepIndex = 1, statevector%numStep

        !- Loop over all the observations
        do headerIndex = 1, numHeaderTile

          if ( myTimeInterpWeight(headerIndex,stepIndex) > 0.0d0 ) then

            !- 2.1 Find the obs position within the analysis grid
            ypos = myYposTile(headerIndex)
            xpos = myXposTile(headerIndex)
            Lat  = myLatTile(headerIndex)
            Lon  = myLonTile(headerIndex)

            !- Make sure we are within bounds
            if ( ypos < real(statevector%myLatBeg,8) .or. &
                 ypos > real(myLatEndP1          ,8) .or. &
                 xpos < real(statevector%myLonBeg,8) .or. &
                 xpos > real(myLonEndP1          ,8) ) then
              write(*,*) 's2c_tl: Obs outside local domain for headerIndex = ', headerIndex
              write(*,*) '        obs lat, lon position               = ',  &
                         Lat*MPC_DEGREES_PER_RADIAN_R8, Lon*MPC_DEGREES_PER_RADIAN_R8
              write(*,*) '        obs x, y     position               = ', xpos, ypos
              write(*,*) '        domain x_start, x_end, y_start, y_end bounds = ',  &
                         statevector%myLonBeg, myLonEndP1, statevector%myLatBeg, myLatEndP1

              ! if obs above or below latitude band, move it to the edge of this latitude band
              if( ypos < real(statevector%myLatBeg,8) ) ypos = real(statevector%myLatBeg,8)
              if( ypos > real(myLatEndP1          ,8) ) ypos = real(myLatEndP1          ,8)

              ! if obs left or right longitude band, move it to the edge of this longitude band
              if( xpos < real(statevector%myLonBeg,8) ) xpos = real(statevector%myLonBeg,8)
              if( xpos > real(myLonEndP1          ,8) ) xpos = real(myLonEndP1          ,8)
              write(*,*) ' new   obs x, y     position               = ', xpos, ypos
            end if

            !- 2.2 Find the lower-left grid point next to the observation
            if ( xpos /= real(myLonEndP1,8) ) then
              ilon = floor(xpos)
            else
              ilon = floor(xpos) - 1
            end if

            if ( ypos /= real(myLatEndP1,8) ) then
              ilat = floor(ypos)
            else
              ilat = floor(ypos) - 1
            end if

            !- 2.3 Compute the 4 weights of the bilinear interpolation
            dldx = xpos - real(ilon,8)
            dldy = ypos - real(ilat,8)

            dlw1 = (1.d0-dldx) * (1.d0-dldy)
            dlw2 =       dldx  * (1.d0-dldy)
            dlw3 = (1.d0-dldx) *       dldy
            dlw4 =       dldx  *       dldy

            !- 2.4 Interpolate the model state to the obs point
            do jlev = 1, numLev
              column_tile(jlev,headerIndex) = column_tile(jlev,headerIndex)  +  &
                                  myTimeInterpWeight(headerIndex,stepIndex) * &
                                  ( dlw1 * fieldsWithHalo(jlev,ilon  ,ilat  ,stepIndex) &
                                  + dlw2 * fieldsWithHalo(jlev,ilon+1,ilat  ,stepIndex) &
                                  + dlw3 * fieldsWithHalo(jlev,ilon  ,ilat+1,stepIndex) &
                                  + dlw4 * fieldsWithHalo(jlev,ilon+1,ilat+1,stepIndex) )
            end do

          end if

        end do ! headerIndex

      end do ! stepIndex

      call tmg_start(102,'INTERP_BARR_TL2')
      if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(102)

      call transpose_tileToColumn(column_tile,column_ptr,numLev)
      deallocate(column_tile)

    end subroutine gd2mvo

    !--------------------------------------------------------------------------
    ! transpose_tileToColumn
    !--------------------------------------------------------------------------
    subroutine transpose_tileToColumn(col_ptr_in,col_ptr_out,numLev)
      implicit none
      real(8), pointer :: col_ptr_in(:,:),col_ptr_out(:,:)
      integer :: numLev

      real(8), allocatable :: col_recv(:,:)
      integer :: headerIndex, headerIndex_local, headerIndex_recv, jlev, nsize, ierr, status

      call tmg_start(39,'INTERP_COMM')

      nsize = abs(numlev*(numHeaderTile-numHeaderColumn))

      if(nsize > 0) then

        if(PE_receiver_mpiglobal(mpi_myid+1) /= -1) then
          ! I am a sender from TILE to COLUMN
          col_ptr_out(:,1:numHeaderColumn) = col_ptr_in(:,1:numHeaderColumn)

          call tmg_start(35,'INTERP_SENDRECV')
          call rpn_comm_send(col_ptr_in(1:numlev,numHeaderColumn+1:numHeaderTile),nsize, &
                             'mpi_double_precision',PE_receiver_mpiglobal(mpi_myid+1), &
                             PE_receiver_mpiglobal(mpi_myid+1)*2000+mpi_myid,  &
                             'GRID',ierr)
          call tmg_stop(35)
        elseif(PE_sender_mpiglobal(mpi_myid+1) /= -1) then
          ! I am a receiver from TILE to COLUMN:
          allocate(col_recv(numLev,nsize))
          call tmg_start(35,'INTERP_SENDRECV')
          call rpn_comm_recv(col_recv(:,:),nsize, &
                             'mpi_double_precision',PE_sender_mpiglobal(mpi_myid+1), &
                             mpi_myid*2000+PE_sender_mpiglobal(mpi_myid+1),  &
                             'GRID',status,ierr)
          call tmg_stop(35)
          headerIndex_local = 0
          headerIndex_recv = 0
          do headerIndex = 1, numHeaderColumn
            if(myPEsourceColumn(headerIndex) == mpi_myid) then
              headerIndex_local = headerIndex_local + 1
              do jlev = 1, numLev
                col_ptr_out(jlev,headerIndex) = col_ptr_in(jlev,headerIndex_local)
              end do
            else
              headerIndex_recv = headerIndex_recv + 1
              do jlev = 1, numLev
                col_ptr_out(jlev,headerIndex) = col_recv(jlev,headerIndex_recv)
              end do
            end if
          end do
          deallocate(col_recv)
        else
          call utl_abort('transpose_tileToColumn: NOT SURE WHAT IS GOING ON, CALL MARK')
        end if

      else

        ! No communication
        if(numHeaderColumn > 0) then
          col_ptr_out(:,1:numHeaderColumn) = col_ptr_in(:,1:numHeaderColumn)
        end if
      end if

      call tmg_stop(39)

    end subroutine transpose_tileToColumn

    !--------------------------------------------------------------------------
    ! UVROT2UV
    !--------------------------------------------------------------------------
    subroutine uvrot2uv (UUvarName,VVvarName,numLev)
      !
      !- uvrot2uv - Transforms tangential (U,V) wind components at observation
      !             locations on GEM rotated frame to the real sphere.
      use WindRotation_mod
      implicit none

      character(len=*), intent(in) :: UUvarName
      character(len=*), intent(in) :: VVvarName
      integer,          intent(in) :: numLev
    
      real(8) :: lat, lon, latrot, lonrot, xpos, ypos
      real(8), pointer :: UUcolumn(:), VVcolumn(:)
      integer :: headerIndex

      !
      !- 1.  Loop over all the observation locations
      !
      do headerIndex = 1, col_getNumCol(column)

        !- 1.1 Extract (rotated) wind profiles
        UUColumn => col_getColumn(column,headerIndex,UUvarName)
        VVColumn => col_getColumn(column,headerIndex,VVvarName)
       
        !- 1.2 Find the latitudes and longitudes
        call col_getLatLon( column, headerIndex,                   & ! IN
                            Lat, Lon, ypos, xpos, LatRot, LonRot )   ! OUT

        !- 1.3 Rotate Winds
        call uvr_RotateWind( UUColumn, VVColumn,       & ! INOUT
                             Lat, Lon, LatRot, LonRot, & ! IN
                             'ToMetWind', numLev )       ! IN

      end do

    end subroutine uvrot2uv

  end subroutine s2c_tl


  !---------------------------------------------------------
  ! Adjoint operator (replaces oda_LT)
  !---------------------------------------------------------
  subroutine s2c_ad(statevector,column,columng,obsSpaceData)
    implicit none
    !
    ! Purpose: Adjoint of horizontal interpolations
    !          Compute adjoint statevector object from adjoint column object
    !
    ! Author:
    ! Revisions:
    !           Yves Rochon and Mike Sitwell, Oct 2016
    !           - Addition of gsv_varExist test for calls to att2phi*
    !           

    type(struct_columnData) :: column,columng
    type(struct_obs)        :: obsSpaceData
    type(struct_gsv)        :: statevector

    integer :: ierr, myLonEndP1, myLatEndP1

    call tmg_start(103,'INTERP_BARR_AD')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(103)

    if(.not.initialized) call s2c_setup(statevector,column,obsSpaceData)

    myLonEndP1 = statevector%myLonEnd + 1
    myLatEndP1 = min(statevector%myLatEnd + 1, statevector%nj)

    !
    !- 2.  Variable conversions
    !

    !- 2.2 Rotated wind to Meteorological wind
    if ( gsv_varExist(statevector,'UU') .and. gsv_varExist(statevector,'VV') .and. &
         statevector%hco%Rotated ) then
       write(*,*) 'uvrot2uvAdj Active' 
       call uvrot2uvAdj('UU', 'VV', col_getNumLev(column,'MM')) ! IN
    end if
  
    !- 2.1 Mass fields (TT,PS,HU) to hydrostatic geopotential
    if (col_getNumLev(columng,'MM') > 1 .and. gsv_varExist(statevector,'TT') .and. gsv_varExist(statevector,'HU') &
        .and. gsv_varExist(statevector,'P0') ) then
       call tmg_start(37,'INTERP_TT2PHI_AD')
       call tt2phi_ad(column,columng)
       call tmg_stop(37)
    end if

    !
    !- 1.  Interpolation to obs location
    !
    if(mpi_myid == 0) write(*,*) 's2c_ad - Adjoint of horizontal interpolation StateVector --> ColumnData'
  
    call gd2mvoad

    !- 1.1 Communicate extra latitude needed for interpolation
    call tmg_start(39,'INTERP_COMM')
    call commLatLonAd(statevector)
    call tmg_stop(39)

  CONTAINS

    !--------------------------------------------------------------------------
    ! GD2MVOAD
    !--------------------------------------------------------------------------
    subroutine gd2mvoad
      !
      ! s/r GD2MVOAD  - Adjoint of the bilinear horizontal interpolation of the 
      !                 model variables in grid-point space to observation locations.
      !
      !    Purpose:  Update the estimate of GD from the gradient components
      !              at the observation points which have been stored in
      !              column
      !
      implicit none
    
      integer :: numLev    
      integer :: jlev, headerIndex, stepIndex, ilon, ilat
      real(8) :: dldx, dldy, dlw1, dlw2, dlw3, dlw4, DInterpWeight
      real(8) :: xpos, ypos, lat, lon, latrot, lonrot
      real(8), pointer :: column_ptr(:,:), column_tile(:,:)

      if (numHeaderColumn > 0) then
        column_ptr => col_getAllColumns(column)
        numLev = size(column_ptr,1)
      else
        numLev = 1
      end if

      if (numHeaderTile > 0) then
        allocate(column_tile(numLev,numHeaderTile))
      else
        allocate(column_tile(numLev,1))
      end if
      column_tile(:,:) = 0.0d0

      call tmg_start(104,'INTERP_BARR_AD2')
      if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(104)

      call transpose_columnToTile(column_tile,column_ptr,numLev)

      ! Note: We assume here the all the obs between the poles and the last grid points
      !       (i.e. outside the grid) have been moved within the grid previously


      do stepIndex = 1, statevector%numStep

        fieldsWithHalo(:,:,:,stepIndex) = 0.0D0
      
        !
        !- 2.  LOOP OVER ALL THE OBSERVATIONS
        !
        do headerIndex = 1, numHeaderTile

          if ( myTimeInterpWeight(headerIndex,stepIndex) > 0.0d0 ) then

            !- 2.1 Find the obs position within the analysis grid
            ypos = myYposTile(headerIndex)
            xpos = myXposTile(headerIndex)
            Lat  = myLatTile(headerIndex)
            Lon  = myLonTile(headerIndex)

            !- Make sure we are within bounds
            if ( ypos < real(statevector%myLatBeg,8) .or. &
                 ypos > real(myLatEndP1          ,8) .or. &
                 xpos < real(statevector%myLonBeg,8) .or. &
                 xpos > real(myLonEndP1          ,8) ) then
              write(*,*) 's2c_ad: Obs outside local domain for job = ', headerIndex
              write(*,*) '  obs x, y position  = ', xpos, ypos
              write(*,*) '  domain x_start, x_end, y_start, y_end bounds = ',  &
                         statevector%myLonBeg, myLonEndP1, statevector%myLatBeg, myLatEndP1

              ! if obs above or below latitude band, move it to the edge of this latitude band
              if( ypos < real(statevector%myLatBeg,8) ) ypos = real(statevector%myLatBeg,8)
              if( ypos > real(myLatEndP1          ,8) ) ypos = real(myLatEndP1          ,8)

              ! abort if obs is to the left or right of the analysis domain
              if( xpos < real(statevector%myLonBeg,8) ) xpos = real(statevector%myLonBeg,8)
              if( xpos > real(myLonEndP1          ,8) ) xpos = real(myLonEndP1          ,8)
            end if

            !- 2.2 Find the lower-left grid point next to the observation
            if ( xpos /= real(myLonEndP1,8) ) then
              ilon = floor(xpos)
            else
              ilon = floor(xpos) - 1
            end if

            if ( ypos /= real(myLatEndP1,8) ) then
              ilat = floor(ypos)
            else
              ilat = floor(ypos) - 1
            end if

            !- 2.3 COMPUTE THE 4 WEIGHTS OF THE BILINEAR INTERPOLATION
            dldx = xpos - real(ilon,8)
            dldy = ypos - real(ilat,8)

            DInterpWeight = myTimeInterpWeight(headerIndex,stepIndex)
          
            dlw1  = DInterpWeight * (1.d0-dldx) * (1.d0-dldy)
            dlw2  = DInterpWeight *       dldx  * (1.d0-dldy)
            dlw3  = DInterpWeight * (1.d0-dldx) *       dldy
            dlw4  = DInterpWeight *       dldx  *       dldy

            !- 2.4 Interpolate the model state to the obs point
            do jlev = 1, numLev
              fieldsWithHalo(jlev,ilon  ,ilat,  stepIndex) = fieldsWithHalo(jlev,ilon  ,ilat,  stepIndex)    &
                                          + dlw1 * column_tile(jlev,headerIndex)
              fieldsWithHalo(jlev,ilon+1,ilat,  stepIndex) = fieldsWithHalo(jlev,ilon+1,ilat,  stepIndex)    &
                                          + dlw2 * column_tile(jlev,headerIndex)
              fieldsWithHalo(jlev,ilon  ,ilat+1,stepIndex) = fieldsWithHalo(jlev,ilon  ,ilat+1,stepIndex)    &
                                          + dlw3 * column_tile(jlev,headerIndex)
              fieldsWithHalo(jlev,ilon+1,ilat+1,stepIndex) = fieldsWithHalo(jlev,ilon+1,ilat+1,stepIndex)    &
                                          + dlw4 * column_tile(jlev,headerIndex)
            end do

          end if

        end do ! headerIndex

      end do ! stepIndex

      deallocate(column_tile)

    end subroutine gd2mvoad

    !--------------------------------------------------------------------------
    ! transpose_columnToTile
    !--------------------------------------------------------------------------
    subroutine transpose_columnToTile(col_ptr_out,col_ptr_in,numLev)
      implicit none
      real(8), pointer :: col_ptr_in(:,:),col_ptr_out(:,:)
      integer :: numLev

      real(8), allocatable :: col_send(:,:)
      integer :: headerIndex, headerIndex_local, headerIndex_send, jlev, nsize, ierr, status

      call tmg_start(39,'INTERP_COMM')

      nsize = abs(numlev*(numHeaderTile-numHeaderColumn))

      if(nsize > 0) then

        if(PE_receiver_mpiglobal(mpi_myid+1) /= -1) then
          ! I am a receiver from COLUMN to TILE
          col_ptr_out(:,1:numHeaderColumn) = col_ptr_in(:,1:numHeaderColumn)

          call tmg_start(35,'INTERP_SENDRECV')
          call rpn_comm_recv(col_ptr_out(1:numlev,numHeaderColumn+1:numHeaderTile),nsize, &
                             'mpi_double_precision',PE_receiver_mpiglobal(mpi_myid+1), &
                             PE_receiver_mpiglobal(mpi_myid+1)*2000+mpi_myid,  &
                             'GRID',status,ierr)
          call tmg_stop(35)


        elseif(PE_sender_mpiglobal(mpi_myid+1) /= -1) then
          ! I am a sender from COLUMN to TILE:
          allocate(col_send(numLev,nsize))
          headerIndex_local = 0
          headerIndex_send = 0
          do headerIndex = 1, numHeaderColumn
            if(myPEsourceColumn(headerIndex) == mpi_myid) then
              headerIndex_local = headerIndex_local + 1
              do jlev = 1, numLev
                col_ptr_out(jlev,headerIndex_local) = col_ptr_in(jlev,headerIndex)
              end do
            else
              headerIndex_send = headerIndex_send + 1
              do jlev = 1, numLev
                col_send(jlev,headerIndex_send) = col_ptr_in(jlev,headerIndex)
              end do
            end if
          end do
          call tmg_start(35,'INTERP_SENDRECV')
          call rpn_comm_send(col_send(:,:),nsize, &
                             'mpi_double_precision',PE_sender_mpiglobal(mpi_myid+1), &
                             mpi_myid*2000+PE_sender_mpiglobal(mpi_myid+1),  &
                             'GRID',ierr)
          call tmg_stop(35)
          deallocate(col_send)
        else
          call utl_abort('transpose_columnToTile: NOT SURE WHAT IS GOING ON, CALL MARK')
        end if

      else

        ! No communication
        if(numHeaderColumn > 0) then
          col_ptr_out(:,1:numHeaderColumn) = col_ptr_in(:,1:numHeaderColumn)
        end if

      end if

      call tmg_stop(39)

    end subroutine transpose_columnToTile

    !--------------------------------------------------------------------------
    ! UVROT2UVADJ
    !--------------------------------------------------------------------------
    subroutine uvrot2uvAdj(UUvarName,VVvarName,numLev)
      !
      !- uvrot2uv - Transforms tangential (U,V) wind components at observation
      !             locations on GEM rotated frame to the real sphere.
      use WindRotation_mod
      implicit none

      character(len=*), intent(in) :: UUvarName
      character(len=*), intent(in) :: VVvarName
      integer,          intent(in) :: numLev
    
      real(8) :: lat, lon, latrot, lonrot, xpos, ypos
      real(8), pointer :: UUcolumn(:), VVcolumn(:)
      integer :: headerIndex

      !
      !- 1.  Loop over all the observation locations
      !
      do headerIndex = 1, col_getNumCol(column)

        !- 1.1 Extract (rotated) wind profiles
        UUColumn => col_getColumn(column,headerIndex,UUvarName)
        VVColumn => col_getColumn(column,headerIndex,VVvarName)
       
        !- 1.2 Find the latitudes and longitudes
        call col_getLatLon( column, headerIndex,                   & ! IN
                            Lat, Lon, ypos, xpos, LatRot, LonRot )   ! OUT

        !- 1.3 Rotate Winds
        call uvr_RotateWindAdj( UUColumn, VVColumn,       & ! INOUT
                                Lat, Lon, LatRot, LonRot, & ! IN
                                'ToMetWind', numLev )       ! IN

      end do

    end subroutine uvrot2uvAdj

  end subroutine s2c_ad


  subroutine commLatLon(statevector_in)
    implicit none
    type(struct_gsv) :: statevector_in
    integer :: nsize, ierr, status, latPerPEhalo, myLatEndP1
    integer :: jlat, jstep, jlev, jlon
    real(8), pointer :: field_ptr(:,:,:,:) 

    ! copy statevector into fieldsWithHalo
    field_ptr => gsv_getField_r8(statevector_in)
!$OMP PARALLEL DO PRIVATE (jlat,jstep,jlev,jlon)    
    do jstep = 1, statevector_in%numStep
      do jlev = 1, statevector_in%nk
        do jlat = statevector_in%myLatBeg, statevector_in%myLatEnd
          do jlon = statevector_in%myLonBeg, statevector_in%myLonEnd
            fieldsWithHalo(jlev,jlon,jlat,jstep) = field_ptr(jlon,jlat,jlev,jstep)
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

    ! ******First send latitude halos

    if(mpi_npey > 1) then  ! only do exchange when more than one mpi task in Y direction

      nsize = statevector_in%lonPerPE * statevector_in%nk * statevector_in%numStep

      ! northern most latitude band
      if(mpi_myidy == (mpi_npey-1)) then
        call rpn_comm_send(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                                            statevector_in%myLatBeg:statevector_in%myLatBeg,:),nsize, &
                           'mpi_double_precision',mpi_myidy-1,mpi_myidy*500+(mpi_myidy-1),  &
                           'NS',ierr)
      end if

      ! all latitude bands not at the north or south poles
      if(mpi_myidy /= 0 .and. mpi_myidy /= (mpi_npey-1)) then
        call rpn_comm_sendrecv(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                                                statevector_in%myLatBeg:statevector_in%myLatBeg,:), &
                               nsize,'mpi_double_precision',mpi_myidy-1,mpi_myidy*500+(mpi_myidy-1), &
                               fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd, &
                                                (statevector_in%myLatEnd+1):(statevector_in%myLatEnd+1),:), &
                               nsize,'mpi_double_precision',mpi_myidy+1,(mpi_myidy+1)*500+mpi_myidy, &
                               'NS',status,ierr)
      end if

      ! southern most latitude band
      if(mpi_myidy == 0) then
        call rpn_comm_recv(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd, &
                                          (statevector_in%myLatEnd+1):(statevector_in%myLatEnd+1),:), &
                           nsize,'mpi_double_precision',mpi_myidy+1,(mpi_myidy+1)*500+mpi_myidy,  &
                           'NS',status,ierr)
      end if

    end if ! mpi_npey > 1


    ! ******Now send longitude halos

    if(mpi_myidy == (mpi_npey-1)) then
      ! northern most latitude band does not have a latitude halo to the north
      latPerPEhalo = statevector_in%latPerPE
      myLatEndP1 = statevector_in%myLatEnd
    else
      ! all others do
      latPerPEhalo = statevector_in%latPerPE + 1
      myLatEndP1 = statevector_in%myLatEnd + 1
    end if

    if(mpi_npex > 1) then  ! only do exchange when more than one mpi task in X direction

      nsize = latPerPEhalo * statevector_in%nk * statevector_in%numStep

      ! eastern most longitude band
      if(mpi_myidx == (mpi_npex-1)) then
        call rpn_comm_send(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonBeg,  &
                                            statevector_in%myLatBeg:myLatEndP1,:),  &
                           nsize, 'mpi_double_precision',mpi_myidx-1,mpi_myidx*500+(mpi_myidx-1), &
                           'EW',ierr)
      end if

      ! all other longitude bands (not first nor last)
      if(mpi_myidx /= 0 .and. mpi_myidx /= (mpi_npex-1)) then
        call rpn_comm_sendrecv(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonBeg,  &
                                                statevector_in%myLatBeg:myLatEndP1,:), &
                               nsize,'mpi_double_precision',mpi_myidx-1,mpi_myidx*500+(mpi_myidx-1), &
                               fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                                statevector_in%myLatBeg:myLatEndP1,:), &
                               nsize,'mpi_double_precision',mpi_myidx+1,(mpi_myidx+1)*500+mpi_myidx, &
                               'EW',status,ierr)
      end if

      ! western most longitude band
      if(mpi_myidx == 0) then
        call rpn_comm_recv(fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                            statevector_in%myLatBeg:myLatEndP1,:),nsize, &
          'mpi_double_precision',mpi_myidx+1,(mpi_myidx+1)*500+mpi_myidx,'EW',status,ierr)
      end if

      ! periodic, so also send the first meridian on myidx=0 to the last meridian on myidx=(npex-1)
      if(mpi_myidx == 0) then
        call rpn_comm_send(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonBeg,  &
                                          statevector_in%myLatBeg:myLatEndP1,:),  &
                           nsize, 'mpi_double_precision', mpi_npex-1, mpi_myidx*500+(mpi_npex-1), &
                           'EW', ierr)
      end if
      if(mpi_myidx == (mpi_npex-1)) then
        call rpn_comm_recv(fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                          statevector_in%myLatBeg:myLatEndP1,:), &
                           nsize, 'mpi_double_precision',0,0*500+mpi_myidx,  &
                           'EW',status,ierr)
      end if

    else ! only one mpi task in X direction, so just copy first meridian to last (plus 1)
      
      fieldsWithHalo(:,statevector_in%myLonEnd+1,statevector_in%myLatBeg:myLatEndP1,:) = &
        fieldsWithHalo(:,1                      ,statevector_in%myLatBeg:myLatEndP1,:)

    end if

  end subroutine commlatlon


  subroutine commLatLonAd(statevector_in)
    implicit none
    type(struct_gsv) :: statevector_in
    integer :: nsize, ierr, status, latPerPEhalo, myLatEndP1
    real*8, allocatable :: latHalo(:,:,:,:)
    real*8, allocatable :: lonHalo(:,:,:,:)
    integer :: jlat, jstep, jlev, jlon
    real(8), pointer :: field_ptr(:,:,:,:) 

    ! ******Adjoint of sending longitude halos

    if(mpi_myidy == (mpi_npey-1)) then
      ! northern most latitude band does not have a latitude halo to the north
      latPerPEhalo = statevector_in%latPerPE
      myLatEndP1 = statevector_in%myLatEnd
    else
      ! all others do
      latPerPEhalo = statevector_in%latPerPE + 1
      myLatEndP1 = statevector_in%myLatEnd + 1
    end if

    if(mpi_npex > 1) then  ! only do adjoint of exchange when more than one mpi task in X direction

      allocate(lonHalo(statevector_in%nk, 1, latPerPEhalo, statevector_in%numStep))

      nsize = latPerPEhalo*statevector_in%nk*statevector_in%numStep

      ! periodic, so also do adjoint of sending the first meridian on myidx=0 to the last meridian on myidx=(npex-1)
      if(mpi_myidx == (mpi_npex-1)) then
        call rpn_comm_send(fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                          statevector_in%myLatBeg:myLatEndP1,:),nsize, &
                           'mpi_double_precision',0,0*500+mpi_myidx,  &
                           'EW',ierr)
      end if
      if(mpi_myidx == 0) then
        call rpn_comm_recv(lonHalo,nsize, &
                           'mpi_double_precision',mpi_npex-1,mpi_myidx*500+(mpi_npex-1),'EW',status,ierr)
      end if

      ! western most longitude band
      if(mpi_myidx == 0) then
        call rpn_comm_send(fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                          statevector_in%myLatBeg:myLatEndP1,:),nsize, &
                           'mpi_double_precision',mpi_myidx+1,(mpi_myidx+1)*500+mpi_myidx,  &
                           'EW',ierr)
      end if

      ! all other longitude bands (not first nor last)
      if(mpi_myidx /= 0 .and. mpi_myidx /= (mpi_npex-1)) then
        call rpn_comm_sendrecv(fieldsWithHalo(:,(statevector_in%myLonEnd+1):(statevector_in%myLonEnd+1), &
                                              statevector_in%myLatBeg:myLatEndP1,:), &
                               nsize,'mpi_double_precision',mpi_myidx+1,(mpi_myidx+1)*500+mpi_myidx, &
                               lonHalo, &
                               nsize,'mpi_double_precision',mpi_myidx-1,mpi_myidx*500+(mpi_myidx-1), &
                               'EW',status,ierr)
      end if

      ! eastern most longitude band
      if(mpi_myidx == (mpi_npex-1)) then
        call rpn_comm_recv(lonHalo,nsize, &
                           'mpi_double_precision',mpi_myidx-1,mpi_myidx*500+(mpi_myidx-1),'EW',status,ierr)
      end if

      ! add the sensitivity from the halo to the in situ sensitivity
      fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonBeg,  &
                     statevector_in%myLatBeg:myLatEndP1,:) = &
        fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonBeg,  &
                       statevector_in%myLatBeg:myLatEndP1,:) + lonHalo(:,:,:,:)

      ! to make sure sensitivity from the halo is not double counted, set to zero
      fieldsWithHalo(:,statevector_in%myLonEnd+1,statevector_in%myLatBeg:myLatEndP1,:) = 0.0d0

      deallocate(lonHalo)

    else ! only one mpi task in X direction, so just adjoint of copying first meridian to last (plus 1)
      
      fieldsWithHalo(:,1                          ,statevector_in%myLatBeg:myLatEndP1,:) =  &
        fieldsWithHalo(:,1                        ,statevector_in%myLatBeg:myLatEndP1,:) +  &
        fieldsWithHalo(:,statevector_in%myLonEnd+1,statevector_in%myLatBeg:myLatEndP1,:)
        
      ! to make sure sensitivity from the halo is not double counted, set to zero
      fieldsWithHalo(:,statevector_in%myLonEnd+1,statevector_in%myLatBeg:myLatEndP1,:) = 0.0d0

    end if

    ! ******Adjoint of sending latitude halos

    if(mpi_npey > 1) then  ! only do exchange when more than one mpi task in Y direction

      allocate(latHalo(statevector_in%nk, statevector_in%lonPerPE, 1, statevector_in%numStep))

      nsize = statevector_in%lonPerPE*statevector_in%nk*statevector_in%numStep

      ! southern most latitude band
      if(mpi_myidy == 0) then
        call rpn_comm_send(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                                          (statevector_in%myLatEnd+1):(statevector_in%myLatEnd+1),:), &
                           nsize,'mpi_double_precision',mpi_myidy+1,(mpi_myidy+1)*500+mpi_myidy,'NS',ierr)
      end if

      ! all latitude bands not at the north or south poles
      if(mpi_myidy /= 0 .and. mpi_myidy /= (mpi_npey-1)) then
        call rpn_comm_sendrecv(fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                                              (statevector_in%myLatEnd+1):(statevector_in%myLatEnd+1),:), &
                               nsize,'mpi_double_precision',mpi_myidy+1,(mpi_myidy+1)*500+mpi_myidy, &
                               latHalo, &
                               nsize,'mpi_double_precision',mpi_myidy-1,mpi_myidy*500+(mpi_myidy-1), &
                               'NS',status,ierr)
      end if

      ! northern most latitude band
      if(mpi_myidy == (mpi_npey-1)) then
        call rpn_comm_recv(latHalo,nsize,  &
                           'mpi_double_precision',mpi_myidy-1,mpi_myidy*500+(mpi_myidy-1), &
                           'NS',status,ierr)
      end if

      ! add the sensitivity from the halo to the in situ sensitivity
      if(mpi_myidy /= 0) then
        fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                       statevector_in%myLatBeg:statevector_in%myLatBeg,:) = &
               fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,  &
                              statevector_in%myLatBeg:statevector_in%myLatBeg,:) + &
               latHalo(:,:,:,:)
      end if

      ! to make sure sensitivity from the halo is not double counted, set to zero
      if(mpi_myidy /= (mpi_npey-1)) then
        fieldsWithHalo(:,statevector_in%myLonBeg:statevector_in%myLonEnd,statevector_in%myLatEnd+1,:) = 0.0d0
      end if

      deallocate(latHalo)

    end if

    ! copy statevector into fieldsWithHalo
    field_ptr => gsv_getField_r8(statevector_in)
!$OMP PARALLEL DO PRIVATE (jlat,jstep,jlev,jlon)    
    do jstep = 1, statevector_in%numStep
      do jlev = 1, statevector_in%nk
        do jlat = statevector_in%myLatBeg, statevector_in%myLatEnd
          do jlon = statevector_in%myLonBeg, statevector_in%myLonEnd
            field_ptr(jlon,jlat,jlev,jstep) = fieldsWithHalo(jlev,jlon,jlat,jstep)
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine commLatLonAd

  !--------------------------------------------------------------------------
  ! S2C_COLUMN_HBILIN  
  !--------------------------------------------------------------------------
  subroutine s2c_column_hbilin(field,vlev,nlong,nlat,nlev,xlong,xlat,plong,plat,vprof,vlevout,nlevout)
  !
  ! Author:  Y. Rochon, Nov 2015 
  !
  ! Purpose: Horizontal bilinear interpolation from a 3D field to a profile at (plong,plat).
  !          Assumes vertical interpolation not needed or already done.
  !
  !          This version can be used with fields that are not part of the background state,
  !          such as climatologies.
  !
  !          This version does not depend in column_data and gridstatevector modules.
  !
  ! Arguments:
  !
  !   Input
  !      
  !     field(nlong,nlat,nlev)  3D field
  !     nlong         number or latitudes
  !     nlat          number of longitudes
  !     nlev          number of vertical levels
  !     xlong         longitudes (radians)
  !     xlat          latitudes (radians)
  !     vlev          vertical levels of input field (in pressure)
  !     plat          target latitude (radian)
  !     plong         target longitude (radians) 
  !     nlevout       Number of target vertical levels
  !     vlevout       Target vertical levels (in pressure)
  !
  !   Output
  !
  !     vprof(nlev)   Profile at (plong,plat) 
  !
  !-------------------------------------------------------------------------------------------  

    implicit none

    integer, intent(in) :: nlong,nlat,nlev,nlevout
    real(8), intent(in) :: field(nlong,nlat,nlev),vlev(nlev),xlong(nlong),xlat(nlat),plong,plat
    real(8), intent(in) :: vlevout(nlevout)
    real(8), intent(out) :: vprof(nlevout)
    
    real(8) :: lnvlev(nlev),lnvlevout(nlevout),plong2
    integer :: ilev,ilon,ilat,i,j

    real(8) :: DLDX, DLDY, DLDP, DLW1, DLW2, DLW3, DLW4

    ! Find near lat/long grid points
    
    plong2 = plong
    if (plong2 < 0.0) plong2 = 2.D0*MPC_PI_R8 + plong2
    do ilon = 2,nlong
       if  (xlong(ilon-1) < xlong(ilon)) then
           if (plong2 >= xlong(ilon-1) .and. plong2 <= xlong(ilon)) exit
       else 
           ! Assumes this is a transition between 360 to 0 (if it exists). Skip over.
       end if
    end do
    ilon = ilon-1
       
    do ilat = 2,nlat
       if (plat <= xlat(ilat)) exit
    end do
    ilat = ilat-1
    
    ! Set lat/long interpolation weights
    
    DLDX = (plong - xlong(ilon))/(xlong(ilon+1)-xlong(ilon))
    DLDY = (plat - xlat(ilat))/(xlat(ilat+1)-xlat(ilat))

    DLW1 = (1.d0-DLDX) * (1.d0-DLDY)
    DLW2 =       DLDX  * (1.d0-DLDY)
    DLW3 = (1.d0-DLDX) *       DLDY
    DLW4 =       DLDX  *       DLDY

    ! Set vertical interpolation weights (assumes pressure vertical coordinate)
    
    lnvlevout(:) = log(vlevout(:))    
    lnvlev(:) = log(vlev(:))    
!    lnvlev(:) = DLW1 * log(vlev(ilon,:,ilat)) &
!               + DLW2 * log(vlev(ilon+1,:,ilat)) &
!               + DLW3 * log(vlev(ilon,:,ilat+1)) &
!               + DLW4 * log(vlev(ilon+1,:,ilat+1)) 
         
    ilev = 1
    do i = 1, nlevout
       do j = ilev, nlev          
          if (lnvlevout(i) < lnvlev(j)) exit    ! assumes both lnvlevout and lnvlev increase with increasing index value
       end do
       ilev = j-1
       if (ilev < 1) then
          ilev = 1
       else if (ilev >= nlev) then
           ilev = nlev-1
       end if
       
       DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
       vprof(i) = DLDP* (DLW1 * field(ilon,ilev,ilat) &
                       + DLW2 * field(ilon+1,ilev,ilat) &
                       + DLW3 * field(ilon,ilev,ilat+1) &
                       + DLW4 * field(ilon+1,ilev,ilat+1)) &
         + (1.d0-DLDP)* (DLW1 * field(ilon,ilev+1,ilat) &
                       + DLW2 * field(ilon+1,ilev+1,ilat) &
                       + DLW3 * field(ilon,ilev+1,ilat+1) &
                       + DLW4 * field(ilon+1,ilev+1,ilat+1))                               
    end do
        
  end subroutine s2c_column_hbilin   

end module stateToColumn_mod