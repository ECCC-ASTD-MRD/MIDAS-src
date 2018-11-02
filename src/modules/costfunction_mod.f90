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
!! MODULE utilities_mod,(prefix="cfn" category='1')
!!
!! *Purpose*: To compute Jo term
!!
!--------------------------------------------------------------------------
module costfunction_mod
  use mpi_mod
  use mpivar_mod
  use obsSpaceData_mod
  use rmatrix_mod
  use rttov_const, only : inst_name, platform_name
  use tovs_nl_mod

  implicit none

  save
  private

  ! public procedures

  public :: cfn_calcJo, cfn_sumJo, cfn_RsqrtInverse

contains

  !--------------------------------------------------------------------------
  ! cfn_RsqrtInverse
  !--------------------------------------------------------------------------
  subroutine cfn_RsqrtInverse( lobsSpaceData, elem_dest_i, elem_src_i )

    implicit none

    type(struct_obs), intent(inout) :: lobsSpaceData
    integer, intent(in)  :: elem_dest_i ! destination index
    integer, intent(in)  :: elem_src_i  ! source index
    !
    ! Applied observation error variances to ROBDATA8(k_src,*)
    ! and store it in the elem_src_s of lobsSpaceData
    !
    integer :: bodyIndex, iass, headerIndex
    integer :: idata, idatend, idatyp, count, ichn
    real(8) :: x( tvs_maxChannelNumber ), y( tvs_maxChannelNumber )
    integer :: list_chan( tvs_maxChannelNumber )

    ! NOTE I tried using openMP on this loop, but it increased the cost from 4sec to 80sec!!!
    do headerIndex =1, obs_numHeader(lobsSpaceData)

      idata   = obs_headElem_i( lobsSpaceData, OBS_RLN, headerIndex )
      idatend = obs_headElem_i( lobsSpaceData, OBS_NLV, headerIndex ) + idata - 1
      idatyp  = obs_headElem_i( lobsSpaceData, OBS_ITY, headerIndex )

      if ( tvs_isIdBurpTovs(idatyp) .and. rmat_lnondiagr ) then

        count = 0

        do bodyIndex = idata, idatend
        
          iass = obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex )

          if (iass == 1 .or. iass == -1 ) then
            ichn = nint( obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex ))
            ichn = max( 0, min( ichn, tvs_maxChannelNumber + 1 ))
            count = count + 1
            list_chan( count ) = ichn
            x( count ) = obs_bodyElem_r( lobsSpaceData, elem_src_i, bodyIndex )
          end if
        
        end do

        if ( count > 0 ) then

          call rmat_sqrtRm1( tvs_lsensor( tvs_ltovsno( headerIndex )), count, x(1:count), y(1:count), list_chan(1:count), tvs_ltovsno(headerIndex) )

          count = 0
          do bodyIndex = idata, idatend

            iass = obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex )
            if ( iass == 1 .or. iass == -1) then
              count = count + 1
              call obs_bodySet_r(lobsSpaceData, elem_dest_i, bodyIndex,y(count))
            end if
          end do
        
        end if

      else

        do bodyIndex = idata, idatend
          iass = obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex )
          if (iass == 1 .or. iass == -1) then
            call obs_bodySet_r( lobsSpaceData, elem_dest_i, bodyIndex, &
                     obs_bodyElem_r( lobsSpaceData, elem_src_i, bodyIndex) / obs_bodyElem_r( lobsSpaceData, OBS_OER, bodyIndex ))
          end if
        end do
       end if

    end do

  end subroutine cfn_RsqrtInverse

  !--------------------------------------------------------------------------
  ! cfn_calcJo
  !--------------------------------------------------------------------------
  subroutine cfn_calcJo(lobsSpaceData)
    implicit none

    ! Compute JO contribution of each assimilated and diagnosed datum
    ! and store the result in OBS_JOBS

    type(struct_obs) :: lobsSpaceData

    integer :: bodyIndex

!$OMP PARALLEL DO PRIVATE(bodyIndex)
    do bodyIndex=1,obs_numbody(lobsSpaceData)

      if ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex) == 1) then
        call obs_bodySet_r( lobsSpaceData, OBS_JOBS, bodyIndex, &
             ( obs_bodyElem_r( lobsSpaceData, OBS_WORK, bodyIndex ) &
             * obs_bodyElem_r( lobsSpaceData, OBS_WORK, bodyIndex ) &
             ) / 2.d0 )
      else
        call obs_bodySet_r(lobsSpaceData,OBS_JOBS,bodyIndex, 0.0d0)
      end if

    end do
!$OMP END PARALLEL DO

  end subroutine cfn_calcJo

  !--------------------------------------------------------------------------
  ! cfn_sumJo
  !--------------------------------------------------------------------------
  subroutine cfn_sumJo( lobsSpaceData, pjo )

    ! Purpose:
    !   Compute the sum of Jo contributions saved in OBS_JOBS
    !   Also, compute contribution of each family of observation (for
    !   diagnostic purposes)

    implicit none

    real(8) :: pjo ! Total observation cost function

    type(struct_obs) :: lobsSpaceData
    integer :: bodyIndex, itvs, isens, headerIndex, idata, idatend

    real(8) :: dljoraob, dljoairep, dljosatwind, dljoscat, dljosurfc, dljotov, dljosst
    real(8) :: dljoprof, dljogpsro, dljogpsztd, dljochm, pjo_1
    real(8) :: dljotov_sensors( tvs_nsensors )
    integer :: ierr

    call tmg_start(81,'SUMJO')

    dljogpsztd = 0.d0
    dljoraob = 0.d0
    dljoairep = 0.d0
    dljosatwind = 0.d0
    dljosurfc = 0.d0
    dljoscat = 0.d0
    dljotov = 0.d0
    dljogpsro = 0.d0
    dljoprof = 0.d0
    dljochm = 0.d0
    dljosst = 0.0d0
    dljotov_sensors(:) = 0.d0

    do bodyIndex = 1, obs_numbody( lobsSpaceData )

      pjo_1 = obs_bodyElem_r( lobsSpaceData, OBS_JOBS, bodyIndex )

      ! total observation cost function
      pjo   = pjo + pjo_1

      ! subcomponents of observation cost function (diagnostic only)
      select case( obs_getFamily( lobsSpaceData, bodyIndex = bodyIndex ))
      case('UA')
        dljoraob    = dljoraob    + pjo_1
      case('AI')
        dljoairep   = dljoairep   + pjo_1
      case('SW')
        dljosatwind = dljosatwind + pjo_1
      case('SF')
        dljosurfc   = dljosurfc   + pjo_1
      case('SC')
        dljoscat    = dljoscat    + pjo_1
      case('TO')
        dljotov     = dljotov     + pjo_1
      case('RO')
        dljogpsro   = dljogpsro   + pjo_1
      case('PR')
        dljoprof    = dljoprof    + pjo_1
      case('GP')
        dljogpsztd  = dljogpsztd  + pjo_1
      case('CH')
        dljochm     = dljochm     + pjo_1
      case('TM')
        dljosst     = dljosst     + pjo_1
      end select
    enddo

    do itvs = 1, tvs_nobtov
      headerIndex = tvs_lobsno( itvs )
      if ( headerIndex > 0 ) then
        idata   = obs_headElem_i( lobsSpaceData, OBS_RLN, headerIndex )
        idatend = obs_headElem_i( lobsSpaceData, OBS_NLV, headerIndex ) + idata - 1
        do bodyIndex = idata, idatend
          pjo_1 = obs_bodyElem_r( lobsSpaceData, OBS_JOBS, bodyIndex)
          isens = tvs_lsensor (itvs)
          dljotov_sensors(isens) =  dljotov_sensors(isens) + pjo_1
        end do
      end if
    end do

    call mpi_allreduce_sumreal8scalar( pjo, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljoraob, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljoairep, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljosatwind, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljosurfc, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljoscat, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljotov, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljogpsro, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljoprof, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljogpsztd, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljochm, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljosst, "GRID" )
    do isens = 1, tvs_nsensors
       call mpi_allreduce_sumreal8scalar( dljotov_sensors(isens), "GRID" )
    end do

    if ( mpi_myid == 0 ) then
      write(*,'(a15,f25.17)') 'Jo(UA)   = ', dljoraob
      write(*,'(a15,f25.17)') 'Jo(AI)   = ', dljoairep
      write(*,'(a15,f25.17)') 'Jo(SF)   = ', dljosurfc
      write(*,'(a15,f25.17)') 'Jo(SC)   = ', dljoscat
      write(*,'(a15,f25.17)') 'Jo(TO)   = ', dljotov
      write(*,'(a15,f25.17)') 'Jo(SW)   = ', dljosatwind
      write(*,'(a15,f25.17)') 'Jo(PR)   = ', dljoprof
      write(*,'(a15,f25.17)') 'Jo(RO)   = ', dljogpsro
      write(*,'(a15,f25.17)') 'Jo(GP)   = ', dljogpsztd
      write(*,'(a15,f25.17)') 'Jo(CH)   = ', dljochm
      write(*,'(a15,f25.17)') 'Jo(TM)   = ', dljosst
      write(*,*) ' '
      if ( tvs_nsensors > 0 ) then
        write(*,'(1x,a)') 'For TOVS decomposition by sensor:'
        write(*,'(1x,a)') '#  plt sat ins    Jo'
        do isens = 1, tvs_nsensors
          write(*,'(i2,1x,a,1x,a,1x,i2,1x,f25.17)') isens, inst_name(tvs_instruments(isens)), &
               platform_name(tvs_platforms(isens)), tvs_satellites(isens), dljotov_sensors(isens)
        end do
        write(*,*) ' '
      end if
    end if

    call tmg_stop(81)

  end subroutine cfn_sumJo

end module costfunction_mod
