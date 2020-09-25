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

module costfunction_mod
  ! MODULE costfunction_mod, (prefix="cfn" category='1. High-level functionality')
  !
  ! :Purpose: To compute Jo term
  !
  use mpi_mod
  use mpivar_mod
  use obsSpaceData_mod
  use rmatrix_mod
  use rttov_const, only : inst_name, platform_name
  use tovs_nl_mod
  use codeprecision_mod
  use MathPhysConstants_mod
  use utilities_mod
  use obserrors_mod

  implicit none

  save
  private

  ! public procedures

  public :: cfn_calcJo, cfn_sumJo, cfn_computeNlTovsJo

contains

  !--------------------------------------------------------------------------
  ! cfn_calcJo
  !--------------------------------------------------------------------------
  subroutine cfn_calcJo(lobsSpaceData)
    !
    !:Purpose: To compute JO contribution of each assimilated and diagnosed
    !          datum, and to store the result in OBS_JOBS
    implicit none

    ! Arguments:
    type(struct_obs) :: lobsSpaceData

    ! Locals:
    integer :: bodyIndex

!$OMP PARALLEL DO PRIVATE(bodyIndex)
    do bodyIndex=1,obs_numbody(lobsSpaceData)

      if ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
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
    !
    !:Purpose: To compute the sum of Jo contributions saved in OBS_JOBS. Also,
    !          to compute contribution of each family of observation (for
    !          diagnostic purposes)
    implicit none

    ! Arguments:
    type(struct_obs) :: lobsSpaceData
    real(8) :: pjo ! Total observation cost function

    ! Locals:
    integer :: bodyIndex, itvs, isens, headerIndex, idata, idatend

    real(8) :: dljoraob, dljoairep, dljosatwind, dljoscat, dljosurfc, dljotov, dljosst, dljoice
    real(8) :: dljoprof, dljogpsro, dljogpsztd, dljochm, pjo_1, dljoaladin, dljohydro, dljoradar
    real(8) :: dljotov_sensors( tvs_nsensors )

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
    dljoaladin = 0.d0
    dljoice = 0.0d0
    dljotov_sensors(:) = 0.d0
    dljohydro = 0.0d0
    dljoradar = 0.0d0

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
      case('AL')
        dljoaladin  = dljoaladin  + pjo_1
      case('GL')
        dljoice     = dljoice     + pjo_1
      case('HY')
        dljohydro   = dljohydro   + pjo_1
      case('RA')
        dljoradar   = dljoradar   + pjo_1
      end select
    enddo

    do itvs = 1, tvs_nobtov
      headerIndex = tvs_headerIndex( itvs )
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
    call mpi_allreduce_sumreal8scalar( dljoaladin,"GRID")
    call mpi_allreduce_sumreal8scalar( dljoice, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljohydro, "GRID" )
    call mpi_allreduce_sumreal8scalar( dljoradar, "GRID" )
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
      write(*,'(a15,f25.17)') 'Jo(AL)   = ', dljoaladin
      write(*,'(a15,f25.17)') 'Jo(GL)   = ', dljoice
      write(*,'(a15,f25.17)') 'Jo(HY)   = ', dljohydro
      write(*,'(a15,f25.17)') 'Jo(RA)   = ', dljoradar
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


  !--------------------------------------------------------------------------
  !  cfn_computeNlTovsJo
  !--------------------------------------------------------------------------
  subroutine cfn_computeNlTovsJo(pjo, obsSpaceData, dest_obs)
    !
    ! *Purpose*: Computation of Jo for TOVS observations and transfer O-F
    !                     to obsSpaceData data column dest_obs
    !
    implicit none

    !Arguments:
    real(8),          intent(out)   :: pjo          ! Computed Tovs observation cost fucntion
    type(struct_obs), intent(inout) :: obsSpaceData! obsSpacaData structure
    integer,          intent(in)    :: dest_obs     ! obsSpaceData body destinationcolumn (ex: OBS_OMP or OBS_OMA)

    ! Locals:
    integer :: sensorIndex, channelIndex, tovsIndex
    
    real(pre_obsReal) :: zdtb
    integer :: idatyp, channelNumber, ichobs_a
    integer :: headerIndex, bodyIndex, count
    real(8) :: obsIn(tvs_maxChannelNumber), obsOut(tvs_maxChannelNumber)
    real(8) :: sigmaObs, dlsum
    integer :: list_chan(tvs_maxChannelNumber)
    real(8) :: list_oer(tvs_maxChannelNumber)

    pjo = 0.d0

    if ( tvs_nobtov == 0) return    ! exit if there are not tovs data

    ! 1.  Computation of (hx - z)/sigma for tovs data only

    dlsum    = 0.d0
    

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! Extract general information for this observation point
      !      ------------------------------------------------------

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) cycle HEADER

      tovsIndex = tvs_tovsIndex(headerIndex)
      if ( tovsIndex == -1 ) cycle HEADER
       
      sensorIndex = tvs_lsensor(tovsIndex)

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      count = 0
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        ! Only consider if flagged for assimilation
        if ( obs_bodyElem_i(obsSpaceData,obs_aSS,bodyIndex) /= obs_assimilated ) cycle BODY                

        channelNumber = nint(obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex))
        channelNumber = max( 0 , min( channelNumber , tvs_maxChannelNumber + 1))
        ichobs_a = channelNumber
        channelNumber = channelNumber - tvs_channelOffset(sensorIndex)
        channelIndex = utl_findArrayIndex(tvs_ichan(:,sensorIndex),tvs_nchan(sensorIndex),channelNumber)
        if ( channelIndex == 0 ) then
          write(*,'(A)') ' cfn_computeNlTovsJo: error with channel number'
          call utl_abort('cfn_computeNlTovsJo')
        end if

        zdtb = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex) - &
             tvs_radiance (tovsIndex) % bt(channelIndex)
        if ( tvs_debug ) then
          write(*,'(a,i4,2f8.2,f6.2)') ' channelNumber,sim,obs,diff= ', &
               channelNumber,  tvs_radiance (tovsIndex) % bt(channelIndex), &
               obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex), -zdtb
        end if
        call obs_bodySet_r(obsSpaceData,dest_obs,bodyIndex, zdtb)

        call oer_computeInflatedStateDepSigmaObs(obsSpaceData, headerIndex, bodyIndex, sensorIndex, dest_obs, beSilent_opt=.false.)

        sigmaObs = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)

        if ( sigmaObs == MPC_missingValue_R8) cycle body

        ! Comment out the modification of Jobs due to varqc for now, since this is probably
        ! only needed for use of nonlinear obs operator in minimization, which is not yet
        ! functional, but this interferes with doing ensemble of analyses (M. Buehner, Dec. 2013)
        !if (.not. min_lvarqc .or. obs_bodyElem_r(obsSpaceData,OBS_POB,bodyIndex).eq.0.0d0) then
        !dlsum =  dlsum &
        !     + (zdtb * zdtb) / (2.d0 * sigmaObs * sigmaObs)
        !else
        !  compute contribution of data with varqc
        !   zgami = obs_bodyElem_r(obsSpaceData,OBS_POB,bodyIndex)
        !   zjon = (zdtb* &
        !           zdtb)/2.d0
        !   zqcarg = zgami + exp(-1.0d0*zjon)
        !   dlsum= dlsum - log(zqcarg/(zgami+1.d0))
        !end if

        count = count + 1
        obsIn(count) = zdtb
        if (.not. rmat_lnondiagr) obsOut(count) = obsIn(count) / sigmaObs
        list_chan(count) = channelNumber
        list_oer(count) = sigmaObs
       
      end do BODY

      if (count >0) then
        if (rmat_lnondiagr) then
          call rmat_RsqrtInverseOneObs(sensorIndex,count,obsIn(1:count),obsOut(1:count),list_chan(1:count),list_oer(1:count),tovsIndex)
        end if
        dlsum =  dlsum + 0.5d0 * dot_product(obsOut(1:count), obsOut(1:count))
      end if

    end do HEADER

    pjo = dlsum

  end subroutine cfn_computeNlTovsJo


end module costfunction_mod
