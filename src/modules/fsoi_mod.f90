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

module fsoi_mod 
  ! MODULE fsoi_mod (prefix='fso' category='1. High-level functionality')
  !
  ! :Purpose: Observation impact (FSOI) library
  !
  use mpi_mod
  use Vgrid_Descriptors
  use gridStateVector_mod
  use MathPhysConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public subroutines and functions
  public :: fso_multEnergyNorm

  contains

  !--------------------------------------------------------------------------
  ! fso_multEnergyNorm
  !--------------------------------------------------------------------------
  subroutine fso_multEnergyNorm(statevector_inout, statevector_ref,  &
                                latMin, latMax, lonMin, lonMax,      &
                                uvNorm,ttNorm,p0Norm,huNorm,tgNorm)
    !
    ! :Purpose: Computes energy norms
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)  :: statevector_inout
    type(struct_gsv), intent(in)     :: statevector_ref
    real(8),          intent(in)     :: latMin, latMax, lonMin, lonMax
    logical,          intent(in)     :: uvNorm, ttNorm, p0Norm, huNorm, tgNorm

    ! Locals:
    integer              :: stepIndex, lonIndex, levIndex, latIndex, lonIndex2, latIndex2, status, nLev_M, nLev_T
    real(8)              :: scaleFactor, scaleFactorConst, scaleFactorLat, scaleFactorLon, scaleFactorLev
    real(8)              :: pfac, tfac, qfac
    real(8)              :: sumScale , sumeu, sumev, sumep, sumet, sumeq
    real(8), pointer     :: field_UU(:,:,:,:), field_VV(:,:,:,:), field_T(:,:,:,:), field_LQ(:,:,:,:)
    real(8), pointer     :: field_Psfc(:,:,:,:), field_TG(:,:,:,:),Psfc_ptr(:,:,:)
    real(8), pointer     :: Press_T(:,:,:) 
    real(8), pointer     :: Press_M(:,:,:)
    real(8), allocatable :: Psfc_ref(:,:)
    real(8), parameter   :: T_r = 280.0D0
    real(8), parameter   :: Psfc_r = 100000.0D0 ! unit Pa

    if (mpi_myid == 0) write(*,*) 'fso_multEnergyNorm: START'
    nullify(Press_T,Press_M)

    ! the factors for TT, HU and Ps (for wind is 1)
    tfac = mpc_cp_dry_air_r8/T_r                                 ! temperature factor (c_p/T_r)
    qfac = mpc_heat_condens_water_r8**2/(mpc_cp_dry_air_r8*T_r)  ! humidity factor ( (l_p*l_p)/(c_p*T_r) )
    pfac = mpc_rgas_dry_air_r8*T_r/(Psfc_r**2)                   ! surface pressure factor (R*T_r/Psfc_r^2)

    if (.not. gsv_isAllocated(statevector_inout)) then
      call utl_abort('fso_multEnergyNorm: gridStateVector_inout not yet allocated')
    end if

    nLev_M = gsv_getNumLev(statevector_inout,'MM')
    nLev_T = gsv_getNumLev(statevector_inout,'TH')

    ! compute 3D log pressure fields
    call gsv_getField(statevector_ref,Psfc_ptr,'P0')
    allocate(Psfc_ref(statevector_inout%lonPerPEmax,statevector_inout%latPerPEmax))
    Psfc_ref(:,:) =  &
                  Psfc_ptr(statevector_inout%myLonBeg:statevector_inout%myLonEnd,  &
                  statevector_inout%myLatBeg:statevector_inout%myLatEnd, 1)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_T,  &
                        levels=Press_T,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_M,  &
                        levels=Press_M,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    ! dlat * dlon
    scaleFactorConst = statevector_inout%hco%dlat*statevector_inout%hco%dlon

    ! for wind components if to include in Norm calculation
    call gsv_getField(statevector_inout,field_UU,'UU')
    call gsv_getField(statevector_inout,field_VV,'VV')
    sumeu = 0.0D0
    sumev = 0.0D0
    sumScale = 0.0D0
    if (uvNorm) then
      do levIndex = 1, nLev_M
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if ( levIndex == nLev_M) then
                scaleFactorLev = Press_M(lonIndex2, latIndex2, nLev_M)-Press_T(lonIndex2, latIndex2, nLev_T-1) 
              else if ( Press_T(lonIndex2, latIndex2, levIndex) < 10000.0D0) then 
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_T(lonIndex2, latIndex2, levIndex+1) -  Press_T(lonIndex2, latIndex2, levIndex)
              end if

              scaleFactor = scaleFactorConst * scaleFactorLat* scaleFactorLon * scaleFactorLev
              sumScale = sumScale + scaleFactor

              sumeu = sumeu + &
                      0.5 * field_UU(lonIndex,latIndex,levIndex,stepIndex) * field_UU(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor
              sumev = sumev + &
                      0.5 * field_VV(lonIndex,latIndex,levIndex,stepIndex) * field_VV(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor

              field_UU(lonIndex,latIndex,levIndex,stepIndex) = &
                   field_UU(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor
              field_VV(lonIndex,latIndex,levIndex,stepIndex) = &
                   field_VV(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor
            end do !lonIndex
          end do !latIndex
        end do ! stepIndex
      end do ! levIndex
      call mpi_allreduce_sumreal8scalar(sumeu,'grid')
      call mpi_allreduce_sumreal8scalar(sumev,'grid')
      call mpi_allreduce_sumreal8scalar(sumScale,'grid')

      sumeu = sumeu/sumScale
      sumev = sumev/sumScale

      field_UU(:,:,:,:) = field_UU(:,:,:,:)/sumScale
      field_VV(:,:,:,:) = field_VV(:,:,:,:)/sumScale
    else
      field_UU(:,:,:,:) = field_UU(:,:,:,:)*0.0D0
      field_VV(:,:,:,:) = field_VV(:,:,:,:)*0.0D0
    end if ! if uvNorm

    if (mpi_myid == 0)  write(*,*) 'energy for UU=', sumeu
    if (mpi_myid == 0)  write(*,*) 'energy for VV=', sumev

    ! for Temperature
    call gsv_getField(statevector_inout,field_T,'TT')
    sumScale = 0.0D0
    sumet = 0.0D0
    if (ttNorm) then
      do levIndex = 1, nLev_T
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if (levIndex == nLev_T) then  !surface
                scaleFactorLev =  Press_T(lonIndex2, latIndex2, nLev_T)-Press_T(lonIndex2, latIndex2, nLev_T-1) 
              else if (levIndex == 1)  then  ! top
                scaleFactorLev = 0.0D0
              else if ( Press_M(lonIndex2, latIndex2, levIndex-1) < 10000.0D0) then 
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_M(lonIndex2, latIndex2, levIndex ) - Press_M(lonIndex2, latIndex2, levIndex-1)
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon * scaleFactorLev
              sumet = sumet + &
                   0.5 * tfac * field_T(lonIndex,latIndex,levIndex,stepIndex) * field_T(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor
              sumScale = sumScale + scaleFactor
              field_T(lonIndex,latIndex,levIndex,stepIndex) = &
                           field_T(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * tfac * scaleFactor
            end do
          end do
        end do ! stepIndex
      end do ! levIndex
      call mpi_allreduce_sumreal8scalar(sumet,'grid')
      call mpi_allreduce_sumreal8scalar(sumScale,'grid')
      sumet = sumet/sumScale
      
      field_T(:,:,:,:) = field_T(:,:,:,:)/sumScale
    else 
      field_T(:,:,:,:) = field_T(:,:,:,:)*0.0D0
    end if ! if ttNorm

    if (mpi_myid == 0)  write(*,*) 'energy for TT=', sumet


    ! humidity (set to zero, for now)
    call gsv_getField(statevector_inout,field_LQ,'HU')
    sumScale = 0.0D0
    sumeq = 0.0D0
    if (huNorm) then
      do levIndex = 1, nLev_T
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            latIndex2 = latIndex - statevector_inout%myLatBeg + 1
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              lonIndex2 = lonIndex - statevector_inout%myLonBeg + 1
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              ! do all thermo levels for which there is a momentum level above and below
              if ( levIndex == nLev_T) then !surface
                scaleFactorLev =  Press_T(lonIndex2, latIndex2, nLev_T) - Press_T(lonIndex2, latIndex2, nLev_T-1) 
              else if (levIndex == 1)  then  ! top
                scaleFactorLev = 0.0D0
              else if ( Press_M(lonIndex2, latIndex2, levIndex-1) < 10000.0D0) then 
                scaleFactorLev = 0.0D0
              else
                scaleFactorLev = Press_M(lonIndex2, latIndex2, levIndex ) - Press_M(lonIndex2, latIndex2, levIndex-1)
              end if

              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon * scaleFactorLev
              sumScale = sumScale + scaleFactor

              sumeq = sumeq + 0.5 * qfac * &
                    field_LQ(lonIndex,latIndex,levIndex,stepIndex) * field_LQ(lonIndex,latIndex,levIndex,stepIndex) * scaleFactor

              field_LQ(lonIndex,latIndex,levIndex,stepIndex) = &
                       field_LQ(lonIndex,latIndex,levIndex,stepIndex) * 0.5 * scaleFactor * qfac * 0.0

            end do
          end do
        end do ! stepIndex
      end do ! latIndex
      call mpi_allreduce_sumreal8scalar(sumScale,'grid')
      field_LQ(:,:,:,:) = field_LQ(:,:,:,:)/sumScale
    else 
      field_LQ(:,:,:,:) = field_LQ(:,:,:,:)*0.0D0
    end if ! if huNorm

    if (mpi_myid == 0)  write(*,*) 'energy for HU=', sumeq

    ! surface pressure
    call gsv_getField(statevector_inout,field_Psfc,'P0')
    sumScale = 0.0D0
    sumep = 0.0
    if (p0Norm) then
      do stepIndex = 1, statevector_inout%numStep
        do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon
              sumScale = sumScale + scaleFactor
              sumep = sumep + 0.5 * pfac * &
                  field_Psfc(lonIndex,latIndex,1,stepIndex) * field_Psfc(lonIndex,latIndex,1,stepIndex) * scaleFactor
              field_Psfc(lonIndex,latIndex,1,stepIndex) = &
              field_Psfc(lonIndex,latIndex,1,stepIndex) * 0.5 * scaleFactor * pfac
          end do
        end do ! latIndex
      end do ! stepIndex

      call mpi_allreduce_sumreal8scalar(sumep,'grid')
      call mpi_allreduce_sumreal8scalar(sumScale,'grid')
      sumep = sumep/sumScale
      field_Psfc(:,:,:,:) =  field_Psfc(:,:,:,:)/sumScale
    else
      field_Psfc(:,:,:,:) =  field_Psfc(:,:,:,:)*0.0D0
    end if ! if p0Norm

    if (mpi_myid == 0)  write(*,*) 'energy for Ps=', sumep


    ! skin temperature (set to zero for now)
    call gsv_getField(statevector_inout,field_TG,'TG')
    sumScale = 0.0D0
    if (tgNorm) then
      do stepIndex = 1, statevector_inout%numStep
        do latIndex = statevector_inout%myLatBeg, statevector_inout%myLatEnd
            ! IF lat is out of the domain where we want to compute the NRJ norm, we put scaleFactorLat = 0.
            if (statevector_inout%hco%lat(latIndex) >= latMin .and. statevector_inout%hco%lat(latIndex) <= latMax) then
              scaleFactorLat = cos(statevector_inout%hco%lat(latIndex))
            else
              scaleFactorLat = 0.0D0
            end if
            do lonIndex = statevector_inout%myLonBeg, statevector_inout%myLonEnd
              ! Similarly, if lon is out of the domain where we want to compute the NRJ norm, we put scaleFactorLon = 0.
              if (statevector_inout%hco%lon(lonIndex) >= lonMin .and. statevector_inout%hco%lon(lonIndex) <= lonMax) then
                scaleFactorLon = 1.0D0
              else
                scaleFactorLon = 0.0D0
              end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLon
              sumScale = sumScale + scaleFactor
              field_TG(lonIndex,latIndex,1,stepIndex) = &
              field_TG(lonIndex,latIndex,1,stepIndex) * 0.5 * scaleFactor * 0.0
          end do
        end do ! latIndex
      end do ! stepIndex
      call mpi_allreduce_sumreal8scalar(sumScale,'grid')
      field_TG(:,:,:,:) = field_TG(:,:,:,:)/sumScale 
    else
      field_TG(:,:,:,:) = field_TG(:,:,:,:)*0.0D0
    end if ! if tgNorm

    if (mpi_myid == 0) write(*,*) 'energy for total=', sumeu + sumev + sumet + sumep
    deallocate(Press_T,Press_M)
    deallocate(Psfc_ref)

    if (mpi_myid == 0) write(*,*) 'fso_multEnergyNorm: END'

  end subroutine fso_multEnergyNorm

end module fsoi_mod
