
module residual_mod
  ! MODULE residual_mod (prefix='res' category='1. High-level functionality')
  !
  ! :Purpose: To compute OMA (= OMP - H dx) and its adjoint
  !
  use obsSpaceData_mod
  implicit none
  save
  private

  ! public procedures
  public :: res_compute , res_computeAd

contains

  subroutine res_compute(obsSpaceData)
    !
    !:Purpose: Computes residual of observation - analysis from Hdx.
    !          Takes as input OBS_WORK (observation increment Hdx) and 
    !          OBS_OMP (innovation) and computes OBS_OMA (observation 
    !          - analysis)
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData 

    ! Locals:
    integer :: index_body

    !$OMP PARALLEL DO PRIVATE(index_body)
    do index_body=1,obs_numbody(obsSpaceData)
      if(obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obs_assimilated) then
        call obs_bodySet_r(obsSpaceData,OBS_OMA,index_body, &
             obs_bodyElem_r(obsSpaceData,OBS_OMP,index_body) &
             -obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body))
      endif
    enddo
    !$OMP END PARALLEL DO

  end subroutine res_compute

  subroutine res_computeAd(obsSpaceData)
    !
    !:Purpose: Adjoint of computing residuals to observations.
    !          OBS_WORK contains input and output.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer :: index_body
    
    !$OMP PARALLEL DO PRIVATE(index_body)
    do index_body=1,obs_numbody(obsSpaceData)
      if(obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obs_assimilated) then
        call obs_bodySet_r(obsSpaceData,OBS_WORK,index_body, &
             -obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body))
      endif
    enddo
    !$OMP END PARALLEL DO

  end subroutine res_computeAd

end module residual_mod
