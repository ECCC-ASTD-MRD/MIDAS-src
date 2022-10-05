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

    type(struct_obs) :: obsSpaceData 
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

    type(struct_obs) :: obsSpaceData
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
