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
!! *Purpose*: Horizontal interpolation of the model variables
!!            in grid-point space to observation points.
!!            Bilinear interpolation from the 4 nearest grid points.
!!
!--------------------------------------------------------------------------
SUBROUTINE BILIN(lcolumn,statevector,lobsSpaceData)
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod
  use gridStateVector_mod
  IMPLICIT NONE
  
  type(struct_columnData) :: lcolumn
  type(struct_gsv) :: statevector
  type(struct_obs) :: lobsSpaceData
  
  INTEGER   JLEV, JK, JK2, JGL, JLON, JOBS
  INTEGER   ILON, ILA, IERR
  
  REAL*8  :: LAT, LON
  REAL*8  :: DLDY, DLW1, DLW2, DLW3, DLW4, DLDX, ypos, xpos
  real(8) :: latrot, lonrot
  
  real*8, allocatable ::zgd(:,:,:)
  
  real*8, pointer :: uu_column(:),vv_column(:),hu_column(:)
  real*8, pointer :: tt_column(:),tr_column(:),ps_column(:),tg_column(:)
  real*8, pointer :: field_ptr(:,:,:), uu_ptr(:,:,:), vv_ptr(:,:,:)
  
  ! Note: We assume here the all the obs between the poles and the last grid points
  !       (i.e. outside the grid) have been moved within the grid by suprep

  allocate(zgd(statevector%ni+1,statevector%nj,statevector%nk))
  
  zgd(:,:,:)=0.0d0
  field_ptr => gsv_getField3D_r8(statevector)
  zgd(1:statevector%ni,1:statevector%nj,1:statevector%nk)= &
       field_ptr(1:statevector%ni,1:statevector%nj,1:statevector%nk)

  !
  !- 1.  EXPAND Field BY REPEATING MERIDIAN 1 into INTO MERIDIAN NI+1
  !
  do jk = 1, statevector%nk
     do jgl = 1, statevector%nj
        zgd(statevector%ni+1,jgl,jk) = zgd( 1,jgl,jk)
     end do
  end do
  
  !
  !- 2.  LOOP OVER ALL THE OBSERVATIONS
  !
  do jobs = 1, col_getNumCol(lcolumn)

     !- 2.1 Find the obs positin within the analysis grid
     call col_getLatLon( lcolumn, jobs,                       & ! IN
                         Lat, Lon, ypos, xpos, LatRot, LonRot ) ! OUT

     !- Make sure we are within bounds
     if ( ypos < 1.d0                           .or. &
          ypos > real(statevector%nj      ,8) .or. &
          xpos < 1.d0                           .or. &
          xpos > real(statevector%ni + 1  ,8) ) then
       write(*,*) 'bilin: Obs outside local domain for job = ', jobs
       write(*,*) '  obs    lat, lon position            = ', Lat*MPC_DEGREES_PER_RADIAN_R8, Lon*MPC_DEGREES_PER_RADIAN_R8
       write(*,*) '  obs    x, y     position            = ', xpos, ypos
       write(*,*) '  domain x_end, y_end bounds          = ', statevector%ni + 1, statevector%nj
       stop
     end if

     !- 2.2 Find the lower-left grid point next to the observation
     if ( xpos /= real(statevector%ni + 1,8) ) then
       ILON = floor(xpos)
     else
       ILON = floor(xpos) - 1
     end if

     if ( ypos /= real(statevector%nj,8) ) then
       ILA = floor(ypos)
     else
       ILA = floor(ypos) - 1
     end if

     !- 2.3 COMPUTE THE 4 WEIGHTS OF THE BILINEAR INTERPOLATION
     dldx = xpos - real(ILON,8)
     dldy = ypos - real(ILA,8)

     dlw1 = (1.d0-dldx) * (1.d0-dldy)
     dlw2 =       dldx  * (1.d0-dldy)
     dlw3 = (1.d0-dldx) *       dldy
     dlw4 =       dldx  *       dldy
     
     !- 2.4 Interpolate the model state to the obs point
     if(col_varExist('UU')) uu_column => col_getColumn(lcolumn,jobs,'UU')
     if(col_varExist('VV')) vv_column => col_getColumn(lcolumn,jobs,'VV')
     if(col_varExist('HU')) hu_column => col_getColumn(lcolumn,jobs,'HU')
     if(col_varExist('TT')) tt_column => col_getColumn(lcolumn,jobs,'TT')
     if(col_varExist('P0')) ps_column => col_getColumn(lcolumn,jobs,'P0')
     if(col_varExist('TG')) tg_column => col_getColumn(lcolumn,jobs,'TG')
     
     do jk = 1, gsv_getNumLev(statevector,'MM')
        if(gsv_varExist(statevector,'UU')) then
           jk2=jk+gsv_getOffsetFromVarName(statevector,'UU')
           uu_column(jk) =   dlw1*zgd(ilon  ,ila,jk2)  &
                           + dlw2*zgd(ilon+1,ila,jk2)  &
                           + dlw3*zgd(ilon  ,ila+1,jk2)  &
                           + dlw4*zgd(ilon+1,ila+1,jk2)
        endif
        if(gsv_varExist(statevector,'VV')) then
           jk2=jk+gsv_getOffsetFromVarName(statevector,'VV')
           vv_column(jk) =   dlw1*zgd(ilon  ,ila,jk2)  &
                           + dlw2*zgd(ilon+1,ila,jk2)  &
                           + dlw3*zgd(ilon  ,ila+1,jk2)  &
                           + dlw4*zgd(ilon+1,ila+1,jk2)
        endif
     enddo
     do jk = 1, gsv_getNumLev(statevector,'TH')
        if(gsv_varExist(statevector,'HU')) then
           jk2=jk+gsv_getOffsetFromVarName(statevector,'HU')
           hu_column(jk) =   dlw1*zgd(ilon  ,ila,jk2)  &
                           + dlw2*zgd(ilon+1,ila,jk2)  &
                           + dlw3*zgd(ilon  ,ila+1,jk2)  &
                           + dlw4*zgd(ilon+1,ila+1,jk2)
        endif
        if(gsv_varExist(statevector,'TT')) then
           jk2=jk+gsv_getOffsetFromVarName(statevector,'TT')
           tt_column(jk) =   dlw1*zgd(ilon  ,ila,jk2)  &
                           + dlw2*zgd(ilon+1,ila,jk2)  &
                           + dlw3*zgd(ilon  ,ila+1,jk2)  &
                           + dlw4*zgd(ilon+1,ila+1,jk2)
        endif
     enddo
     if(gsv_varExist(statevector,'P0')) then
        jk2=1+gsv_getOffsetFromVarName(statevector,'P0')
        ps_column(1) =   dlw1*zgd(ilon  ,ila,jk2)  &
                       + dlw2*zgd(ilon+1,ila,jk2)  &
                       + dlw3*zgd(ilon  ,ila+1,jk2)  &
                       + dlw4*zgd(ilon+1,ila+1,jk2)
     endif
     if(gsv_varExist(statevector,'TG')) then
        jk2=1+gsv_getOffsetFromVarName(statevector,'TG')
        tg_column(1) =   dlw1*zgd(ilon  ,ila,jk2)  &
                       + dlw2*zgd(ilon+1,ila,jk2)  &
                       + dlw3*zgd(ilon  ,ila+1,jk2)  &
                       + dlw4*zgd(ilon+1,ila+1,jk2)
     endif
  end do

  deallocate(zgd)
  
END SUBROUTINE BILIN