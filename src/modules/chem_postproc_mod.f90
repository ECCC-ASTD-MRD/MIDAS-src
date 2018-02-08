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
!! MODULE chem_postproc_mod
!!
!! *Purpose*: Provides post-analysis routines.
!!
!! @author Mike Sitwell and Yves Rochon (ARQI/AQRD)
!!
!! Public routines:
!!v       - "chm_transform_final_increments" for any transformations or boundary
!!v          values to apply to final increments.
!--------------------------------------------------------------------------
module chem_postproc_mod

  use utilities_mod
  use chem_setup_mod, only: chm_var_maxnumber, chm_setup_get_str, &
                            chm_setup_get_int, chm_setup_get_float
  use gridStateVector_mod
  use verticalCoord_mod
  use HorizontalCoord_mod
  use varNameList_mod
  use bmatrixchem_mod, only: bchm_getBgSigma
 
  implicit none
  private

! public procedures
! -----------------

  public :: chm_transform_final_increments

contains

!--------------------------- MODEL and ANALYSIS SPACE POST-PROC --------------------
!-----------------------------------------------------------------------------------

!--------------------------------------------------------------------------
!! *Purpose*: Apply any required adjustments and/or transformations to increments
!!            prior to saving in rebm file.
!!
!! @author M. Sitwell Sept 2015
!!
!! InOut
!!
!!v     statevector_increment     statevector for the increment (with same resolution as the background)
!!
!! Comments:
!!v    - Possible adjustments/transformations:
!!v        - Convert to dx when the analysis increments are in a different form (e.g. dlnx).
!!v        - Checks for negative analysis values and if any negative values found will
!!v          modify the increment so that the analysis  is reset to low_cutoff*background
!!v          at these points. This check is done with fields of the same resolution as
!!v          the increment. Correspondingly impose upper limits.
!!
!! Revisions:
!!v       Y. Rochon, Nov 2015
!!v         - Added upper bound.
!!v         - Use of col_varExist.
!--------------------------------------------------------------------------
  subroutine chm_transform_final_increments(statevector_increment,hco_anl,vco_anl)

    implicit none
    
    type(struct_gsv), intent(inout) :: statevector_increment
    type(struct_hco), pointer, intent(inout) :: hco_anl
    type(struct_vco), pointer, intent(inout) :: vco_anl

    type(struct_gsv) :: statevector_trial
    integer :: jvar,unit,ier
    character(len=4) :: varName

    integer, external :: fclos

    write(*,*) "chm_transform_final_increments: Beginning final increment transformations."

    call gsv_readTrials(hco_anl,vco_anl,statevector_trial)

    do jvar = 1, vnl_numvarmax
      if (gsv_varExist(varName=vnl_varNameList(jvar))) then 
         if (vnl_varKindFromVarname(vnl_varNameList(jvar)).eq.'CH') then

            varName=vnl_varNameList(jvar)

            ! Check that the background and increment fields have the same index boundaries

            if (statevector_trial%myLonBeg.ne.statevector_increment%myLonBeg .or. statevector_trial%myLonEnd.ne.statevector_increment%myLonEnd .or. &
                statevector_trial%myLatBeg.ne.statevector_increment%myLatBeg .or. statevector_trial%myLatEnd.ne.statevector_increment%myLatEnd .or. &
                statevector_trial%varNumLev(vnl_varListIndex(varName)).ne.statevector_increment%varNumLev(vnl_varListIndex(varName)) .or. &
                statevector_trial%numStep.ne.statevector_increment%numStep) then

               call utl_open_asciifile(chm_setup_get_str('message'),unit)

               write(unit,'(A)') "chm_transform_final_increments: Warning - background and increment field index boundaries not equal."
               write(unit,'(A)') "   index      background     increment"
               write(unit,'(A,6X,I5,10X,I5)') 'myLonBeg ',statevector_trial%myLonBeg,statevector_increment%myLonBeg
               write(unit,'(A,6X,I5,10X,I5)') 'myLonEnd ',statevector_trial%myLonEnd,statevector_increment%myLonEnd
               write(unit,'(A,6X,I5,10X,I5)') 'myLatBeg ',statevector_trial%myLatBeg,statevector_increment%myLatBeg
               write(unit,'(A,6X,I5,10X,I5)') 'myLatEnd ',statevector_trial%myLatEnd,statevector_increment%myLatEnd
               write(unit,'(A,6X,I5,10X,I5)') 'varNumLev',statevector_trial%varNumLev(vnl_varListIndex(varName)),statevector_increment%varNumLev(vnl_varListIndex(varName))
               write(unit,'(A,6X,I5,10X,I5)') 'numStep  ',statevector_trial%numStep,statevector_increment%numStep

               ier=fclos(unit)

            end if

            ! Analysis increment transformations to dx if needed.
           
            call chm_apply_transform(statevector_trial, varName, l_reverse=.true., mode=2, increment_opt=statevector_increment)

            ! Apply min and max boundaries if needed.

            call chm_apply_bounds(statevector_trial,statevector_increment,varName)
             
         end if
      end if
    end do
    
    call gsv_deallocate(statevector_trial)

    write(*,*) "chm_transform_final_increments: Ending final increment transformations."

  end subroutine chm_transform_final_increments

!--------------------------------------------------------------------------
!! *Purpose*: Apply transform (or its inverse) of 4D background field or increment field.
!!
!! @author Y. Rochon, Nov 2015
!!
!! Input
!!
!!v     varName          Field name (nomvar)
!!v     l_reverse        Reverse/inverse transformation
!!v     mode             selected sub-transformation type (defined for each transformation
!!v                      given in chm_setup_get_int('transform',*)
!!
!! InOut
!!
!!v     background       statevector for the background
!!v     increment        statevector for the increment (with same resolution as the background)
!!
!! Revisions:
!!v        M. Sitwell, April 2016
!!v          - Added input integer mode for selection of transform sub-type.
!!v          - Create new sub-function log_matrix from previous code for handling
!!v            negative values before taking the log. 
!--------------------------------------------------------------------------
  subroutine chm_apply_transform(background,varName,l_reverse,mode,increment_opt)

    implicit none

    type(struct_gsv), intent(inout) :: background
    type(struct_gsv), intent(inout), optional :: increment_opt
    character(len=*), intent(in) :: varName
    integer, intent(in) :: mode
    logical, intent(in) :: l_reverse

    real(8), pointer :: background_field(:,:,:,:)        
    real(8), pointer :: increment_field(:,:,:,:)   

    integer :: lon1,lon2,lev1,lev2,lat1,lat2,step1,step2
    integer :: iconstituent_id

    iconstituent_id = vnl_varnumFromVarname(varName)
    if (iconstituent_id.lt.0.or.iconstituent_id.gt.chm_var_maxnumber()) return
    
    select case(chm_setup_get_int('transform',iconstituent_id))
    case(0)
       write(*,'(A,A)') "chm_apply_transform: No transform to be applied for field ",trim(varName)
       return
    case(1)
    
       write(*,'(A,I2,A,A)') "chm_apply_transform: Applying logarithm transform with mode = ",mode," for field ",trim(varName)

       ! Transform lnx to/from x or dlnx to/from dx
       
       if (.not.present(increment_opt).and.(mode.eq.1.or.mode.eq.2)) then
          write(*,'(A,I2)') "chm_apply_transform: increment must be provided for mode = ",mode
          call utl_abort("chm_apply_transform")
       end if
       
       ! Get the pointers to the fields
       background_field => gsv_getField_r8(background,varName)
       if (present(increment_opt)) increment_field => gsv_getField_r8(increment_opt,varName)

       if (.not.l_reverse) then

          ! Forward transformation
          
          ! Get lat,lon,time,height index range, assumed the same for the background and increment
          lon1 = background%myLonBeg
          lon2 = background%myLonEnd
          lat1 = background%myLatBeg
          lat2 = background%myLatEnd
          lev1 = 1
          lev2 = background%varNumLev(vnl_varListIndex(varName))
          step1 = 1
          step2 = background%numStep

          select case(mode)
          case(0)
             background_field = log_field(background_field,lon1,lon2,lev1,lev2,lat1,lat2,step1,step2)
          case(1)
             background_field = log_field(background_field+increment_field,lon1,lon2,lev1,lev2,lat1,lat2,step1,step2)
          case(2)
             increment_field = log_field(background_field+increment_field,lon1,lon2,lev1,lev2,lat1,lat2,step1,step2) &
                             - log_field(background_field,lon1,lon2,lev1,lev2,lat1,lat2,step1,step2)           
          end select

       else
    
          ! Reverse transformation
       
          select case(mode)
          case(0)
             background_field = exp(background_field)
          case(1)
             background_field = exp(background_field + increment_field)
          case(2)
             increment_field = exp(background_field + increment_field) - exp(background_field)
          end select 
             
       end if

    case default
       write(*,'(A,I2,A)') "chm_apply_transform: transform selection ", &
               chm_setup_get_int('transform',iconstituent_id)," not currenly defined." 
       call utl_abort("chm_apply_transform")
    end select
    
  contains

    function log_field(field,lon1,lon2,lev1,lev2,lat1,lat2,step1,step2)
    !
    ! Author: Y. Rochon, Nov 2015 (made into function by M. Sitwell)
    !
    ! Purpose: Helper function for taking the log of a field that might contain negative values.
    !          Places where negative values occur will be set as the log of the minimum positive
    !          value along the longitudinal dimension.
    !
    !-----------------------------------------------------------------------------------

      implicit none

      integer, intent(in) :: lon1,lon2,lev1,lev2,lat1,lat2,step1,step2
      real(8), intent(in) :: field(lon1:lon2,lev1:lev2,lat1:lat2,step1:step2)
      real(8) :: log_field(lon1:lon2,lev1:lev2,lat1:lat2,step1:step2)

      integer :: jstep,jlon,jlev,jlat
      real(8) :: valmin
      
      do jstep=step1,step2
         do jlev=lev1,lev2
            do jlat=lat1,lat2
                
               valmin = minval(field(lon1:lon2,jlat,jlev,jstep), mask=field(lon1:lon2,jlat,jlev,jstep).gt.0.)
               if (valmin.gt.1.E30) valmin=1.E-20

               do jlon=lon1,lon2
                  if (field(jlon,jlev,jlat,jstep).gt.0) then
                     log_field(jlon,jlat,jlev,jstep) = log(field(jlon,jlat,jlev,jstep))
                  else
                     log_field(jlon,jlat,jlev,jstep) = log(valmin)
                  end if
               end do

            end do
         end do
      end do

    end function log_field

  end subroutine chm_apply_transform

!--------------------------------------------------------------------------
!! *Purpose*: Checks for negative analysis values and if any negative values found will
!!            modify the increment so that the analysis  is reset to low_cutoff*background
!!            at these points. This check is done with fields of the same resolution as
!!            the increment. Also applies upper bound.
!!
!! @author M. Sitwell Sept 2015
!!
!! Input
!!
!!v    varName        Field name (nomvar)
!!
!! InOut
!!
!!v    background       statevector for the background
!!v    increment        statevector for the increment (with same resolution as the background)
!!
!! Revisions:
!!v       Y. Rochon, Nov 2015, Mar 2017
!!v          - Modified to chm_apply_bounds
!--------------------------------------------------------------------------
  subroutine chm_apply_bounds(background,increment,varName)

    implicit none
    
    type(struct_gsv), intent(inout) :: background,increment
    character(len=*), intent(in) :: varName

    real(8), pointer :: background_field(:,:,:,:),increment_field(:,:,:,:)
    real(8) :: bkgrnd,inc,new_inc,refval,bg_sigma
    integer :: iconstituent_id,count,ier,unit,jvar
    integer :: jstep,jlon,jlev,jlat,lon1,lon2,lat1,lat2,step1,step2,lev1,lev2

    integer, external :: fclos
    
    iconstituent_id = vnl_varnumFromVarname(varName)
    if (iconstituent_id.lt.0.or.iconstituent_id.gt.chm_var_maxnumber()) return

    ! Open output file
    call utl_open_asciifile(chm_setup_get_str('message'),unit)

    ! Get lat,lon,time,height index range, assumed the same for the background and increment
    lon1 = background%myLonBeg
    lon2 = background%myLonEnd
    lat1 = background%myLatBeg
    lat2 = background%myLatEnd
    lev1 = 1
    lev2 = background%varNumLev(vnl_varListIndex(varName))
    step1 = 1
    step2 = background%numStep

    ! Get the pointers to the fields
    background_field => gsv_getField_r8(background,varName)
    increment_field => gsv_getField_r8(increment,varName)

    refval=chm_setup_get_float('sigma_cutoff',iconstituent_id)
    if (refval.gt.0.0D0.and.refval.lt.0.11) then

       ! Find the field index within the background error covariance matrix

       jvar=0
       ier=-1
       do jstep = 1, vnl_numvarmax   
          if(gsv_varExist(varName=vnl_varNameList(jstep))) then
             if (vnl_varKindFromVarname(vnl_varNameList(jstep)) == 'CH') then
                jvar=jvar+1
                if (varName.eq.vnl_varNameList(jstep)) then
                   ier=0
                   exit
                end if
             end if
          end if
       end do
       if (ier.ne.0) call utl_abort('chm_apply_bounds: Field name not recognized -' // trim(varName))

       ! Check if increment size and background are less than the background error std. dev. scaled by refval.
       ! If so, set increment to zero.

       count = 0

       do jstep=step1,step2
          do jlev=lev1,lev2
             do jlat=lat1,lat2
                do jlon=lon1,lon2

                   bkgrnd = background_field(jlon,jlat,jlev,jstep)
                   inc = increment_field(jlon,jlat,jlev,jstep)
                   bg_sigma = bchm_getBgSigma(jlon,jlat,jlev,jvar)

                   if (abs(inc).lt.refval*bg_sigma .and. bkgrnd.lt.refval*bg_sigma) then
                      increment_field(jlon,jlat,jlev,jstep) = 0.0
                      count = count+1
                   end if

                end do
             end do
          end do
       end do

       if (count.gt.0) then
          write(unit,*)
          write(unit,'(A,I7,A,A,A)') "Total of ",count," for small increments set to zero for field ",trim(varName),"."
          write(unit,*)
       end if


    end if

    refval=chm_setup_get_float('low_cutoff',iconstituent_id)
    if (refval.ge.0.) then

       ! Check if both background and analysis are negative. If so set the increment to force an analysis of zero at these locations

       count = 0

       do jstep=step1,step2
          do jlev=lev1,lev2
             do jlat=lat1,lat2
                do jlon=lon1,lon2
                   
                   bkgrnd = background_field(jlon,jlat,jlev,jstep)
                   inc = increment_field(jlon,jlat,jlev,jstep)

                   if (bkgrnd.lt.0.0 .and. bkgrnd+inc.lt.0.0) then

                      if (count.eq.0) then
                         write(unit,'(A,A,A)') "chm_apply_bounds: Negative background and analysis values were found for field ",trim(varName),"."
                         write(unit,'(A)') "Modifying the increment so the analysis is instead equal to zero at these locations."       
                         write(unit,'(A)') "JLON JLEV JLAT TSTEP   Background   Init. incr."
                      end if

                      write(unit,'(4(I4,1X),2G12.2)') jlon,jlev,jlat,jstep,bkgrnd,inc

                      increment_field(jlon,jlat,jlev,jstep) = -bkgrnd
                      count = count+1

                   end if

                end do
             end do
          end do
       end do

       if (count.gt.0) then
          write(unit,*)
          write(unit,'(A,I7,A,A,A)') "Total of ",count," negative background and analysis values were found for field ",trim(varName),"."
          write(unit,*)
       end if

            
       ! Check if the analysis is below the trial field cutoff. If so set analysis to this cut-off value at these locations.

       count = 0

       do jstep=step1,step2
          do jlev=lev1,lev2
             do jlat=lat1,lat2
                do jlon=lon1,lon2

                   bkgrnd = background_field(jlon,jlat,jlev,jstep)
                   inc = increment_field(jlon,jlat,jlev,jstep)

                   if (bkgrnd+inc.lt.refval*bkgrnd .and. bkgrnd.ge.0.) then

                      if (count.eq.0) then
                         write(unit,'(A,F4.1,A,A,A)') "chm_apply_bounds: Analysis values were found below the cut-off of ", &
                              100*chm_setup_get_float('low_cutoff',iconstituent_id),"% of the trial field for field ",trim(varName),"."
                          write(unit,'(A)') "Modifying the increment so the analysis is instead equal to the lower bound cut-off value at these locations."
                          write(unit,'(A)') "JLON JLEV JLAT TSTEP   Background  Init. incr.  New incr."
                      end if

                      new_inc = (refval-1.)*bkgrnd

                      write(unit,'(4(I4,1X),3G12.2)') jlon,jlev,jlat,jstep,bkgrnd,inc,new_inc

                      increment_field(jlon,jlat,jlev,jstep) = new_inc
                      count = count+1

                   end if

                end do
             end do
          end do
       end do

       if (count.gt.0) then
          write(unit,*)
          write(unit,'(A,I7,A,A,A)') "Total of ",count," analysis values found below the cut-off for field ",trim(varName),"."
          write(unit,*)
       end if

    end if

    refval=chm_setup_get_float('high_cutoff',iconstituent_id)
    if (refval.ge.0.) then

       ! Check if the analysis is above the imposed upper bound. If so set analysis to this cut-off value at these locations.

       count = 0

       do jstep=step1,step2
          do jlev=lev1,lev2
             do jlat=lat1,lat2
                do jlon=lon1,lon2

                   bkgrnd = background_field(jlon,jlat,jlev,jstep)
                   inc = increment_field(jlon,jlat,jlev,jstep)

                   if (bkgrnd+inc.gt.refval*bkgrnd .and. bkgrnd.ge.0.) then

                      if (count.eq.0) then
                         write(unit,'(A,A,F4.1,A,A,A)') "chm_apply_bounds: Analysis values were found above the cut-off of ", &
                              100*chm_setup_get_float('high_cutoff',iconstituent_id)," times the trial field for field ",trim(varName),"."
                         write(unit,'(A)') "Modifying the increment so the analysis is instead equal to the upper bound cut-off value at these locations."
                         write(unit,'(A)') "JLON JLEV JLAT TSTEP   Background  Init.-incr.  New-incr."
                      end if

                      new_inc = (refval-1.)*bkgrnd

                      write(unit,'(4(I4,1X),3G12.2)') jlon,jlev,jlat,jstep,bkgrnd,inc,new_inc

                      increment_field(jlon,jlat,jlev,jstep) = new_inc
                      count = count+1

                   end if

                end do
             end do
          end do
       end do

       if (count.gt.0) then
          write(unit,*)
          write(unit,'(A,I7,A,A,A)') "Total of ",count," analysis values found above the cut-off for field ",trim(varName),"."
          write(unit,*)
       end if

    end if

    ier=fclos(unit)

  end subroutine chm_apply_bounds

end module chem_postproc_mod
