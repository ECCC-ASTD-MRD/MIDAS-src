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
!-------------------------------------- LICENCE end --------------------------------------

module oceanBackground_mod
  ! MODULE oceanBackground_mod (prefix='obgd' category='1. High-level functionality')
  !
  ! :Purpose: storage for ocean background related subroutines
  !
  use midasMpi_mod
  use utilities_mod
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use codePrecision_mod  
  
  implicit none

  save
  private

  ! Public functions/subroutines
  public :: obgd_computeSSTrial
 
  ! External functions
  integer, external :: fnom, fclos  

  contains

  !----------------------------------------------------------------------------------------
  ! obgd_computeSSTrial
  !----------------------------------------------------------------------------------------
  subroutine obgd_computeSSTrial(hco, vco, trialDateStamp, analysisDateStamp, &
                                 nmonthsClim, datestampClim, alphaClim, etiket)
    !
    !: Purpose: to compute SST background   
    !           xb(t) = (xa(t-1) - xclim(t-1))*alpha + xclim(t)         
    
    implicit none

    ! Arguments:
    type(struct_hco) , intent(in), pointer :: hco               ! horizontal grid
    type(struct_vco) , intent(in), pointer :: vco               ! vertical grid
    integer          , intent(in)          :: trialDateStamp    ! trial datestamp
    integer          , intent(in)          :: analysisDateStamp ! datestamp for last analysis 
    integer          , intent(in)          :: nmonthsClim       ! number of climatological fields (= 12)
    integer          , intent(in)          :: datestampClim(:)  ! datestamps of input climatology fields
    real(4)          , intent(in)          :: alphaClim         ! scalling factor to relax towards climatology
    character(len=10), intent(in)          :: etiket            ! etiket from namelist and for trial
    
    ! Locals:
    type(struct_gsv) :: stateVector
    real(4), pointer :: stateVector_ptr(:, :, :)
    integer          :: lonIndex, latIndex
    real(8)          :: climatology_m1(hco % ni, hco % nj), climatology(hco % ni, hco % nj)
    
    write(*,*) 'obgd_computeSSTrial: starting...'  
      
    ! get SST analysis
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 4, datestamp_opt = analysisDateStamp, &
                      mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVector, './analysis', ' ','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector, stateVector_ptr)
    
    call obgd_getClimatology(analysisDateStamp, hco, vco, nmonthsClim, datestampClim, climatology_m1)
    call obgd_getClimatology(trialDateStamp   , hco, vco, nmonthsClim, datestampClim, climatology)
    
    ! compute background state
    ! xb(t) = (xa(t-1) - xclim(t-1))*alpha + xclim(t)
    do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd 
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        stateVector_ptr(lonIndex, latIndex, 1) = (stateVector_ptr(lonIndex, latIndex, 1) - &
                                                 climatology_m1(lonIndex, latIndex)) * &
                                                 alphaClim + climatology(lonIndex, latIndex)
      end do
    end do      	 
  
    ! modify dateStamp (from analysis) with trial dateStamp
    call gsv_modifyDate( stateVector, trialDateStamp, modifyDateOrigin_opt = .true. )
    
    ! save trial field
    call gio_writeToFile(stateVector, './trial', etiket, typvar_opt = 'P@')

    call gsv_deallocate(stateVector)

  end subroutine obgd_computeSSTrial
  
  !----------------------------------------------------------------------------------------
  ! obgd_getClimatology
  !----------------------------------------------------------------------------------------
  subroutine obgd_getClimatology(dateStamp, hco, vco, nmonthsClim, datestampClim, output)
    !
    !: Purpose: 1) to read SST climatological fields from a std file
    !           2) to interpolate the field in time using the current day (t) in current month (m)    
    !           SST(t) = SST_clim(m) + (t-1)/(ndays-1) * (SST_clim(m+1) - SST_clim(m)),
    !           where ndays is a number of days in current month        
    
    implicit none

    ! Arguments:
    integer         , intent(in)          :: dateStamp        ! date stamp for the current day
    type(struct_hco), intent(in), pointer :: hco              ! horizontal grid
    type(struct_vco), intent(in), pointer :: vco              ! vertical grid
    integer         , intent(in)          :: nmonthsClim      ! number of records in the climatology file 
    integer         , intent(in)          :: datestampClim(:) ! datestamps in the climatology file
    real(8)         , intent(inout)       :: output(:,:)      ! interpolated SST field from climatology
  
    ! locals
    integer          :: hour, day, month, yyyy, ndays, nextMonth
    type(struct_gsv) :: stateVector, stateVector_nextMonth
    real(4), pointer :: clim_ptr(:, :, :), clim_nextMonth_ptr(:, :, :)
    integer          :: lonIndex, latIndex
   
    call tim_dateStampToYYYYMMDDHH(dateStamp, hour, day, month, ndays, yyyy)
    write(*,'(a,3i5,a,i12,a)') 'obgd_getClimatology: interpolating climatology for day/month/year (datestamp): ', &
    day, month, yyyy, '(', datestamp, ')'

    if (month == nmonthsClim) then
      nextMonth = 1
    else
      nextMonth = month + 1
    end if 
     
    ! get climatology, current month
    write(*,*) 'obgd_getClimatology: reading climatology, month: ', month, ', datestamp: ', datestampClim(month) 
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = datestampClim(month), mpi_local_opt = .true., &
                      varNames_opt = (/'TM'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVector, './climatology', ' ',' ', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector, clim_ptr)
    
    ! get climatology, next month
    write(*,*) 'obgd_getClimatology: reading climatology, month: ', nextMonth, ', datestamp: ', datestampClim(nextMonth) 
    call gsv_allocate(stateVector_nextMonth, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = datestampClim(nextMonth), mpi_local_opt = .true., &
                      varNames_opt = (/'TM'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVector_nextMonth, './climatology', ' ',' ', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_nextMonth, clim_nextMonth_ptr)

    do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd 
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        output(lonIndex, latIndex) = clim_ptr(lonIndex, latIndex, 1) + real(day - 1) / real(ndays - 1) * &
                                     (clim_nextMonth_ptr(lonIndex, latIndex, 1) - clim_ptr(lonIndex, latIndex, 1))	 
      end do
    end do
    
    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVector_nextMonth)
  
  end subroutine obgd_getClimatology 
  
end module oceanBackground_mod  
