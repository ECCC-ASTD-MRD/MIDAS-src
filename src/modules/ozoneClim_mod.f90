
module ozoneClim_mod
  ! MODULE ozoneClim_mod (prefix='ozo' category='5. Observation operators')
  !
  ! :Purpose: Climatological ozone (1998)
  !
  use obsSpaceData_mod
  use presProfileOperators_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures  
  public :: OZO_GET_PROFILE, OZO_READ_CLIMATOLOGY


  ! Number of latitudes and vertical levels in climatology file
  INTEGER, PARAMETER    :: NLATO3=19, NLEVO3=28

  ! Climatological ozone field (ppmv) and total ozone 
  REAL                  :: FOZO_r4(NLATO3,NLEVO3), TOTOZO_r4(NLATO3,12)

  ! Pressure height of climatology file vertical levels (mb)
  REAL(8)               :: PO3(NLEVO3)

  DATA PO3 /    0.010D0, 0.015D0, 0.022D0, 0.032D0, 0.046D0, 0.068D0, 0.100D0,   &
       0.150D0, 0.200D0, 0.300D0, 0.500D0, 1.000D0, 2.000D0, 3.000D0,   &
       5.000D0, 7.000D0, 10.00D0, 20.00D0, 30.00D0, 50.00D0, 70.00D0,   &
       100.0D0, 150.0D0, 200.0D0, 300.0D0, 500.0D0, 700.0D0, 1000.D0 / 

contains

  !--------------------------------------------------------------------------
  ! ozo_get_profile
  !--------------------------------------------------------------------------
  subroutine ozo_get_profile(o3p,zlat,plev,nlev,nprf)
    !
    !:Purpose: Get ozone profile from climatology interpolated to desired P levels
    !
    implicit none

    ! Arguments:
    integer ,intent(in) :: nlev            ! NUMBER OF VERTICAL LEVELS
    integer ,intent(in) :: nprf            ! NUMBER OF PROFILES
    REAL(8),intent(in)  :: ZLAT(NPRF)      ! ARRAY OF LATITUDE (-90S TO 90N)
    REAL(8),intent(in)  :: PLEV(NLEV,NPRF) ! PRESSURE LEVELS (HPA)
    REAL(8),intent(out) :: O3P(NLEV,NPRF)  ! OZONE PROFILES (PPMV)

    ! Locals:
    INTEGER   :: JN, K, NUMLAT
    REAL(8)   :: QO3B(NLEVO3,NPRF)
    REAL(8)   :: PRO3(NLEVO3,NPRF)

    !* assign default qgas values if need be

    DO JN = 1, NPRF
       NUMLAT = NINT( (ZLAT(JN)+90.D0) / (180.D0/(REAL(NLATO3-1,8))) ) + 1
       DO K = 1, NLEVO3
          QO3B(K,JN) = FOZO_r4(NUMLAT,K)
       END DO
    END DO

    !* interpolation of field QO3B at NLEVO3 levels of height PO3mbb
    !* into field O3P at NLEV levels of height PLEV

    FORALL(K=1:NLEVO3) PRO3(K,:) = PO3(K)

    CALL ppo_LINTV(pro3,qo3b,nlevo3,nprf,nlev,plev,O3P)

  end subroutine ozo_get_profile

  !--------------------------------------------------------------------------
  ! ozo_read_climatology
  !--------------------------------------------------------------------------
  subroutine ozo_read_climatology(datestamp,nlat_opt,nlev_opt,press_opt,ozone_opt)
    !
    !:Purpose: READ OZONE CLIMATOLOGICAL FIELDS
    !
    implicit none
    
    ! Arguments:
    integer            :: datestamp            ! Datestamp
    integer, intent(out), optional :: nlat_opt ! Number of latitudes
    integer, intent(out), optional :: nlev_opt ! Number of vertical levels
    real(8), allocatable, intent(out), optional :: ozone_opt(:,:) ! Ozone field
    real(8), allocatable, intent(out), optional :: press_opt(:)   ! Pressure levels

    ! Locals:
    INTEGER            :: IJOUR,ITIME,IMONTH,IJ,IER
    CHARACTER(len=100) :: CFILE
    INTEGER            :: NIOZO,NJOZO,NKOZO
    INTEGER, EXTERNAL  :: FNOM,FSTOUV,FSTLIR,FSTFRM,FCLOS,NEWDATE
    integer            :: IOZTEST
    integer            :: iv1,iv2,iv3,iv4,iv5,iv6

    ier = newdate(datestamp,ijour,itime,-3)

    IJ= IJOUR/100
    IMONTH = IJ - (IJ/100)*100

    ioztest=0

    CFILE='ozoneclim98'
    IV1=FNOM(IOZTEST,CFILE,'RND+R/O',0)
    IV2=FSTOUV(IOZTEST,'RND')
    IV3=FSTLIR(FOZO_r4,IOZTEST,NIOZO,NJOZO,NKOZO,-1,' ',-1,-1,IMONTH,' ','O3')
    IV4=FSTLIR(TOTOZO_r4,IOZTEST,NIOZO,NJOZO,NKOZO,-1,' ',-1,-1,-1,' ','TO')
    IV5=FSTFRM(IOZTEST)
    IV6=FCLOS(IOZTEST)

    if(iv1.lt.0.or.iv2.lt.0.or.iv3.lt.0.or.iv4.lt.0.or.iv5.lt.0.or.iv6.lt.0) then
       write(*,*) 'LES IV DE OZO_READ_CLIMATOLOGY ',iv1,iv2,iv3,iv4,iv5,iv6
       write(*,*) 'THESE NUMBERS SHOULD NOT BE NEGATIVE'
       write(*,*) 'datestamp,ijour,itime,imonth = ',datestamp,ijour,itime,imonth
       call utl_abort('Problem with file in ozo_read_climatology (ozoneclim_mod)')
    endif

    if (present(nlat_opt)) then
       nlat_opt=NLATO3
       nlev_opt=NLEVO3
       allocate(press_opt(nlevo3),ozone_opt(nlato3,nlevo3))
       press_opt(1:nlevo3)=PO3(1:nlevo3)
       ozone_opt(1:nlato3,1:nlevo3)=FOZO_r4(1:nlato3,1:nlevo3)
    endif
    
  end subroutine OZO_READ_CLIMATOLOGY

end module ozoneClim_mod
