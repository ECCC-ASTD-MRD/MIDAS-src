
module codePrecision_mod
  ! MODULE codePrecision_mod (prefix='pre' category='8. Low-level utilities and constants')
  !
  !:Purpose: A module to specify the precision, mostly for floating
  !          point variables
  !
  use mpi_f08, only: MPI_Datatype, MPI_Type_get_name, MPI_MAX_OBJECT_NAME ! this is the Fortran 2008 MPI library module
  use midasMpi_mod
  implicit none
  save
  public

  !
  ! Precision for columns in obsSpaceData - default is real8
  !
#if !defined(CODEPRECISION_OBS_REAL_SINGLE) && !defined(CODEPRECISION_OBS_REAL_DOUBLE)
#define CODEPRECISION_OBS_REAL_DOUBLE
#endif

#ifdef CODEPRECISION_OBS_REAL_DOUBLE
  integer,                   parameter :: pre_obsReal       = selected_real_kind(15)
  type(MPI_Datatype),        parameter :: pre_obsMpiReal    = mmpi_real8
#endif

#ifdef CODEPRECISION_OBS_REAL_SINGLE
  integer,                   parameter :: pre_obsReal       = selected_real_kind(6)
  type(MPI_Datatype),        parameter :: pre_obsMpiReal    = mmpi_real4
#endif

  !
  ! Precision for calculation of analysis increment in variational analysis - default is real8
  !
#if !defined(CODEPRECISION_INCR_REAL_SINGLE) && !defined(CODEPRECISION_INCR_REAL_DOUBLE)
#define CODEPRECISION_INCR_REAL_DOUBLE
#endif

#ifdef CODEPRECISION_INCR_REAL_DOUBLE
  integer, parameter :: pre_incrReal = selected_real_kind(15)
#endif

#ifdef CODEPRECISION_INCR_REAL_SINGLE
  integer, parameter :: pre_incrReal = selected_real_kind(6)
#endif

  !
  ! Precision for mpi transposes in (global) spectral transform - default is real8
  !
#if !defined(CODEPRECISION_SPECTRANS_REAL_SINGLE) && !defined(CODEPRECISION_SPECTRANS_REAL_DOUBLE)
#define CODEPRECISION_SPECTRANS_REAL_DOUBLE
#endif

#ifdef CODEPRECISION_SPECTRANS_REAL_DOUBLE
  integer,                   parameter :: pre_specTransReal       = selected_real_kind(15)
  type(MPI_Datatype),        parameter :: pre_specTransMpiType    = mmpi_real8
#endif

#ifdef CODEPRECISION_SPECTRANS_REAL_SINGLE
  integer,                   parameter :: pre_specTransReal       = selected_real_kind(6)
  type(MPI_Datatype),        parameter :: pre_specTransMpiType    = mmpi_real4
#endif

contains

  subroutine pre_printPrecisions
    !
    !:Purpose: To print precision parameters in the listing.
    !
    implicit none

    ! Locals:
    integer :: nameLength
    character(len=MPI_MAX_OBJECT_NAME) :: mpiTypeName

    write(*,*)
    write(*,*) " <<<<<< Code precision parameters >>>>>>"
    write(*,"(A36, I2)")  "         pre_obsReal= ",    pre_obsReal
    call MPI_Type_get_name(pre_obsMpiReal,       mpiTypeName, nameLength)
    write(*,"(A36, A20)") "      pre_obsMpiReal= ",    mpiTypeName
    write(*,"(A36, I2)")  "        pre_incrReal= ",    pre_incrReal
    write(*,"(A36, I2)")  "   pre_specTransReal= ",    pre_specTransReal
    call MPI_Type_get_name(pre_specTransMpiType, mpiTypeName, nameLength)
    write(*,*)  "             pre_specTransMpiType= ", mpiTypeName
    write(*,*)

  end subroutine pre_printPrecisions

end module codePrecision_mod
