
module randomNumber_mod
  ! MODULE randomNumber_mod (prefix='rng' category='8. Low-level utilities and constants')
  !
  !:Purpose: A Gaussian random number generator (RNG) module. The actual calculation
  !          is performed by functions in an external library.
  !
  use ISO_C_BINDING
  implicit none

  include 'randomfunctions.inc'

  save
  private

  ! Public procedures
  public :: rng_setup, rng_gaussian, rng_uniform

  logical :: initialized = .false.

  type(RANDOM_STREAM) :: randomStream

contains

  !--------------------------------------------------------------------------
  ! rng_Setup
  !--------------------------------------------------------------------------
  subroutine rng_setup(seed, beSilent_opt)
    !
    !:Purpose: Initialize the random number generator with a supplied seed.
    !
    implicit none

    ! Arguments:
    integer,           intent(in) :: seed
    logical, optional, intent(in) :: beSilent_opt

    ! Locals:
    integer, dimension(1) :: seeds
    type(RANDOM_STREAM)   :: null_stream
    logical               :: beSilent

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if (initialized .and. .not.beSilent) then
      write(*,*) 'rng_setup: WARNING: you are re-initializing the module!!!'
    end if

    seeds(1) = seed

    null_stream = RANDOM_STREAM(C_NULL_PTR)

    ! 'seeds' is an array of dimension 1
    call Ran_R250_new_stream(randomStream, null_stream, seeds, size(seeds))
    call RanSetSeed_gaussian_stream(randomStream, seeds, size(seeds))

    initialized = .true.

    write(*,*) 'rng_setup: done using seed = ', seed
 
  end subroutine rng_setup
  
  !--------------------------------------------------------------------------
  ! rng_Gaussian
  !--------------------------------------------------------------------------
  function rng_gaussian() result(randomNumberGaussian)
    !
    !:Purpose: Returns a normally distributed deviate
    !          with zero mean and unit variance
    !
    implicit none

    ! Result:
    real(8) :: randomNumberGaussian

    randomNumberGaussian = DRan_gaussian_stream(randomStream)
  end function rng_gaussian
  
  !--------------------------------------------------------------------------
  ! random
  !--------------------------------------------------------------------------
  function rng_uniform() result(randomNumberUniform)
    !
    !:Purpose: Returns a random deviate between 0.0 and 1.0.
    !
    implicit none

    ! Result:
    real(8) :: randomNumberUniform
    
    randomNumberUniform = DRan_generic_stream(randomStream)
    
  end function rng_uniform
  
end module randomNumber_mod
