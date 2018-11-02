MODULE EarthConstants_mod
   !
   ! MODULE earthConstants_mod (prefix='' category='8')
   !
   ! The following constants should ultimately be taken from module
   ! modgps02wgs84const OR module modgps06gravity.  They have been placed here
   ! as an intermediate step.
   real(8), parameter :: ROMEGA  = 7.29200000000000000000D-05
   real(8), parameter :: RG      = 9.80616000000000000000D+00
   real(8), parameter :: RA      = 6371229.00000000000000D+00
   real(8), parameter :: R1SA    = 1.56955588945241177170D-07
   real(8), parameter :: RV      = 4.61524993308387870310D+02
   real(8), parameter :: rayt    = 0.6371220000000d+07
   real(8), parameter :: grav    = 0.9806160000000d+01

end MODULE EarthConstants_mod
