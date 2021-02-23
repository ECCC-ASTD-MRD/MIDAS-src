MODULE EarthConstants_mod
   ! MODULE EarthConstants_mod (prefix='' category='8. Global constants and interfaces')
   !
   ! The following constants should ultimately be taken from module
   ! modgps02wgs84const OR module modgps06gravity.  They have been placed here
   ! as an intermediate step.
   ! double precision
   integer, parameter  :: dp=kind(0.d0)  
   real(8), parameter  :: ROMEGA   = 7.29200000000000000000D-05
   real(8), parameter  :: RG       = 9.80616000000000000000D+00
   real(8), parameter  :: RA       = 6371229.00000000000000D+00
   ! Radius of sphere of equal area (WSG_R2 is the same value from gps_mod.f90)
   real(dp), parameter :: earth_r2 = 6371007.1809_dp    
   real(8), parameter  :: R1SA     = 1.56955588945241177170D-07
   real(8), parameter  :: RV       = 4.61524993308387870310D+02
   real(8), parameter  :: rayt     = 0.6371220000000d+07
   real(8), parameter  :: grav     = 0.9806160000000d+01

end MODULE EarthConstants_mod
