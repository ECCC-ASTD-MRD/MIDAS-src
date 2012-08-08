! This file centralizes compile-time (pre-processing) switches.  In order for
! these switches to take effect in a particular source-code file, the present
! file must be included therein, thus:
!              #include "compileswitches.inc"
!
! GENERAL SWITCHES:
! If one of the general switches is defined, it specifies for which system the 
! files will be compiled.  All other switches are automatically set to match.  If
! none of these switches is set, then the last if-branch must be used to set the
! subordinate switches individually.  More than one general switch should not be
! defined simultaneously.  It is intended that one general switch be #define'd
! just before #include'ing this file.  These are the choices:

!#define      ENKF   - Ensemble-Kalman-Filter assimilation
!#define      VAR3D  - 3-dimensional variational assimilation


! SUBORDINATE SWITCHES
#ifdef ENKF
    #define OBS_REAL 4

#elif VAR3D
    #define OBS_REAL 8

#else
    ! General switch has not been set.  Set the switches individually here.
    #define OBS_REAL 2
#endif 
