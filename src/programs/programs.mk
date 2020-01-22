
#===| DOUBLE PRECISION |=========================
PGM_DBL := var thinning seaIce randomPert obsSelection obsImpact oMinusF\
			letkf genCoeff ensembleH ensManip diagHBHt diagBmatrix\
			calcStats bgckMW advector adjointTest addIncrement\
			ensPostProcess 
# TODO: automate that ^

ABS_DBL := $(addsuffix .Abs,$(PGM_DBL)) 

#--------------------------------------
var.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
		rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
thinning.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random 

#--------------------------------------
seaIce.Abs: LIBAPPL = netcdff burp_module $(VGRID_LIBNAME) irc\
		$(MPILIB) random

#--------------------------------------
randomPert.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
obsSelection.Abs: LIBAPPL = rttov_coef_io rttov_hdf\
		rttov_parallel  rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB)\
		f90sqlite udfsqlite random

#--------------------------------------
obsImpact.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
oMinusF.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
letkf.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
genCoeff.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) $(MPILIB) random

#--------------------------------------
ensembleH.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
 		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
ensPostProcess.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
 		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
ensManip.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
	  	rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
 		$(VGRID_LIBNAME) irc $(MPILIB)

#--------------------------------------
diagHBHt.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
diagBmatrix.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
calcStats.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB)

#--------------------------------------
bgckMW.Abs: LIBAPPL = burp_module irc $(MPILIB)

#--------------------------------------
advector.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB)

#--------------------------------------
adjointTest.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
addIncrement.Abs: LIBAPPL = rttov_coef_io rttov_hdf\
		rttov_parallel  rttov_main rttov_emis_atlas rttov_other\
		$(HDF5_LIBS) f90sqlite udfsqlite burp_module $(VGRID_LIBNAME) irc\
		$(MPILIB) random

#--------------------------------------
letkf.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel rttov_main\
		rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module $(VGRID_LIBNAME)\
		irc $(MPILIB) f90sqlite udfsqlite random

#===| SINGLE PRECISION |=========================
PGM_SMPL := prepcma obsIO
ABS_SMPL := $(addsuffix .Abs,$(PGM_SMPL))

#--------------------------------------
prepcma.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) f90sqlite\
		udfsqlite burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
obsIO.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite  

# vim: set noexpandtab noautoindent:
