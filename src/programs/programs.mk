#--------------------------------------
adjointTest.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
advector.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
analysisErrorOI.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
calcStats.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
diagBmatrix.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
diagHBHt.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
ensembleH.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
ensDiagnostics.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
ensManip.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
                rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module\
                $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
ensPostProcess.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
                rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module\
                $(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
extractBmatrixFor1Dvar.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
								rttov_main rttov_emis_atlas rttov_other $(HDF5_LIBS) burp_module\
	 							$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
genCoeff.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) $(MPILIB) random

#--------------------------------------
letkf.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module $(VGRID_LIBNAME)\
		irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
obsImpact.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
obsSelection.Abs: LIBAPPL = rttov_coef_io rttov_hdf\
		rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB)\
		f90sqlite udfsqlite random

#--------------------------------------
oMinusF.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) burp_module\
		$(VGRID_LIBNAME) irc $(MPILIB) f90sqlite udfsqlite random

#--------------------------------------
prepcma.Abs: LIBAPPL = rttov_coef_io rttov_hdf rttov_parallel\
		rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main $(HDF5_LIBS) f90sqlite\
		udfsqlite burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
pseudoSSTobs.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
randomPert.Abs: LIBAPPL = $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
SSTbias.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
SSTtrial.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
thinning.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io\
		rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
		$(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random 

#--------------------------------------
var.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
var1D.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

#--------------------------------------
pseudoSSTobs.Abs: LIBAPPL = f90sqlite udfsqlite rttov_coef_io rttov_hdf\
                rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main\
                $(HDF5_LIBS) burp_module $(VGRID_LIBNAME) irc $(MPILIB) random

# vim: set noexpandtab noautoindent nolist:
