#-----------------------------------------------------------------------
# This is the makefile for bfield_library_jdl
#	4/18/2011 - 2014 JDL
#
# After updating the compiler options, if necessary, run "make".
# "make clean" will remove all .o and .mod file from the build directory
# "make cleaner" will also remove .a files.
#
# .mod, .o, and .a files are made in current directory
#
#-----------------------------------------------------------------------

# The following flags should work for the pgf90 compiler
#FCOMP = pgf90
#FFLAGS = -fast

# The following flags should work for the gfortran compiler
FCOMP = gfortran
FFLAGS = -O3
FFLAGS_DEBUG = -g -Wall -frange-check -fbounds-check 

USE_M3DC1 = 1  # Set to 1 to compile with m3dc1

SOURCEPATH = ../src/

# Output
LIBFILE = libbjdl.a

# External library paths and libraries for m3dc1
ifeq ($(strip $(USE_M3DC1)),1)
  DEFINES = -DHAVE_M3DC1
  LIBDIR_m3dc1 = $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
  LIBS_m3dc1 = -lfusionio -lm3dc1_fortran -lm3dc1
  INCLUDE_m3dc1 = $(M3DC1_INSTALL_DIR)/include/_$(M3DC1_ARCH)
endif

RANLIB = ranlib

#-----------------------------------------------------------------------
# Dependencies for the library that use standard recipe
OBJS90 = \
bfield_module.o \
bspline90_22.o \
DIIID_routines_mod.o \
g3df_math_routines_mod.o \
g3d_module.o \
kind_mod.o \
NSTX_routines_mod.o \
numerics_module.o \
phys_const.o \
rmp_module.o \
screening_module.o 

# objects that require precompilation
OBJS90_PRE = \
fieldline_follow_mod.o \
util_routines.o
MODS90_PRE = \
fieldline_follow_mod.mod \
util_routines.mod

# Dependencies that have their own recipes
ifeq ($(strip $(USE_M3DC1)),1)
  OBJS90_SPC = M3DC1_routines_mod.o
endif

# Modules to be compiled which have identically named .f90 files
MODS90 = \
g3d_module.mod \
numerics_module.mod \
g3df_math_routines_mod.mod  \
bfield_module.mod \
screening_module.mod \
phys_const.mod \
kind_mod.mod 

#-----------------------------------------------------------------------
.PHONY: debug
.default: $(LIBFILE)

$(LIBFILE): $(OBJS90) $(OBJS90_SPC) $(OBJS90_PRE)
	@echo "Creating "$(LIBFILE)" in directory "$(PWD)
	$(AR) rc $@ $^ #$(OBJS90) $(OBJS90_SPC)
	-$(RANLIB) $@


debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(LIBFILE)

#-----------------------------------------------------------------------

$(MODS90) : %.mod : $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

$(OBJS90) : %.o :  $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

$(OBJS90_PRE) : %.o : $(SOURCEPATH)%.F90
	$(FCOMP) $(FFLAGS) $(DEFINES) -c $<

$(MODS90_PRE) : %.mod : $(SOURCEPATH)%.F90
	$(FCOMP) $(FFLAGS) $(DEFINES) -c $<

M3DC1_routines_mod.o: %.o : $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -L$(LIBDIR_m3dc1) $(LIBS_m3dc1) -I$(INCLUDE_m3dc1) -c $<

#-----------------------------------------------------------------------
bfield_module.o bfield_module.mod:   kind_mod.mod
DIIID_routines_mod.o:   kind_mod.mod phys_const.mod
fieldline_follow_mod.o fieldline_follow_mod.mod: kind_mod.mod g3d_module.mod rmpcoil_module.mod screening_module.mod bfield_module.mod gfile_var_pass.mod phys_const.mod
g3df_math_routines_mod.o g3df_math_routines_mod.mod:   kind_mod.mod
g3d_module.o g3d_module.mod:   kind_mod.mod bspline.mod g3df_math_routines_mod.mod
NSTX_routines_mod.o: kind_mod.mod phys_const.mod
numerics_module.o numerics_module.mod: kind_mod.mod phys_const.mod
phys_const.o phys_const.mod : kind_mod.mod
rmp_module.o  : kind_mod.mod diiid_routines_mod.mod nstx_routines_mod.mod bfield_module.mod
screening_module.o screening_module.mod : kind_mod.mod bspline.mod phys_const.mod math_geo_module.mod
util_routines.o util_routines.mod : kind_mod.mod fieldline_follow_mod.mod gfile_var_pass.mod math_geo_module.mod phys_const.mod g3d_module.mod

ifeq ($(strip $(USE_M3DC1)),1)
  fieldline_follow_mod.o fieldline_follow_mod.mod: m3dc1_routines_mod.mod
  util_routines.o util_routines.mod : m3dc1_routines_mod.mod
  M3DC1_routines_mod.o  : kind_mod.mod $(INCLUDE_m3dc1)/fusion_io.mod
  # These modules are inside of files without the same name
  m3dc1_routines_mod.mod: M3DC1_routines_mod.o
endif


# These modules are inside of files without the same name
rmpcoil_module.mod: rmp_module.o
gfile_var_pass.mod: g3d_module.o
bspline.mod numeric.mod: bspline90_22.o
math_geo_module.mod: numerics_module.o
nstx_routines_mod.mod: NSTX_routines_mod.o
diiid_routines_mod.mod: DIIID_routines_mod.o


#-----------------------------------------------------------------------
.PHONY: clean cleaner
clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."
cleaner:
	-rm *.o *.mod *.a
	@echo "Cleaned object, library and mod files from "$(PWD)"."



