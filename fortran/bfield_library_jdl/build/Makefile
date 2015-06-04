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

SOURCEPATH = ../src/

# Output
LIBFILE = libbjdl.a

# External library paths and libraries for m3dc1
LIBDIR_m = -L$(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
LIBS_m = -lfusionio -lm3dc1_fortran -lm3dc1
INCLUDE_m = -I$(M3DC1_INSTALL_DIR)/include/_$(M3DC1_ARCH)

RANLIB = ranlib

#-----------------------------------------------------------------------
# Dependencies for the libray that use default rule
OBJS90 = g3d_module.o numerics_module.o fieldline_follow_mod.o\
g3df_math_routines_mod.o bspline90_22.o NSTX_routines_mod.o bfield_module.o \
rmp_module.o screening_module.o DIIID_routines_mod.o util_routines.o

# Dependencies that have their own rules
OBJS90_SPC = M3DC1_routines_mod.o

# Modules to be compiled which have identically named .f90 files
MODS90 = g3d_module.mod numerics_module.mod fieldline_follow_mod.mod\
g3df_math_routines_mod.mod  bfield_module.mod util_routines.mod \
screening_module.mod 

#-----------------------------------------------------------------------
$(LIBFILE): $(OBJS90) $(OBJS90_SPC)
	@echo "Creating "$(LIBFILE)" in directory "$(PWD)
	$(AR) rc $@ $(OBJS90) $(OBJS90_SPC)
	-$(RANLIB) $@

.PHONY: debug
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(LIBFILE)


#-----------------------------------------------------------------------

# The following are static pattern rules that provide
# instructions for creating object files from the fortran
# source files.  The older suffix rule would have been ".f.o:"

# The $< at the end is a macro for the current source file.

$(MODS90) : %.mod :  $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

$(OBJS90) : %.o:  $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

M3DC1_routines_mod.o: %.o : $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) $(LIBDIR_m) $(LIBS_m) $(INCLUDE_m) -c $<

#-----------------------------------------------------------------------
# A list of module dependencies ensures that module information
# specific to a particular source is available.
g3d_module.o: kind_mod.mod g3df_math_routines_mod.mod bspline.mod 
fieldline_follow_mod.o: rmpcoil_module.mod bfield_module.mod kind_mod.mod screening_module.mod m3dc1_routines_mod.mod
rmpcoil_module.mod: diiid_routines_mod.mod nstx_routines_mod.mod bfield_module.mod
m3dc1_routines_mod.mod: M3DC1_routines_mod.o
kind_mod.mod: numerics_module.o
bspline.mod: bspline90_22.o
rmpcoil_module.mod: rmp_module.o
nstx_routines_mod.mod: NSTX_routines_mod.o
diiid_routines_mod.mod: DIIID_routines_mod.o
util_routines.mod: kind_mod.mod fieldline_follow_mod.mod
#-----------------------------------------------------------------------
# The dash in front of a command means to continue even if that
# command returns a nonzero error code.
.PHONY: clean cleaner
clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."
cleaner:
	-rm *.o *.mod *.a
	@echo "Cleaned object, library and mod files from "$(PWD)"."



