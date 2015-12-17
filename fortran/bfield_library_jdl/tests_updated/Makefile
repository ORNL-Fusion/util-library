#-----------------------------------------------------------------------
# This is the makefile for bfield_library_jdl tests
#	4/18/2011 - 2015 JDL
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
FC = gfortran
FF = -O3
FF_DEBUG = -g -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace

USE_M3DC1 = 1  # Set to 1 to compile with m3dc1

SOURCEPATH = ./

# Location of bfield library and .mod files
LIBBJDL_DIR = $(HOME)/util-library/fortran/bfield_library_jdl/build

# Libraries needed to build the executable
LIBDIR = -L$(LIBBJDL_DIR) 
LIBS = -lbjdl -llapack

# External library paths and libraries for m3dc1
ifeq ($(strip $(USE_M3DC1)),1)
  DEFINES = -DHAVE_M3DC1
  LIBDIR_m3dc1 = $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
  LIBS_m3dc1 = -lfusionio -lm3dc1_fortran -lm3dc1 -lstdc++ -lhdf5
#  LIBS_m3dc1_full = libfusionio.a libm3dc1_fortran.a libm3dc1.a
  INCLUDE_m3dc1 = $(M3DC1_INSTALL_DIR)/include/_$(M3DC1_ARCH)
endif


#-----------------------------------------------------------------------
# Dependencies for the library that use standard recipe
OBJS90 = 

# objects that require precompilation
OBJS90_PRE = 
MODS90_PRE = 

# Modules to be compiled which have identically named .f90 files
MODS90 = 
#-----------------------------------------------------------------------
.PHONY: debug all_tests
.default: all_tests

all_tests: test_int_two_lines.x test_int_curve_curve.x test_find_xpt_gonly.x


test_int_two_lines.x: test_int_two_lines.o $(LIBBJDL_DIR)/libbjdl.a 
	@echo "Creating "$@
	${FC} -o $@ $(FF) $^ $(LIBDIR) $(LIBS)

test_int_two_lines.o : test_int_two_lines.f90 
	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

test_int_curve_curve.x: test_int_curve_curve.o $(LIBBJDL_DIR)/libbjdl.a 
	@echo "Creating "$@
	${FC} -o $@ $(FF) $^ $(LIBDIR) $(LIBS)

test_int_curve_curve.o : test_int_curve_curve.f90 
	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

test_find_xpt_gonly.x: %.x :%.o $(LIBBJDL_DIR)/libbjdl.a 
	@echo "Creating "$@
	${FC} -o $@ $(FF) $^ $(LIBDIR) $(LIBS) -L$(LIBDIR_m3dc1) $(LIBS_m3dc1) -I$(INCLUDE_m3dc1)

test_find_xpt_gonly.o : %.o : %.f90 
	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

debug: FF = $(FF_DEBUG)
debug: all_tests

#-----------------------------------------------------------------------
.PHONY: clean cleaner
clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."
cleaner:
	-rm *.o *.mod *.x
	@echo "Cleaned object, library and mod files from "$(PWD)"."


