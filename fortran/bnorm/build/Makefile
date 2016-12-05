#-----------------------------------------------------------------------
# This is the makefile for bnorm
#
#  - Intended to be run from the 'build' directory
#
# targets:
#  [default] -  Makes executable
#   debug    -  Makes exec with debug options
#   clean    -  Removes *.mod, *.o from build directory
#   cleaner  -  Same as clean AND removes  executables from bin directory
#   
# .mod, .o, and .a files are made in current directory
#
#	7/14/2011 - current JDL
#-----------------------------------------------------------------------

SRC_DIR = ../src/
BIN_DIR = ../bin/

# Location of bfield library and .mod files
LIBBJDL_DIR = $(HOME)/util-library/fortran/bfield_library_jdl/lib
LIBBJDL_INCDIR = $(HOME)/util-library/fortran/bfield_library_jdl/include

# Location of m3dc1 libraries
LIBM3DC1_DIR = $(FIO_INSTALL_DIR)/lib

# The following flags should work for the gfortran compiler
FC = gfortran
FF = -O3 
FF_DEBUG = -g -Wall -frange-check -fbounds-check

EXEC  = $(BIN_DIR)bnorm.x

#-----------------------------------------------------------------------

# Libraries needed to build the executable
LIBDIR = -L$(LIBBJDL_DIR) -L$(LIBM3DC1_DIR)
LIBS = -lbjdl -lfusionio -lm3dc1 -lz -lhdf5 -lstdc++ -llapack

# Objects needed to build the executable that use the default rule
OBJS90 =  \
bnorm_driver.o \
bnorm_routines.o

# Modules to be compiled which have identically named .f90 files
MODS90 = \
bnorm_bfield_mod.mod

#-----------------------------------------------------------------------
all: $(EXEC)

$(EXEC): $(OBJS90) $(LIBBJDL_DIR)/libbjdl.a 
	@echo "Creating "$@" in directory "$(BIN_DIR)
	${FC} -o $@ $(FF) $^ $(LIBDIR) $(LIBS)

.PHONY: debug
debug: FF=$(FF_DEBUG)
debug: all

#-----------------------------------------------------------------------
# Static patterns

$(OBJS90) : %.o:  $(SRC_DIR)%.f90
	$(FC) $(FF) -c -I$(LIBBJDL_INCDIR) $(LIBDIR) $(LIBS) $<

#!$(MODS90) : %.mod:  $(SRC_DIR)%.f90
#	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

# dependencies
bnorm_driver.o : bnorm_bfield_mod.mod
bnorm_bfield_mod.mod : bnorm_routines.o

.PHONY: clean cleaner
clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)

cleaner:
	-rm *.o *.mod *.a $(EXEC)
	@echo "Cleaned object, mod, library and executable files from "$(PWD)
