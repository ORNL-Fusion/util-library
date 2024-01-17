cmake_minimum_required(VERSION 2.8)
project ( bjdl_lib )
enable_language (Fortran)

#
# To compile without M3D-C1 set variable USE_M3DC1 to 0 below
#   -- M3D-C1 library variable FIO_INSTALL_DIR assumed to be set up properly
#
# Similar logic for XDR routines
#
# This is now automated by 
# 1) cd build
# 2) ./cmake_setup.sh [debug]
#
# Debug and release flags are set below.
#
# 1) cd build
# 2) cmake -DCMAKE_BUILD_TYPE=Release .. 
# 3) make 
#
# To see the commands during compilation use "make VERBOSE=1"
#
# To compile with debug options:
#
# 2) cmake -DCMAKE_BUILD_TYPE=Debug ..
# 


# M3D-C1 Options
set (USE_M3DC1 0)
set ( M3DC1_INCLUDE_PATH $ENV{FIO_INSTALL_DIR}/include)

# XDR options
set (USE_FXDR 0)

# Compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
  # gfortran
  Message ("Using gfortran flags")
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG "-g -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -Og -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace")  
  #-ffpe-trap=zero,underflow,overflow,invalid)
else ()
  Message ("Did not recognize compiler -- no flags set!")
endif ()

# Set up directories and options
message ("CMAKE_CURRENT_SOURCE_DIR:" ${CMAKE_CURRENT_SOURCE_DIR})
set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set (INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set (CMAKE_Fortran_MODULE_DIRECTORY ${INC_DIR})
include_directories(${INC_DIR})
file( MAKE_DIRECTORY ${INC_DIR} )

# Handle M3D-C1 settings
if (USE_M3DC1 MATCHES 1)
  set (EXTRA_FILES ${SRC_DIR}/M3DC1_routines_mod.f90)
  message ("Using M3DC1 files from FIO_INSTALL_DIR:" $ENV{FIO_INSTALL_DIR})
  add_definitions (-DHAVE_M3DC1)
  include_directories (${M3DC1_INCLUDE_PATH})  
endif ()

# Handle FXDR settings
if (USE_FXDR MATCHES 1)
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrinit.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrint.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrclose.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrreal.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrdouble.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrmat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrdmat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrimat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrstring.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrreal64.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrewind.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrmat64.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrshort.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/fxdr_module.f90")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/xdr_routines_mod.f90")
  message ("Using FXDR")
  add_definitions (-DHAVE_FXDR)
endif ()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR})

add_library( bjdl STATIC
  ${SRC_DIR}/bfield_module.f90
  ${SRC_DIR}/bspline90_22.f90
  ${SRC_DIR}/DIIID_routines_mod.f90
  ${SRC_DIR}/fieldline_follow_mod.F90
  ${SRC_DIR}/g3d_module.f90
  ${SRC_DIR}/kind_mod.f90
  ${SRC_DIR}/NSTX_routines_mod.f90
  ${SRC_DIR}/numerics_module.f90
  ${SRC_DIR}/phys_const.f90
  ${SRC_DIR}/rmp_module.f90
  ${SRC_DIR}/screening_module.f90
  ${SRC_DIR}/util_routines.F90
  ${SRC_DIR}/IPEC_routines_mod.f90
  ${SRC_DIR}/XPAND_routines_mod.f90
  ${SRC_DIR}/biotsavart_module.f90
  ${SRC_DIR}/coil_typedef.f90
  ${SRC_DIR}/bnorm_routines.f90
  ${SRC_DIR}/setup_bfield.f90
  ${SRC_DIR}/VMEC_routines_mod.f90
  ${SRC_DIR}/bgrid_routines_mod.f90
  ${SRC_DIR}/LIBSTELL_LIGHT/tolower.f90
  ${SRC_DIR}/LIBSTELL_LIGHT/stel_kinds.f
  ${SRC_DIR}/LIBSTELL_LIGHT/safe_open_mod.f90
  ${SRC_DIR}/LIBSTELL_LIGHT/biotsavart.f
  ${SRC_DIR}/LIBSTELL_LIGHT/bsc_T.f
  ${EXTRA_FILES}
  )




# Installation
#install( DIRECTORY DESTINATION ${LIB_DIR} )
#install( TARGETS bjdl
#  LIBRARY DESTINATION ${LIB_DIR}
#  ARCHIVE DESTINATION ${LIB_DIR})
