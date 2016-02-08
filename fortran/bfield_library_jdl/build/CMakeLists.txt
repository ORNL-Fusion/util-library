cmake_minimum_required(VERSION 3.3)
project ( bjdl_lib )
enable_language (Fortran)

set (USE_M3DC1 0)   # set to zero to turn off compilation with M3D-C1
#option (USE_M3DC1 "Build libary against M3D-C1?" OFF)
set ( M3DC1_LIB_PATH $ENV{FIO_INSTALL_DIR}/lib)
set ( M3DC1_INCLUDE_PATH $ENV{FIO_INSTALL_DIR}/include)


if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif ()

if (USE_M3DC1 MATCHES 1)
  set (EXTRA_FILES ../src/M3DC1_routines_mod.f90)
  message ("Using M3DC1 files from FIO_INSTALL_DIR:" $ENV{FIO_INSTALL_DIR})
  find_library( M3DC1_F NAMES fusionio PATHS ${M3DC1_LIB_PATH} NO_DEFAULT_PATH)
  find_library( M3DC1_F NAMES m3dc1 PATHS ${M3DC1_LIB_PATH} NO_DEFAULT_PATH)
  add_definitions (-DHAVE_M3DC1)
  include_directories (${M3DC1_INCLUDE_PATH})  
endif ()

set_source_files_properties(
  ../src/fieldline_follow_mod.F90
  ../src/util_routines.F90
  )
add_library( bjdl STATIC
  ../src/bfield_module.f90
  ../src/bspline90_22.f90
  ../src/DIIID_routines_mod.f90
  ../src/fieldline_follow_mod.F90
  ../src/g3d_module.f90
  ../src/kind_mod.f90
  ../src/NSTX_routines_mod.f90
  ../src/numerics_module.f90
  ../src/phys_const.f90
  ../src/rmp_module.f90
  ../src/screening_module.f90
  ../src/util_routines.F90
  ${EXTRA_FILES}
  )

#file (TO_NATIVE_PATH ${M3DC1_INCLUDE_PATH} M3DC1_INCLUDE_PATH)
#find_file (M3DC1_I m3dc1_


install( DIRECTORY DESTINATION ${CMAKE_SOURCE_DIR}/../lib )
install( TARGETS bjdl
  LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../lib
  ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/../lib)



#TARGET_LINK_LIBRARIES( bjdl 

#set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}../bin)

#add_subdirectory (../src)
