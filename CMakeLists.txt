###################################################
#
# CMakeLists.txt for pets3code
# $Id: $
#
###################################################
# Project information
project(pets3code)
set(VERSION_MAJOR "3")
set(VERSION_MINOR "4")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8.8)

#######################################################
#
# Options
#
#######################################################

option(LINK_WITH_MKL "Link with MKL" OFF)
option(LINK_WITH_MKL_SHARED "Link with MKL shared libraries" OFF)
option(ENABLE_SHARED "Enable shared libraries" FALSE)

###################################################
#
# Initialization
#
###################################################

if (APPLE)
  add_definitions(-D__OSX)
else (APPLE)
  if (LINUX)
    add_definitions(-D__LINUX)
  endif (LINUX)
endif (APPLE)

#set(BUILD_WITH_SHARED_RUNTIME FALSE)
if (WIN32)
  # this is to prevent scimake to remove the /MD compiler option
  set(ENABLE_SHARED TRUE)
endif()

include(${CMAKE_SOURCE_DIR}/scimake/SciInit.cmake)

set(PESTCODE_VERSION "${VERSION}-${PROJECT_REV}")

###################################################
#
# Find packages and set includes
#
###################################################


find_package(SciMdsplus)
if (NOT MDSPLUS_FOUND) 
  add_definitions(-D__NOMDSPLUS)
endif (NOT MDSPLUS_FOUND)

find_package(SciDoxygen)
# requires hdf5 version >= 1.8, TODO: check for this
find_package(SciHdf5)
if (HDF5_FOUND)
  # Depending on build options hdf5 may require libz
  find_package(SciZ) 
endif ()
# netcdf AND netcdff 
find_package(SciNetcdff REQUIRED)

if(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED)
    find_package(SciMkl REQUIRED)
else(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED)
  find_package(SciBlasLapack REQUIRED)
endif(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL PGI)
  # fix CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES on hopper 
  # getting eh_util.c:(.text+0x227): undefined reference to `__zceh_uncaught_exception'
  set (fortran_link_libraries "")
  set (c_link_libraries "")
  foreach (lib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    message(STATUS ".... looking at ${lib}")
    if ("${lib}" STREQUAL stdc++)
      message(STATUS "==> removing ${lib} from the list of Fortran implicit libraries")
    else ()
      set (fortran_link_libraries ${fortran_link_libraries} ${lib})
    endif ()
  endforeach ()
  foreach (lib ${CMAKE_C_IMPLICIT_LINK_LIBRARIES})
    if (NOT "${lib}" STREQUAL "stdc++")
      set (c_link_libraries ${fortran_link_libraries} ${lib})
    else ()
      message(STATUS "==> removing ${lib} from the list of C implicit libraries")
    endif ()
  endforeach ()
  set (CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ${fortran_link_libraries})
  set (CMAKE_C_IMPLICIT_LINK_LIBRARIES ${c_link_libraries})
  message(STATUS "==> CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES = ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  message(STATUS "==> CMAKE_CXX_IMPLICIT_LINK_LIBRARIES = ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}")
  message(STATUS "==> CMAKE_C_IMPLICIT_LINK_LIBRARIES = ${CMAKE_C_IMPLICIT_LINK_LIBRARIES}")
endif ()

###################################################
#
# CTest
#
###################################################

enable_testing()

###################################################
#
# Set libraries to include
#
###################################################

set(math_libs)
if(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED)
  set (math_libs ${Mkl_LIBRARIES})
else(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED)
  if(ENABLE_SHARED)
    set (math_libs ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  else(ENABLE_SHARED)
    set (math_libs ${NETCDF_LIBRARIES} ${LAPACK_STLIBS} ${BLAS_STLIBS})
  endif(ENABLE_SHARED)
endif(LINK_WITH_MKL OR LINK_WITH_MKL_SHARED) 

set(h5_libs)
if(HDF5_FOUND)
  if (ENABLE_SHARED)
    set(h5_libs ${Hdf5_LIBRARIES} ${Z_LIBRARIES})
  else(ENABLE_SHARED)
    set(h5_libs ${Hdf5_STLIBS} ${Z_STLIBS})
  endif(ENABLE_SHARED)
  if (UNIX)
    # Required for hdf5 version 1.8.11 and greater
    set(h5_libs ${h5_libs} ${CMAKE_DL_LIBS}) 
  endif (UNIX)
endif(HDF5_FOUND)

###################################################
#
# Subdirectories
#
###################################################

add_subdirectory(generic_dummy)
add_subdirectory(portlib)
add_subdirectory(vaxonly)
add_subdirectory(lsode_linpack)
add_subdirectory(lsode)
add_subdirectory(comput)
add_subdirectory(smlib)
add_subdirectory(fluxav)
add_subdirectory(r8bloat)
add_subdirectory(mdstransp)
add_subdirectory(esc)
add_subdirectory(ezcdf)
add_subdirectory(nscrunch)
add_subdirectory(pspline)
add_subdirectory(geqdsk_mds)
add_subdirectory(xplasma2)
add_subdirectory(old_xplasma)
add_subdirectory(i2mex)
add_subdirectory(pest3)

###################################################
#
# Configure files
#
###################################################

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

