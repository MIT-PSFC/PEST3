# - FindSciGPerfTools: Module to find include directories and libraries
#   for GPerfTools. This module was implemented as there is no stock
#   CMake module for GPerfTools.
#
# This module can be included in CMake builds in find_package:
#   find_package(SciGPerfTools REQUIRED)
#
# This module will define the following variables:
#  HAVE_GPERFTOOLS         = Whether have the GPerfTools library
#  GPerfTools_INCLUDE_DIRS = Location of GPerfTools includes
#  GPerfTools_LIBRARY_DIRS = Location of GPerfTools libraries
#  GPerfTools_LIBRARIES    = Required libraries
#  GPerfTools_STLIBS       = Location of GPerfTools static library

######################################################################
#
# SciFindGPerfTools: find includes and libraries for GPerfTools.
#
# $Rev$ $Date$
#
# Copyright 2016-2017, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
######################################################################
set(SUPRA_SEARCH_PATH ${SUPRA_SEARCH_PATH})

if (WIN32)
  set(GPERFTOOLS_LIB_PREFIX "")
else (WIN32)
  set(GPERFTOOLS_LIB_PREFIX "lib")
endif (WIN32)

if (WIN32)
  set(GPERFTOOLS_LIB_SUFFIX "lib")
else (WIN32)
  set(GPERFTOOLS_LIB_SUFFIX "a")
endif (WIN32)

if (NOT DEFINED GPerfTools_SEARCH)
  if (DEFINED GPerfTools_FIND_VERSION)
    set(GPerfTools_SEARCH "gperftools${GPerfTools_FIND_VERSION}")
  else ()
    set(GPerfTools_SEARCH "gperftools")
  endif ()
endif ()
if (NOT DEFINED GPerfTools_SEARCH_HEADERS)
  set(GPerfTools_SEARCH_HEADERS "gperftools/tcmalloc.h")
endif ()
if (NOT DEFINED GPerfTools_SEARCH_LIBS)
  set(GPerfTools_SEARCH_LIBS "${GPERFTOOLS_LIB_PREFIX}tcmalloc_minimal.${GPERFTOOLS_LIB_SUFFIX}")
endif ()

SciFindPackage(PACKAGE "GPerfTools"
              INSTALL_DIR ${GPerfTools_SEARCH}
              HEADERS ${GPerfTools_SEARCH_HEADERS}
              LIBRARIES ${GPerfTools_SEARCH_LIBS}
              )

if (GPERFTOOLS_FOUND)
  message(STATUS "Found GperfTools")
  set(HAVE_GPERFTOOLS 1 CACHE BOOL "Whether have the GPerfTools library")
else ()
  message(STATUS "Did not find GPerfTools. Use -DGPERFTOOLS_DIR to specify the installation directory.")
  if (SciGPerfTools_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()

