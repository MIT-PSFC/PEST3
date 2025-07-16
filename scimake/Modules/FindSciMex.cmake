# - FindSciMex: Module to find include directories and
#   libraries for Mex.
#
# Module usage:
#   find_package(SciMex ...)
#
# May need to be changed to use SciFindPackge()

########################################################################
#
# FindSciMex
#
# $Rev$ $Date$
#
# Copyright 2012-2017, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
########################################################################

find_file(MEX NAMES mex.bat mex HINTS ENV PATH)

if (MEX)
  set(MEX_FOUND TRUE)
endif ()

if (MEX_FOUND)
  if (NOT MEX_FIND_QUIETLY)
    message(STATUS "Found MEX: ${MEX}")
  endif ()
  set(HAVE_MEX 1 CACHE BOOL "Whether have MEX")
else ()
   if (SciMex_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find MEX")
   endif ()
endif ()

