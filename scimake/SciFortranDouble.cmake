######################################################################
#
# : Compute ntcc specific options
#
# $Rev$ $Date$
#
# Copyright 2012-2017, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
######################################################################

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FC_DOUBLE_FLAGS}")
# message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

