# - FindSciSowing: Module to find sowing
#
# Module usage:
#   find_package(SciSowing ...)
#
# This module will define the following variables:
#  SOWING_FOUND         = Whether Sowing was found
#  Sowing_sowing    = Path to sowing executables

######################################################################
#
# SciSowing: Find Sowing
#
# $Rev$ $Date$
#
# Copyright 2012-2017, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
######################################################################

message("")
message("--------- FindSciSowing looking for sowing ---------")

SciFindPackage(PACKAGE Sowing
  PROGRAMS "bfort"
)

if (SOWING_FOUND)
      message(STATUS "Sowing_sowing found.")
      message(STATUS "Sowing_sowing = ${Sowing_sowing}")
      get_filename_component(sdir ${Sowing_PROGRAMS}/../.. REALPATH)
      set(Sowing_ROOT_DIR ${sdir})
      SciPrintVar(Sowing_ROOT_DIR)
endif ()

