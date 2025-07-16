# - FindOpenCascade: Module to find include directories and
#   libraries for Opencascade Community Edition
#
# Module usage:
#   find_package(OpenCascade ...)
#
# This module will define the following variables:
#  HAVE_OPENCASCADE, OPENCASCADE_FOUND = Whether libraries and includes are found
#  OpenCascade_INCLUDE_DIRS    = Location of OpenCascade includes
#  OpenCascade_LIBRARY_DIRS    = Location of OpenCascade libraries
#  OpenCascade_LIBRARIES       = Required libraries

######################################################################
#
# FindOpenCascade: find includes and libraries for OPENCASCADE
#
# $Rev$ $Date$
#
# Copyright 2012-2017, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
######################################################################

# Uses of these libs found from CMakeLists.txt in OPENCASCADE, and
# in doxygen documentation, which shows
#    Module FoundationClasses
#    Module ModelingData
#    Module ModelingAlgorithms
#    Module Visualization
#    Module ApplicationFramework
#    Module DataExchange
#    Module Draw
#
# Analyze dependencies using otool -L on OS X for more discreteness.
# Below is a layered list, top to bottom, left to right.
# ToDo: define the SEARCHHDRS
#
# Data exchange
# TKVRML

message (STATUS "")
message (STATUS "--------- FindSciOpenCascade seeking OpenCascade -----------")

set(OpenCascadeXdeIges_SEARCHLIBS TKXDEIGES)
set(OpenCascadeXdeStep_SEARCHLIBS TKXDESTEP)
# The libs below were required by TKXDEIGES
set(OpenCascadeXde_SEARCHLIBS TKCDF TKV3d TKService TKHLR TKOffset TKXmlXCAF TKXmlXCAF TKXmlL TKXml TKBinL TKBinXCAF TKBin)
set(OpenCascadeXde_SEARCHLIBS ${OpenCascadeXde_SEARCHLIBS} TKXCAF TKXmlXCAF TKBinXCAF TKCAF TKTObj TKLCAF TKVCAF)
# These disappeared in 7.1?
# set(OpenCascadeXde_SEARCHLIBS TKPLCAF PTKernel TKPShape TKShapeSchema TKPCAF TKStdLSchema)
# set(OpenCascadeXde_SEARCHLIBS ${OpenCascadeXde_SEARCHLIBS} TKXCAFSchema)
# Mesh contains triangulation
set(OpenCascadeMesh_SEARCHLIBS TKXMesh TKMesh)
set(OpenCascadeMesh_SEARCHHDRS XBRepMesh.hxx) # contains triangulation
set(OpenCascadeIges_SEARCHLIBS TKIGES)
set(OpenCascadeIges_SEARCHHDRS IGESFile_Read.hxx)
# IGES dependends on AdvAlgo
set(OpenCascadeAdvAlgo_SEARCHLIBS TKFillet TKBool TKPrim TKBO)
set(OpenCascadeStep_SEARCHLIBS TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase)
set(OpenCascadeStep_SEARCHHDRS STEPControl_Reader.hxx)
# STEP and IGES depend on this, but not STL
set(OpenCascadeIoBase_SEARCHLIBS TKXSBase)
set(OpenCascadeStl_SEARCHLIBS TKSTL)
set(OpenCascadeAlgo_SEARCHLIBS TKFeat TKShHealing TKTopAlgo TKGeomAlgo)
set(OpenCascadeModelData_SEARCHLIBS TKBRep TKG3d TKG2d TKGeomBase)
# AdvTools gone as of OPENCASCADE-0.17
# set(OpenCascadeTools_SEARCHLIBS TKMath TKAdvTools)
set(OpenCascadeTools_SEARCHLIBS TKMath)
set(OpenCascadeKernel_SEARCHLIBS TKernel)

# All the components
set(SciOpenCascade_ALL_COMPONENTS XdeIges XdeStep Xde Mesh Iges AdvAlgo Step IoBase Stl Algo ModelData Tools Kernel)

foreach (comp ${SciOpenCascade_FIND_COMPONENTS})
  set(OpenCascade${comp}_FIND TRUE)
endforeach ()

message(STATUS "Looking for components, ${SciOpenCascade_FIND_COMPONENTS}.")

# Enforce dependencies
if (OpenCascadeXdeIges_FIND)
  set(OpenCascadeIges_FIND TRUE)
  set(OpenCascadeXde_FIND TRUE)
endif ()
if (OpenCascadeXdeStep_FIND)
  set(OpenCascadeStep_FIND TRUE)
  set(OpenCascadeXde_FIND TRUE)
endif ()
if (OpenCascadeMesh_FIND)
  set(OpenCascadeBrep_FIND TRUE)
endif ()
if (OpenCascadeIges_FIND)
  set(OpenCascadeAdvAlgo_FIND TRUE)
  set(OpenCascadeIoBase_FIND TRUE)
endif ()
if (OpenCascadeStep_FIND)
  set(OpenCascadeIoBase_FIND TRUE)
endif ()
if (OpenCascadeIoBase_FIND OR OpenCascadeStl_FIND)
  set(OpenCascadeAlgo_FIND TRUE)
endif ()
if (OpenCascadeAlgo_FIND)
  set(OpenCascadeModelData_FIND TRUE)
endif ()
if (OpenCascadeModelData_FIND)
  set(OpenCascadeTools_FIND TRUE)
endif ()
if (OpenCascadeTools_FIND)
  set(OpenCascadeKernel_FIND TRUE)
endif ()

# Set the libraries
set(OpenCascade_SEARCHLIBS)
set(OpenCascade_comps)
foreach (pkg XdeIges XdeStep Xde Mesh Iges AdvAlgo Step IoBase Stl Algo ModelData Tools Kernel)
  if (DEBUG_CMAKE)
    message(STATUS "OpenCascade${pkg}_FIND = ${OpenCascade${pkg}_FIND}.")
  endif ()
  if (OpenCascade${pkg}_FIND)
    set(OpenCascade_comps ${OpenCascade_comps} ${pkg})
    set(OpenCascade_SEARCHLIBS ${OpenCascade_SEARCHLIBS} ${OpenCascade${pkg}_SEARCHLIBS})
  endif ()
endforeach ()
message(STATUS "After dependencies, looking for components, ${OpenCascade_comps}.")
message(STATUS "OpenCascade_SEARCHLIBS = ${OpenCascade_SEARCHLIBS}.")

# Worry about data exchange later

# To Do: Set variables for each group individually

# Set library subdirs
if (WIN32)
  set(incsubdir inc)
  if (${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(libsubdir win64)
  else ()
    set(libsubdir win32)
  endif ()
  if (CXX_VERSION MATCHES "^18\\.")
   set(libsubdir ${libsubdir}/vc12)
  endif ()
else ()
  set(incsubdir include/opencascade)
  set(libsubdir)
endif ()

# Only sersh build exists

# All the components
set(SEARCH_RESULTS PROGRAMS FILES INCLUDE_DIRS MODULE_DIRS LIBFLAGS LIBRARY_DIRS LIBRARY_NAMES LIBRARIES STLIBS)
if (WIN32)
  set(SEARCH_RESULTS ${SEARCH_RESULTS} DLLS)
endif ()
foreach (res ${SEARCH_RESULTS})
  set(OpenCascade_${res})
endforeach ()
set(OPENCASCADE_FOUND TRUE)
# Set the installation search directory for opencascade with no component suffix
if (USE_OPENCASCADE_SHARED)
  if (USE_PYC_LIBS)
    set(instdirs opencascade-pycsh opencascade-sersh)
  else ()
    set(instdirs opencascade-sersh opencascade-pycsh)
  endif ()
else ()
  SciGetInstSubdirs(opencascade instdirs)
endif ()

foreach (comp ${SciOpenCascade_ALL_COMPONENTS})
  if (OpenCascade${comp}_FIND)
    set(OpenCascade${comp}_ROOT_DIR ${OpenCascade_ROOT_DIR})
    SciFindPackage(PACKAGE OpenCascade${comp}
      INSTALL_DIRS "${instdirs}"
      HEADERS "${OpenCascade${comp}_SEARCHHDRS}"
      INCLUDE_SUBDIRS ${incsubdir}
      LIBRARIES "${OpenCascade${comp}_SEARCHLIBS}"
      LIBRARY_SUBDIRS "${libsubdir}/lib"
      PROGRAM_SUBDIRS "${libsubdir}/bin"
      FIND_QUIETLY
    )
    foreach (res ${SEARCH_RESULTS})
      set(OpenCascade_${res} ${OpenCascade_${res}} ${OpenCascade${comp}_${res}})
      set(OpenCascade${comp}_${res}
        ${OpenCascade${comp}_${res}}
        CACHE STRING "List of all ${res} for ${OpenCascade_${res}}"
      )
      endforeach ()
    string(TOUPPER OpenCascade${comp} pkguc)
    if (NOT ${pkguc}_FOUND)
      message(WARNING "${pkguc}_FOUND = ${${pkguc}_FOUND}.")
      set(OPENCASCADE_FOUND FALSE)
    endif ()
  endif ()
endforeach ()
foreach (res ${SEARCH_RESULTS})
  if (OpenCascade_${res})
    list(REMOVE_DUPLICATES OpenCascade_${res})
  endif ()
endforeach ()

find_library(OpenCascade_PLUGINS
  NAMES FWOSPlugin
  PATHS ${OpenCascade_LIBRARY_DIRS}
  NO_DEFAULT_PATH)

if (OPENCASCADE_FOUND)
  # message(STATUS "Found OpenCascade.")
  set(HAVE_OPENCASCADE 1 CACHE BOOL "Whether have OpenCascade library")
  SciPrintCMakeResults(OpenCascade)
else ()
  message(STATUS "Did not find OpenCascade.  Use -DOPENCASCADE_ROOT_DIR to specify the installation directory.")
  if (OpenCascade_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()

message (STATUS "--------- FindSciOpenCascade done with OpenCascade -----------")
message (STATUS "")

