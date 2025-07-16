# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
ezcdf.f90
ezcdf.hh
ezcdf_attrib.f90
ezcdf_exe.link
ezcdf_genget.f90
ezcdf_genput.f90
ezcdf_inqvar.f90
ezcdf_lib.public
ezcdf_opncls.f90
eztest.F90
eztest2.cc
handle_err.f90
)

add_library(ezcdf ${LIB_FILES})
