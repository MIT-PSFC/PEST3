# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
isize_get_lsode.c
lsode_r8.F
srcom.f
xerrwv.f
)

add_library(lsode ${LIB_FILES})
