# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
cgefa.f
cgesl.f
dgbfa.f
dgbsl.f
dgefa.f
dgesl.f
sgbfa.f
sgbsl.f
sgefa.f
sgesl.f
zgefa.f
zgesl.f
)

add_library(lsode_linpack ${LIB_FILES})
