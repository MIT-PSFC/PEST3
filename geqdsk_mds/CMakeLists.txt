# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${CMAKE_BINARY_DIR}/ezcdf
${CMAKE_BINARY_DIR}/pspline
)

set(LIB_FILES
cgeq_paths.f90
geqdsk_aux.f90
geqdsk_index1_sub.f90
geqdsk_index_mod.f90
geqdsk_lunmsg.f90
geqdsk_mod.F90
geqdsk_subs.f90
)

add_library(geqdsk_mds ${LIB_FILES})

