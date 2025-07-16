# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${CMAKE_BINARY_DIR}/ezcdf
${CMAKE_BINARY_DIR}/pspline
)

if (MDSPLUS_FOUND)
  include_directories(${Mdsplus_mdsdescrip_h_INCLUDE_DIR})
endif (MDSPLUS_FOUND)

set(LIB_FILES
f_mdscacheconnect.F90
find_image_symbol.F
idescr_short.F
lmds_errstat.f90
mds__copywnull.f
mds__leniftrim.f
mds__short.c
mds_add_node.F
mds_add_tag.F
mds_addblanks.f
mds_close.F
mds_create_node.F
mds_delete_node.F
mds_descr_dum1.c
mds_dummy.c
mds_errstr.F
mds_open.F
mds_put.F
mds_put1d.F
mds_put2d.F
mds_put3d.F
mds_remove_tag.F
mds_setserver.F
mds_value.F
mdscacheconnect.c
mdserr.f
mdsfconnect.c
mdssetsocket_dum1.c
)

add_library(mdstransp ${LIB_FILES})

