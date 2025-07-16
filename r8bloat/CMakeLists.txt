# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
r8bloat_lib.csh
r8bloat_names.sub
r8bloata0.F
r8bloats0.F
r8fpolar.f
r8lokjac.F
r8ntk.f
r8tkbcon.F
r8tkbmmc.F
r8tkbmry.f
r8tkbndr.f
r8tkbparm.inc
r8tkbrad.f
r8tkbry2.f
r8tkbrys.F
r8tksbcon.f
)

add_library(r8bloat ${LIB_FILES})
