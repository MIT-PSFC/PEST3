# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
copyy.f
ddrop.f
fileps.f
filers.f
filfcn.f
filfn6.f
filter.f
filtr6.f
filtr90.F
ierfcn.f
pctran.f
r4_qksmooth.f90
r8_pctran.f
r8_qksmooth.f90
r8ddrop.f
r8fileps.f
r8filers.f
r8filfn6.f
r8filter.f
r8filtr6.f
r8filtr90.f
r8ierfcn.f
r8wxcint.f
r8wxfint.f
wxcint.f
wxfint.f
)

add_library(smlib ${LIB_FILES})
