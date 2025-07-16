# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
Esc2Any.c
bln.f
ecB.c
ecD.c
ecE.c
ecM.c
ecT.c
esc.h
escA.c
escB.c
escC.c
escD.c
escF.c
escG.c
escH.c
escIPC.c
escM.c
escO.c
escP.c
escQ.c
escS.c
escShm.c
escSt.c
escZ.c
esc_local.h
esc_test.f90
escget.c
escget.c.3
zkhF.f
zkhL.c
zkhS.c
)

add_library(esc ${LIB_FILES})
