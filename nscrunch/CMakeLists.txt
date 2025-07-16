# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
)

set(LIB_FILES
eq_arc_findarc.f90
eq_arc_fmidplane.f90
eq_arc_scrunch.f90
eq_arc_scrunch_mod.F90
scraborter.f90
scradphase.f90
scramplituddes.F90
scrclear.f90
scrcptimi.f90
scrdescur.f90
scrdesdbg.F90
scrdesreg.F90
screvolvedes.f90
scrfixaraydes.f90
scrfree.f90
scrfsetup.f90
scrfunct1des.F90
scrgctnorm.f90
scrgetangldes.F90
scringuessdes.F90
scriscrunch.f90
scrmemory.f90
scrmfalsi.f90
scrmomintrp.f90
scrorderdes.f90
scrphasead.f90
scrpreset.f90
scrprintdes.f90
scrrestartdes.F90
scrscrmod.f90
scrscrunchme.f90
scrsetlun.f90
scrsincos.f90
scrunch_cleanup.f90
scrunch_rz.f90
scrxceptcon.f90
scryceptcon.f90
)

add_library(nscrunch ${LIB_FILES})

