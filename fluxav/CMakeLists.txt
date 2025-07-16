# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
dqng_frhoth.F
dqng_fth.F
fluxav.f90
fluxav_balloc.f90
fluxav_brk_add.f90
fluxav_brk_add0.f90
fluxav_brk_init.f90
fluxav_opts.f90
)

add_library(fluxav ${LIB_FILES})
