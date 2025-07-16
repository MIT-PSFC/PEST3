# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${CMAKE_BINARY_DIR}/ezcdf
${CMAKE_BINARY_DIR}/pspline
${CMAKE_BINARY_DIR}/geqdsk_mds
)

set(LIB_FILES
cont_mod.f90
contour.f90
delta.f90
freeqbe_mod.f90
gausspoints.f90
getboozer.f90
getdelstar.f90
getg.f90
getj.f90
getp.f90
getphi.f90
getq.f90
gets.f90
getx.f90
getz.f90
glasser.f90
graphics.f90
i2mex_decimate.f90
i2mex_error.f90
i2mex_free.f90
i2mex_init.f90
i2mex_mod.f90
i2mex_theta_orient.f90
input.f90
input_geqdsk.f90
integral.f90
integrate1d.f90
interp1d.f90
metric.f90
mex2eqs.f90
mex2gtc.f90
original.f90
output.f90
parameters.f90
poloidalremap.f90
ratsurf.f90
refine.f90
save.f90
scaleg.f90
scaleq.f90
toaxis.f90
)

add_library(i2mex ${LIB_FILES})

