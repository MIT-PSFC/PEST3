# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${CMAKE_SOURCE_DIR}/ezcdf
${CMAKE_BINARY_DIR}/ezcdf
${CMAKE_BINARY_DIR}/i2mex
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
addvac.f90
aleg.f90
alva.f90
beta.f90
bigsol.f90
blkcpy.f90
blockt.f90
boundy.f90
bounknk.f90
bshift.f90
bsstep.f90
bzio.f90
cald.f90
cardmo.f90
cbxlu.f90
cdlhxv.f90
checkmetric.f90
chkinv.f90
com.f90
compmetric.f90
conalr.f90
concom.f90
conorm.f90
decomc.f90
decomp.f90
def.f90
defglo.f90
deltap.f90
delwbb.f90
delwbb2.f90
delwrb.f90
delwrr.f90
dsubi.f90
eigch.f90
eigen.f90
eigeng.f90
eigens.f90
empty.f90
ent23.f90
ent34.f90
epslon.f90
errmes.f90
exec.f90
exsolx.f90
fetchm.f90
ff2prp.f90
fft2.f90
fitter.f90
free.f90
free_mem.f90
frobe3.f90
funint.f90
fxfu.f90
gammln.f90
gelg.f90
get.f90
ggjcoef.f90
green.f90
htribk.f90
htridi.f90
hypdrv.f90
hypgeo.f90
hypgeo_info.f90
hypser.f90
indint.f90
inglo.f90
init.f90
integral.f90
inverc.f90
kernel.f90
labrt.f90
laeigen.f90
loadw.f90
mapper.hh
matelt.f90
matwrt.f90
memory.f90
metrealloc.f90
metric.f90
mmid.f90
nodes.f90
nr.f90
nrtype.f90
nrutil.f90
odeint.f90
output.f90
pest3.hh
pest3x.cpp
pstcheckmercier.f90
psttoaxis1.f90
pythag.f90
pzextr.f90
ratsu2.f90
remap2.f90
repet2.f90
reset.f90
resolv.f90
saver.f90
search.f90
set.f90
setlnk.f90
shape.f90
shifft.f90
solve.f90
solvec.f90
spl1d1.f90
spl1d2.f90
storrb.f90
storrr.f90
strhs2.f90
thefft.f90
thfft2.f90
tidyer.f90
tql2.f90
tqlrat.f90
uf.f90
ufpr.f90
vacda2.f90
vacds2.f90
vacdsk.f90
vacuum.f90
vax.f90
vbx.f90
vekit.f90
vf.f90
vfpr.f90
wallck.f90
wresult.f90
writem.f90
zeror.f90
zerot.f90
)

set(HEADER_FILES
CUCCCC.inc
mapper.hh
pest3.hh
)

add_library(pest3 ${LIB_FILES})
add_executable(pest3x pest3x.cpp)
target_link_libraries(pest3x pest3 
i2mex old_xplasma xplasma2 geqdsk_mds pspline esc nscrunch ezcdf mdstransp r8bloat fluxav smlib comput lsode lsode_linpack vaxonly portlib generic_dummy
${Netcdf_netcdff_LIBRARY} ${Netcdf_netcdf_LIBRARY} ${Mdsplus_LIBRARIES} ${math_libs}
)


install(FILES ${HEADER_FILES} DESTINATION INCLUDE)


