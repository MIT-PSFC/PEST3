# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
${Netcdf_INCLUDE_DIRS}
)

set(LIB_FILES
arisyr.f
augdyr.f
build_date.F
checkid.f
ckreserv.f
d3dyr.f
eastyr.f
fldlim.f
get_tokyr.f90
iagtrim.f
ilnurd.f
iscmp0.f
iscmp1.f
isolver_config.f90
iteryr.f90
jetyr.f
kstryr.f
lslcps2.F
lslcps2s.F
ltxyr.f
lxasc.f
mastyr.f
mkbeastid.f
mkfrid.f
mkixrid.F
nstxyr.f
numfld.f
pbinrd.F
pbinwr.F
str_cbc.f
tftryr.f
trmds_parse.f
trmflnam.F
txtblok_in.f
txtblok_out.f
ubi3dcod.F
ubi3ncod.F
ubufio.F
udbsub.f
ufopen.F
vmsieee.F
vprotec.F
xfopen.F
)

add_library(vaxonly ${LIB_FILES})
