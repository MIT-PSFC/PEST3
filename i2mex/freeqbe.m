%
% m file plot content of mapdsk.cdf file
% a pletzer oct 13 2000
%
%clear all

addpath /usr/local/mexcdf/netcdf
ncstartup

file = 'nstx_07.ncd'
file = input('enter file name (nstx_07.ncd): ','s')

nc = netcdf(file, 'nowrite');

description = nc.description(:);
allvars = var(nc);                                 % Get variable data.
alldims = dim(nc);                                 % Get the dimensions.
allatts = att(nc);                                 % Get all attributes.

sprintf('there are %d variables in file %s ', length(allvars), file)

npsi = nc{'npsi'}(:);
nr = nc{'nr'}(:);
nz = nc{'nz'}(:);
ulength_m = nc{'ulength_m'}(:);
upsi_Wb__rad = nc{'upsi_Wb__rad'}(:);
up_Pa = nc{'up_Pa'}(:);
uRB_mT = nc{'uRB_mT'}(:);
psi0b=nc{'psi0b'}(:);
rzaxis = nc{'rzaxis'}(:);
psi = nc{'psi'}(:);
p = nc{'p'}(:);
RB = nc{'RB'}(:);
r = nc{'r'}(:);
z = nc{'z'}(:);
psi2 = nc{'psi2'}(:);

contour(r, z, psi2, 41), axis('image')
colorbar('vert')




nc = close(nc);                                      % Close the file.


