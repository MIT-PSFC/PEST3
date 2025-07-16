%
% m file plot content of mapdsk.cdf file
% a pletzer oct 13 2000
%
%clear all

addpath /usr/local/mexcdf/netcdf
ncstartup

file = 'map01.cdf'
%file = input('enter file name: ','s')

nc = netcdf(file, 'nowrite');

description = nc.description(:);
allvars = var(nc);                                 % Get variable data.
alldims = dim(nc);                                 % Get the dimensions.
allatts = att(nc);                                 % Get all attributes.

sprintf('there are %d variables in file %s ', length(allvars), file)

sprintf('get 0-d data...')
title_nc = nc{'title'}(:);
date_nc = nc{'time'}(:);
sprintf('equilibrium generated on %s ', date_nc)
sprintf('%s ', title_nc)
nr = nc{'nr'}(:);
nz = nc{'nz'}(:);
mth = nc{'mth'}(:);
nosurf = nc{'nosurf'}(:);
sprintf('equilibrium dimensions nr * nz = %d * %d ', nr, nz)
sprintf('global parameters:')
b0 = nc{'B0'}(:);
ip = nc{'Ip'}(:);
beta = nc{'Beta'}(:);
betastar = nc{'BetaStar'}(:);
betan = nc{'BetaN'}(:);
li = nc{'li'}(:);
ppf = nc{'PPF'}(:);
sprintf('B0=%10.4f Ip=%10.4f li=%10.4f', b0, ip, li)
sprintf('Beta=%10.4f Beta*=%10.4f BetaN=%10.4f PPF=%10.4f ', beta, betastar, betan, ppf)

sprintf('get 1-d data...')
psi = nc{'psibar'}(:);
q = nc{'q'}(:);
sprintf('get 1-d data...')

figure(1)
plot(psi, q), title('q(\psi)')

sprintf('get (psi, the) data...')
xflux = nc{'x'}(:);
zflux = nc{'z'}(:);
xflux = xflux(1:nosurf, 1:mth+1);
zflux = zflux(1:nosurf, 1:mth+1);
psiflux = psi * ones(1, mth+1);

sprintf('get (R,Z) data...')
x = nc{'xcoord'}(:);
z = nc{'zcoord'}(:);
psixz = nc{'psixz'}(:);
bx = nc{'Bx'}(:);
bz = nc{'Bz'}(:);
bphi = nc{'Bphi'}(:);

figure(2)
subplot(2,2,1); pcolor(x, z, psixz), title('\psi')
axis('image'), colorbar('vert'), shading('interp'),
hold on, contour(xflux, zflux, psiflux, 21, 'w'), 
plot(xflux(nosurf,:), zflux(nosurf,:), 'k'), hold off
subplot(2,2,2); pcolor(x, z, bx), title('B_x'), 
axis('image'), colorbar('vert'), shading('interp')
hold on, contour(x, z, bx, 21, 'w'), 
plot(xflux(nosurf,:), zflux(nosurf,:), 'k'), hold off
subplot(2,2,3); pcolor(x, z, bz), title('B_z'),
axis('image'), colorbar('vert'), shading('interp')
hold on, contour(x, z, bz, 21, 'w'), 
plot(xflux(nosurf,:), zflux(nosurf,:), 'k'), hold off
subplot(2,2,4); pcolor(x, z, bphi), title('B_\phi'), 
axis('image'), colorbar('vert'), shading('interp')
hold on, contour(x, z, bphi, 21, 'w'), 
plot(xflux(nosurf,:), zflux(nosurf,:), 'k'), hold off




nc = close(nc);                                      % Close the file.


