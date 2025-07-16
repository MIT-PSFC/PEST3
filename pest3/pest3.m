%
% m file plot eigensolution form the pest3.cdf file
% a pletzer oct 19 2000
%
clear all

addpath /usr/local/mexcdf/netcdf
ncstartup

file = 'pest3.cdf'
nc = netcdf(file, 'nowrite');

description = nc.description(:);
allvars = var(nc);                                 % Get variable data.
alldims = dim(nc);                                 % Get the dimensions.
allatts = att(nc);                                 % Get all attributes.

pa = nc{'pa'}(:);
qa = nc{'qa'}(:);
ga = nc{'ga'}(:);
dr = nc{'dr'}(:);
di = nc{'di'}(:);

cmatch =  nc{'cmatch'}(:);
dprim_re =  nc{'dprim_re'}(:);
dprim_im =  nc{'dprim_im'}(:);
gprim_re =  nc{'gprim_re'}(:);
gprim_im =  nc{'gprim_im'}(:);
psisin   =  nc{'psisin'}(:);
xmu   =  nc{'xmu'}(:);

sprintf('psi_s^(2 mu) Delta-prime = %10.5f + i*%10.5f', dprim_re*cmatch*psisin^(2*xmu), dprim_im*cmatch*psisin^(2*xmu))
sprintf('psi_s^(2 mu) Gamma-prime = %10.5f + i*%10.5f', gprim_re*cmatch*psisin^(2*xmu), gprim_im*cmatch*psisin^(2*xmu))

psinod = nc{'psinod'}(:);
psinew = nc{'psinew'}(:);
xa = nc{'xa'}(:);
za = nc{'za'}(:);
x1frbo_re = nc{'x1frbo_re'}(:);
x1frbo_im = nc{'x1frbo_im'}(:);
xisolo_re = nc{'xisolo_re'}(:);
xisolo_im = nc{'xisolo_im'}(:);
xisole_re = nc{'xisole_re'}(:);
xisole_im = nc{'xisole_im'}(:);


nfourier = size(x1frbo_re, 1);
mf = -(nfourier-1)/2:(nfourier-1)/2;
[ns, nt1] = size(xa);
t = linspace(0, 2*pi, nt1);
cosmt = cos(t'*mf)';
sinmt = sin(t'*mf)';

large_solution = x1frbo_re'*cosmt - x1frbo_im'*sinmt;

small_solution = 0.5*( ...
dprim_re*interp1(psinod, xisolo_re', psinew)*cosmt - ...
dprim_im*interp1(psinod, xisolo_im', psinew)*sinmt + ...
gprim_re*interp1(psinod, xisole_re', psinew)*cosmt - ...
gprim_im*interp1(psinod, xisole_im', psinew)*sinmt );

total_solution = large_solution+small_solution;

figure(1)
s = psinew/max(psinew);
subplot(2,2,1), plot(s, pa), title('pressure vs \psi/\psi_a [\mu_0 Pa]')
subplot(2,2,2), plot(s, qa), title('q vs \psi/\psi_a')
subplot(2,2,3), plot(s, ga), title('g vs \psi/\psi_a [T m]')
subplot(2,2,4), plot(s, di, s, dr, '--'), title('D_I(-) and D_R(--) vs \psi/\psi_a'), axis([0 1 -1.25 1.])

figure(2)
pcolor(xa, za, total_solution)
cmax = max(max(total_solution))/5; cmin = -cmax;
axis('image'), shading('interp'), caxis([cmin cmax])
hold on
contour(xa, za, total_solution, linspace(cmin, cmax, 11), 'k')
title('\xi . \nabla \psi')
xlabel('X'), ylabel('Z')


nc = close(nc);                                      % Close the file.


