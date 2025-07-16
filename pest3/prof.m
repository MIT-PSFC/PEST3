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

nsin = size(xmu);
for is=1:nsin
	sprintf('psi_s^(2 mu) Delta-prime = %10.5f + i*%10.5f', ...
		dprim_re(is)*cmatch(is)*psisin(is)^(2*xmu(is)), ...
		dprim_im(is)*cmatch(is)*psisin(is)^(2*xmu(is)))
	sprintf('psi_s^(2 mu) Gamma-prime = %10.5f + i*%10.5f', ...
		gprim_re(is)*cmatch(is)*psisin(is)^(2*xmu(is)), ...
		gprim_im(is)*cmatch(is)*psisin(is)^(2*xmu(is)))
end

psinod = nc{'psinod'}(:);
psinew = nc{'psinew'}(:);
xa = nc{'xa'}(:);
za = nc{'za'}(:);

figure(1)
s = psinew/max(psinew);
subplot(2,2,1), plot(s, pa), title('pressure vs \psi/\psi_a [\mu_0 Pa]')
subplot(2,2,2), plot(s, qa), title('q vs \psi/\psi_a')
subplot(2,2,3), plot(s, ga), title('g vs \psi/\psi_a [T m]')
subplot(2,2,4), plot(s, di, s, dr, '--'), title('D_I(-) and D_R(--) vs \psi/\psi_a'), axis([0 1 -1.25 1.])


nc = close(nc);                                      % Close the file.


