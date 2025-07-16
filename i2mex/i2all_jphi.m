clear all

addpath /usr/local/mexcdf/netcdf
ncstartup
nc = netcdf('MAPDSK.CDF', 'nowrite');

mth = nc{'mth'}(:);
nt1 = mth + 1;
ns = nc{'nosurf'}(:);

pp = nc{'ppa'}(:);
g = nc{'ga'}(:);
gp = nc{'gpa'}(:);
psi = nc{'PsiBig'}(:)/(2*pi);

x = nc{'xa'}(:);
z = nc{'za'}(:);

s = (psi-psi(1)) / (psi(ns)-psi(1));

figure(3)
plot(s, g-g(1))
title('g-g1')


jsubphi = - ((g .* gp) *ones(1, nt1)) - x.^2 .* (pp*ones(1,nt1));
figure(4) 
plot(x(:,2), jsubphi(:,2)), title(' covar J_\phi vs R (MAPDSK.CDF) ')

