%
% m file plot content of mapdsk.cdf file
% a pletzer oct 13 2000
%
%clear all

addpath /usr/local/mexcdf/netcdf
ncstartup

file = 'mapdsk.cdf'
file = input('enter file name: ','s')

nc = netcdf(file, 'nowrite');

description = nc.description(:);
allvars = var(nc);                                 % Get variable data.
alldims = dim(nc);                                 % Get the dimensions.
allatts = att(nc);                                 % Get all attributes.

sprintf('there are %d variables in file %s ', length(allvars), file)

sprintf('get 0-d data...')
title_nc = nc{'title'}(:);
date_nc = nc{'date'}(:);
sprintf('equilibrium generated on %s ', date_nc)
sprintf('%s ', title_nc)
mth = nc{'mth'}(:);
nosurf = nc{'nosurf'}(:);
sprintf('equilibrium dimensions mth * nosurf = %d * %d ', mth, nosurf)
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
psibig = nc{'PsiBig'}(:);
psi = psibig/(2*pi);
p = nc{'pa'}(:);
pp = nc{'ppa'}(:);
q = nc{'qa'}(:);
qp = nc{'qpa'}(:);
g = nc{'ga'}(:);
gp = nc{'gpa'}(:);
p = nc{'pa'}(:);

figure(1)
subplot(3,2,1); plot(psi, p), title('p(\psi)') 
subplot(3,2,2); plot(psi, pp), title('p\prime(\psi)')
subplot(3,2,3); plot(psi, q), title('q(\psi)')
subplot(3,2,4); plot(psi, qp), title('q\prime(\psi)')
subplot(3,2,5); plot(psi, g), title('g(\psi)')
subplot(3,2,6); plot(psi, gp), title('g\prime(\psi)')

jdotb = nc{'JdotBOverBSquare'}(:);
alpha = nc{'ballooning_alpha'}(:);
shear = nc{'local_shear_s'}(:);
rho   = nc{'surface_averaged_radius'}(:);
di   = nc{'Di'}(:);
dr   = nc{'Dr'}(:);
figure(2)
psimax = max(psi)
subplot(3,2,1); plot(psi, shear), title('s(\psi)') 
subplot(3,2,2); plot(psi, alpha), title('\alpha(\psi)')
subplot(3,2,3); plot(psi, rho), title('\rho(\psi)')
subplot(3,2,4); plot(psi, jdotb), title('<j.B>/<B^2>(\psi)')
subplot(3,2,5); plot(psi/psimax, di), title('D_I(\psi)'), axis([0,1,-5,5])
subplot(3,2,6); plot(psi/psimax, dr), title('D_R(\psi)'), axis([0,1,-5,5])


sprintf('get 2-d data...')
x = nc{'xa'}(:);
z = nc{'za'}(:);
jac = nc{'xjacob'}(:);
jacp = nc{'xjprym'}(:);
delta = nc{'delta'}(:);
qdelp = nc{'qdelp'}(:);
grpssq = nc{'grpssq'}(:);
grpsth = nc{'grpsth'}(:);
gptdth = nc{'gptdth'}(:);
gpsdth = nc{'gpsdth'}(:);
gserror = nc{'gserror'}(:);

xp = nc{'xpsi'}(:);
xt = nc{'xpth'}(:);
zp = nc{'zpsi'}(:);
zt = nc{'zpth'}(:);

psimax = max(psi);

figure(3)
subplot(2,2,1); pcolor(x, z, jac), title('J')
axis('image'), colorbar('vert'), shading('flat')
maxjp = 4*max(max(jac))/psimax;
subplot(2,2,2); pcolor(x, z, jacp), title('J\prime'), caxis([0 maxjp])
axis('image'), colorbar('vert'), shading('flat')
maxd = 2*pi/mth;
subplot(2,2,3); pcolor(x, z, delta), title('\delta'), caxis([-maxd maxd])
axis('image'), colorbar('vert'), shading('flat')
maxqdelp = max(q)*maxd*nosurf/psimax;
subplot(2,2,4); pcolor(x, z, qdelp), title('(q \delta)\prime'), 
caxis([-maxqdelp maxqdelp])
axis('image'), colorbar('vert'), shading('flat')


dx = 2*max(max(z));
figure(4)
subplot(2,2,1); pcolor(x, z, grpssq), title('|grad \psi|^2'), 
axis('image'), colorbar('vert'), shading('flat')
maxmet = psimax*2*pi/dx^2;
subplot(2,2,2); pcolor(x, z, grpsth), title('grad \psi . grad \theta'), 
caxis([-maxmet maxmet]), axis('image'), colorbar('vert'), shading('flat')
maxmet = 4*psimax*2*pi/(2*pi*dx^2);
subplot(2,2,3); pcolor(x, z, gptdth), title('(grad \psi . grad \theta)_\theta'), 
caxis([-maxmet maxmet]), axis('image'), colorbar('vert'), shading('flat')
maxmet = 4*psimax/dx^2;
subplot(2,2,4); pcolor(x, z, gpsdth), title('(|grad \psi|^2)_\theta'), 
caxis([-maxmet maxmet])
axis('image'), colorbar('vert'), shading('flat')


figure(5)
% plot the mesh
nlevels = 21; nrays = 32;
is = round(nosurf/nlevels); it = round(mth/nrays);
hold on
for i=nosurf:-is:1, plot(x(i,:),z(i,:),'b-'), end
for i=1:it:mth, plot(x(:,i),z(:,i), 'c-'), end
axis('image')
hold off

figure(6)
subplot(2,2,1); pcolor(x, z, xp), title('\partial_\psi X'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,xp,21,'k')
subplot(2,2,2); pcolor(x, z, xt), title('\partial_\theta X'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,xt,21,'k')
subplot(2,2,3); pcolor(x, z, zp), title('\partial_\psi Z'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,zp,21,'k')
subplot(2,2,4); pcolor(x, z, zt), title('\partial_\theta Z'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,zt,21,'k')

figure(7)
subplot(2,2,1); pcolor(x, z, xp), title('\partial_\psi X'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,xp,21,'k')
subplot(2,2,2); pcolor(x, z, jac), title('J'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,jac,21,'w')
subplot(2,2,3); pcolor(x, z, jacp), title('\partial_\psi J'), axis('image')
shading('flat'), colorbar('vert'), hold on, 
contour(x,z,jacp,linspace(-300,300.,21),'w')
subplot(2,2,4); pcolor(x, z, grpssq), title('|grad \psi|^2'), axis('image')
shading('flat'), colorbar('vert'), hold on, contour(x,z,grpssq,21,'w')

figure(8)
pcolor(x, z, log(abs(gserror))), shading('flat'),
title('Log rel. GS error'), axis('image')
caxis([-10 0]), colorbar('vert'), hold on 
contour(x, z, log(abs(gserror)), 11, 'w'), hold off


nc = close(nc);                                      % Close the file.


