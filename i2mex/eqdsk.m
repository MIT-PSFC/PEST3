% read eqdsk.cdf
%
% a. pletzer july 26 1999
%
% 
addpath /u/pletzer/mexcdf/netcdf

ncstartup

file=input('enter NetCDF file ')

nc = netcdf(file, 'nowrite');
description = nc.description(:)
variables = var(nc);                                 % Get variable data.
for i = 1:length(variables)
   sprintf('%d  %s', i, name(variables{i}))
   %disp([name(variables{i}) ' =']), disp(' ')
   %disp(variables{i}(:))
end

k = input('x data number? ')
x = variables{k}(:);
k = input('y data number? ')
y = variables{k}(:);
k = input('z data number? ')
z = variables{k}(:);

surf(x,y,z)

nc = close(nc);                                      % Close the file.

