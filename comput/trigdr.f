C---------------------------------------
C  these are single precision REAL trig routines
C  sindr,cosdr,tandr -- arguments in degrees
C  asindr,acosdr,atandr -- inverse routines -- result in degrees
C
C  these are intended as (single precision only) substitutes for the
C  non-standard f77 intrinsic functions sind,cosd,tand,asind,acosd,atand
C  the latter not being available on all machines.
C
C  dmc 23 June 1999
C
C---------------------------------------
      real function sindr(zdeg)
      real zdeg
c
c  sine, arg in degrees, single precision
c
      data zpi/3.14159265/
 
      zrad=zdeg*zpi/180.0
      sindr=sin(zrad)
 
      return
      end
C---------------------------------------
      real function cosdr(zdeg)
      real zdeg
c
c  cosine, arg in degrees, single precision
c
      data zpi/3.14159265/
 
      zrad=zdeg*zpi/180.0
      cosdr=cos(zrad)
 
      return
      end
C---------------------------------------
      real function tandr(zdeg)
      real zdeg
c
c  tangent, arg in degrees, single precision
c
      data zpi/3.14159265/
 
      zrad=zdeg*zpi/180.0
      tandr=tan(zrad)
 
      return
      end
C---------------------------------------
      real function asindr(zarg)
c
c  arcsine, returned in degrees, single precision
c
      real zarg
 
      data zpi/3.14159265/
 
      zrad = asin(zarg)
      zdeg = zrad*180.0/zpi
 
      asindr = zdeg
      return
      end
C---------------------------------------
      real function acosdr(zarg)
c
c  arc-cosine, returned in degrees, single precision
c
      real zarg
 
      data zpi/3.14159265/
 
      zrad = acos(zarg)
      zdeg = zrad*180.0/zpi
 
      acosdr = zdeg
      return
      end
C---------------------------------------
      real function atandr(zarg)
c
c  arctangent, returned in degrees, single precision
c
      real zarg
 
      data zpi/3.14159265/
 
      zrad = atan(zarg)
      zdeg = zrad*180.0/zpi
 
      atandr = zdeg
      return
      end
