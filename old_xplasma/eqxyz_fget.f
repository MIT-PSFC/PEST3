      subroutine eqxyz_fget(ivec,xx,yy,zz,ztol,ifcn,fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results not returned.
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) target coordinates
c
      real*8 ztol                       ! rel. tolerance for inv. map
c
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(ivec)                 ! the function values returned
      integer nregion(ivec)             ! region code
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zrho(ivec),zchi(ivec),zphi(ivec)
c
c--------------------------
c
      call eqxyzmap_fget(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,ifcn,
     >   fval,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eqxyz_fgrad(ivec,xx,yy,zz,ztol,ifcn,
     >   fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c  get cartesian gradient of f (df/dx,df/dy,df/dz) @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c  inverse map results not returned.
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) target coordinates
c
      real*8 ztol                       ! rel. tolerance for inv. map
c
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(ivec)                 ! the function values returned
      real*8 fgrad(3,ivec)              ! the function gradients returned
      integer nregion(ivec)             ! region code
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zrho(ivec),zchi(ivec),zphi(ivec)
c
c--------------------------
c
      call eqxyzmap_fgrad(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,ifcn,
     >   fval,fgrad,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqxyzmap_fget(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   ifcn,fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results are returned.
c
c  method:
C  transform to (R,Z,phi), get fcn value in this space
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for (x,y,z)->flux coords
C
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer ifcn                      ! function number
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,phi) (flux coords)
      REAL*8 fval(ivec)                 ! function value
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      REAL*8 zR(ivec),zdum
C
C--------------------------
C  transform to cylindric coords
      call eq_rcyl(ivec,xx,yy,zR,zphi)
C
C  get values & vectors
C
      call eqrzmap_fget(ivec,zR,zZ,zphi,zrho,zchi,ztol,ifcn,
     >   fval,nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqxyzmap_fgrad(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   ifcn,fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c  get cartesian gradient of f (df/dx,df/dy,df/dz) @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c  inverse map results not returned.
c
C  transform to (R,Z,phi), get fcn value and gradient in this space,
C  then transform back to (x,y,z)
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for (x,y,z)->flux coords
C
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer ifcn                      ! function number
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,phi) (flux coords)
      REAL*8 fval(ivec)                 ! function value
      REAL*8 fgrad(3,ivec)              ! cartesian gradient of f
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      REAL*8 zR(ivec),zdum
C
C--------------------------
C  transform to cylindric coords
      call eq_rcyl(ivec,xx,yy,zR,zphi)
C
C  get values & vectors
C
      call eqrzmap_fgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,ifcn,
     >   fval,fgrad,nregion,ierr)
C
C  map (R,Z,phi) vectors to (x,y,z) system
C
      if(ierr.eq.0) then
         call eqi_rzxvec(ivec,zR,zZ,zphi,0,zdum,zdum,1,ivec,fval,fgrad)
      endif
C
      return
      end
C--------------------------------------------------------------------
