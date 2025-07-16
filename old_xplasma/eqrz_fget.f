      subroutine eqrz_fget(ivec,zR,zZ,zphi,ztol,ifcn,fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (R,Z,phi)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results not returned.
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zR(ivec),zZ(ivec),zphi(ivec) ! (R,Z,phi) target coordinates
c
      real*8 ztol                       ! rel. tolerance for inv. map
c
c (R,Z,phi)->(rho,chi,phi) by Newton's method
c          (rho,cho) maps to (R',Z'); Newton iteration stops when
c                          max(|R-R'|,|Z-Z'|).le.ztol*Raxis
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(ivec)                 ! the function values returned
      integer nregion(ivec)             ! region code, 1= inside core
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zrho(ivec),zchi(ivec)
c
c--------------------------
c
      call eqrzmap_fget(ivec,zR,zZ,zphi,zrho,zchi,ztol,ifcn,
     >   fval,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eqrz_fgrad(ivec,zR,zZ,zphi,ztol,ifcn,
     >   fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (R,Z,phi)
c  get cylindric coord.gradient of f (df/dR,df/dZ,(1/R)*df/dphi) @ (R,Z,phi)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c  inverse map results not returned.
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zR(ivec),zZ(ivec),zphi(ivec) ! (R,Z,phi) target coordinates
c
      real*8 ztol                       ! rel. tolerance for inv. map
c
c (R,Z,phi)->(rho,chi,phi) by Newton's method
c          (rho,cho) maps to (R',Z'); Newton iteration stops when
c                          max(|R-R'|,|Z-Z'|).le.ztol*Raxis
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(ivec)                 ! the function values returned
      real*8 fgrad(3,ivec)              ! the function gradients returned
      integer nregion(ivec)             ! region code, 1=inside core
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zrho(ivec),zchi(ivec)
c
c--------------------------
c
      call eqrzmap_fgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,ifcn,
     >   fval,fgrad,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_fget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   ifcn,fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (R,Z,phi)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results *are* returned.
c
c  method:
C  transform to (R,Z,phi), get fcn value in this space
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! (R,Z,phi) location
      REAL*8 ztol                       ! rel. tol. for (R,Z,phi)->flux coords
C
c (R,Z,phi)->(rho,chi,phi) by Newton's method
c          (rho,cho) maps to (R',Z'); Newton iteration stops when
c                          max(|R-R'|,|Z-Z'|).le.ztol*Raxis
c
      integer ifcn                      ! function number
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec)      ! (rho,chi,[phi]) (flux coords)
      REAL*8 fval(ivec)                 ! function value
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      REAL*8 zdum
      REAL*8 fgrad(3,ivec)              ! (R,Z,phi) gradient of f
C
C--------------------------
C
      nlist=1
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,0,
     >   0,zdum,zdum,nlist,ifcn,ivec,fval,fgrad,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_fgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   ifcn,fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (R,Z,phi)
c  get cylindric coord.gradient of f (df/dR,df/dZ,(1/R)*df/dphi) @ (R,Z,phi)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c  inverse map results *are* returned.
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! (R,Z,phi) location
      REAL*8 ztol                       ! rel. tol. for (R,Z,phi)->flux coords
C
c (R,Z,phi)->(rho,chi,phi) by Newton's method
c          (rho,cho) maps to (R',Z'); Newton iteration stops when
c                          max(|R-R'|,|Z-Z'|).le.ztol*Raxis
c
      integer ifcn                      ! function number
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec)      ! (rho,chi,phi) (flux coords)
      REAL*8 fval(ivec)                 ! function value
      REAL*8 fgrad(3,ivec)              ! cartesian gradient of f
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      REAL*8 zdum
C
C--------------------------
C
      nlist=1
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,1,
     >   0,zdum,zdum,nlist,ifcn,ivec,fval,fgrad,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
