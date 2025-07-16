      subroutine eqxyz_ffget(ivec,xx,yy,zz,ztol,nlist,ifcns,
     >   ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (x,y,z)
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
      integer nlist                     ! number of functions
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! vector dimension -- output arrays
c
c output:
      real*8 fvals(ivecd,nlist)         ! the function values returned
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
      call eqxyzmap_ffget(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   nlist,ifcns,
     >   ivecd,fvals,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eqxyz_ffgrad(ivec,xx,yy,zz,ztol,nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
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
      integer nlist                     ! no. of fcns in list
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! vector dimension -- output arrays
c
c output:
      real*8 fvals(ivecd,nlist)         ! the function values returned
      real*8 fgrads(3,ivecd,nlist)      ! the function gradients returned
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
      call eqxyzmap_ffgrad(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqxyzmap_ffget(ivec,xx,yy,zz,zrho,zchi,zphi,
     >   ztol,nlist,ifcns,
     >   ivecd,fvals,nregion,ierr)
C
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results are returned.
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer nlist                     ! number of functions
      integer ifcns(nlist)              ! function id numbers
c
      integer ivecd                     ! output vector dimension, .ge.ivec
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,phi) (flux coords)
      REAL*8 fvals(ivecd,nlist)         ! function values
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
      call eqrzmap_ffget(ivec,zR,zZ,zphi,zrho,zchi,ztol,nlist,ifcns,
     >   ivecd,fvals,nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqxyzmap_ffgrad(ivec,xx,yy,zz,zrho,zchi,zphi,
     >   ztol,nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (x,y,z)
c  get cartesian gradient of f (df/dx,df/dy,df/dz) @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c  inverse map results are returned.
C
C  transform to (R,Z,phi), get fcn values and gradients in this space,
C  then transform back to (x,y,z)
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
c (x,y,z)->(rho,chi,phi) by Newton's method
c          (rho,cho,phi) maps to (x',y',z'); Newton iteration stops when
c                          max(|x-x'|,|y-y'|,|z-z'|).le.ztol*Raxis
c
      integer nlist                     ! number of functions
      integer ifcns(nlist)              ! function id numbers
c
      integer ivecd                     ! output vector dimension, .ge.ivec
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,phi) (flux coords)
      REAL*8 fvals(ivecd,nlist)         ! function values
      REAL*8 fgrads(3,ivecd,nlist)      ! cartesian gradients of functions
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
      call eqrzmap_ffgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
C
C  map (R,Z,phi) vectors to (x,y,z) system
C
      if(ierr.eq.0) then
         call eqi_rzxvec(ivec,zR,zZ,zphi,0,zdum,zdum,
     >      nlist,ivecd,fvals,fgrads)
      endif
C
      return
      end
C--------------------------------------------------------------------
