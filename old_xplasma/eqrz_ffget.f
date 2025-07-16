      subroutine eqrz_ffget(ivec,zR,zZ,zphi,ztol,
     >   nlist,ifcns,ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
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
      integer nlist                     ! no. of functions to evaluate
      integer ifcns(nlist)              ! list of function ids
c
      integer ivecd                     ! output vector dimension
c
c output:
      real*8 fvals(ivecd,nlist)         ! the function values returned
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
      call eqrzmap_ffget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   nlist,ifcns,ivecd,fvals,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eqrz_ffgrad(ivec,zR,zZ,zphi,ztol,
     >   nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
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
      integer nlist                     ! no. of functions to evaluate
      integer ifcns(nlist)              ! list of function ids
c
      integer ivecd                     ! output vector dimension
c
c output:
      real*8 fvals(ivecd,nlist)         ! the function values returned
      real*8 fgrads(3,ivecd,nlist)      ! the function gradients returned
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
      call eqrzmap_ffgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_ffget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   nlist,ifcns,ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
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
      integer nlist                     ! no. of functions to evaluate
      integer ifcns(nlist)              ! list of function ids
c
      integer ivecd                     ! output vector dimension
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec)      ! (rho,chi,[phi]) (flux coords)
      REAL*8 fvals(ivecd,nlist)         ! function value
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer i
      REAL*8 zdum
      REAL*8 fgrads(3,ivec,nlist)       ! (R,Z,phi) gradients of functions
C
C--------------------------
C
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,0,
     >   0,zdum,zdum,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_ffgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
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
      integer nlist                     ! no. of functions to evaluate
      integer ifcns(nlist)              ! list of function ids
c
      integer ivecd                     ! output vector dimension
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec)      ! (rho,chi,phi) (flux coords)
      REAL*8 fvals(ivecd,nlist)         ! function value
      REAL*8 fgrads(3,ivecd,nlist)      ! cartesian gradient of f
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer i
      REAL*8 zdum
C
C--------------------------
C
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,1,
     >   0,zdum,zdum,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
