      subroutine eqrz_bffget(ivec,zR,zZ,zphi,ztol,
     >   bvec,nlist,ifcns,ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
c  get (BR,BZ,Bphi) (BR,BZ,Bphi) @ (R,Z,phi)
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
      real*8 bvec(3,ivec)               ! the B field values returned
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
      call eqrzmap_bffget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,nlist,ifcns,ivecd,fvals,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eqrz_bffgrad(ivec,zR,zZ,zphi,ztol,
     >   bvec,gbtensr,nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
c  get cylindric coord.gradient of f (df/dR,df/dZ,(1/R)*df/dphi) @ (R,Z,phi)
c  get indicated (BR,BZ,Bphi) @ (R,Z,phi)
c  get [d/dR,d/dZ,(1/R)d/dphi) of each B field coordinate also.
c      gbtensr(1:3,j,k)=(d/dR,d/dZ,(1/R)d/dphi)(bvec(j,k))
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
      real*8 bvec(3,ivec)               ! the field values returned
      REAL*8 gbtensr(3,3,ivec)          ! gbtensor(*,j)=d[bvec(j)]/d[R,Z,phi]
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
      call eqrzmap_bffgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,gbtensr,nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_bffget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,nlist,ifcns,ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
c  get indicated (BR,BZ,Bphi) @ (R,Z,phi)
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
      real*8 bvec(3,ivec)               ! the field values returned
      REAL*8 fvals(ivecd,nlist)         ! function value
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      REAL*8 gbtensr(3,3,ivec)          ! dummy gbtensor
      REAL*8 fgrads(3,ivecd,max(1,nlist)) ! dummy gradients array
C
C--------------------------
C
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,0,
     >   1,bvec,gbtensr,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqrzmap_bffgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,gbtensr,nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (R,Z,phi)
c  get cylindric coord.gradient of f (df/dR,df/dZ,(1/R)*df/dphi) @ (R,Z,phi)
c  get indicated (BR,BZ,Bphi) @ (R,Z,phi)
c  get [d/dR,d/dZ,(1/R)d/dphi) of each B field coordinate also.
c      gbtensr(1:3,j,k)=(d/dR,d/dZ,(1/R)d/dphi)(bvec(j,k))
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
      real*8 bvec(3,ivec)               ! the field values returned
      REAL*8 gbtensr(3,3,ivec)          ! gbtensor(*,j)=d[bvec(j)]/d[R,Z,phi]
      REAL*8 fvals(ivecd,nlist)         ! function value
      REAL*8 fgrads(3,ivecd,nlist)      ! cartesian gradient of f
C
      integer nregion(ivec)             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C--------------------------
C
      call eqi_getall(ivec,.TRUE.,zR,zZ,zphi,1,zrho,zchi,ztol,1,
     >   1,bvec,gbtensr,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
