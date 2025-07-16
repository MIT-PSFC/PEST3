      subroutine eq_fget(ivec,zrho,zchi,zphi,ifcn,
     >   fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(abs(ivec))                 ! the function values returned
      integer nregion(abs(ivec))             ! region code, 1= inside core
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zR(abs(ivec)),zZ(abs(ivec))
c
c--------------------------
c
      call eqmap_fget(ivec,zrho,zchi,zphi,zR,zZ,ifcn,
     >   fval,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eq_fgrad(ivec,zrho,zchi,zphi,irzmode,ifcn,
     >   fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c  get gradient of f
c     if irzmode=1:  (df/dR,df/dZ,(1/R)*df/dphi) @ (rho,chi,phi)
c     if irzmode=0:  (df/dsperp,df/dschi,(1/R)*df/dphi) @ (rho,chi,phi)
c          where sperp = dl in direction normal to flux surface,
c                schi  = dl in chi direction, normal to phi & tangent
c                           to flux surface
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer irzmode                   ! =1: grad. vector in cyl coords.
c
      integer ifcn                      ! xplasma fcn id, function wanted.
c
c output:
      real*8 fval(abs(ivec))                 ! the function values returned
      real*8 fgrad(3,abs(ivec))              ! the function gradients returned
      integer nregion(abs(ivec))             ! region code, 1=inside core
c
      integer ierr                      ! completion code, 0=OK
c
c--------------------------
c
      real*8 zR(abs(ivec)),zZ(abs(ivec))
c
c--------------------------
c
      call eqmap_fgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,ifcn,
     >   fval,fgrad,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_fget(ivec,zrho,zchi,zphi,zR,zZ,
     >   ifcn,fval,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
c  also return (R,Z) of target points
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
C
      integer ifcn                      ! function number
C
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
      REAL*8 fval(abs(ivec))                 ! function value
C
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      integer imap
      REAL*8 zdum,ztol
      REAL*8 fgrad(3,abs(ivec))              ! (R,Z,phi) gradient of f
      logical :: ccwflag
      integer :: iveca
C
C--------------------------
C
      nlist=1
      ztol=0.0d0
      imap=2
      iveca=abs(ivec)
      ccwflag = ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,0,
     >   0,zdum,zdum,nlist,ifcn,iveca,fval,fgrad,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_fgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   ifcn,fval,fgrad,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated function f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c  get gradient of f
c     if irzmode=1:  (df/dR,df/dZ,(1/R)*df/dphi) @ (rho,chi,phi)
c     if irzmode=0:  (df/dsperp,df/dschi,(1/R)*df/dphi) @ (rho,chi,phi)
c          where sperp = dl in direction normal to flux surface,
c                schi  = dl in chi direction, normal to phi & tangent
c                           to flux surface
c
c  also return (R,Z) of target points
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer irzmode                   ! =1: grad. vector in cyl coords.
C
      integer ifcn                      ! function number
C
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
      REAL*8 fval(abs(ivec))                 ! function values returned
      REAL*8 fgrad(3,abs(ivec))              ! function gradients returned
C
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      integer imap
      REAL*8 zdum,ztol
      logical :: ccwflag
      integer :: iveca,lunerr
C
C--------------------------
C
      if((irzmode.lt.0).or.(irzmode.gt.1)) then
         call eq_get_lunerr(lunerr)
         write(lunerr,*) ' ?eqmap_fget:  invalid irzmode=',irzmode
         ierr=1
      endif
c
      if(ierr.eq.1) then
         nregion=0
         fval=0.0d0
         return
      endif
c
      nlist=1
      imap=3-irzmode
      ztol=0.0d0
      iveca=abs(ivec)
      ccwflag = ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,1,
     >   0,zdum,zdum,nlist,ifcn,iveca,fval,fgrad,
     >   nregion,ierr)
C
      return
      end
