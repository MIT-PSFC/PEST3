      subroutine eq_ffget(ivec,zrho,zchi,zphi,nlist,ifcns,
     >   ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer nlist                     ! no. of fcns in list
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! output vector dimension
c
c output:
      real*8 fvals(ivecd,nlist)         ! the functions values returned
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
      call eqmap_ffget(ivec,zrho,zchi,zphi,zR,zZ,nlist,ifcns,ivecd,
     >   fvals,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eq_ffgrad(ivec,zrho,zchi,zphi,irzmode,nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (rho,chi,phi)
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
      integer nlist                     ! no. of fcns wanted
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! output vector dimension
c
c output:
      real*8 fvals(ivecd,nlist)         ! the functions values returned
      real*8 fgrads(3,ivecd,nlist)      ! the functions gradients returned
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
      call eqmap_ffgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,nlist,ifcns,
     >   ivecd,fvals,fgrads,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_ffget(ivec,zrho,zchi,zphi,zR,zZ,
     >   nlist,ifcns,ivecd,fvals,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
c  also return (R,Z) of target points
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer nlist                     ! no. of fcns in list
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! output vector dimension
C
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
      REAL*8 fvals(ivecd,nlist)         ! functions value
C
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer imap
      REAL*8 zdum,ztol
      REAL*8 fgrads(3,abs(ivec),nlist)       ! (R,Z,phi) gradient of f
      integer iveca
      logical ccwflag
C
C--------------------------
C
      ztol=0.0d0
      imap=2
      iveca=abs(ivec)
      ccwflag = ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,0,
     >   0,zdum,zdum,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_ffgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   nlist,ifcns,ivecd,fvals,fgrads,nregion,ierr)
c
c  vectorized interpolation routine.
c  get indicated functions f @ (rho,chi,phi)
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
c
      integer nlist                     ! no. of fcns in list
      integer ifcns(nlist)              ! xplasma fcn ids, functions wanted.
c
      integer ivecd                     ! output vector dimension
C
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
      REAL*8 fvals(ivecd,nlist)         ! functions values returned
      REAL*8 fgrads(3,ivecd,nlist)      ! functions gradients returned
C
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer imap
      REAL*8 zdum,ztol
      integer iveca,lunerr
      logical ccwflag
C
C--------------------------
C
      if((irzmode.lt.0).or.(irzmode.gt.1)) then
         call eq_get_lunerr(lunerr)
         write(lunerr,*) ' ?eqmap_ffget:  invalid irzmode=',irzmode
         ierr=1
      endif
c
      if(ierr.eq.1) then
         nregion=0
         return
      endif
c
      imap=3-irzmode
      ztol=0.0d0
      iveca=abs(ivec)
      ccwflag = ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,1,
     >   0,zdum,zdum,nlist,ifcns,ivecd,fvals,fgrads,
     >   nregion,ierr)
C
      return
      end
