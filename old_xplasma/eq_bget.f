      subroutine eq_bget(ivec,zrho,zchi,zphi,irzmode,
     >   bvec,nregion,ierr)
c
c  vectorized interpolation routine.
c  get B field @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
c    if irzmode=1:  BR,BZ,Bphi
c    if irzmode=0:  Brho,Bchi,Bphi -- Brho in dir. normal to flux surface,
c                                     should be zero.
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer irzmode                   ! =1: B vector in cyl coords.
c
c output:
      real*8 bvec(3,abs(ivec))               ! the B field components
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
      call eqmap_bget(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   bvec,nregion,ierr)
c
      return
      end
c---------------------------------------------------------------------
      subroutine eq_bgrad(ivec,zrho,zchi,zphi,irzmode,
     >   bvec,gbtensr,nregion,ierr)
c
c  vectorized interpolation routine.
c  get B field @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c  get gradient of B
c     if irzmode=1:  (df/dR,df/dZ,(1/R)*df/dphi) of (BR,BZ,Bphi)
c     if irzmode=0:  (df/dsperp,df/dschi,(1/R)*df/dphi) of (Brho,Bchi,Bphi)
c          where sperp = dl in direction normal to flux surface,
c                schi  = dl in chi direction, normal to phi & tangent
c                           to flux surface
c
      implicit NONE
c input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer irzmode                   ! =1: B & grad. vector in cyl coords.
c
c output:
      real*8 bvec(3,abs(ivec))               ! the B field components
      real*8 gbtensr(3,3,abs(ivec))          ! gbtensr(1:3,j,iv)=grad(Bvec(j,iv))
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
      call eqmap_bgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   bvec,gbtensr,nregion,ierr)
c
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_bget(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   bvec,nregion,ierr)
c
c  vectorized interpolation routine.
c  get B field @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c
c    if irzmode=1:  BR,BZ,Bphi
c    if irzmode=0:  Brho,Bchi,Bphi -- Brho in dir. normal to flux surface,
c                                     should be zero.
c
c  also return (R,Z) of target points
c
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi) target coords
c
      integer irzmode                   ! =1: B vector in cyl coords.
C
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
      REAL*8 bvec(3,abs(ivec))               ! B field vector
C
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      integer imap
      REAL*8 zdumv(abs(ivec),1),zdumgv(3,abs(ivec),1)
      REAL*8 ztol
      REAL*8 fgrad(3,abs(ivec))              ! (R,Z,phi) gradient of f
      real*8 gbtensr(3,3,abs(ivec))
      integer iveca,lunerr
      logical :: ccwflag
C
C--------------------------
C
      if((irzmode.lt.0).or.(irzmode.gt.1)) then
         call eq_get_lunerr(lunerr)
         write(lunerr,*) ' ?eqmap_bget:  invalid irzmode=',irzmode
         ierr=1
      endif
c
      if(ierr.eq.1) then
         nregion=0
         bvec=0.0d0
         return
      endif
c
      nlist=0
      imap=3-irzmode
      ztol=0.0d0
      iveca=abs(ivec)
      ccwflag=ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,0,
     >   1,bvec,gbtensr,nlist,(/ 0 /),iveca,zdumv,zdumgv,
     >   nregion,ierr)
C
      return
      end
C--------------------------------------------------------------------
      subroutine eqmap_bgrad(ivec,zrho,zchi,zphi,zR,zZ,irzmode,
     >   bvec,gbtensr,nregion,ierr)
c
c  vectorized interpolation routine.
c  get B field @ (rho,chi,phi)
c        ivec < 0 indicates reversed chi
c  get gradient of B
c     if irzmode=1:  (df/dR,df/dZ,(1/R)*df/dphi) of (BR,BZ,Bphi)
c     if irzmode=0:  (df/dsperp,df/dschi,(1/R)*df/dphi) of (Brho,Bchi,Bphi)
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
C  output:
      real*8 zR(abs(ivec)),zZ(abs(ivec))          ! (R,Z) of target points
C
      real*8 bvec(3,abs(ivec))               ! the B field components
      real*8 gbtensr(3,3,abs(ivec))          ! gbtensr(1:3,j,iv)=grad(Bvec(j,iv))
      integer nregion(abs(ivec))             ! region code for (R,Z,phi) 1=inside
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      integer nlist
      integer imap
      REAL*8 zdumv(abs(ivec),1),zdumgv(3,abs(ivec),1)
      REAL*8 ztol
      integer iveca,lunerr
      logical :: ccwflag
C
C--------------------------
C
      if((irzmode.lt.0).or.(irzmode.gt.1)) then
         call eq_get_lunerr(lunerr)
         write(lunerr,*) ' ?eqmap_bget:  invalid irzmode=',irzmode
         ierr=1
      endif
c
      if(ierr.eq.1) then
         nregion=0
         return
      endif
c
      nlist=0
      imap=3-irzmode
      ztol=0.0d0
      iveca=abs(ivec)
      ccwflag=ivec.gt.0
c
      call eqi_getall(iveca,ccwflag,zR,zZ,zphi,imap,zrho,zchi,ztol,1,
     >   1,bvec,gbtensr,nlist,(/ 0 /),iveca,zdumv,zdumgv,
     >   nregion,ierr)
c
      return
      end
