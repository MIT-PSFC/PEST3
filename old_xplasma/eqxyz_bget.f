      subroutine eqxyz_bget(ivec,xx,yy,zz,ztol,bvec,nregion,ierr)
C
c  vectorized interpolation routine.
c  get (Bx,By,Bz) vector @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results not returned.
C
C  transform to (R,Z,phi), get BR,BZ,Bphi
C  then transform back to Bx,By,Bz
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
C  output:
      REAL*8 bvec(3,ivec)               ! B vector (Bx,By,Bz)
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,[phi]) (flux coords)
      real*8 zR(ivec),zdum,btmp(3),zcosphi,zsinphi
      integer iv
C
C--------------------------
C
      call eqxyzmap_bget(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   bvec,nregion,ierr)
C
      return
      end
C----------------------------------------------------------------------
      subroutine eqxyzmap_bget(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   bvec,nregion,ierr)
c
c  vectorized interpolation routine.
c  get (Bx,By,Bz) vector @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results are returned.
C
C  transform to (R,Z,phi), get BR,BZ,Bphi
C  then transform back to Bx,By,Bz
C
      implicit NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,[phi]) (flux coords)
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
C  output:
      REAL*8 bvec(3,ivec)               ! B vector (Bx,By,Bz)
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
c
C-----------------------------------------
c
      integer iv
      real*8 zR(ivec),zcosphi,zsinphi,btmp(3)
c
c-----------------------------------------
c
C  transform to cylindric coords
      call eq_rcyl(ivec,xx,yy,zR,zphi)
C
C  get vectors
      call eqrzmap_bget(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,nregion,ierr)
C
C  map (R,Z,phi) vectors to (x,y,z) system
C
      do iv=1,ivec
         zcosphi=cos(zphi(iv))
         zsinphi=sin(zphi(iv))
C
C  back transform B field
C
         btmp(3)=bvec(2,iv)
         btmp(1)=bvec(1,iv)*zcosphi -bvec(3,iv)*zsinphi
         btmp(2)=bvec(1,iv)*zsinphi +bvec(3,iv)*zcosphi
C
         bvec(1:3,iv)=btmp
      enddo
C
      return
      end
C---------------------------------------------------------------------
      subroutine eqxyz_bgrad(ivec,xx,yy,zz,ztol,bvec,gbtensr,nregion,
     >   ierr)
C
c  vectorized interpolation routine.
c  get (Bx,By,Bz) vector @ (x,y,z)
c
c      since this involves an inverse map to (rho,chi,phi)
c      a relative accuracy tolerance is required.
c
c      inverse map results not returned.
C
C  transform to (R,Z,phi), get BR,BZ,Bphi
C  then transform back to Bx,By,Bz
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
C  output:
      REAL*8 bvec(3,ivec)               ! B vector (Bx,By,Bz)
      REAL*8 gbtensr(3,3,ivec)          ! gbtensor(i,j)=d[bvec(j)]/d[x,y,z]
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,[phi]) (flux coords)
      real*8 zR(ivec),zdum,btmp(3),zcosphi,zsinphi
      integer iv
C
C--------------------------
C
      call eqxyzmap_bgrad(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   bvec,gbtensr,nregion,ierr)
C
      return
      end
C---------------------------------------------------------------------
      subroutine eqxyzmap_bgrad(ivec,xx,yy,zz,zrho,zchi,zphi,ztol,
     >   bvec,gbtensr,nregion,ierr)
C
C  **Vectorized**
C  get vector B(x,y,z) and tensor grad(Bx)/grad(By)/grad(Bz)
C
C  transform to (R,Z,phi), get BR,BZ,Bphi and the gradient in this space,
C  then transform back to (x,y,z)
C
      IMPLICIT NONE
c
c  input:
      integer ivec                      ! vector dimensioning
      REAL*8 xx(ivec),yy(ivec),zz(ivec) ! (x,y,z) location
      REAL*8 ztol                       ! rel. tol. for x,y,z -> flux coords
C
C  output:
      REAL*8 zrho(ivec),zchi(ivec),zphi(ivec) ! (rho,chi,[phi]) (flux coords)
      REAL*8 bvec(3,ivec)               ! B vector (Bx,By,Bz)
      REAL*8 gbtensr(3,3,ivec)          ! gbtensor(i,j)=d[bvec(j)]/d[x,y,z]
C
      integer nregion(ivec)             ! region code for (x,y,z) (cf eq_xinv).
      integer ierr                      ! completion code, 0=OK
C
C--------------------------
C
      real*8 zR(ivec),zdum
C
C--------------------------
C  transform to cylindric coords
      call eq_rcyl(ivec,xx,yy,zR,zphi)
C
C  get vectors
      call eqrzmap_bgrad(ivec,zR,zZ,zphi,zrho,zchi,ztol,
     >   bvec,gbtensr,nregion,ierr)
C
C  map (R,Z,phi) vectors to (x,y,z) system
C
      if(ierr.eq.0) then
         call eqi_rzxvec(ivec,zR,zZ,zphi,1,bvec,gbtensr,
     >      0,ivec,zdum,zdum)
      endif
C
      return
      end
