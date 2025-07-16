      subroutine r8bsmoo(zprof,zoption)
C
C  smooth radial profile -- triangular smooth
C  (sum of cubics i.e. definite integrals of products of linear pieces).
C
C     *** smoothing half-width (in xi):  namelist control DXBSMOO ***
C
C  weighting fcn about each point x0 is:
C
C   w(x)=
C     zero for |x-x0| >= dxbsmoo
C     (1/dxbsmoo)-|x-x0|/(dxbsmoo*dxbsmoo) for |x-x0| <= dxbsmoo
C
C  this makes a symmetrical triangle of area 1 centered on x0 and with base
C  2*dxbsmoo wide, height 1/dxbsmoo.
C
C  the smoothed function is treated as a sequence of linear segments.
C
C  f(x) = f(j)-x(j)*S + x*S    ... for x btw x(j) and x(j+1)
C                                  and S = (f(j+1)-f(j))/(x(j+1)-x(j))
C
C  piecewise integration of w(x)*f(x) (quadratic pieces) btw data pts
C  yields a sum of cubics evaluated as definite integrals.
C
C  curve will be "normalized" via znorm factor, prior to smooth, then
C  "un-normalized" on the way back out.  Allows user to control the "units"
C  of the space in which the smooth operation is carried out
C
C  zprof(...) -- profile to be smoothed
C  zoption(...) -- if zoption(lcentr).gt.0.0 then
C                    ...copy whole array into znorm
C                    ...use dxbsmoo as smoothing parameter
C                  if zoption(lcentr).le.0.0 then
C                    ...set znorm to 1.0
C                    ...set iopt = -zoption(lcentr)
C                       ...iopt = 0 means standard bdy conditions
C                       ...iopt = 1 means:  use reflection to fix point
C                          at outer bdy, position lep1.
C                    ...use zoption(lcp1) as smoothing parameter
C
C
C
C----------
C
      use r8bsmoo_mod
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER iopt,j,inzp1,in3zp1,iz
      INTEGER ii,icen,ism,iamdone
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zxsmoo,zedge,zdpt,zdiff,zdbxi,zdbx2,zxcen
      REAL*8 zx1b,zx2b,zf1b,zf2b,zx1a,zx2a,zf1a,zf2a,zfac,zdelx
      REAL*8 zsumx,zsmx2,zslope,zcns,za,zb,zc,zans
!============
      REAL*8 zprof(MJ)                    ! in:  unsmoothed, out:  smoothed
      REAL*8 zoption(MJ)                  ! controls
C
      REAL*8 znorm(MJ)                    ! presmooth normalization
C
      REAL*8 zwork(MJ)                    ! workspace
      REAL*8 zworkx(3*MJ)
      REAL*8 zworkf(3*MJ)
C
C-----------------------------
C check the smoothing parameter
C
      if(zoption(lcentr).gt.0.0D0) then
C beam code output smooth
         if(dxbsmoo.le.0.0D0) return
         dxbsmoo=max(0.01D0,min(0.2D0,dxbsmoo))
         zxsmoo=dxbsmoo
         iopt=0
         do j=1,mj
            znorm(j)=zoption(j)
         enddo
      else
C NCLASS code input smooth
         iopt=-zoption(lcentr)
         zxsmoo=zoption(lcp1)
         if(zxsmoo.le.0.0D0) return
         zxsmoo=max(0.01D0,min(0.2D0,zxsmoo))
         do j=1,mj
            znorm(j)=1.0D0
         enddo
      endif
C
C----------------
C
C  boundary conditions:
C    extrapolation by symmetric reflection at center
C    flat extrapolation at edge
C
      inzp1=nzones+1
      in3zp1=3*nzones+1
      do iz=1,nzones
         j=iz+lcentr-1
         zworkf(iz+nzones)=zprof(j)/znorm(j)
         zworkf(inzp1-iz) =zprof(j)/znorm(j)     ! reflection
         zworkf(in3zp1-iz)=zprof(ledge)/znorm(ledge)  ! flat extrap.
         zworkx(iz+nzones)=xi(j,2)
         zworkx(inzp1-iz)=-xi(j,2)
         zworkx(in3zp1-iz)=2*xi(lep1,1)-xi(j,2)
      enddo
C
      if(iopt.eq.1) then
         zedge=zprof(lep1)/znorm(lep1)
         zworkf(2*nzones+1)=zedge
C  anti-symmetric reflection beyond edge point; this causes the
C  symmetric smoothing operator to fix the edge point.
         do ii=2,nzones
            j=ledge+2-ii
            zdpt=zprof(j)/znorm(j)
            zdiff=zdpt-zedge
            zworkf(2*nzones+ii)=zedge-zdiff
         enddo
      endif
C
C  working from zworkf & zworkx compute smoothed results and store
C  back in zprof
C
      zdbxi=1.0D0/zxsmoo
      zdbx2=zdbxi*zdbxi
C
      do iz=1,nzones
         j=iz+lcentr-1
         zprof(j)=0.0D0  ! clear the sum
         icen=iz+nzones
         zxcen=zworkx(icen)
         zx1b=0.0D0
         zx2b=0.0D0
         zf1b=zworkf(icen)
         zf2b=zf1b
         do ism=1,nzones
            zx1a=zx1b
            zx2a=zx2b
            zf1a=zf1b
            zf2a=zf2b
            zx1b=zxcen-zworkx(icen-ism) ! reflected, for convenience
            zx2b=zworkx(icen+ism)-zxcen
            zf1b=zworkf(icen-ism)
            zf2b=zworkf(icen+ism)
            iamdone=0
C
C  deal with crossing left edge of weighting triangle
            if(zx1b.ge.zxsmoo) then
               iamdone=iamdone+1
               zfac=(zxsmoo-zx1a)/(zx1b-zx1a)
               zf1b=zf1a+(zf1b-zf1a)*zfac
               zx1b=zxsmoo
            endif
C
C  deal with crossing right edge of weighting triangle
            if(zx2b.ge.zxsmoo) then
               iamdone=iamdone+1
               zfac=(zxsmoo-zx2a)/(zx2b-zx2a)
               zf2b=zf2a+(zf2b-zf2a)*zfac
               zx2b=zxsmoo
            endif
C
C  add in left piece
            if(zx1a.lt.zxsmoo) then
               zdelx=zx1b-zx1a
               zsumx=zx1b+zx1a
               zsmx2=zx1b*zx1b+zx1b*zx1a+zx1a*zx1a
               zslope=(zf1b-zf1a)/zdelx
               zcns=(zf1a-zx1a*zslope)
               za=zcns*zdbxi
               zb=0.5D0*(zslope*zdbxi-zcns*zdbx2)
               zc=-0.33333333333333333D0*zslope*zdbx2
               zans=zdelx*(za+ zsumx*zb +zsmx2*zc)
               zprof(j)=zprof(j)+zans
            endif
C
C  add in right piece
            if(zx2a.lt.zxsmoo) then
               zdelx=zx2b-zx2a
               zsumx=zx2b+zx2a
               zsmx2=zx2b*zx2b+zx2b*zx2a+zx2a*zx2a
               zslope=(zf2b-zf2a)/zdelx
               zcns=(zf2a-zx2a*zslope)
               za=zcns*zdbxi
               zb=0.5D0*(zslope*zdbxi-zcns*zdbx2)
               zc=-0.33333333333333333D0*zslope*zdbx2
               zans=zdelx*(za+ zsumx*zb +zsmx2*zc)
               zprof(j)=zprof(j)+zans
            endif
C
C  test if done with this triangle
            if(iamdone.eq.2) go to 100
         enddo                          ! ism, triangle piece loop
 100     continue
         zprof(j)=zprof(j)*znorm(j)     ! un-normalize
      enddo                             ! iz, profile zone loop
C
C  all done
C
      return
      end
