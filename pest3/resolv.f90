      SUBROUTINE pstresolv(g,l,n,om,on)
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER L
      INTEGER NCX
      INTEGER NCY
      INTEGER NCZ
      INTEGER NCX1
      INTEGER NCY1
      INTEGER NCZ1
      INTEGER NCXY
      INTEGER NCX2
      INTEGER NCY2
      INTEGER NCZ2
      INTEGER NZ
      INTEGER K
      INTEGER KBAR
      INTEGER NCYLIM
      INTEGER IK
      INTEGER IKBAR
      INTEGER J
      INTEGER NCXLIM
      INTEGER JBAR
      INTEGER IJ
      INTEGER IJBAR
      INTEGER IBAR
      INTEGER I
      INTEGER IS


 COMPLEX*16, DIMENSION(*) 	 :: om
 COMPLEX*16, DIMENSION(*) 	 :: on
 COMPLEX*16	 	 	 :: u
 COMPLEX*16	 	 	 :: v
 COMPLEX*16	 	 	 :: t
 COMPLEX*16, DIMENSION(*) 	 :: g
 LOGICAL	 	 	 :: kcent
 LOGICAL	 	 	 :: jcent
 INTEGER, DIMENSION(3) 	 :: n
      ncx=n(1)
      ncy=n(2)
      ncz=n(3)
      ncx1=ncx-1
      ncy1=ncy-1
      ncz1=ncz-1
      ncxy=ncx*ncy
      ncx2=ncx/2
      ncy2=ncy/2
      ncz2=ncz/2
      nz=0
      if(l-1) 10000,10000,20000
10000 continue
      do 10 k=nz,ncz2
      kcent=k == nz.or.k.eq.ncz2
      kbar=ncz-k
      ncylim=ncy1
      if(.not.kcent) go to 111
      kbar=k
      ncylim=ncy2
  111 continue
      ik=1+ncxy*k
      ikbar=1+ncxy*kbar
      do 20 j=nz,ncylim
      jcent=j == nz.or.j.eq.ncy2
      ncxlim=ncx1
      jbar=ncy-j
      if(.not.jcent) go to 222
      jbar=j
      if(kcent) ncxlim=ncx2
  222 continue
      jcent=jcent.and.kcent
      ij=ik+ncx*j
      ijbar=ikbar+ncx*jbar
!   is=0 case
      u=g(ij)+conjg(g(ijbar))
      v=g(ij)-conjg(g(ijbar))
      v=cmplx(aimag(v),-real(v))
!   equivalence to multication by om(1)=(0,-1)
      g(ij)=u+v
      g(ijbar)=conjg(u-v)
      if(jcent) g(ij)=cmplx(real(u+v),real(u-v))
      ibar=ijbar+ncx-1
      i=ij+1
      do 30 is=1,ncxlim
      u=g(i)+conjg(g(ibar))
      t=om(is+1)
      v=g(i)-conjg(g(ibar))
      v=cmplx(real(v)*real(t)-aimag(v)*aimag(t),real(v)*aimag(t) &
   +real(t)*aimag(v))
      g(i)=u+v
      i=i+1
      g(ibar)=conjg(u-v)
      ibar=ibar-1
   30 continue
   20 continue
   10 continue
      return
20000 continue
      do 11 k=nz,ncz2
      kcent=k == nz.or.k.eq.ncz2
      ncylim=ncy1
      kbar=ncz-k
      if(.not.kcent) go to 333
      kbar=k
      ncylim=ncy2
  333 continue
      ik=1+ncxy*k
      ikbar=1+ncxy*kbar
      do 21 j=nz,ncylim
      jcent=j == nz.or.j.eq.ncy2
      jbar=ncy-j
      ncxlim=ncx1
      if(.not.jcent) go to 444
      jbar=j
      if(kcent) ncxlim=ncx2
  444 continue
      jcent=jcent.and.kcent
      ij=ik+ncx*j
      ijbar=ikbar+ncx*jbar
!   is=0 case
      u=g(ij)+conjg(g(ijbar))
      v=g(ij)-conjg(g(ijbar))
      v=cmplx(-aimag(v),real(v))
!   equivalence to multiplication by on(1)=(0,1)
      g(ij)=u+v
      g(ijbar)=conjg(u-v)
      if(jcent) g(ij)=cmplx(real(u-v),real(u+v))* 0.5_r8 
      i=ij+1
      ibar=ijbar+ncx-1
      do 31 is=1,ncxlim
      u=g(i)+conjg(g(ibar))
      t=on(is+1)
      v=g(i)-conjg(g(ibar))
      v=cmplx(real(v)*real(t)-aimag(v)*aimag(t),real(v)*aimag(t) &
   +real(t)*aimag(v))
      g(i)=u+v
      i=i+1
      g(ibar)=conjg(u-v)
      ibar=ibar-1
   31 continue
   21 continue
   11 continue
      return
      end

