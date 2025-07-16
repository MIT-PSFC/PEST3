!....... changed fast fourier transform routines to eliminate the implicit
!........ use of memory through entry statements. jan 7, 1988 jm
!........ requires change in calling program
!........
!........   CALL pstfft2(a,ifset)  ---> CALL pstfft2(a,ifset,m,inv,s,om,on,iferr)
!........
!***************************************************************************
      SUBROUTINE pstff2prp(m,inv,s,om,on,iferr)
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IFERR
      INTEGER N1
      INTEGER N2
      INTEGER N3
      real*8 PI
      INTEGER NNT
      INTEGER NN
      real*8 ARGUM
      real*8 CARGUM
      real*8 WW
      real*8 VV
      INTEGER K
      INTEGER KC
      INTEGER NTV2
      real*8 THETA
      INTEGER JSTEP
      INTEGER JDIF
      INTEGER L
      INTEGER JSTEP2
      INTEGER JC1
      INTEGER JLAST
      INTEGER J
      INTEGER JC
      INTEGER JD
      INTEGER MTLEXP
      INTEGER LM1EXP
      INTEGER JJ


      equivalence (n1,n(1)),(n2,n(2)),(n3,n(3))
!
!     the following program computes the sin and inv tables.
!
 REAL*8, DIMENSION(*) 	 :: s
 REAL*8, DIMENSION(*) 	 :: om
 REAL*8, DIMENSION(*) 	 :: on
 INTEGER, DIMENSION(*) 	 :: inv
 INTEGER, DIMENSION(3) 	 :: n
 INTEGER, DIMENSION(*) 	 :: m
 INTEGER	 	 	 :: mt
 INTEGER	 	 	 :: nt
 COMMON /c2f90mt/ mt 
 COMMON /c2f90nt/ nt 
 COMMON /c2f90n/ n 
  900 mt=max(m(1),m(2),m(3)) -2
      mt = max(2,mt)
  904 if (mt-18) 906,906,905
  905 iferr=1
      return
  906 iferr=0
      n(1)=2**m(1)
      pi= 3.1415926536_r8 
      nnt=n1+1
      nn=n1*2+2
      argum=pi/float(n1)
      cargum= 0._r8 
      ww= 1._r8 
      vv= 0._r8 
! array(k,2)=-conjg(array(k,1))
      do 35 k=1,nnt,2
      kc=nn-k
!  om(k)=real(carg)
!  on(k+1)=aimag(carg)
!  om(k+1)=on(k+1)
!  on(k)=-om(k)
      om(k+1)=-ww
      on(k)=vv
      om(k)=on(k)
      on(k+1)=-om(k+1)
      if(k == 1 .OR. k == nnt) go to 34
      on(kc)=om(k)
      om(kc)=on(k)
      on(kc+1)=om(k+1)
      om(kc+1)=on(k+1)
   34 cargum=cargum+argum
      ww=cos(cargum)
      vv=-sin(cargum)
   35 continue
!  complex omc(1),onc(1)
!  equivalence (om(1),omc(1)),(on(1),onc(1))
!  icomp=cmplx(0, -1._r8 ) 
!  do 666 i=1,n1
!  omc(i)=icomp*omc(i)
!  666 onc(i)=icomp*onc(i)
      n(2)=2**m(2)
      n(3)=2**m(3)
      nt=2**mt
      ntv2=nt/2
!
!     set up sin table
!     theta=pie/2**(l+1) for l=1
  910 theta= .7853981634_r8 
!
!     jstep=2**(mt-l+1) for l=1
      jstep=nt
!
!     jdif=2**(mt-l) for l=1
      jdif=ntv2
      s(jdif)=sin(theta)
      do 950 l=2,mt
      theta=theta/ 2._r8 
      jstep2=jstep
      jstep=jdif
      jdif=jstep/2
      s(jdif)=sin(theta)
      jc1=nt-jdif
      s(jc1)=cos(theta)
      jlast=nt-jstep2
      if(jlast - jstep) 950,920,920
  920 do 940 j=jstep,jlast,jstep
      jc=nt-j
      jd=j+jdif
  940 s(jd)=s(j)*s(jc1)+s(jdif)*s(jc)
  950 continue
!
!     set up inv(j) table
!
  960 mtlexp=ntv2
!
!     mtlexp=2**(mt-l). for l=1
      lm1exp=1
!
!     lm1exp=2**(l-1). for l=1
      inv(1)=0
      do 980 l=1,mt
      inv(lm1exp+1) = mtlexp
      do 970 j=2,lm1exp
      jj=j+lm1exp
  970 inv(jj)=inv(j)+mtlexp
      mtlexp=mtlexp/2
  980 lm1exp=lm1exp*2
      return
      end



