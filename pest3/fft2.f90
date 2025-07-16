      SUBROUTINE pstfft2(a,ifset,m,inv,s,om,on,iferr)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IFSET
      INTEGER IFERR
      INTEGER N1
      INTEGER N2
      INTEGER N3
      INTEGER MTT
      real*8 ROOT2
      INTEGER NX
      real*8 FN
      INTEGER I
      INTEGER ID
      INTEGER IL
      INTEGER IL1
      INTEGER MI
      INTEGER IDIF
      INTEGER KBIT
      INTEGER MEV
      INTEGER KL
      INTEGER KLAST
      INTEGER K
      INTEGER KD
      real*8 T
      INTEGER LFIRST
      INTEGER JLAST
      INTEGER L
      INTEGER JJDIF
      INTEGER K1
      INTEGER K2
      INTEGER K3
      real*8 R
      INTEGER JJ
      INTEGER ILAST
      real*8 AWR
      real*8 AWI
      INTEGER J
      INTEGER IC
      INTEGER I2
      INTEGER I2C
      INTEGER I2CC
      INTEGER I3
      INTEGER I3C
      INTEGER I3CC
      INTEGER I3CCC
      INTEGER NTSQ
      INTEGER M3MT
      INTEGER IGO3
      INTEGER N3VNT
      INTEGER MINN3
      INTEGER NTVN3
      INTEGER JJD3
      INTEGER M2MT
      INTEGER IGO2
      INTEGER N2VNT
      INTEGER MINN2
      INTEGER NTVN2
      INTEGER JJD2
      INTEGER M1MT
      INTEGER IGO1
      INTEGER N1VNT
      INTEGER MINN1
      INTEGER NTVN1
      INTEGER JJD1
      INTEGER JJ3
      INTEGER JPP3
      INTEGER IPP3
      INTEGER JP3
      INTEGER IP3
      INTEGER JJ2
      INTEGER JPP2
      INTEGER IPP2
      INTEGER JP2
      INTEGER IP2
      INTEGER JJ1
      INTEGER JPP1
      INTEGER IPP1
      INTEGER JP1
      INTEGER IP1

      equivalence (n1,n(1)),(n2,n(2)),(n3,n(3))
 REAL*8, DIMENSION(*) 	 :: a
 REAL*8, DIMENSION(*) 	 :: s
 REAL*8, DIMENSION(*) 	 :: om
 REAL*8, DIMENSION(*) 	 :: on
 REAL*8, DIMENSION(2) 	 :: w
 REAL*8, DIMENSION(2) 	 :: w2
 REAL*8, DIMENSION(2) 	 :: w3
 INTEGER, DIMENSION(*) 	 :: m
 INTEGER, DIMENSION(*) 	 :: inv
 INTEGER, DIMENSION(3) 	 :: np
 INTEGER, DIMENSION(3) 	 :: n
 INTEGER	 	 	 :: mt
 INTEGER	 	 	 :: nt
 COMMON /c2f90mt/ mt 
 COMMON /c2f90nt/ nt 
 COMMON /c2f90n/ n 
   12 mtt=max(m(1),m(2),m(3)) -2
      root2 = sqrt( 2._r8 ) 
      if (mtt-mt ) 14,14,13
   13 iferr=1
      return
   14 iferr=0
   16 if(ifset) 18,18,20
   18 nx= n1*n2*n3
      fn=4*nx
      do 19 i = 1,nx
      a(2*i-1) = a(2*i-1)/fn
   19 a(2*i) = -a(2*i)/fn
      go to 21
   20 CALL pstresolv(a,2,n,om,on)
   21 continue
      np(1)=n1*2
      np(2)= np(1)*n2
      np(3)=np(2)*n3
      do 250 id=1,3
      il = np(3)-np(id)
      il1 = il+1
      mi = m(id)
      if (mi)250,250,30
   30 idif=np(id)
      kbit=np(id)
      mev = 2*(mi/2)
      if (mi - mev )60,60,40
!
!     m is odd. do l=1 case
   40 kbit=kbit/2
      kl=kbit-2
      do 50 i=1,il1,idif
      klast=kl+i
      do 50 k=i,klast,2
      kd=k+kbit
!
!     do one step with l=1,j=0
!     a(k)=a(k)+a(kd)
!     a(kd)=a(k)-a(kd)
!
      t=a(kd)
      a(kd)=a(k)-t
      a(k)=a(k)+t
      t=a(kd+1)
      a(kd+1)=a(k+1)-t
   50 a(k+1)=a(k+1)+t
      if (mi - 1)250,250,52
   52 lfirst =3
!
!     def - jlast = 2**(l-2) -1
      jlast=1
      go to 70
!
!     m is even
   60 lfirst = 2
      jlast=0
   70 do 240 l=lfirst,mi,2
      jjdif=kbit
      kbit=kbit/4
      kl=kbit-2
!
!     do for j=0
      do 80 i=1,il1,idif
      klast=i+kl
      do 80 k=i,klast,2
      k1=k+kbit
      k2=k1+kbit
      k3=k2+kbit
!
!     do two steps with j=0
!     a(k)=a(k)+a(k2)
!     a(k2)=a(k)-a(k2)
!     a(k1)=a(k1)+a(k3)
!     a(k3)=a(k1)-a(k3)
!
!     a(k)=a(k)+a(k1)
!     a(k1)=a(k)-a(k1)
!     a(k2)=a(k2)+a(k3)*i
!     a(k3)=a(k2)-a(k3)*i
!
      t=a(k2)
      a(k2)=a(k)-t
      a(k)=a(k)+t
      t=a(k2+1)
      a(k2+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
!
      t=a(k3)
      a(k3)=a(k1)-t
      a(k1)=a(k1)+t
      t=a(k3+1)
      a(k3+1)=a(k1+1)-t
      a(k1+1)=a(k1+1)+t
!
      t=a(k1)
      a(k1)=a(k)-t
      a(k)=a(k)+t
      t=a(k1+1)
      a(k1+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
!
      r=-a(k3+1)
      t = a(k3)
      a(k3)=a(k2)-r
      a(k2)=a(k2)+r
      a(k3+1)=a(k2+1)-t
   80 a(k2+1)=a(k2+1)+t
      if (jlast) 235,235,82
   82 jj=jjdif   +1
!
!     do for j=1
      ilast= il +jj
      do 85 i = jj,ilast,idif
      klast = kl+i
      do 85 k=i,klast,2
      k1 = k+kbit
      k2 = k1+kbit
      k3 = k2+kbit
!
!     letting w=(1+i)/root2,w3=(-1+i)/root2,w2=i,
!     a(k)=a(k)+a(k2)*i
!     a(k2)=a(k)-a(k2)*i
!     a(k1)=a(k1)*w+a(k3)*w3
!     a(k3)=a(k1)*w-a(k3)*w3
!
!     a(k)=a(k)+a(k1)
!     a(k1)=a(k)-a(k1)
!     a(k2)=a(k2)+a(k3)*i
!     a(k3)=a(k2)-a(k3)*i
!
      r =-a(k2+1)
      t = a(k2)
      a(k2) = a(k)-r
      a(k) = a(k)+r
      a(k2+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
!
      awr=a(k1)-a(k1+1)
      awi = a(k1+1)+a(k1)
      r=-a(k3)-a(k3+1)
      t=a(k3)-a(k3+1)
      a(k3)=(awr-r)/root2
      a(k3+1)=(awi-t)/root2
      a(k1)=(awr+r)/root2
      a(k1+1)=(awi+t)/root2
      t= a(k1)
      a(k1)=a(k)-t
      a(k)=a(k)+t
      t=a(k1+1)
      a(k1+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
      r=-a(k3+1)
      t=a(k3)
      a(k3)=a(k2)-r
      a(k2)=a(k2)+r
      a(k3+1)=a(k2+1)-t
   85 a(k2+1)=a(k2+1)+t
      if(jlast-1) 235,235,90
   90 jj= jj + jjdif
!
!     now do the remaining j"s
      do 230 j=2,jlast
!
!     fetch w"s
!     def- w=w**inv(j), w2=w**2, w3=w**3
   96 i=inv(j+1)
   98 ic=nt-i
      w(1)=s(ic)
      w(2)=s(i)
      i2=2*i
      i2c=nt-i2
      if(i2c)120,110,100
!
!     2*i is in first quadrant
  100 w2(1)=s(i2c)
      w2(2)=s(i2)
      go to 130
  110 w2(1)= 0._r8 
      w2(2)= 1._r8 
      go to 130
!
!     2*i is in second quadrant
  120 i2cc = i2c+nt
      i2c=-i2c
      w2(1)=-s(i2c)
      w2(2)=s(i2cc)
  130 i3=i+i2
      i3c=nt-i3
      if(i3c)160,150,140
!
!     i3 in first quadrant
  140 w3(1)=s(i3c)
      w3(2)=s(i3)
      go to 200
  150 w3(1)= 0._r8 
      w3(2)= 1._r8 
      go to 200
!
  160 i3cc=i3c+nt
      if(i3cc)190,180,170
!
!     i3 in second quadrant
  170 i3c=-i3c
      w3(1)=-s(i3c)
      w3(2)=s(i3cc)
      go to 200
  180 w3(1)=- 1._r8 
      w3(2)= 0._r8 
      go to 200
!
!     3*i in third quadrant
  190 i3ccc=nt+i3cc
      i3cc = -i3cc
      w3(1)=-s(i3ccc)
      w3(2)=-s(i3cc)
  200 ilast=il+jj
      do 220 i=jj,ilast,idif
      klast=kl+i
      do 220 k=i,klast,2
      k1=k+kbit
      k2=k1+kbit
      k3=k2+kbit
!
!     do two steps with j not 0
!     a(k)=a(k)+a(k2)*w2
!     a(k2)=a(k)-a(k2)*w2
!     a(k1)=a(k1)*w+a(k3)*w3
!     a(k3)=a(k1)*w-a(k3)*w3
!
!     a(k)=a(k)+a(k1)
!     a(k1)=a(k)-a(k1)
!     a(k2)=a(k2)+a(k3)*i
!     a(k3)=a(k2)-a(k3)*i
!
      r=a(k2)*w2(1)-a(k2+1)*w2(2)
      t=a(k2)*w2(2)+a(k2+1)*w2(1)
      a(k2)=a(k)-r
      a(k)=a(k)+r
      a(k2+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
!
      r=a(k3)*w3(1)-a(k3+1)*w3(2)
      t=a(k3)*w3(2)+a(k3+1)*w3(1)
      awr=a(k1)*w(1)-a(k1+1)*w(2)
      awi=a(k1)*w(2)+a(k1+1)*w(1)
      a(k3)=awr-r
      a(k3+1)=awi-t
      a(k1)=awr+r
      a(k1+1)=awi+t
      t=a(k1)
      a(k1)=a(k)-t
      a(k)=a(k)+t
      t=a(k1+1)
      a(k1+1)=a(k+1)-t
      a(k+1)=a(k+1)+t
      r=-a(k3+1)
      t=a(k3)
      a(k3)=a(k2)-r
      a(k2)=a(k2)+r
      a(k3+1)=a(k2+1)-t
  220 a(k2+1)=a(k2+1)+t
!     end of i and k loops
!
  230 jj=jjdif+jj
!     end of j-loop
!
  235 jlast=4*jlast+3
  240 continue
!     end of  l  loop
!
  250 continue
!     end of  id  loop
!
!     we now have the complex fourier sums but their addresses are
!     bit-reversed.  the following routine puts them in order
      ntsq=nt*nt
      m3mt=m(3)-mt
  350 if(m3mt) 370,360,360
!
!     m3 gr. or eq. mt
  360 igo3=1
      n3vnt=n3/nt
      minn3=nt
      go to 380
!
!     m3 less than mt
  370 igo3=2
      n3vnt=1
      ntvn3=nt/n3
      minn3=n3
  380 jjd3 = ntsq/n3
      m2mt=m(2)-mt
  450 if (m2mt)470,460,460
!
!     m2 gr. or eq. mt
  460 igo2=1
      n2vnt=n2/nt
      minn2=nt
      go to 480
!
!     m2 less than mt
  470 igo2 = 2
      n2vnt=1
      ntvn2=nt/n2
      minn2=n2
  480 jjd2=ntsq/n2
      m1mt=m(1)-mt
  550 if(m1mt)570,560,560
!
!     m1 gr. or eq. mt
  560 igo1=1
      n1vnt=n1/nt
      minn1=nt
      go to 580
!
!     m1 less than mt
  570 igo1=2
      n1vnt=1
      ntvn1=nt/n1
      minn1=n1
  580 jjd1=ntsq/n1
  600 jj3=1
      j=1
      do 880 jpp3=1,n3vnt
      ipp3=inv(jj3)
      do 870 jp3=1,minn3
      go to (610,620),igo3
  610 ip3=inv(jp3)*n3vnt
      go to 630
  620 ip3=inv(jp3)/ntvn3
  630 i3=(ipp3+ip3)*n2
  700 jj2=1
      do 870 jpp2=1,n2vnt
      ipp2=inv(jj2)+i3
      do 860 jp2=1,minn2
      go to (710,720),igo2
  710 ip2=inv(jp2)*n2vnt
      go to 730
  720 ip2=inv(jp2)/ntvn2
  730 i2=(ipp2+ip2)*n1
  800 jj1=1
      do 860 jpp1=1,n1vnt
      ipp1=inv(jj1)+i2
      do 850 jp1=1,minn1
         go to (810,820),igo1
810      ip1=inv(jp1)*n1vnt
         go to 830
820      ip1=inv(jp1)/ntvn1
830      i=2*(ipp1+ip1)+1
         if (j-i) 840,850,850
840      t=a(i)
         a(i)=a(j)
         a(j)=t
         t=a(i+1)
         a(i+1)=a(j+1)
         a(j+1)=t
  850 j=j+2
  860 jj1=jj1+jjd1
!     end of jpp1 and jp2
!
  870 jj2=jj2+jjd2
!     end of jpp2 and jp3 loops
!
  880 jj3 = jj3+jjd3
!     end of jpp3 loop
!
  890 if(ifset)891,895,895
  891 do  i = 1,nx
        a(2*i) = -a(2*i)
      end do
      CALL pstresolv(a,1,n,om,on)
  895 return
      end


