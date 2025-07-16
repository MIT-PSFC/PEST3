      subroutine pstspl1d1(n,x,f,w,iop,ij,a,b,c)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER IJ
      real*8 ZZ
      real*8 OZ
      real*8 TZ
      real*8 SZ
      INTEGER K
      INTEGER I
      INTEGER M
      INTEGER J1
      INTEGER J2
      real*8 CON
      real*8 DON
      real*8 E
      INTEGER K1
      INTEGER K2
      INTEGER J
      INTEGER K3
      INTEGER K4
      real*8 C1
      real*8 C2
      real*8 A1
      real*8 A2
      real*8 A3
      real*8 A4
      real*8 A5
      real*8 A6
      real*8 B1
      real*8 B2
      real*8 B3
      real*8 B4
      real*8 B5
      real*8 B6
      INTEGER L1
      INTEGER L2
      INTEGER L3
      real*8 BOB
      real*8 BILL
      INTEGER MK
      INTEGER ML
      real*8 D1
      real*8 D2

!     lcm (spl1d1)
!     where n= number of points in the interpolation
!           x= origin of table of independent variable
!           f= origin of table of dependent variable
!           w= an array of dimension n which contains the calculated
!              second derivatives upon return
!           iop= an array of dimension 2 which contains combinations of
!                the integers 1 thru 5 used to specify the boundary
!                conditions
!           ij= spacing in the f and w tables
!           a,b,c= arrays of dimension n used for temporary storage
!
!    comm made real in order to hold 8 characters in ibm-land
 CHARACTER*48 	 :: comm 
      data  comm  &
        /'pstspl1d1 n less than  4._r8  results incorrect.        '/
      data zz,oz,tz,sz/ 0.0e0_r8 , 1.0e0_r8 , 3.0e0_r8 ,  6.0e0_r8 / 
 REAL*8, DIMENSION(*) 	 :: x
 REAL*8, DIMENSION(2) 	 :: f
 REAL*8, DIMENSION(2) 	 :: w
 REAL*8, DIMENSION(2) 	 :: a
 REAL*8, DIMENSION(2) 	 :: b
 REAL*8, DIMENSION(2) 	 :: c
 INTEGER, DIMENSION(2) 	 :: iop
      k=n-1
      a(2)=-(x(2)-x(1))/sz
      b(2)=(x(3)-x(1))/tz
      w(ij+1)=(f(2*ij+1)-f(ij+1))/(x(3)-x(2))-(f(ij+1)-f(1)) &
/(x(2)-x(1))
      if (n-3)3,4,3
    3 do 10 i=3,k
      m=(i-1)*ij+1
      j1=m+ij
      j2=m-ij
      con=(x(i+1)-x(i-1))/tz
      don=(x(i)-x(i-1))/sz
      b(i)=con-(don**2)/b(i-1)
      e=(f(j1)-f(m))/(x(i+1)-x(i))-(f(m)-f(j2))/ &
(x(i)-x(i-1))
      w(m)=e-(don*w(j2))/b(i-1)
   10 a(i)=-(don*a(i-1))/b(i-1)
    4 k1=(n-2)*ij+1
      c(n-1)=-((x(n)-x(n-1))/sz)/b(n-1)
      w(k1)=w(k1)/b(n-1)
      a(n-1)=a(n-1)/b(n-1)
      k2=k-1
      if (n-3)7,8,7
    7 do 20 i=2,k2
      j=n-i
      con=(x(j+1)-x(j))/sz
      a(j)=(a(j)-con*a(j+1))/b(j)
      c(j)=-(con*c(j+1))/b(j)
      k3=(j-1)*ij+1
      m=k3+ij
   20 w(k3)=(w(k3)-con*w(m))/b(j)
    8 k4=(n-1)*ij+1
      if (iop(1)-5) 201,200,201
  201 c1=w(1)
      if (iop(2)-5) 203,202,203
  203 c2=w(k4)
      go to 205
  200 if (n-4)300,302,302
  302 a1=x(1)-x(2)
      a2=x(1)-x(3)
      a3=x(1)-x(4)
      a4=x(2)-x(3)
      a5=x(2)-x(4)
      a6=x(3)-x(4)
      w(1)=f(1)*(oz/a1+oz/a2+oz/a3)-a2*a3*f(ij+1)/(a1*a4*a5)+ &
 a1*a3*f(2*ij+1)/(a2*a4*a6 )-a1*a2*f(3*ij+1)/(a3*a5*a6)
      go to 201
  202 if (n-4)300,303,303
  303 b1=x(n)-x(n-3)
      b2=x(n)-x(n-2)
      b3=x(n)-x(n-1)
      b4=x(n-1)-x(n-3)
      b5=x(n-1)-x(n-2)
      b6=x(n-2)-x(n-3)
      l1=k4-ij
      l2=l1-ij
      l3=l2-ij
      w(k4)=-b2*b3*f(l3)/(b6*b4*b1)+b1*b3*f(l2)/(b6*b5*b2) &
 -b1*b2*f(l1)/(b4*b5*b3)+f(k4)*(oz/b1+oz/b2+oz/b3)
      go to 203
!  cdc compiler apparently permits transfer into the range of a do-loop
!   corrected for ibm  compilers (br)
! 205 do 50 i=1,k
 205          i    =    1
 2051 continue
      m=(i-1)*ij+1
      go to 60
   70 if (i-1)80,50,80
   80 w(1)=w(1)-bob*w(m)
      w(k4)=w(k4)-bill*w(m)
      a(1)=a(1)-bob*a(i)
      a(n)=a(n)-bill*a(i)
      c(1)=c(1)-bob*c(i)
      c(n)=c(n)-bill*c(i)
!     see   note at statement label 250
   50 continue
      i=i+1
      if ( i  <=  k )   go to 2051
!  50 continue
      go to 100
   60 mk=iop(1)
      go to (62,64,66,68,66),mk
   62 if (i-1)71,63,71
   63 a(1)=-oz
      c(1)=zz
      go to 500
   71 bob=zz
      go to 500
   64 if (i-1)73,76,73
   76 a(1)=-oz
      c(1)=zz
      w(1)=zz
      go to 500
   73 if (i-2)81,81,82
   81 bob=-c1
      go to 500
   82 bob=zz
      go to 500
   66 if (i-1)83,84,83
   84 a(1)=-(x(2)-x(1))/tz
      c(1)=zz
      w(1)=-c1+(f(ij+1)-f(1))/(x(2)-x(1))
      go to 500
   83 if (i-2)85,85,86
   85 bob=(x(2)-x(1))/sz
      go to 500
   86 bob=zz
      go to 500
   68 if (i-1)87,88,87
   88 a(1)=-oz
      c(1)=oz
      w(1)=zz
      go to 500
   87 bob=zz
  500 ml=iop(2)
      go to (120,130,140,150,140),ml
  120 if (i-1)121,122,121
  122 a(n)=zz
      c(n)=-oz
      go to 70
  121 bill=zz
      go to 70
  130 if (i-1)131,132,131
  132 a(n)=zz
      c(n)=-oz
      w(k4)=zz
      go to 70
  131 if (i-k)134,133,134
  133 bill=-c2
      go to 70
  134 bill=zz
      go to 70
  140 if (i-1)141,142,141
  142 a(n)=zz
      c(n)=(x(n-1)-x(n))/tz
      w(k4)=c2-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  141 if (i-k)143,144,143
  144 bill=(x(n)-x(n-1))/sz
      go to 70
  143 bill=zz
      go to 70
  150 if (i-1)151,152,151
  152 a(n)=zz
      c(n)=(x(n-1)+x(1)-x(n)-x(2))/tz
      w(k4)=(f(ij+1)-f(1))/(x(2)-x(1))-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  151 if (i-2)153,154,153
  154 bill=(x(2)-x(1))/sz
      go to 70
  153 if (i-k)155,156,155
  156 bill=(x(n)-x(n-1))/sz
      go to 70
  155 bill=zz
      go to 70
  100 con=a(1)*c(n)-c(1)*a(n)
      d1=-w(1)
      d2=-w(k4)
      w(1)=(d1*c(n)-c(1)*d2)/con
      w(k4)=(a(1)*d2-d1*a(n))/con
      do 110 i=2,k
      m=(i-1)*ij+1
  110 w(m)=w(m)+a(i)*w(1)+c(i)*w(k4)
      go to 305
  300 CALL pstlabrt(1,comm,1)
  305 return
      end


