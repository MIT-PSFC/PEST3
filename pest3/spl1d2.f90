      subroutine pstspl1d2(n,x,f,w,ij,y,tab)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER N
      INTEGER IJ
      real*8 Y
      real*8 WZ
      real*8 SZ
      INTEGER I
      INTEGER MFLAG
      INTEGER MI
      INTEGER K1
      real*8 FLK
      real*8 A
      real*8 B
      real*8 C

      data wz,sz/ 2.0e0_r8 , 6.0e0_r8 / 
!     where n= number of points in the interpolation
!           x= origin of table of the independent variable
!           f= origin of table of the dependent variable
!           w= origin of table of second derivatives as calculated by
!              spl1d1
!           ij= spacing in the tables f and w
!           y= the point at which interpolation is desired
!           tab= an array of dimension 3 which contains the function
!                value, first derivative, and second derivative at y
!
!
!     locate y in the x table
!
 REAL*8, DIMENSION(3) 	 :: x
 REAL*8, DIMENSION(3) 	 :: f
 REAL*8, DIMENSION(3) 	 :: w
 REAL*8, DIMENSION(3) 	 :: tab
      if(y-x(1))10,10,20
   10 i=1
      go to 30
   20 if(y-x(n))15,40,40
   40 i=n-1
      go to 30
   15 call pstsearch(y,x,n,i,mflag)
   30 mi=(i-1)*ij+1
      k1=mi+ij
      flk=x(i+1)-x(i)
!
!     calculate f(y)
!
      a=(w(mi)*(x(i+1)-y)**3+w(k1)*(y-x(i))**3)/(sz*flk)
      b=(f(k1)/flk-w(k1)*flk/sz)*(y-x(i))
      c=(f(mi)/flk-flk*w(mi)/sz)*(x(i+1)-y)
      tab(1)=a+b+c
!
!     calculate the first derivative at y
!
      a=(w(k1)*(y-x(i))**2-w(mi)*(x(i+1)-y)**2)/(wz*flk)
      b=(f(k1)-f(mi))/flk
      c=flk*(w(mi)-w(k1))/sz
      tab(2)=a+b+c
!
!     calculate the second derivative at y
!
      tab(3)=(w(mi)*(x(i+1)-y)+w(k1)*(y-x(i)))/flk
      return
      end


