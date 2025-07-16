      subroutine pstsearch(xbar,x,n,i,mflag)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 XBAR
      INTEGER N
      INTEGER I
      INTEGER MFLAG
      INTEGER IXBAR
      INTEGER IX1
      INTEGER IXN
      INTEGER K
      INTEGER J
      INTEGER L
      real*8 A
      real*8 B

!
!      dec-10 version    rg  1979._r8 
 CHARACTER*64            :: com1
 REAL*8, DIMENSION(*) 	 :: x
      data com1/'pstsearch xbar s outside range of table    '/
      mflag=0
      i=n
      if(xbar == x(n))return
      i=1
      if(n <= 1) return
      ixbar=sign(1.0_r8,xbar)
      ix1=sign(1.0_r8,x(1))
      ixn=sign(1.0_r8,x(n))
      do 5 k=1,n
      j=i+i
      if(j >= n) go to 6
  5   i=j
 6    k=i
      mflag = 1
      do 115  l=n,2,-1
      a=x(l-1)
      b=x(l)
      if(sign(1.0_r8,a)-sign(1.0_r8,b)) 7,113,8
 113   if(a-b)7,115,8
 115    continue
  7    j=1
      if(ixbar < ix1 .OR. (ixbar == ix1 .AND. xbar < x(1)) .OR. ixbar &
.gt.ixn.or.(ixbar.eq.ixn.and.xbar.gt.x(n))) go to 16
      go to 10
  8    j=2
      if(ixbar < ixn .OR. (ixbar == ixn .AND. xbar < x(n)) .OR. ixbar &
.gt.ix1.or.(ixbar.eq.ix1.and.xbar.gt.x(1))) go to 16
   10 k=k/2
       a=x(i)
      go to (11,20),j
  11   if(ixbar-sign(1.0_r8,a)) 111,1111,2111
 1111 if(xbar-a)111,14,2111
 2111 b=x(i+1)
      if(ixbar-sign(1.0_r8,b)) 2112,2113,12
 2113   if(xbar >= b) go to 12
 2112   return
 111  i = i-k
      go to 13
   12 i = i+k
   13 if (i < n) go to 10
      k=k/2
      go to 111
   14 mflag=0
      return
   16 CALL pstlabrt(1,com1,1)
      mflag=2
      return
  20   if(ixbar-sign(1.0_r8,a) ) 2120,2121,111
 2121  if(xbar-a) 2120,14,111
 2120 b=x(i+1)
       if(ixbar-sign(1.0_r8,b)) 12,2122,2112
 2122  if(xbar-b) 12,12,2112
 8888    return
      end


