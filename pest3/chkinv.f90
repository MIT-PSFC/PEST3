!
!......................................................................
      SUBROUTINE pstchkinv(msm)
!......................................................................
!      routine to verify the accuracy of the inversion. this is done by
!      iterating using the diff of the rhs and the solution * a for the
!      rhs. ie. rhs(new) = rhs(old) - a * x(old)
 USE pstcom

 USE combla

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER MSM
      INTEGER IS
      INTEGER MF
      INTEGER ML
      INTEGER M3
      INTEGER NLONG
      INTEGER IGIVUP
      INTEGER I
      INTEGER NADRES
      INTEGER IP
      INTEGER I1
      INTEGER I2
      INTEGER I2P
      INTEGER I1P
      INTEGER ICOLST
      INTEGER ICOLND
      INTEGER IROW
      INTEGER IVEC
      INTEGER ICOL
      INTEGER ICOL1


      double precision sum
!
 COMPLEX*16, DIMENSION(nad) 	 :: a1
 COMPLEX*16, DIMENSION(nfmg,nfmgc) 	 :: a2
 COMPLEX*16, DIMENSION(nmtg) 	 :: xsav
 REAL*8, DIMENSION(nsg,nsg) 	 :: dgampr
 REAL*8, DIMENSION(nsg,nsg) 	 :: ddelpr
      is = 0
      mf = jsub(1)
      ml = 2*mf
      m3 = 3*mf
      nlong = ml * ( ml+1 ) /2
      igivup = 1
      do 60 i = 1, mp
!      get a
      nadres = (i-1) * nlong + 1
      a1(1:nlong) = wpot(nadres:nadres+nlong-1)
!
!      form a2 and symmetrise
!
      ip = 0
      do i1 = 1, mf
         do i2 = i1, ml
            i2p = i2 + mf
            ip = ip + 1
            a2(i1,i2p) = a1(ip)
         end do
      end do
!
!      fill in lower triangle
!
      do i2 = 2, mf
         i2p = i2 + mf
         do i1 = 1, i2-1
            i1p = i1 + mf
            a2(i2,i1p) = conjg( a2(i1,i2p) )
         end do
      end do
!
!      form a * x - b
      icolst = 1
      if(i == 1)icolst = mf + 1
      icolnd = m3
      if(i == mp)icolnd = ml
!
      do irow = 1, mf
         sum =  0._r8 
         ivec = (i-2) * mf
         if(i == 1) ivec = 0
         do icol = icolst, icolnd
            ivec = ivec + 1
            sum = sum + a2(irow,icol) * ut(ivec)
         end do
         is = is + 1
         xsav(is) = ut(is)
         xt(is) = xt(is) - sum
      end do
!
!      now shift block for use again
      if (i == mp) go to 60

      do irow=1, mf
         do icol = 1, mf
            icol1 = ml + icol
            a2(icol,irow) = a2(irow,icol1)
         end do
      end do
!
   60 continue
!      now find the new x
      CALL pstcbxlu(1)
      CALL pstcdlhxv
!      write output
!      do 70 iw = 1, mp
!      iw1 = (iw-1) * mf +1
!      iw2 = iw1 + mf -1
!      write(outmod,1000)(ut(k),k=iw1,iw2)
!      write(outmod,1000)(xsav(k),k=iw1,iw2)
!   70 continue
!
      if(msm  /=  nosing) return
 1000 format(10(10e12.5))
      return
 7000 CALL psterrmes(outpst,'chkinv',1)
      end


