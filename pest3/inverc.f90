!-----------------------------------------------------------------------
! aplet  19.4_r8 .94      SUBROUTINE pstinverc ( a,ldim,l )
      SUBROUTINE pstinverc ( a,ldim,l , ierr)
!-----------------------------------------------------------------------
!     lcm(invert)
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER LDIM
      INTEGER L
      INTEGER IERR
      real*8 DET
      INTEGER I
      INTEGER IDET
      INTEGER J
      INTEGER J1
      INTEGER J2



 COMPLEX*16, DIMENSION(ldim,ldim) 	 :: a
 COMPLEX*16, DIMENSION(ldim) 	 :: b
 COMPLEX*16, DIMENSION(ldim,ldim) 	 :: ao
 COMPLEX*16, DIMENSION(ldim,ldim) 	 :: as
 COMPLEX*16	 	 	 :: s
 INTEGER, DIMENSION(ldim) 	 :: ip
!
      ao(1:L,1:L) = a(1:L,1:L)
!
      CALL pstdecomc(l,ldim,a,ip )
!ccccccccccccccc compute determinant  rg 12/03/79
      det =  0.0_r8 
      do i = 1, l
         det = det +log10( real( abs( a(i,i) ) ) )
      end do
      idet = INT( det )
      det = 10 ** ( det-idet )
      s = ( 1._r8 , 0._r8 ) 
      do i = 1 , l
         s = s * a(i,i)/ abs(a(i,i))
      end do
      s = s * ip(l)
      if ( real(s) <=   0.0_r8  )  go to 99
      if ( ip(l)  ==  0 )  go to 99
      do 10 i = 1, l
         do 5 j = 1, l
            b(j) = ( 0._r8 , 0._r8 ) 
            if ( i ==  j )  b(j) = ( 1._r8 , 0._r8 ) 
5        continue
         call pstsolvec( l,ldim,a,b,ip  )
         do j = 1, l
           as(j,i) = b(j)
         end do
10   continue
   a(1:L,1:L) = as(1:L,1:L)
!
      return
!cccccccccccccc  error exit  rg 01/12/79
   99 continue
      ierr = 1
!
      write(itty,8001) s,det,idet
      write(outmod,8001) s,det,idet
!
! aplet  19.4_r8 . 94._r8 ..
!
      write(itty, 8002) 
      write(outmod,8002)
      do j1 = 1, l
      write(itty, 8003) (real(ao(j1,j2)), j2=1,l)
      write(outmod,8003) (real(ao(j1,j2)), j2=1,l)
      end do
      write(itty, 8004)
      write(outmod,8004)
      do j1 = 1, l
      write(itty, 8003) (aimag(ao(j1,j2)), j2=1,l)
      write(outmod,8003) (aimag(ao(j1,j2)), j2=1,l)
      end do
!
!
 8001 format(1x," determinant: s = ",2f4.1," det = ",f8.5, &
  " idet = ",i5," in inverc le zero ")
 8002 format(1x," re a = ")
 8003 format(9(1x,e15.8))
 8004 format(1x," im a = ")
! aplet  19.4_r8 .94      stop
      return
!ccccccccccccccc
      end



