!------------------------------------------------------
      subroutine pststorrb
!------------------------------------------------------
! quadrature of basis function times rhs...
!
! a. pletzer aug.90
!.......................................................................
!
 USE pstcom

 USE l22com

 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 THIRD
      real*8 TWOTH
      real*8 FOURTH
      INTEGER M1L
      INTEGER I
      INTEGER M11
      real*8 W
      INTEGER IS
      real*8 DPS
      INTEGER ILEFT
      INTEGER IRIGH
      real*8 WL
      real*8 WR
      INTEGER J1
      real*8 pstUF
      real*8 pstUFPR


!
!
!        some constants....
!
 COMPLEX*16	 	 	 :: phis1
 COMPLEX*16	 	 	 :: phis2
 INTEGER	 	 	 :: surfl
 INTEGER	 	 	 :: surfc
 INTEGER	 	 	 :: surfr
 INTEGER	 	 	 :: surfcp
      third = one / three
      twoth = two * third
      fourth = four * third
!
! fix surface number limits...
!
      surfl = 1
      m1l = m1-2

      if (m1l >= 1) then
         do i= 1, m1l
            surfl = surfl + msub(i)
         end do
      end if
!
      surfc = 1
      if ( m1  >   1 ) surfc = surfl + msub(m1-1)
!
! aplet  13.7.98      surfr = surfc + msub(m1)
! aplet  13.7.98       if ( m1  ==  1 )  surfr = surfl + msub(1)
      surfr = mp
      if ( m1 < mp ) surfr = surfc + msub(m1)
!
      scl  = surfl
      scr  = surfr
! aplet  13.7.98      if ( m1  ==  mp )   scr = surfc
      m11  = m1
!
      w0l1ou(m1,1:jmax1,ms) =  0.0_r8 
      w0l1eu(m1,1:jmax1,ms) =  0.0_r8 
!
!     calculation of the psi integrals.
!
      do    38 surf = scl, scr
         !
         !        quadrature weight
         !        for simpson"s rule pquad =  1./ 3.
         !
         w  = fourth
         !aplet  7.5.93         w  =  2._r8 *( 1._r8  - pquad)
         !aplet  7.5.93
         is  = surf - scl
         if ( (is/2)*2  ==  is )   w = twoth
         !aplet  7.5.93         if ( (is/2)*2  ==  is )   w =  2._r8 * pquad
         !aplet  7.5.93
         if ( surf  ==  scl   .OR.   surf  ==  scr )  w = third
         !aplet  7.5.93    if ( surf  ==  scl   .OR.   surf  ==  scr )  w = pquad
         !aplet  7.5.93
         dps = psinew(surf+1) - psinew(surf)
         if ( m1  ==  mp  .AND.  surf  ==  scr ) &
              dps = psinew(surf) - psinew(surf-1)
         if ( surf  ==  surfr ) dps = psinew(surfr)-psinew(surfr-1)
         w = w * dps
         psiv = psinew(surf)
         !
         if( m1 >  2 ) bit= 1.0E-12_r8  
         !aplet  20.12.95      if( .not.(lcub) .OR. (m1 == 1) ) then
         if( surf == scl ) psiv = psiv + bit
         if( surf == scr ) psiv = psiv - bit
         !plaet  20.12.95      end if
         !
         !
         if ( surf == surfc ) then      
            ileft = max( surf-1, 1)
            irigh = min( surf+1, 2*(mp-1)+1)
            wl = ( psinew(surf ) - psinew(ileft) )/ 3.0_r8 
            wr = ( psinew(irigh) - psinew(surf ) )/ 3.0_r8 
            w = wl + wr
         end if
         !
         do 37 j1 = 1,jmax1
            !
            phis1 = - pstuf(m11)*alxi1o(surf,j1) &
                 + pstuf(m11)*qdchio(surf,j1) &
                 + pstufpr(m11)*ddchio(surf,j1)
            phis2 = - pstuf(m11)*alxi1e(surf,j1) &
                 + pstuf(m11)*qdchie(surf,j1) &
                 + pstufpr(m11)*ddchie(surf,j1)
            w0l1ou(m1,j1,ms) = w0l1ou(m1,j1,ms) + w*phis1
            w0l1eu(m1,j1,ms) = w0l1eu(m1,j1,ms) + w*phis2
            !
37          continue
            !   ------ end of l loop -----
            !
38          continue
            !   ------ end of surf loop ------
!
      return
      end


