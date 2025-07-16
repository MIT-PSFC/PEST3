!----------------------------------------------------------------------
      subroutine pststorrr
!
!.......................................................................
!                  this subroutine evaluates and stores values of
!                  certain integrals over psi and theta   needed
!                  in constructing the delta-w and kinetic energy matrix
!                  elements.
!                  since these elements are tridiagonal in the m"s
!                  (finite element number), the diagonal, super-diagonal
!                  and sub-diagonal blocks in m are calculated only.
!                  fourier transforms are included up to order
!                  max ( lmax(1)-lmin(1), lmax(1)+lmin(1) ).
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
      INTEGER IPPS
      INTEGER M11
      INTEGER M12
      real*8 W
      INTEGER IS
      real*8 DPS
      real*8 WL
      real*8 PHIS1L
      real*8 pstUF
      real*8 PHIS2L
      real*8 pstUFPR
      real*8 PHIS3L
      real*8 WR
      real*8 PHIS1R
      real*8 PHIS2R
      real*8 PHIS3R
      real*8 PHIS1
      real*8 PHIS2
      real*8 PHIS3


!
!        some constants....
 INTEGER	 	 	 :: surfl
 INTEGER	 	 	 :: surfc
 INTEGER	 	 	 :: surfr
 INTEGER	 	 	 :: surfcp
      third = one / three
      twoth = two * third
      fourth = four * third
!      4.1_r8 .2.1   fix surface number limits, and unless doing first row
!               shift (old) t rights into (new) t lefts.
!
      surfl = 1
      m1l = m1-2
      if ( m1l  <  1 )  go to 6
      do 5 i = 1, m1l
      surfl = surfl + msub(i)
    5 continue
    6 continue
      surfc = 1
      if ( m1  >   1 ) surfc = surfl + msub(m1-1)
      surfr = surfc + msub(m1)
!         fix for first element...
      if ( m1  ==  1 )  surfl = 1
      if ( m1  >   2 )   go to    10
!         first element is only half tent...
      if ( m1  >   1 )  go to    9
!        special for m1=1
      scl= 1
      scr = scl + msub(1)
      go to    20
    9 continue
!         special for m1=2
      scl = msub(1) + 2
      scr = scl + msub(2) - 1
      go to    20
   10 continue
      call pstshifft
!        last element is only half tent.
        if ( m1  ==  mp )  go to    22
      scl = msub(m1-1) + 2
      scr = scl + msub(m1) - 1
!
!      4.1_r8 .2.2   calculate and store all fft's needed on surfaces to the
!               right of the element center.(except for last element)
   20 continue
      surfcp = surfc + 1
!      if (( m1 == mend)  .AND.  (m1 /= mp))  scr=scr+1
      CALL pstzerot ( scl,scr )
      do 21 ipps = scl, scr
      scount = ipps
      surf = surfl + scount - 1
!!$  26.2_r8 .97     call thfft2
!call thefft
call pstthfft2
!
   21 continue
!
   22 continue
!               zero off elements first.
!
                CALL pstzeror
!
      do    39 i = 1,3
         if ( m1  ==  mp  .AND.   i  == 1 )  go to    39
         if ( i  >=   2 )   go to    31
!
!................................ lower band ...........................
!
         scl  = surfc
         scr  = surfr
         m11  = m1 + 1
         m12  = m1
         go to    33
   31    if ( i  >   2 )   go to    32
!................................ diagonal .............................
         scl  = surfl
         scr  = surfr
         if ( m1  ==  mp )   scr = surfc
         m11  = m1
         m12  = m1
         go to    33
   32    continue
!................................ upper band ...........................
         if ( m1  ==  mp  .AND.  i  ==  3 )   go to    39
         scl  = surfc
         scr  = surfr
         m11  = m1
         m12  = m1 + 1
   33    continue
!
!     calculation of the psi integrals.
!
      do    38 surf = scl, scr
!
!        quadrature weight 
!        for simpson"s rule pquad =  1._r8 / 3._r8 
!
!         w  = fourth
         w  =  2._r8 *( 1._r8  - pquad)
!aplet  7.5_r8 .93
         is  = surf - scl
!         if ( (is/2)*2  ==  is )   w = twoth
         if ( (is/2)*2  ==  is )   w =  2._r8 * pquad
!aplet  7.5_r8 .93
!         if ( surf  ==  scl   .OR.   surf == scr )  w = third
          if ( surf  ==  scl   .OR.   surf == scr )  w = pquad
!aplet  7.5_r8 .93
      dps = psinew(surf+1) - psinew(surf)
      if ( m1  ==  mp  .AND.  surf  ==  scr ) &
                  dps = psinew(surf) - psinew(surf-1)
      if ( surf  ==  surfr ) dps = psinew(surfr)-psinew(surfr-1)
      w = w * dps
      psiv = psinew(surf)
!++++++++++++++
      if ( m1  >   2 ) bit =  1.0E-12_r8  
!aplet  20.12_r8 .95      if ( .not.(m1 == 1)  .AND.  lcub  )  go to 8000
         if ( surf  ==  scl )   psiv = psiv + bit
         if ( surf  ==  scr )   psiv = psiv - bit
!aplet  20.12_r8 .95  8000 continue
!++++++++++++++
         scount  = surf - surfl + 1
!
!        make adjustments for step functions.
      if ( m1  ==  1 )  go to    34
      if ( m1  ==  mp )  go to    34
      if ( i  /=  2 )  go to    34
!
!        diagonal element-- check to see if approaching the middle of
!        the two finite element cells...if so add 1/2 contribution from
!        each side.
!
      if ( surf  /=  surfc )  go to    34
!
      psiv = psiv - bit
      wl = ( psinew(surfc)-psinew(surfc-1) ) /  3.0_r8 
      phis1l =  pstuf(m11) * pstuf(m12)
      phis2l =  pstuf(m11) * pstufpr(m12)
      phis3l =  pstufpr(m11) * pstufpr(m12)
      psiv = psiv + two*bit
      wr = ( psinew(surfc+1)-psinew(surfc) ) /  3.0_r8 
      phis1r =  pstuf(m11) * pstuf(m12)
      phis2r =  pstuf(m11) * pstufpr(m12)
      phis3r =  pstufpr(m11) * pstufpr(m12)
      phis1 = wl*phis1l+wr*phis1r
      phis2 = wl*phis2l+wr*phis2r
      phis3 = wl*phis3l+wr*phis3r
!
      go to    35
   34 continue
      phis1 = w* pstuf(m11) * pstuf(m12)
      phis2 = w* pstuf(m11) * pstufpr(m12)
      phis3 = w* pstufpr(m11) * pstufpr(m12)
!
   35 continue
!
!      run over all l-values...
!
      rw1(i,1:jmax1,1:jmax1) = rw1(i,1:jmax1,1:jmax1) + tw(scount,1:jmax1,1:jmax1) * phis1
      ry2(i,1:jmax1,1:jmax1) = ry2(i,1:jmax1,1:jmax1) + ty(scount,1:jmax1,1:jmax1) * phis2
      rz3(i,1:jmax1,1:jmax1) = rz3(i,1:jmax1,1:jmax1) + tz(scount,1:jmax1,1:jmax1) * phis3
!
      rk1(i,1:jmax1,1:jmax1) = rk1(i,1:jmax1,1:jmax1) + td(scount,1:jmax1,1:jmax1) * phis1
!
   38 continue
!   ------ end of surf loop ------
   39 continue
!   ----- end of i loop -----
!
      return
      end


