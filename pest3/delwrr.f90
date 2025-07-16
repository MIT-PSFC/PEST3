!         construct and store delta-w in the plasma.
!.......................................................................
      SUBROUTINE pstdelwrr
!
!      this SUBROUTINE pstcomputes the matrix elements of delta-w
!      between regular finite elements.
!
!      ray grimm, ppl, sept, 1978.
!.............................................................
!
 USE pstcom

 USE l22com

 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JLOW1
      INTEGER JLOW2
      INTEGER JUPP1
      INTEGER JUPP2
      INTEGER JMAX3
      INTEGER M11
      INTEGER JDIF1
      INTEGER JDIF2
      INTEGER J1
      INTEGER JM1T
      INTEGER J2
      INTEGER JM2T
      INTEGER KMAX
      INTEGER JM1
      INTEGER JM2
      INTEGER IS


! check allocation status

if ( .not. allocated(wreg) ) print *,'*** wreg not allocated'
if ( .not. allocated(wkin) ) print *,'*** wkin not allocated'
jlow1 =  lbound( wreg,1 )
jlow2 = lbound( wreg,2 )
jupp1 = ubound( wreg,1 )
jupp2 = ubound( wreg,2 )
if( jlow1 /= 1  .OR.  jlow2 /= 1 ) print *,'*** lower bound =',jlow1,jlow2,' for wreg'
if( jupp1 /= nkb  .OR.  jupp2 /= nkc ) print *,'*** upper  bound =',jupp1,jupp2,' for wreg'
jlow1 =  lbound( wkin,1 )
jlow2 = lbound( wkin,2 )
jupp1 = ubound( wkin,1 )
jupp2 = ubound( wkin,2 )
if( jlow1 /= 1  .OR.  jlow2 /= 1 ) print *,'*** lower bound =',jlow1,jlow2,' for wkin'
if( jupp1 /= nkb  .OR.  jupp2 /= nkc ) print *,'*** upper  bound =',jupp1,jupp2,' for wkin'


!
!           zero out matrices.
!
       wreg = ( 0.0_r8 , 0.0_r8 ) 
       wkin = ( 0.0_r8 , 0.0_r8 ) 
 !
!      run out to next singular surface...
!
!
      do 90 m1 = mbeg, mend
!
!        limit of fourier modes.
!
      jmax1  = lmax(1) - lmin(1) + 1
      jmax3 = ibas*jmax1
!
      llp = jmax1
         m11 = m1
         if ( m1  <  mp  )   m11 = m1 + 1
!
!      do integrals for matrix elements around node m1...
!
         call pststorrr
!
      if ( m1  ==  1 ) go to 14
      jdif1 =   jsub(mp)
      if ( m1  <  mp ) jdif1 = jtot(m1+1) - jtot(m1)
      jdif2 = jdif1
      if ( m1  <  mp ) jdif2 = jtot(m1+2) - jtot(m1)
         wreg(1:jdif1,1:jdif2) = ( 0.0_r8 , 0.0_r8 )  
         wkin(1:jdif1,1:jdif2) = ( 0.0_r8 , 0.0_r8 ) 
   14 continue
!
         do 80 m2 = m1,m11
!
              im1 = 2
              if ( m2  >   m1 )   im1 = 3
              im2 = 2
              if ( m2  >   m1 )   im2 = 1
!
      jtot1 = jtot(1)
      jtot2 = jtot(m2) - jtot(m1) + 1
!
!      run over all fourier modes...
!
      do j1 = 1, jmax3
         jm1t = jtot1 + ibas* (j1-1)
         do j2 = 1, jmax3
            jm2t = jtot2 + ibas* (j2-1)
            !
            kmax = 1
            !
            jm1 = jm1t
            jm2 = jm2t
            !
            if ( m1  ==  1  .AND.  lcirc )   go to 52
            if ( m1  ==  1   .AND.   fast )  go to 52
            !
            !----- 1 1 elt -----
            wreg(jm1,jm2) = rw1(im1,j1,j2) + ry2(im2,j1,j2) &
                 +conjg(ry2(im1,j2,j1)) + rz3(im1,j1,j2)
            wkin(jm1,jm2) = rk1(im1,j1,j2)

            go to 55
            !
52          continue
            !
            !      scrap element at origin for circular,cylindrical cases..
            !      also fill extra rows/columns for singular elements...
            !
            wreg(jm1,jm2) = ( 0.0_r8 , 0.0_r8 ) 
            if ( jm1  ==  jm2 )  wreg(jm1,jm2) = ( 1.0_r8 , 0.0_r8 ) 
            wkin(jm1,jm2) = ( 0.0_r8 , 0.0_r8 ) 
            if ( jm1  ==  jm2 )  wkin(jm1,jm2) = ( 1.0_r8 , 0.0_r8 ) 
            !
55          continue
            wreg(jm1,jm2) = twopi2*wreg(jm1,jm2)
            wkin(jm1,jm2) = twopi2*wkin(jm1,jm2)
            !
         end do
      end do
!
   80         continue
!
!        write out to disk...
!
      CALL pstblockt
!
!
   90 continue
!
      return
 7000 CALL psterrmes(outpst, 'delwrr', is)
!
      end


