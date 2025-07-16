! 
!     purpose      calculate and store the theta integrals on psi
!                  surface number surf. in temporary storage scount
!                  points to the appropriate surface lying between the
!                  centre and r.h.s of element support.
!.......................................................................
!      ray grimm, jan, 1980._r8  j manickam feb. 1982 a pletzer 1993
!
! use tarnsformed operator form (aplet 18/8/95)
!.......................................................................
!
      subroutine pstthfft2
!.......................................................................
!
 USE pstcom

 USE l21com

 USE l22com

 USE newmet

 USE r33com

 USE comggf

 USE comfft
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER ISET
      INTEGER MFFT
      INTEGER MREM
      INTEGER IFR
      real*8 F
      real*8 G
      real*8 G2
      real*8 GGP
      real*8 Q
      real*8 QSQ
      real*8 QN
      real*8 QP
      real*8 PS
      real*8 PP
      real*8 SIG
      real*8 SIGP
      INTEGER LM1
      INTEGER LM2
      INTEGER JINC
      INTEGER J1
      INTEGER J2
      INTEGER L121
      INTEGER L1
      INTEGER L2
      INTEGER L12
      INTEGER LC1
      INTEGER JJ1
      INTEGER JJ2
      INTEGER IERR
      real*8 ALNQ1
      real*8 ALNQ2
      INTEGER KK1
      INTEGER JCOL
      INTEGER JTH
      INTEGER JL
      INTEGER I1
      INTEGER I2
      INTEGER I3
      INTEGER I4
      real*8 Z1
      real*8 Z2
      real*8 Z3
      real*8 Z4
      INTEGER LGIVUP



!
      save lfft
!...................
!
!  
      data iset /0/
!
!        for each t construct an array of theta values on the magnetic
!        surface, then integrate over theta with appropriate sin or cos.
!
!      set up for fft's on this surface...
!
 COMPLEX*16, DIMENSION(nths) 	 :: work
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: ttw
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: tty
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: ttz
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: ttd
 COMPLEX*16	 	 	 :: sum
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: workc
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: workg
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zidentity
 INTEGER, DIMENSION(nfn) :: IPIV
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: workf
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zt1
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zt2
 COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zt4
 REAL*8, DIMENSION(nths2) 	 :: zf
 REAL*8, DIMENSION(nfn,nfn) 	 :: zr
 REAL*8, DIMENSION(nfn,nfn) 	 :: zi
 INTEGER	 	 	 :: t
 INTEGER, DIMENSION(3) 	 :: lfft
!
!
      iset = iset + 1
      if ( iset  >   1  .AND.  surf /= 1 )  go to 5
      mfft = 1
      mrem = mth
    2 continue
      mrem = mrem / 2
      if ( mrem  ==  1 )  go to 3
      mfft = mfft + 1
      go to 2
    3 continue
      if ( jmax2  <=  mth/2-1 )  go to 4
!      error in use of fft...
      write(outmod,9000)
 9000 format( 1x, " ****error in calling fft...too large llp ****")
      stop 'ERROR in thfft: too many Fourier modes'
!      ......
!
    4 continue
      lfft(1) = mfft - 1
      lfft(2) = 0
      lfft(3) = 0
      CALL pstff2prp( lfft,invfft,sfft,omfft,onfft,ifr )
!
    5 continue
!
      f = fa(surf)
      g = r*ga(surf)
      gp = r*gpa(surf)
      g2 = g*g
      ggp = g * gp
      q = qa(surf)
      qsq = q*q
      qn = n * q
      qp = qpa(surf)
      ps = pa(surf)
      pp = ppa(surf)
      sig = sigavr(surf)
      sigp = sigavp(surf)

      lm1 = lmin(1) -1
      lm2 = - (jmax2-1)/2 - 1
      if ( jmax2  ==  1 )  lm2 = lm1
      jinc = lm1 - lm2

!
! common matrices...
!
  zf(1:mth1) = xjacob(1:mth1,surf)/xsq(1:mth1,surf)
!
       CALL pstloadw(mth1, zf, work)
       CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
       do j1 = 1,jmax2
       zt1(j1,j1) = real(work(1))
       do j2 = 1,j1-1
       l121 = abs(j1-j2) + 1
       zt1(j1,j2) = work(l121)
       zt1(j2,j1) = conjg( zt1(j1,j2) )
       end do
       end do
!
 zf(1:mth1) = xsq(1:mth1,surf) / (xjacob(1:mth1,surf)*grpssq(1:mth1,surf))
!
       CALL pstloadw(mth1, zf, work)
       CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
       do j1 = 1,jmax2
       zt2(j1,j1) = real( work(1) )
       do j2 = 1,j1-1
       l121 = abs(j1-j2) + 1
       zt2(j1,j2) = work( l121 )
       zt2(j2,j1) = conjg( zt2(j1,j2) )
       end do
       end do
!
 zf(1:mth1) =  1._r8  / grpssq(1:mth1,surf)
!
       CALL pstloadw(mth1, zf, work)
       CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
       do j1 = 1,jmax2
       zt4(j1,j1) = real( work(1) )
       do j2 = 1,j1-1
       l121 = abs(j1-j2) + 1
       zt4(j1,j2) = work( l121 )
       zt4(j2,j1) = conjg( zt4(j1,j2) )
       end do
       end do
!
!      start with gee-matrix..
!
       workg(1:jmax2,1:jmax2) =  0.0_r8 
!
! e1 + e2
!
 zf(1:mth1) = xjacob(1:mth1,surf)/xsq(1:mth1,surf) &
    + qsq * xsq(1:mth1,surf)* &
 (f * xjacob(1:mth1,surf)/xsq(1:mth1,surf) -  1.0_r8   )**2 &
 /(xjacob(1:mth1,surf)*grpssq(1:mth1,surf))
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         workg(j1,j1) = ensq * real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            l12 = abs(l1-l2)
            lc1 = l12 + 1
            workg(j1,j2) = ensq * work(lc1)
            workg(j2,j1) = conjg( workg(j1,j2) )
         end do
      end do
!
! f
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         do j2 = 1, jmax2
            l2 = j2 + lm2
            workg(j1,j2) = workg(j1,j2) + l1 * l2 * zt2(j1,j2)
         end do
      end do
!
! v
!
 zf(1:mth1) = (f-xsq(1:mth1,surf)/xjacob(1:mth1,surf))/grpssq(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         workg(j1,j1) = workg(j1,j1) + qn*(l1+l1)*real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            lc1 = abs(l1-l2) + 1
            workg(j1,j2) = workg(j1,j2) + qn*(l1+l2)*work(lc1)
            workg(j2,j1) = conjg(workg(j1,j2))
         end do
      end do
!
      do j1 = 1, jmax1
         jj1 = j1 + jinc
         do j2 = 1, jmax1
            jj2 = j2 + jinc
            gm1(j1,j2) = workg(jj1,jj2)
         end do
      end do
!
      ierr = 0
!!$      CALL pstinverc ( workg,nfn,jmax2, ierr)
 ! invert workg 

     zidentity = (0.0_r8, 0.0_r8)
     do j1=1, nfn
        zidentity(j1, j1) = (1.0_r8, 0.0_r8)
     enddo
     call ZGESV( nfn, nfn, workg, nfn, IPIV, zidentity, nfn, ierr)
     workg = zidentity ! load solution
     if(ierr/=0) print*,'ERROR ',ierr,' occurred while attempting to invert workg'
!
      if( ierr /=  0 ) go to 99
!
!aplet  20.5_r8 .94 check inversion accuracy...
!     if(check1) then
!     zerr =  0.0_r8 
!     do j1 = 1,jmax1
!     jj1 = j1 + jinc
!     do j2 = 1,jmax1
!     jj2 = j2 + jinc
!     zr(j1,j2) =  0.0_r8 
!     zi(j1,j2) =  0.0_r8 
!
!     do j3 = 1, jmax1
!     jj3 = j3 + jinc
!     zr(j1,j2) = zr(j1,j2) + real(gm1(j1,j3))* real(workg(jj3,jj2))
!    .                      -aimag(gm1(j1,j3))*aimag(workg(jj3,jj2))
!     zi(j1,j2) = zi(j1,j2) + real(gm1(j1,j3))*aimag(workg(jj3,jj2))
!    .                      +aimag(gm1(j1,j3))* real(workg(jj3,jj2))
!     end do
!     if( j1 == j2 ) zr(j1,j2) = zr(j1,j2) -  1.0_r8 
!     zerr = zerr + sqrt( zr(j1,j2)**2 + zi(j1,j2)**2 )
!     end do
!     end do
!
!     if( zerr >   1.E-6_r8  ) then
!     write(outmod,*)'***thfft2: inaccurate inversion of g^-1'
!    .              ,' et m1 = ', m1, ' error = ', zerr
!     write(itty  ,*)'***thfft2: inaccurate inversion of g^-1'
!    .              ,' et m1 = ', m1, ' error = ', zerr
!     end if
!     end if
!caplet  20.5_r8 .94
!
! aplet  19.4_r8 .94
!
      geem(1:jmax1,1:jmax1) = workg(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
!
! start with w terms not involving q...
!
      ttw(1:jmax2,1:jmax2) =  0.0_r8 
!
! j
!
 zf(1:mth1) =  1._r8  / ( xjacob(1:mth1,surf)*grpssq(1:mth1,surf) )
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         ttw(j1,j1) = alnq1*alnq1*real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            alnq2 = l2 - qn
            l12 = abs(l1-l2)
            lc1 = l12 + 1
            ttw(j1,j2) = alnq1*alnq2 * work(lc1)
            ttw(j2,j1) = conjg( ttw(j1,j2) )
         end do
      end do
!
! m2 + m4
!
 zf(1:mth1) = ( xsq(1:mth1,surf) * pp + ggp ) *                  &
  ( (xjprym(1:mth1,surf)-xjacob(1:mth1,surf)*xsqdps(1:mth1,surf) &
        /xsq(1:mth1,surf)) /xsq(1:mth1,surf)                     &
       - xjacob(1:mth1,surf)*                                    &
 ( xsq(1:mth1,surf) * pp + ggp )                                 &
 /(xsq(1:mth1,surf)*grpssq(1:mth1,surf))   )
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         ttw(j1,j1) = ttw(j1,j1) + real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            lc1 = abs(l1-l2) + 1
            ttw(j1,j2) = ttw(j1,j2) + work(lc1)
            ttw(j2,j1) = ttw(j2,j1) + conjg( work(lc1) )  
         end do
      end do
      
!
! m1
!
 zf(1:mth1) = pp * xjacob(1:mth1,surf)*xsqdps(1:mth1,surf) /xsq(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         ttw(j1,j1) = ttw(j1,j1) + real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            lc1 = abs(l1-l2) + 1
            ttw(j1,j2) = ttw(j1,j2) + work(lc1)
            ttw(j2,j1) = ttw(j2,j1) + conjg( work(lc1) )
         end do
      end do
!
! m3
!
 zf(1:mth1) =  &
       ( xsq(1:mth1,surf)*pp + ggp )*grpsth(1:mth1,surf)*xjacob(1:mth1,surf) &
       / ( xsq(1:mth1,surf)*grpssq(1:mth1,surf) )
!
      CALL pstloadw(mth1, zf, work)     
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         do j2 = 1,j1-1
            l2 = j2 + lm2
            l12 = abs( l1-l2 )
            lc1 = l12 + 1
            ttw(j1,j2) = ttw(j1,j2) + ( 0._r8 , 1._r8 ) *(l1-l2)* work(lc1)
            ttw(j2,j1) = ttw(j2,j1) + ( 0._r8 , 1._r8 ) *(l2-l1)*  &
                 conjg( work(lc1) )
         end do
      end do
!
   45 continue
!
! start to work on q...
!
!      initialize ttq...
!
   ttq(1:jmax2,1:jmax2) = ( 0.0_r8 , 0.0_r8 ) 
!
! a1
!
 zf(1:mth1) = grpsth(1:mth1,surf) / grpssq(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         ttq(j1,j1) = ( 0._r8 , 1._r8 ) * alnq1**2 * real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            alnq2 = l2 - qn
            l12 = l1-l2
            lc1 = abs(l12) + 1
            ttq(j1,j2) = ( 0._r8 , 1._r8 ) * alnq1* alnq2 * work(lc1)
            ttq(j2,j1) = ( 0._r8 , 1._r8 ) * alnq1* alnq2 *conjg( work(lc1) )
         end do
      end do
!
! d1
!
 zf(1:mth1) = xjacob(1:mth1,surf)*grpsth(1:mth1,surf)/ &
  (xsq(1:mth1,surf)*grpssq(1:mth1,surf))
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         ttq(j1,j1) = ttq(j1,j1) + ( 0._r8 , 1._r8 ) *n*g*alnq1*real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            alnq2 = l2 - qn
            l12 = l1 - l2
            lc1 = abs(l12) + 1
            ttq(j1,j2) = ttq(j1,j2) + ( 0._r8 , 1._r8 ) *n*g*work(lc1)*alnq2
            ttq(j2,j1) = ttq(j2,j1) + ( 0._r8 , 1._r8 ) *n*g*alnq1* &
                 conjg( work(lc1) )
         end do
      end do
!
! d2
!
  zf(1:mth1) =  qdelp(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do    j1 = 1, jmax2
         l1 = j1 + lm2  
         alnq1 = l1 - qn
         ttq(j1,j1) = ttq(j1,j1) + ( 0._r8 , 1._r8 ) *n*alnq1*real( work(1) )
         do    j2 = 1,j1-1
            l2 = j2 + lm2
            alnq2 = l2 - qn
            l12 = l1 - l2
            lc1 = abs(l1-l2) + 1
            ttq(j1,j2) = ttq(j1,j2) + ( 0._r8 , 1._r8 ) *n*work(lc1)*alnq2
            ttq(j2,j1) = ttq(j2,j1) + ( 0._r8 , 1._r8 ) *n*alnq1* &
                 conjg( work(lc1) )
         end do
      end do
!
! b1
!
 zf(1:mth1) = ( xsq(1:mth1,surf)*pp + ggp ) / grpssq(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         ttq(j1,j1) = ttq(j1,j1) + alnq1*real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            alnq2 = l2 - qn
            lc1 = abs(l1-l2) + 1
            ttq(j1,j2) = ttq(j1,j2) + alnq1*work(lc1)
            ttq(j2,j1) = ttq(j2,j1) + conjg( work(lc1) )*alnq2
         end do
      end do
!
! c1
!
 zf(1:mth1) = xjacob(1:mth1,surf)*(xsq(1:mth1,surf) *pp + ggp)/ &
  (grpssq(1:mth1,surf)*xsq(1:mth1,surf))
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         ttq(j1,j1) = ttq(j1,j1) + n*g*real( work(1) )
         do j2 = 1,j1-1
            lc1 = abs(j1-j2) + 1
            ttq(j1,j2) = ttq(j1,j2) + n*g*work(lc1)
            ttq(j2,j1) = ttq(j2,j1) + n*g*conjg( work(lc1) ) 
         end do
      end do
!
!
! c2
!
      ttq(1:jmax2,1:jmax2) = ttq(1:jmax2,1:jmax2) + &
         n*gp*zt1(1:jmax2,1:jmax2) 
!
      do 70 j1 = 1, jmax2
      ttq(j1,j1) = ttq(j1,j1) - n*qp
   70 continue
!
      if (n /=  0._r8   .AND.  ltransf ) then
         !-------
         !
         ! corrections to transformed operators...
         !
         ! a2
         !
         do j1 = 1, jmax2
            l1 = j1 + lm2
            alnq1 = l1 - qn
            do j2 = 1, jmax2
               l2 = j2 + lm2
               alnq2 = l2 - qn
               ttq(j1,j2) = ttq(j1,j2) + alnq1* alnq2 * &
                    sig*zt2(j1,j2)/n
            end do
         end do
         !
         ! d3
         !
         do j1 = 1, jmax2
            do j2 = 1, jmax2
               l2 = j2 + lm2
               alnq2 = l2 - qn
               ttq(j1,j2) = ttq(j1,j2) + g*sig*zt4(j1,j2)*alnq2
            end do
         end do
         !
         ! b1 bis
         !
         do j1 = 1, jmax2
            l1 = j1 + lm2
            alnq1 = l1 - qn
            do j2 = 1, jmax2
               ttq(j1,j2) = ttq(j1,j2) + alnq1*g*sig*zt4(j1,j2)
            end do
         end do
         !
         ! c2 bis
         !
         do j1 = 1,jmax2
            do j2 = 1,jmax2
               ttq(j1,j2) = ttq(j1,j2) + n*sig*zt1(j1,j2)
            end do
         end do
         !
         ! c3
         !
         zf(1:mth1) = xjacob(1:mth1,surf)/(grpssq(1:mth1,surf)*xsq(1:mth1,surf))
         !
         CALL pstloadw(mth1, zf, work)
         CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
         !
         do j1 = 1, jmax2
            ttq(j1,j1) = ttq(j1,j1) + n*sig*g2*real( work(1) )
            do j2 = 1,j1-1
               lc1 = abs(j1-j2) + 1
               ttq(j1,j2) = ttq(j1,j2) + n*sig*g2*work(lc1)
               ttq(j2,j1) = ttq(j2,j1) + n*sig*g2*conjg( work(lc1) )
            end do
         end do
         !
         ! corrections to m matrix
         !
         do j1 = 1,jmax1
            jj1 = j1 + jinc
            l1 = jj1 + lm2
            alnq1 = l1 - qn
            ttw(jj1,jj1) = ttw(jj1,jj1) - sig*qp +  &
                 alnq1*sigp/n
            do j2 = 1,jmax1
               jj2 = j2 + jinc
               ttw(jj1,jj2) = ttw(jj1,jj2) + sig*sig*gm1(j1,j2)/ensq &
                    - sig*( conjg(ttq(jj2,jj1)) + ttq(jj1,jj2) )/n
            end do
         end do
         !
         !-----------
      end if
!
   73 continue
!
!      finish off w-matrix...
!
    workc(1:jmax2,1:jmax2) =  &
         matmul( workg(1:jmax2,1:jmax2), ttq(1:jmax2,1:jmax2) )
    workf(1:jmax2,1:jmax2) =  &
 matmul( conjg(transpose(ttq(1:jmax2,1:jmax2))), workc(1:jmax2,1:jmax2) )
!
         do j1 = 1, jmax2
            do j2 = 1, jmax2
               ttr(j1,j2) = -ttw(j1,j2)
               ttw(j1,j2) = ttw(j1,j2) + workf(j1,j2)
            end do
         end do
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         do j2 = 1, jmax2
            sum = ( 0.0_r8 , 0.0_r8 ) 
            do kk1 = 1, jmax2
               sum = sum + workg(j1,kk1) * ttq(kk1,j2)
            end do
            tty(j1,j2) = alnq1 * sum
         end do
      end do
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         alnq1 = l1 - qn
         do j2 = 1, jmax2
            l2 = j2 + lm2
            alnq2 = l2 - qn
            ttz(j1,j2) = alnq1 * alnq2 * workg(j1,j2)
         end do
      end do
!
!      fix kinetic energy parts...
!
!!$  zf(1:mth1) = rhoa(surf) * xjacob(1:mth1,surf)
!
! take constant rho
!
  zf(1:mth1) = xjacob(1:mth1,surf)
!
      CALL pstloadw(mth1, zf, work)
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do j1 = 1, jmax2
         l1 = j1 + lm2
         ttd(j1,j1) = real( work(1) )
         do j2 = 1,j1-1
            l2 = j2 + lm2
            lc1 = abs(l1-l2) + 1
            ttd(j1,j2) =  work(lc1)
            ttd(j2,j1) = conjg( ttd(j1,j2) )
         end do
      end do
!
      tw(scount,1:jmax1,1:jmax1) = ttw(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      ty(scount,1:jmax1,1:jmax1) = tty(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      tz(scount,1:jmax1,1:jmax1) = ttz(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      td(scount,1:jmax1,1:jmax1) = ttd(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
!
      return
!
! aplet  19.4_r8 . 94._r8 ..
!
! bad inversion
!
 99   continue
!
      write(itty  ,*)' operator g^-1 not positive definite'
      write(outmod,*)' operator g^-1 not positive definite'
      write(itty  ,*)' surf = ', surf,' psinew = ', psinew(surf)
      write(outmod,*)' surf = ', surf,' psinew = ', psinew(surf)
      write(itty  ,*)' check sign of jacobian xjacob = '
      write(outmod,*)' check sign of jacobian xjacob = '
      jcol = 4
      jth = mth/jcol
      do jl = 1, jth
         t = jcol*(jl-1)
         write(itty  ,8000) ( xjacob(t+j2,surf), j2 = 1,jcol )
         write(outmod,8000) ( xjacob(t+j2,surf), j2 = 1,jcol )
      end do
      write(itty  ,*)' delta_theta = '
      write(outmod,*)' delta_theta = '
      do jl = 1, mth/4
         t = jcol*(jl-1) 
         i1 = t + 1
         i2 = t + 2
         i3 = t + 3
         i4 = t + 4
         z1 = f * xjacob(i1,surf)/xsq(i1,surf) -  1.0_r8 
         z2 = f * xjacob(i2,surf)/xsq(i2,surf) -  1.0_r8 
         z3 = f * xjacob(i3,surf)/xsq(i3,surf) -  1.0_r8 
         z4 = f * xjacob(i4,surf)/xsq(i4,surf) -  1.0_r8 
         write(itty  ,8000) z1,z2,z3,z4
         write(outmod,8000) z1,z2,z3,z4
      end do
 8000 format(4(1x,e15.8))
!
      stop ' ERROR in thfft'
!
 7000 CALL psterrmes ( outpst,'pstthefft',lgivup )
      end



