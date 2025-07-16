! 
!     purpose      calculate and store the theta integrals on psi
!                  surface number surf. in temporary storage scount
!                  points to the appropriate surface lying between the
!                  centre and r.h.s of element support.
!.......................................................................
!      ray grimm, jan, 1980._r8  j manickam feb. 1982 a pletzer 1993
!
! use transformed operator form (aplet 18/8/95)
!.......................................................................
!
      subroutine pstthefft
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

 INTEGER :: lgivup, jj1, jj2, j1 ,j2, l1 ,l2, lc1, lc2, kk1, ierr, j, l12 &
            , jinc, lm1, lm2, mfft, mrem, iset, ifr, lsgn
 REAL*8    :: coeff, delth, alnq1, alnq2, xj, sig, sigp, pp, ps, qp, qn, qsq &
            , g2, g, f, ggp, q, r2ggp
!
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
 REAL*8, DIMENSION(nths) 	 :: work

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
      write(itty,9000)
 9000 format( 1x, " ****error in calling fft...too large llp ****")
      stop 'ERROR too many Fourier modes in pstthefft'
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
      r2ggp = g * gp
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

work =  0.0_r8 
ttw =  0.0_r8 
ttz =  0.0_r8 
tty =  0.0_r8 
ttq =  0.0_r8 
ttd =  0.0_r8 

!
!      begin working on w-matrix with parts not involving gee..
!
  zf(1:mth1) =  1.0_r8  / ( xjacob(1:mth1,surf)*grpssq(1:mth1,surf) )
  CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 12 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      ttw(j1,j1) = alnq1*alnq1 * work(1)
      do 12 j2 = 1, j1-1
      l2 = j2 + lm2
      alnq2 = l2 - qn
      l12 = abs(l1-l2)
      lc1 = 2*l12 + 1
      lc2 = lc1 + 1
      ttw(j1,j2) = alnq1*alnq2 * (work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2))
      ttw(j2,j1) = alnq1*alnq2 * (work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2))
   12 continue
!
      do 20 t = 1, mth1
      xj = xsq(t,surf) * pp + r2ggp
      zf(t)=xj*( (xjprym(t,surf)-xjacob(t,surf)*xsqdps(t,surf) &
        /xsq(t,surf)) /xsq(t,surf) &
               - xjacob(t,surf)*xj/(xsq(t,surf)*grpssq(t,surf))   )
   20 continue
CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 22 j1 = 1, jmax2
      l1 = j1 + lm2
      ttw(j1,j1) = ttw(j1,j1) + work(1)
      do 22 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
      ttw(j1,j2) = ttw(j1,j2) + (work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2))
      ttw(j2,j1) = ttw(j2,j1) + (work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2))
   22 continue

      do 30 t = 1, mth1
   30 zf(t) = pp * xjacob(t,surf)*xsqdps(t,surf) /xsq(t,surf)
CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 32 j1 = 1, jmax2
      l1 = j1 + lm2
      ttw(j1,j1) = ttw(j1,j1) + work(1)
      do 32 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
      ttw(j1,j2) = ttw(j1,j2) + (work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2))
      ttw(j2,j1) = ttw(j2,j1) + (work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2))
   32 continue
!
      do 40 t = 1, mth1
      zf(t) = ( xsq(t,surf)*pp + r2ggp )*grpsth(t,surf)*xjacob(t,surf) &
         / ( xsq(t,surf)*grpssq(t,surf) )
   40 continue
CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 42 j1 = 1, jmax2
      l1 = j1 + lm2
      do 42 j2 = 1, j1-1
      l2 = j2 + lm2
      l12 = abs( l1-l2 )
      lc1 = 2*l12 + 1
      lc2 = lc1 + 1
!!$      ttw(j1,j2) = ttw(j1,j2) - l12*work(lc2)
      ttw(j1,j2) = ttw(j1,j2) + l12*(-work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1))
      ttw(j2,j1) = ttw(j2,j1) - l12*( work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1))
   42 continue
!
!
   45 continue
!
!      now compute gee-matrix..
!
      do 50 t = 1, mth1
      delth = f * xjacob(t,surf)/xsq(t,surf) -  1.0_r8 
      zf(t) = xjacob(t,surf)/xsq(t,surf) &
    + qsq * xsq(t,surf)*delth**2/(xjacob(t,surf)*grpssq(t,surf))
   50 continue
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 52 j1 = 1, jmax2
      l1 = j1 + lm2
      workg(j1,j1) = ensq * work(1)
      do 52 j2 = 1, j1-1
      l2 = j2 + lm2
      l12 = abs(l1-l2)
      lc1 = 2*l12 + 1
      lc2 = lc1 + 1
      workg(j1,j2) = ensq * (work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2))
      workg(j2,j1) = ensq * (work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2))
   52 continue
!
      do 53 t = 1, mth1
   53 zf(t) = xsq(t,surf) / ( xjacob(t,surf)*grpssq(t,surf) )
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 54 j1 = 1, jmax2
      l1 = j1 + lm2
      workf(j1,j1) = work(1)
      workg(j1,j1) = workg(j1,j1) + l1*l1*workf(j1,j1)
      do 54 j2 = 1, j1-1
      l2 = j2 + lm2
      l12 = abs(l1-l2)
      lc1 = 2*l12 + 1
      lc2 = lc1 + 1
      workf(j1,j2) = work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2)
      workf(j2,j1) = work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2)
      workg(j1,j2) = workg(j1,j2) + l1*l2*workf(j1,j2)
      workg(j2,j1) = workg(j2,j1) + l2*l1*workf(j2,j1) 
   54 continue
!
      do 55 t = 1, mth1
   55 zf(t) = f * xjacob(t,surf) / xsq(t,surf) -  1.0_r8 
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
!!$      do 56 j1 = 1, jmax2
!!$      l1 = j1 + lm2
!!$      workc(j1,j1) = work(1)
!!$      do 56 j2 = 1, j-1
!!$      l2 = j2 + lm2
!!$      lc1 = 2*abs(l1-l2) + 1
!!$      lc2 = lc1 + 1
!!$      workc(j1,j2) = work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2)
!!$      workc(j2,j1) = work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2)
!!$   56 continue
      do 56 j1 = 1, jmax2
      l1 = j1 + lm2
      workc(j1,j1) = work(1)
      do 56 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
      workc(j1,j2) = work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2)
      workc(j2,j1) = work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2)
   56 continue
!
      do 560 t = 1, mth1
  560 zf(t) = (f-xsq(t,surf)/xjacob(t,surf))/grpssq(t,surf)
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 57 j1 = 1, jmax2
      l1 = j1 + lm2
      workg(j1,j1) = workg(j1,j1) + qn*(l1+l1)*work(1)
      do 57 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
      workg(j1,j2) = workg(j1,j2) + qn*(l1+l2)*(work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2))
      workg(j2,j1) = workg(j2,j1) + qn*(l2+l1)*(work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2))
   57 continue
!
      do 571 j1 = 1, jmax1
      jj1 = j1 + jinc
      do 571 j2 = 1, jmax1
      jj2 = j2 + jinc
      gm1(j1,j2) = workg(jj1,jj2)
  571 continue
  572 continue
!
 ierr = 0
!!$     CALL pstinverc ( workg,nfn,jmax2,ierr )
 
 ! invert workg 

     zidentity = (0.0_r8, 0.0_r8)
     do j1=1, nfn
        zidentity(j1, j1) = (1.0_r8, 0.0_r8)
     enddo
     call ZGESV( nfn, nfn, workg, nfn, IPIV, zidentity, nfn, ierr)
     workg = zidentity ! load solution
     if(ierr/=0) print*,'ERROR ',ierr,' occurred while attempting to invert workg'
!
! aplet nov90
!
      do 90 j1 = 1,jmax1
      jj1 = j1 + jinc
      do 90 j2 = 1,jmax1
      jj2 = j2 + jinc
      geem(j1,j2) = workg(jj1,jj2)
 90   continue
!
      do 60 t = 1, mth1
   60 zf(t) = grpsth(t,surf) / grpssq(t,surf)
CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
!      initialize ttq...
!
      do 63 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      ttq(j1,j1) = - alnq1* l1* ( 0._r8 , 1._r8 ) *work(1)
!!$      do 63 j2 = 1, jmax2
do 63 j2 = 1, j1-1
      l2 = j2 + lm2
      alnq2 = l2 - qn
      l12 = l1-l2
!!$lsgn = 1
!!$if( l12 > 0 ) lsgn = -1
!!$if( l12 == 0 ) lsgn = 0
      lc1 = 2*abs(l12) + 1
      lc2 = lc1 + 1
!!$      ttq(j1,j2) = alnq2* l1*lsgn*work(lc2)
      ttq(j1,j2) = alnq2* l1*(-work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
      ttq(j2,j1) = alnq1* l2*( work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
   63 continue
!
      do 630 t = 1, mth1
  630 zf(t) = grpsth(t,surf)*(f*xjacob(t,surf)/xsq(t,surf)- 1.0_r8 ) &
              / grpssq(t,surf)
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 631 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      ttq(j1,j1) = ttq(j1,j1) - ( 0._r8 , 1._r8 ) *qn*alnq1*work(1)

      do 631 j2 = 1, j1-1
      l2 = j2 + lm2
      alnq2 = l2 - qn
      l12 = l1 - l2
!!$      if ( l12  ==  0 )  go to 631
!!$      lsgn = 1
!!$      if ( l12  >   0 )  lsgn = -1
      lc1 = 2 * abs(l12) + 1
      lc2 = lc1 + 1
!!$      ttq(j1,j2) = ttq(j1,j2) + lsgn*qn*alnq2*work(lc2)
      ttq(j1,j2) = ttq(j1,j2) + qn*alnq2*(-work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
      ttq(j2,j1) = ttq(j2,j1) + qn*alnq2*( work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
  631 continue
!
      do 64 t = 1, mth1
   64 zf(t) = ( xsq(t,surf)*pp + r2ggp ) / grpssq(t,surf)
CALL pstloadw(mth1, zf, work)

      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 66 j1 = 1, jmax2
      l1 = j1 + lm2
      ttq(j1,j1) = ttq(j1,j1) + l1* work(1)  &
                   + n*gp*workc(j1,j1) / f

      do 66 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
!!$      ttq(j1,j2) = ttq(j1,j2) + l1*work(lc1) + n*gp*workc(j1,j2) / f
      ttq(j1,j2) = ttq(j1,j2) + l1*( work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2) ) &
                   + n*gp*workc(j1,j2) / f
      ttq(j2,j1) = ttq(j2,j1) + l2*( work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2) ) &
                   + n*gp*workc(j2,j1) / f

   66 continue
!
!
      do 660 t = 1, mth1
  660 zf(t) = (xsq(t,surf)*pp+r2ggp)*(f*xjacob(t,surf) &
      /xsq(t,surf)- 1.0_r8 )  / grpssq(t,surf)
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 661 j1 = 1, jmax2
      ttq(j1,j1) = ttq(j1,j1) + qn*work(1)
      do 661 j2 = 1, j1-1
      lc1 = 2*abs(j1-j2) + 1
      lc2 = lc1 + 1
!!$      ttq(j1,j2) = ttq(j1,j2) + qn*work(lc1)
      ttq(j1,j2) = ttq(j1,j2) + qn*( work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2) )
      ttq(j2,j1) = ttq(j2,j1) + qn*( work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2) )
  661 continue
!
      do 68 t = 1, mth1
   68 zf(t) = qdelp(t,surf)
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 69 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      ttq(j1,j1) = ttq(j1,j1) - n*alnq1*( 0._r8 , 1._r8 ) *work(1)

      do 69 j2 = 1, j1-1
      l2 = j2 + lm2
      alnq2 = l2 - qn
      l12 = l1 - l2
      if ( l12  ==  0 )  go to 69
      lsgn = 1
      if ( l12  >   0 )  lsgn = -1
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
!!$      ttq(j1,j2) = ttq(j1,j2) + lsgn*n*alnq2*work(lc2)
      ttq(j1,j2) = ttq(j1,j2) + n*alnq2*(-work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
      ttq(j2,j1) = ttq(j2,j1) + n*alnq1*( work(lc2) - ( 0._r8 , 1._r8 ) *work(lc1) )
   69 continue
!
      coeff = n * ( gp/f-qp )
      do 70 j1 = 1, jmax2
      ttq(j1,j1) = ttq(j1,j1) + coeff
   70 continue
!
!
!      finish off w-matrix...
!
!!$      do 75 j1 = 1, jmax2
!!$      do 74 j  = 1, jmax2
!!$   74 workc(j1,j) =  0.0_r8 
!!$      do 75 kk1 = 1, jmax2
!!$      do 75 j2  = 1, jmax2
!!$   75 workc(j1,j2) = workc(j1,j2) + workg(j1,kk1) * ttq(kk1,j2)

workc(1:jmax2,1:jmax2) = matmul( workg(1:jmax2,1:jmax2), ttq(1:jmax2,1:jmax2) )

!!$      do 77 j1 = 1, jmax2
!!$      do 76 j  = 1, jmax2
!!$   76 workf(j1,j)  =  0.0_r8 
!!$      do 77 kk1 = 1, jmax2
!!$      do 77 j2 = 1, jmax2
!!$   77 workf(j1,j2) = workf(j1,j2) + ttq(kk1,j1) * workc(kk1,j2)

workf(1:jmax2,1:jmax2) = matmul( conjg(transpose(ttq(1:jmax2,1:jmax2))), &
         workc(1:jmax2,1:jmax2) )

      ttr(1:jmax2,1:jmax2) = -ttw(1:jmax2,1:jmax2)
      ttw(1:jmax2,1:jmax2) = ttw(1:jmax2,1:jmax2) + workf(1:jmax2,1:jmax2)
!
      do 86 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      do 86 j2 = 1, jmax2
      sum =  0.0_r8 
      do 84 kk1 = 1, jmax2
   84 sum = sum + workg(j1,kk1) * ttq(kk1,j2)
      tty(j1,j2) = alnq1 * sum
   86 continue
!
      do 88 j1 = 1, jmax2
      l1 = j1 + lm2
      alnq1 = l1 - qn
      do 88 j2 = 1, jmax2
      l2 = j2 + lm2
      alnq2 = l2 - qn
      ttz(j1,j2) = alnq1 * alnq2 * workg(j1,j2)
   88 continue
!
!      fix kinetic energy parts...
!
      do 110 t = 1, mth1
110 work(t) = xjacob(t,surf)
CALL pstloadw(mth1, zf, work)
!
      CALL pstfft2 ( work,-1,lfft,invfft,sfft,omfft,onfft,ifr )
!
      do 112 j1 = 1, jmax2
      l1 = j1 + lm2
      ttd(j1,j1) =  work(1)
      do 112 j2 = 1, j1-1
      l2 = j2 + lm2
      lc1 = 2*abs(l1-l2) + 1
      lc2 = lc1 + 1
      ttd(j1,j2) =  work(lc1) - ( 0._r8 , 1._r8 ) *work(lc2)
      ttd(j2,j1) =  work(lc1) + ( 0._r8 , 1._r8 ) *work(lc2)
  112 continue
!
      tw(scount,1:jmax1,1:jmax1) = ttw(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      ty(scount,1:jmax1,1:jmax1) = tty(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      tz(scount,1:jmax1,1:jmax1) = ttz(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
      td(scount,1:jmax1,1:jmax1) = ttd(1+jinc:jmax1+jinc,1+jinc:jmax1+jinc)
!
!
      return
 7000 CALL psterrmes ( outpst,'pstthefft',lgivup )
      end



