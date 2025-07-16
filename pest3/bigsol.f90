subroutine pstBigSolution

  ! compute the big solution arrays

  USE pstcom

  USE l22com

  USE newmet

  USE r33com

  USE comggf

  USE comfft

  implicit none

  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  INTEGER MS1
  INTEGER NRAT
  INTEGER I
  INTEGER JRATL
  INTEGER JRATR
  real*8 ZDPSI
  real*8 ZDPSI2
  real*8 ZDIFF
  INTEGER LM1
  INTEGER LM2
  INTEGER JINC
  INTEGER LS
  INTEGER JS
  INTEGER JJS
  real*8 QP
  real*8 QPN
  real*8 QPP
  real*8 QPPN
  real*8 ZGM1S
  real*8 ZRS
  real*8 ZMUS
  real*8 PB
  real*8 PB1
  real*8 PB2
  real*8 PS
  real*8 PS1
  real*8 PS2
  real*8 ZNORM
  real*8 ZNORM2
  INTEGER JJMAX1
  INTEGER JINC1
  INTEGER JJ1
  INTEGER JJ2
  INTEGER J1
  INTEGER J2
  INTEGER J
  INTEGER JMIN
  INTEGER JMAX
  INTEGER L1
  real*8 ZPBK
  real*8 ZER11
  real*8 ZER12
  real*8 ZER21
  real*8 ZER22
  real*8 ZTOT
  INTEGER J0
  real*8 ZD2
  real*8 ZD1
  real*8 ZD12
  real*8 ZD22
  real*8 ZD110
  real*8 ZD210
  real*8 ZPSING
  INTEGER JSIDE
  INTEGER JSUMIN
  INTEGER JSUMAX
  real*8 ZPARIT
  real*8 ZDLAYB
  real*8 ZDLAY2
  real*8 ZXS
  real*8 ZAXS
  real*8 ZLNX
  real*8 ZD
  real*8 ZARG
  real*8 ZY
  real*8 ZHPLUS
  real*8 ZDHPLS
  real*8 ALNQS
  real*8 ALNQ1
  INTEGER JS1
  INTEGER JS2
  INTEGER JINCR
  real*8 ZDLA10
  real*8 ZNEW
  real*8 ZOLD
  INTEGER JSUR1
  real*8 ZNEWOL
  real*8 ZLN
  real*8 ZPOWER
  real*8 ZEXPAO
  real*8 ZEXPAE
  real*8 ZEXPGO
  real*8 ZEXPGE
  real*8 ZERRAO
  real*8 ZERRAE
  real*8 ZERRGO
  real*8 ZERRGE
  real*8 ZXNEW
  real*8 ZXOLD
  real*8 ZABSX
  real*8 ZDX
  INTEGER J3
  INTEGER JLMIN
  INTEGER JLMAX



  !
  COMPLEX*16	 	 	 :: zc
  COMPLEX*16	 	 	 :: zc1
  COMPLEX*16	 	 	 :: zc2
  COMPLEX*16	 	 	 :: zxibig
  COMPLEX*16	 	 	 :: zchbig
  COMPLEX*16	 	 	 :: zxsmal
  COMPLEX*16	 	 	 :: zcsmal
  COMPLEX*16	 	 	 :: zxadjb
  COMPLEX*16	 	 	 :: zcadjb
  COMPLEX*16	 	 	 :: zxadjs
  COMPLEX*16	 	 	 :: zcadjs
  COMPLEX*16	 	 	 :: z11
  COMPLEX*16	 	 	 :: z12
  COMPLEX*16	 	 	 :: z21
  COMPLEX*16	 	 	 :: z22
  COMPLEX*16	 	 	 :: zhchi
  COMPLEX*16	 	 	 :: zhxi
  COMPLEX*16	 	 	 :: zdhchi
  COMPLEX*16	 	 	 :: zdhxi
  COMPLEX*16	 	 	 :: zqs
  COMPLEX*16	 	 	 :: zqss
  COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zrpr
  COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zqpr
  COMPLEX*16, DIMENSION(nfn,nfn) 	 :: zgm1pr
  COMPLEX*16, DIMENSION(2,2) 	 :: zainv
  COMPLEX*16, DIMENSION(2,2) 	 :: zaa
  COMPLEX*16, DIMENSION(2,2) 	 :: zaaa
  COMPLEX*16, DIMENSION(2) 	 :: zay
  COMPLEX*16, DIMENSION(nfn) 	 :: zaoold
  COMPLEX*16, DIMENSION(nfn) 	 :: zaeold
  COMPLEX*16, DIMENSION(nfn) 	 :: zgoold
  COMPLEX*16, DIMENSION(nfn) 	 :: zgeold
  COMPLEX*16, DIMENSION(nfn) 	 :: zaonew
  COMPLEX*16, DIMENSION(nfn) 	 :: zaenew
  COMPLEX*16, DIMENSION(nfn) 	 :: zgonew
  COMPLEX*16, DIMENSION(nfn) 	 :: zgenew
  COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: qm
  COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: rm
  COMPLEX*16, DIMENSION(2,nfn,nfn) 	 :: gmm1
  REAL*8, DIMENSION(nfm) 	 :: zenr
  REAL*8, DIMENSION(nfm) 	 :: zeni

!
! no of rational surface...
!
      nrat = 1
      do i = 1, ms
         nrat = nrat + ( msing(i+1)-msing(i) )*msub(i)
      end do
      nrat = nrat + msub(ms)/2
      !
      ! no of previous rational surface
      !
      jratl = 1
      if ( ms >  1 ) then
         do i = 1,ms-1
            jratl = jratl + ( msing(i+1)-msing(i) )*msub(i)
         end do
         jratl = jratl + msub(ms-1)/2
      end if
      !
      ! no of following rational surface
      !
      jratr = nusurf
      if ( ms < nosing ) then
         jratr = 1
         do i = 1,ms+1
            jratr = jratr + ( msing(i+1)-msing(i) )*msub(i)
         end do
         jratr = jratr - msub(ms)/2
      end if
!

      lm1 = lmin(1) - 1
      jmax1 = lmax(1) - lmin(1) + 1
      lm2 = -(jmax2-1)/2 - 1
      if ( jmax2  ==  1 )  lm2 = lm1
      jinc = lm1 - lm2
      ls = n*qa(nrat)+ 0.1_r8 
      js = ls - lmin(1) + 1
      jjs = js + jinc


  ! compute source term. the values of the source function is stored
  ! in
  ! arrays alxi1e, alxi1o, qdchie, qdchio, ddchie and ddchio, for big
  ! prescribed frobenius solution of even (e) and odd (o) parity.
  !
  pb = alphb(ms)
  ps = alphs(ms)
  zmus = (ps - pb)/ 2.0_r8 
  zpbk = pb +  2.0_r8 
  if( zmus <=   1.0_r8  ) zpbk = pb +  1.0_r8 
  !
  x1frbo(1:nusurf,1:jmax1,ms) =  0.0_r8 
  x1frbe(1:nusurf,1:jmax1,ms) =  0.0_r8 
  x1fpro(1:nusurf,1:jmax1,ms) =  0.0_r8 
  x1fpre(1:nusurf,1:jmax1,ms) =  0.0_r8 
  c1frbo(1:nusurf,1:jmax1,ms) =  0.0_r8 
  c1frbe(1:nusurf,1:jmax1,ms) =  0.0_r8 
  c1fpro(1:nusurf,1:jmax1,ms) =  0.0_r8 
  c1fpre(1:nusurf,1:jmax1,ms) =  0.0_r8 
  alxi1o(1:nusurf,1:jmax1)    =  0.0_r8 
  alxi1e(1:nusurf,1:jmax1)    =  0.0_r8 
  gm1dco(1:nusurf,1:jmax1)    =  0.0_r8 
  gm1dce(1:nusurf,1:jmax1)    =  0.0_r8 
  qdchio(1:nusurf,1:jmax1)    =  0.0_r8 
  qdchie(1:nusurf,1:jmax1)    =  0.0_r8 
  ddchio(1:nusurf,1:jmax1)    =  0.0_r8 
  ddchie(1:nusurf,1:jmax1)    =  0.0_r8 
  !
  zd2 = psinod( msing(ms+2) ) - psinod( msing(ms+1) )
  zd1 = psinod( msing(ms+1) ) - psinod( msing(ms) )
  !
  zd1 = abs(dlayb) * zd1
  zd2 = abs(dlayb) * zd2
  !
  if ( lsymhi ) then
     !
     ! symmetric shape function
     !
     if ( zd1 < zd2 ) then
        zd2 = zd1
     else
        zd1 = zd2
     end if
     !
  end if
  !
  !aplet  20.8_r8 .93      zd12 = zd1/ 2._r8 
  !aplet  20.8_r8 .93      zd22 = zd2/ 2._r8 
  zd12 = zd1/ 3._r8 
  zd22 = zd2/ 3._r8 
  !aplet  20.8_r8 .93
  zd110 = zd1/ 10._r8 
  zd210 = zd2/ 10._r8 
  !
  ! position of rational surface...
  !
  zpsing = psinew(nrat)
  !
  do jside = 1,2
     !
     ! left-hand side of rational surface...
     !
     jsumin = jratl
     jsumax = nrat - 1
     zparit = - 1._r8 
     zdlayb = zd1
     zdlay2 = zd12
     !
     ! right-hand side of rational surface...
     !
     if( jside == 2) then
        jsumin = nrat + 1
        jsumax = jratr
        zparit = +  1._r8 
        zdlayb = zd2
        zdlay2 = zd22
     end if
     !
     ! run over nodes...
     !
     do surf = jsumin,jsumax
        !
        zxs = psinew(surf) - zpsing
        zaxs = abs( zxs )
        zlnx = log( zaxs )
        zd = zdlayb - zdlay2
        zarg = pye*( zaxs - zdlay2 )/zd
        zy = zxs/zdlayb
        !
        zhplus =  1.0_r8 
        if(zaxs >  zdlay2) zhplus = ( cos(zarg) +  1.0_r8  )/ 2.0_r8 
        if(zaxs >  zdlayb) zhplus =  0.0_r8 
        !aplet  25.8_r8 .93      zhplus =  0.0_r8 
        !aplet  25.8_r8 .93      if(zaxs < zdlayb)
        !aplet  25.8_r8 .93     .  zhplus = exp( - zy**4/( 1._r8  - zy**2) )
        !aplet  16.4_r8 .93    
        !
        zdhpls =  0.0_r8 
        !aplet  20.8_r8 .93      if( (zaxs >  zdlay2) .AND. (zaxs < zdlayb) )
        !aplet  20.8_r8 .93     .     zdhpls = - pye*sin(zarg)*zparit/zdlayb
        if( (zaxs >  zdlay2) .AND. (zaxs < zdlayb) ) &
             zdhpls = - pye*sin(zarg)*zparit/( 2._r8 *zd)
        !aplet  25.8_r8 .93
        !aplet  25.8_r8 .93      if(zaxs < zdlayb)
        !aplet  25.8_r8 .93     .  zdhpls = - zhplus * ( 4._r8 *zy**3 -  2._r8 *zy**5)
        !aplet  25.8_r8 .93     .         /( zdlayb*     ( 1._r8  - zy**2)**2 )
        !aplet  16.4_r8 .93
        !
        qpn = n*qpa(surf)
        alnqs = ls - n*qa(surf)
        !
        do 71 j1 = 1,jmax1
           jj1 = j1 + jinc
           l1 = j1 + lm1
           alnq1 = l1 - n*qa(surf)
           !
           zhchi = zhplus*(chi0(j1,ms)*zxs*zaxs**(pb- 1._r8 ) &
                + chlog(j1,ms)*zlnx + chcon(j1,ms) &
                + chi1(j1,ms)*zxs**2 *zaxs**(pb- 1._r8 ) &
                + chi2(j1,ms)*zxs**3 *zaxs**(pb- 1._r8 ) &
                )
           !
           zhchi = zhchi + zhplus * ( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * ( csml0(j1,ms)*zxs + csml1(j1,ms)*zxs**2 ) &
                + xsmnus(ms)*zaxs**ps &
                * ( csml0(j1,ms)     + csml1(j1,ms)*zxs    ) &
                )
           !
           zhxi  = zhplus*( xix0(j1,ms)*zxs*zaxs**(pb- 1._r8 ) &
                + xilog(j1,ms)*zlnx + xicon(j1,ms) &
                + xix1(j1,ms)*zxs**2 *zaxs**(pb- 1._r8 ) &
                + xix2(j1,ms)*zxs**3 *zaxs**(pb- 1._r8 ) &
                )
           !
           zhxi  = zhxi + zhplus * ( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * ( xsml0(j1,ms)*zxs + xsml1(j1,ms)*zxs**2 ) &
                + xsmnus(ms)*zaxs**ps &
                * ( xsml0(j1,ms)     + xsml1(j1,ms)*zxs    ) &
                )
           !
           zdhchi = zdhpls*(chi0(j1,ms)*zxs*zaxs**(pb- 1._r8 ) &
                + chlog(j1,ms)*zlnx + chcon(j1,ms) &
                + chi1(j1,ms)*zxs**2 *zaxs**(pb- 1._r8 ) &
                + chi2(j1,ms)*zxs**3 *zaxs**(pb- 1._r8 ) &
                ) &
                + zhplus*(chi0(j1,ms)*pb*zaxs**(pb- 1._r8 ) &
                + chlog(j1,ms)/zxs &
                + chi1(j1,ms)*(pb+ 1._r8 )*zxs*zaxs**(pb- 1._r8 ) &
                + chi2(j1,ms)*(pb+ 2._r8 )*zxs**2 *zaxs**(pb- 1._r8 ) &
                )
           zdhchi = zdhchi + zdhpls*( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * ( csml0(j1,ms)*zxs + csml1(j1,ms)*zxs**2 ) &
                + xsmnus(ms)*zaxs**ps &
                * ( csml0(j1,ms)     + csml1(j1,ms)*zxs    ) &
                ) &
                + zhplus*( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * (ps*csml0(j1,ms)     + (ps+ 1._r8 )*csml1(j1,ms)*zxs   ) &
                + xsmnus(ms)*zaxs**(ps- 2._r8 ) &
                * (ps*csml0(j1,ms)*zxs + (ps+ 1._r8 )*csml1(j1,ms)*zxs**2) &
                )
           !
           zdhxi  = zdhpls*(xix0(j1,ms)*zxs*zaxs**(pb- 1._r8 ) &
                + xilog(j1,ms)*zlnx + xicon(j1,ms) &
                + xix1(j1,ms)*zxs**2 *zaxs**(pb- 1._r8 ) &
                + xix2(j1,ms)*zxs**3 *zaxs**(pb- 1._r8 ) &
                ) &
                + zhplus*(xix0(j1,ms)*pb*zaxs**(pb- 1._r8 ) &
                + xilog(j1,ms)/zxs &
                + xix1(j1,ms)*(pb+ 1._r8 )*zxs*zaxs**(pb- 1._r8 ) &
                + xix2(j1,ms)*(pb+ 2._r8 )*zxs**2 *zaxs**(pb- 1._r8 ) &
                )
           !
           zdhxi  = zdhxi  + zdhpls*( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * ( xsml0(j1,ms)*zxs + xsml1(j1,ms)*zxs**2 ) &
                + xsmnus(ms)*zaxs**ps &
                * ( xsml0(j1,ms)     + xsml1(j1,ms)*zxs    ) &
                ) &
                + zhplus*( &
                xsplus(ms)*zaxs**(ps- 1._r8 ) &
                * (ps*xsml0(j1,ms)     + (ps+ 1._r8 )*xsml1(j1,ms)*zxs   ) &
                + xsmnus(ms)*zaxs**(ps- 2._r8 ) &
                * (ps*xsml0(j1,ms)*zxs + (ps+ 1._r8 )*xsml1(j1,ms)*zxs**2) &
                )
           !
           x1frbo(surf,j1,ms) = zhxi
           x1fpro(surf,j1,ms) = zdhxi
           c1frbo(surf,j1,ms) = zhchi
           c1fpro(surf,j1,ms) = zdhchi
           !
           x1frbe(surf,j1,ms) = zparit*x1frbo(surf,j1,ms)
           x1fpre(surf,j1,ms) = zparit*x1fpro(surf,j1,ms)
           c1frbe(surf,j1,ms) = zparit*c1frbo(surf,j1,ms)
           c1fpre(surf,j1,ms) = zparit*c1fpro(surf,j1,ms)
           !
71         continue
           !
        end do
     end do
     !
     ! at rational surface
     !
     surf = nrat
     x1frbo(surf,1:jmax1,ms) =  0.0_r8 
     x1fpro(surf,1:jmax1,ms) =  0.0_r8 
     c1frbo(surf,1:jmax1,ms) =  0.0_r8 
     c1fpro(surf,1:jmax1,ms) =  0.0_r8 
     x1frbe(surf,1:jmax1,ms) =  0.0_r8 
     x1fpre(surf,1:jmax1,ms) =  0.0_r8 
     c1frbe(surf,1:jmax1,ms) =  0.0_r8 
     c1fpre(surf,1:jmax1,ms) =  0.0_r8 
     !
     ! compute alxi1o, alxi1e, gm1dco, gm1dc1, qdchio,
     ! qdchie, ddchio and ddchie...
     !
     !
     !
     if( lchkap ) then
        write(outmod,201)
201     format(3x,'psi',6x,'q',11x,'ttq',22x,'tr',21x,'gm1')
     end if
     !
     do surf = jratl,jratr
        !
!!$  26.2_r8 .97        call thfft2
        !call thefft
        call pstthfft2
        !
        do j1 = 1,jmax1
           jj1 = j1 + jinc
           l1 = j1 + lm1
           alnq1 = l1 - n*qa(surf)
           qpn = n*qpa(surf)
           !
           alxi1o(surf,j1) = - alnq1* c1fpro(surf,j1,ms) &
                + qpn* c1frbo(surf,j1,ms)
           alxi1e(surf,j1) = - alnq1* c1fpre(surf,j1,ms) &
                + qpn* c1frbe(surf,j1,ms)
           gm1dco(surf,j1) = - alnq1* x1fpro(surf,j1,ms)
           gm1dce(surf,j1) = - alnq1* x1fpre(surf,j1,ms)
           !
           do 80 j2 = 1,jmax1
              jj2 = j2 + jinc
              !
              alxi1o(surf,j1) = alxi1o(surf,j1) &
                   + conjg(ttq(jj2,jj1))*c1frbo(surf,j2,ms) &
                   - ttr(jj1,jj2)*x1frbo(surf,j2,ms)
              alxi1e(surf,j1) = alxi1e(surf,j1) &
                   + conjg(ttq(jj2,jj1))*c1frbe(surf,j2,ms) &
                   - ttr(jj1,jj2)*x1frbe(surf,j2,ms)
              gm1dco(surf,j1) =  gm1dco(surf,j1) &
                   - ttq(jj1,jj2)*x1frbo(surf,j2,ms) &
                   + gm1(j1,j2)* c1frbo(surf,j2,ms)
              gm1dce(surf,j1) =  gm1dce(surf,j1) &
                   - ttq(jj1,jj2)*x1frbe(surf,j2,ms) &
                   + gm1(j1,j2)* c1frbe(surf,j2,ms)
              !
80            continue
              !
              !
              !     if ( lchkap ) then
              !     write(outmod,200) psinew(surf),qa(surf)
              !    .     ,(ttq(jj1,jj2),jj2=1+jinc,jmax1+jinc)
              !    .     ,(ttr(jj1,jj2),jj2=1+jinc,jmax1+jinc)
              !    .     ,(gm1(j1,j2),j2=1,jmax1)
              !     end if
200           format(f8.7,18(f8.2))
              !
           end do
        end do
        !
        ! correct behaviour near singular surface to remedy
        ! small inaccuracies of equilibrium quantities...
        !
        !--------------------------
        if ( .not. lsmoth ) go to 107
        !
        !............................
        if( .not. lzerop ) then
           !
           do jside = 1,2
              js1 = jratl
              js2 = nrat - 1
              jincr = +1
              zdla10 = zd110
              if(jside == 2) then
                 js1 = jratr
                 js2 = nrat + 1
                 jincr = -1
                 zdla10 = zd210
              end if
              znew =  1.E+10_r8 
              do surf = js1,js2,jincr
                 zold = znew
                 znew = abs( psinew(surf) - psinew(nrat) )
                 do j1 = 1,jmax1
                    if( (znew < zdla10) .AND. (znew >   1.E-10_r8 ) &
                         .and.(zold >   1.E-10_r8 ) ) then
                       jsur1 = surf - jincr
                       znewol = znew/zold
                       zln = log( znewol )
                       zpower = znewol**zpbk
                       zexpao = log( &
                            abs( alxi1o(surf,j1)/ alxi1o(jsur1,j1) ) &
                            ) / zln
                       zexpae = log( &
                            abs( alxi1e(surf,j1)/ alxi1e(jsur1,j1) ) &
                            ) / zln
                       zexpgo = log( &
                            abs( gm1dco(surf,j1)/ gm1dco(jsur1,j1) ) &
                            ) / zln
                       zexpge = log( &
                            abs( gm1dce(surf,j1)/ gm1dce(jsur1,j1) ) &
                            ) / zln
                       zerrao = abs( zexpao - zpbk )
                       zerrae = abs( zexpae - zpbk )
                       zerrgo = abs( zexpgo - zpbk )
                       zerrge = abs( zexpge - zpbk )
                       if( zerrao >    1._r8  ) &
                            alxi1o(surf,j1) = alxi1o(jsur1,j1)*zpower
                       if( zerrae >    1._r8  ) &
                            alxi1e(surf,j1) = alxi1e(jsur1,j1)*zpower
                       if( zerrgo >    1._r8  ) &
                            gm1dco(surf,j1) = gm1dco(jsur1,j1)*zpower
                       if( zerrge >    1._r8  ) &
                            gm1dce(surf,j1) = gm1dce(jsur1,j1)*zpower
                    end if
                 end do
              end do
           end do
           !.............................
           ! zero pressure correction...
           !
        else
           !
           do jside = 1,2
              js1 = jratl
              js2 = nrat - 1
              jincr = +1
              zdla10 = zd110
              if(jside == 2) then
                 js1 = jratr
                 js2 = nrat + 1
                 jincr = -1
                 zdla10 = zd210
              end if
              !
              zxnew =  1.E+10_r8 
              do j1 = 1, jmax1
                 zaonew(j1) =  1.E+10_r8 
                 zaenew(j1) =  1.E+10_r8 
                 zgonew(j1) =  1.E+10_r8 
                 zgenew(j1) =  1.E+10_r8 
              end do
              !       
              do surf = js1,js2,jincr
                 zxold = zxnew
                 zxnew = psinew(surf) - psinew(nrat)
                 zabsx = abs( zxnew )
                 zdx   = zxnew - zxold
                 zlnx  = log( zabsx )
                 !
                 do 88 j1 = 1, jmax1
                    zpower =  1.0_r8 
                    if( j1 == js ) zpower = zxnew
                    ! 
                    zaoold(j1) = zaonew(j1)
                    zaeold(j1) = zaenew(j1)
                    zgoold(j1) = zgonew(j1)
                    zgeold(j1) = zgenew(j1)
                    zaonew(j1) = alxi1o(surf,j1) / (zlnx*zpower)
                    zaenew(j1) = alxi1e(surf,j1) / (zlnx*zpower)
                    zgonew(j1) = gm1dco(surf,j1) / (zlnx*zpower)
                    zgenew(j1) = gm1dce(surf,j1) / (zlnx*zpower)
                    !
                    if( zabsx < zdla10 ) then
                       !
                       zerrao = abs( psinew(nrat)*(zaonew(j1)-zaoold(j1)) &
                            /(zaonew(j1)*zdx) )
                       zerrae = abs( psinew(nrat)*(zaenew(j1)-zaeold(j1)) &
                            /(zaenew(j1)*zdx) )
                       zerrgo = abs( psinew(nrat)*(zgonew(j1)-zgoold(j1)) &
                            /(zgonew(j1)*zdx) )
                       zerrge = abs( psinew(nrat)*(zgenew(j1)-zgeold(j1)) &
                            /(zgenew(j1)*zdx) )
                       !
                       if( zerrao >    10._r8  ) then
                          alxi1o(surf,j1) = zaoold(j1)*zpower*zlnx
                          zaonew(j1) = zaoold(j1)
                       end if
                       if( zerrae >    10._r8  ) then
                          alxi1e(surf,j1) = zaeold(j1)*zpower*zlnx
                          zaenew(j1) = zaeold(j1)
                       end if
                       if( zerrgo >    10._r8  ) then
                          gm1dco(surf,j1) = zgoold(j1)*zpower*zlnx
                          zgonew(j1) = zgoold(j1)
                       end if
                       if( zerrge >    10._r8  ) then
                          gm1dce(surf,j1) = zgeold(j1)*zpower*zlnx
                          zgenew(j1) = zgeold(j1)
                       end if
                       !
                    end if
                    !
88                  continue
                    !
                 end do
              end do
              !
           end if
           !
107        continue
           !----------------------------
           !
           do surf = jratl,jratr
              !
!!$   26.2_r8 .97      call thfft2
              !call thefft
              call pstthfft2
              !
              do j1 = 1,jmax1
                 qdchio(surf,j1) =  0.0_r8 
                 qdchie(surf,j1) =  0.0_r8 
                 ddchio(surf,j1) =  0.0_r8 
                 ddchie(surf,j1) =  0.0_r8 
                 jj1 = j1 + jinc
                 l1 = j1 + lm1
                 alnq1 = l1 - n*qa(surf)
                 !
                 do j2 = 1,jmax1
                    jj2 = j2 + jinc
                    ddchio(surf,j1) = ddchio(surf,j1) &
                         + alnq1*geem(j1,j2)*gm1dco(surf,j2)
                    ddchie(surf,j1) = ddchie(surf,j1) &
                         + alnq1*geem(j1,j2)*gm1dce(surf,j2)
                    !
                    do j3 = 1,jmax1
                       qdchio(surf,j1) = qdchio(surf,j1) &
                            + conjg(ttq(jj2,jj1))*geem(j2,j3)*gm1dco(surf,j3)
                       qdchie(surf,j1) = qdchie(surf,j1) &
                            + conjg(ttq(jj2,jj1))*geem(j2,j3)*gm1dce(surf,j3)
                       !
                    end do
                 end do
              end do
           end do
           !
72         continue
           !
           if (lchkap) then
              !
              !
              jmin = max(1,js-4)
              jmax = max(jmax1,js+4)
              jlmin = jmin + lm1
              jlmax = jmax + lm1
              !
              !     write(itty,1051)
              write(outmod,1051)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),( real(x1frbo(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1052)
              write(outmod,1052)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),(aimag(x1frbo(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1067)
              write(outmod,1067)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),( real(x1fpro(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1068)
              write(outmod,1068)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),(aimag(x1fpro(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1059)
              write(outmod,1059)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),( real(c1frbo(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1054)
              write(outmod,1054)
              do i = jratl,jratr
                 write(outmod,1053) psinew(i),(aimag(c1frbo(i,j1,ms)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1055)
              write(outmod,1055)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),( real(alxi1o(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1056)
              write(outmod,1056)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),(aimag(alxi1o(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1057)
              write(outmod,1057)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),( real(gm1dco(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1058)
              write(outmod,1058)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),(aimag(gm1dco(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1063)
              write(outmod,1063)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),( real(qdchio(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1064)
              write(outmod,1064)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),(aimag(gm1dco(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1065)
              write(outmod,1065)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),( real(ddchio(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1066)
              write(outmod,1066)
              do i = jratl,jratr
                 write(outmod,1007) psinew(i),(aimag(ddchio(i,j1)), &
                      j1 = 1,jmax1)
              end do
              !
              !     write(itty,1069)
              write(outmod,1069)
              do i = jratl,jratr
                 zenr(j1) = real( &
                      conjg( x1frbo(i,j1,ms) )*alxi1o(i,j1) &
                      - conjg( x1frbo(i,j1,ms) )*qdchio(i,j1) &
                      - conjg( x1fpro(i,j1,ms) )*ddchio(i,j1) &
                      )
                 write(outmod,1007) psinew(i),( zenr(j1), &
                      j1 = 1,jmax1)
              end do
              !     write(itty,1070)
              write(outmod,1070)
              do i = jratl,jratr
                 zeni(j1) =aimag( &
                      conjg( x1frbo(i,j1,ms) )*alxi1o(i,j1) &
                      - conjg( x1frbo(i,j1,ms) )*qdchio(i,j1) &
                      - conjg( x1fpro(i,j1,ms) )*ddchio(i,j1) &
                      )
                 write(outmod,1007) psinew(i),( zeni(j1), &
                      j1 = 1,jmax1)
              end do
              !
106           continue
              !
           end if

           !
1007       format(f10.6,9(1x,e12.5))
1053       format(f10.6,9(1x,e12.5))
1055       format(20x,'re alxi1o')
1051       format(20x,'re x1frbo')
1059       format(20x,'re c1frbo')
1057       format(20x,'re gm1dco')
1052       format(20x,'im x1frbo')
1054       format(20x,'im c1frbo')
1056       format(20x,'im alxi1o')
1058       format(20x,'im gm1dco')
1061       format(20x,'re x1frbo'/'psinew ',9(6x,i2,5x))
1062       format(20x,'im x1frbo'/'psinew ',9(6x,i2,5x))
1063       format(20x,'re qdchio')
1064       format(20x,'im qdchio')
1065       format(20x,'re ddchio')
1066       format(20x,'im ddchio')
1067       format(20x,'re x1fpro')
1068       format(20x,'im x1fpro')
1069       format(20x,'re prescribed energy density')
1070       format(20x,'im prescribed energy density')
           !

      end subroutine 
