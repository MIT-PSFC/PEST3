c======================================================================
c.... High-n Ballooning Stability Code
      subroutine blninit(mt,gtmax,Rext)
      include 'balmc.i'
      integer mt
      real*8 gtmax,Rext
      call defolt()
      almin = 1.e20
      mth=mt
      r=Rext
      thpi=gtmax

      call blnaddr(p,pp,q,qp,g,gp,f,fp,grpssq,xsq,grthsq,grpsth,gptdth
     &     ,xsqdps,gpsdth,xsqdth,xjacob,xjprym,delta,qdelp)

      mth1 = mth+1
      mth2 = mth1+1
      dth = twopi/mth
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 * r2
      return
      end
c======================================================================
      subroutine bln2(a,dideal,gl)
      include 'balmc.i'
      real*8 a,dideal,gl
      integer ik
      aesc	=a
      call blnsurfce()
      call funint
      call varymo
      call perfun
      call mercer
      dideal=di
      if(lmerc) then
         gl=-1.
         return
      endif
      do ik=kmin,kmax
         call varyk(ik)
         call search()
         if(ik.eq.kmin) then
            gl=almda
         else
            if(gl.gt.almda) gl=almda
         endif
         if(lcontm.and.lecon) return
      enddo
      return
      end
c======================================================================
      subroutine bln2z(ifail,a,dideal,gl)
      include 'balmc.i'
      integer ifail
      real*8 a,dideal,gl
      integer ik
      aesc	=a
      call blnsurfce()
      call funint
      call varymo
      call perfun
      call mercer
      dideal=di
      if(lmerc) then
         gl=-1.
         return
      endif
      do ik=kmin,kmax
         call varyk(ik)
         call searchz(ifail)
         if(ik.eq.kmin) then
            gl=almda
         else
            if(gl.gt.almda) gl=almda
         endif
         if(lcontm.and.lecon) return
      enddo
      return
      end
c======================================================================
      subroutine bln3()
      include 'balmc.i'
      end
c======================================================================
      subroutine defolt
      include 'balmc.i'
      zero = 0.0
      one = 1.0
      two = 2.0
      four = 4.0
      pi = dacos(-one)
      pi2 = pi * pi
      twopi = two * pi
      octpi = four * twopi
      s = one
      lconv = .false.
      lconv1 = .false.
      lchek = .false.
      ltrak = .false.
      lecon = .false.
      lp2nrm = .true.
      thpi = 16.0
      thwin = 4.0
      almda0 = 0.0
      ntry = 10
      error = 1.0e-3
      thet0 = 0.
      kmin = 1
      kmax = 11
      delk = 0.05
      iotty = 6
      return
      end
c======================================================================
      subroutine blnsurfce()
c     
c     select surface, and initialize for k-scan
c     
      include 'balmc.i'
      th0=thet0
      almda=almda0
      lmerc=.false.
      lcontm=.false.
      return
      end
c======================================================================
      subroutine varyk(ik)
c     
c     select k-value
c     
      include 'balmc.i'
      integer ik
      tk=twopi*delk*(ik-1)
c     
c     store k/twopi
c     
      return
      end
c======================================================================
      subroutine funint
c     
c     1.4   function initialization.
c.......................................................................
c     
c     redefine g from mapping so q = g/f
c     
      include 'balmc.i'
      g = r*g
      gp = r*gp
      return
      end
c======================================================================
      subroutine varymo
c     2.1   variation of parameters.
      include 'balmc.i'
      real*8 gref,gpref,s2,z1,z2
      integer i
      iref = 0
      if(iref.gt.0) go to 10
      gref = g
      gpref = gp
      iref = iref + 1
      if(dabs(s-one).lt. 1.0e-6) go to 50
 10   continue
c     
c     ..... now use the s scaling parameter...
c     
      s2 = s*s
      g = dsqrt(gref*gref-(one-s2)*r2)
      gp = gref * gpref / g
      z1 = g/gref
      z2 = (gref*gp - gpref*g) / (f*gref)
      do i=1,nths
         qdelp(i) = z1*qdelp(i) + z2*delta(i)
      enddo
 50   continue
      q = g / f
      qp = ( f*gp - g*fp ) / ( f*f )
      f2 = f*f
      g2 = g*g
      q2 = q*q
      return
      end
c     
c======================================================================
      subroutine perfun
c     
c     store periodic parts of alpha, beta, and gamma coefficients
c     
      include 'balmc.i'
      integer i
      real*8 zd,zgpgp,zxxbb,zbb,zfac,zj,zrj,zkaps,zsigma,z1,zb1,zg1
      do 100 i = 1, mth1
         zd = delta(i)
         zgpgp = grpssq(i)
         zxxbb = zgpgp + g2
         zbb = zxxbb/xsq(i)
         zfac = 1./zbb
         zj = xjacob(i)
         zrj = f*zj/xsq(i)
         zkaps = 0.5*g/zj/zxxbb*(xsqdth(i)-gpsdth(i)/zbb)
         zsigma = - gp - g*pp/zbb
         alfa2(i) = zfac*qp*qp*zgpgp
         z1 = q*grpsth(i)*zrj/zgpgp + qdelp(i)
         alfa1(i) = 2.*zfac*qp*zgpgp*z1
         alfa0(i) = zfac*( 1./xsq(i) + g2/xsq(i)/zgpgp
     .        +zgpgp*z1*z1 )
         zb1 = qp*zsigma
         beta1(i) = zb1
         beta0(i) = pp*zrj*grpsth(i)/f/zbb
         zg1 = xjprym(i)/zj*zgpgp + xsqdps(i)/xsq(i)*g2
         gama0(i) = pp/zxxbb*(g*gp + xsq(i)*pp - zg1)
     .        + 2.0*pp*zkaps*qdelp(i)
         gama0d(i) = 2.0 * pp * z1 * zkaps
 100  continue
      return
      end
c======================================================================
      subroutine mercer
c     
c     this subroutine evaluates various surface averages, and
c     calculates mercier di
      include 'balmc.i'
      integer mp,nord,i
      real*8  zyp(6),zdi(6),work(6)
      real*8 range,thmax,h,thet1,zxsq,za,zb,zj,zjp,zbsq,dinew,ditrp4
     &     ,ditrp6,almerct
      mp = 1
      ldi = .true.
c     
c     ...... placate compiler
c     
      work(1) = 0.
      zyp(1) = 0.
      range = 2.0 * pi
      istep = ( range + 0.0001 ) / (two*dth)
      istep = ( (istep+1)/2 ) * 2
      thmax = istep*two*dth*mp
      h = two*dth*mp
      nord = 6
c     .... initialize the integration ...
      thet1 = 0
      do 1 i=1,6
         zdi(i) = 0.0
         avc(i) = zero
 1    continue
c     
c     ...... perform integrations
c     
c     
c......try trapezoidal rule for the integration.
c     
      do 100 i = 1, mth
         zxsq = xsq(i)
         za = alfa2(i)
         zb = beta1(i)
         zj = xjacob(i)
         zdi(1) = zdi(1) + zj*zb/za
         zdi(2) = zdi(2) + zj/za
         zdi(3) = zdi(3) + zj*zb*zb/za
         zdi(4) = zdi(4) + zj*gama0(i) + zb
         zdi(5) = zdi(5) + zj
c.... 
c.... zdi(4) may be written as zdi(6) ...
c     
         zjp = xjprym(i)
         zbsq = ( g2 + grpssq(i) ) / zxsq
         zdi(6) = zdi(6) + pp * ( pp*zj/zbsq - zjp ) - gp*qp
 100  continue
      call deq(thet1,h,nord,avc,zyp,work,istep)
      ldi = .false.
      do 2 i = 1,6
         avc(i) = avc(i) / twopi
         zdi(i) = zdi(i) * dth / twopi
 2    continue
      di = -0.25 - avc(1)*( 1.0 + avc(1) )
     .     + avc(2)*( avc(3) + avc(4) )
      dinew = -0.25 - avc(1)*( 1.0 + avc(1) )
     .     + avc(2)*( avc(3) + avc(6) )
c     almerc = 1e-6
      almerc = - di / avc(2) / avc(5)
      ditrp4 = -0.25 - zdi(1)*( 1.0 + zdi(1) )
     .     + zdi(2)*( zdi(3) + zdi(4) )
      ditrp6 = -0.25 - zdi(1)*( 1.0 + zdi(1) )
     .     + zdi(2)*( zdi(3) + zdi(6) )
c     almerct = 1e-6
      almerct= - ditrp6 / zdi(2) / zdi(5)
c     
c.......use ditrp6, almerct for di and almerc...
c     
      di = ditrp6
      almerc = almerct
c     
c     store q, di, etc.
c     
      if(di.lt.0.) then
         return
      else
         lmerc = .true.
         return
      endif
      end
c======================================================================
      subroutine search()
c     
c     find eigenvalues by shooting method, and store results
c     
c.......avrqnd, ald, etc. are for converting kappa-t, kappa-v to kappa-s
c......and kappa-psi.
c     
      include 'balmc.i'
      real*8 yl,yr,ylp,yrp,zxmax,zymax,zth0,dyp,djump,enorm,err,etest
     &     ,alrq,al1,al2,al3,ald,erl,tkham,zx
      integer ith0,itry,lcont,ist
      integer j

      yl =1.0
      yr = 1.0
      ylp = 0.
      yrp = 0.
c     
c     ...... convert th0 to radians and round off ...
c     
      ith0 = th0*twopi/dth + .0001
      th0 = ith0*dth
      lconv1 = .true.
      do 100 itry = 1, ntry
c     
c     ...... initialize averages ...
c     
         avrqd=0.
         avrqn1=0.
         avrqn2=0.
         avrqn3=0.
         avrqnd = 0.0
         nodes = 0
         j=1
         call intgrt (j,yl,ylp)
         zxmax = xmax
         zymax = ymax
         j=-1
         call intgrt (j,yr,yrp )
         zth0 = th0/twopi - 1e-4
         if( ltrak ) th0 = zxmax
         if(ymax.gt.zymax .and. ltrak) th0 = xmax
         dyp = yrp*yl - ylp*yr
         djump = alfa*dyp/(yr*yl + 1.e-10)
         enorm = yr*yl + 1.e-10
         err = dyp / enorm
         almda = almda - djump/avrqd
c     
c     ...... eliminate continuum solutions
         if(almda.gt.almerc) almda=almerc
         if( (almda.eq.almerc) .and. (itry.gt.1) ) lcontm = .true.
         if( lcontm ) go to 31
         etest = dabs(err)
         if ( etest .lt. error ) go to 31
 100  continue
      lconv1 = .false.
 31   continue
      alrq = (avrqn1 + avrqn2 + avrqn3)/avrqd
      al1 = avrqn1/avrqd
      al2 = avrqn2/avrqd
      al3 = avrqn3/avrqd
      ald = avrqnd/avrqd
      erl = abs(almda-alrq)
      tkham = tk/twopi
      th0 = th0/twopi
c     
c     ...... store lambda  ......
c     
      if((erl/(abs(almerc)+1.e-3) .gt. 1.e-1) .and. lconv1
     .     .and. (almda.ne.almerc)) then
         write(6,38 )'a=',aesc,' q=',q
 38      format(a,1pe10.3,a,1pe10.3," large numerical error")
      endif
      lcont = 0
      if(almda.eq.almerc) then
         lcont = 1
      endif
      if(.not.lconv1) then
         write(6,36)'a=',aesc,' q=',q
 36      format(a,1pe10.3,a,1pe10.3,"failure to converge")
      endif
c     
c     ...... most unstable mode ? ......
c     
      if(almda.ge.almin) go to 150
      almin = almda
      if(.not.lchek) go to 150
      do 126 ist=1,mth1
         zx = thst(ist)
         thstm(ist) = zx
         if(zx.lt.zth0) eigstm(ist) = eigst(ist)/yl
         if(zx.gt.zth0) eigstm(ist) = eigst(ist)/yr
         if(zx.lt.zth0) wstm(ist) = wst(ist)/yl/yl
         if(zx.gt.zth0) wstm(ist) = wst(ist)/yr/yr
 126  continue
 150  continue
      return
      end
c======================================================================
      subroutine searchz(ifail)
c     
c     find eigenvalues by shooting method, and store results
c     
c.......avrqnd, ald, etc. are for converting kappa-t, kappa-v to kappa-s
c......and kappa-psi.
c     
      include 'balmc.i'
      integer ifail
      real*8 yl,yr,ylp,yrp,zxmax,zymax,zth0,dyp,djump,enorm,err,etest
     &     ,alrq,al1,al2,al3,ald,erl,tkham,zx
      integer ith0,itry,lcont,ist
      integer j

      ifail=0
      yl =1.0
      yr = 1.0
      ylp = 0.
      yrp = 0.
c     
c     ...... convert th0 to radians and round off ...
c     
      ith0 = th0*twopi/dth + .0001
      th0 = ith0*dth
      lconv1 = .true.
      do 100 itry = 1, ntry
c     
c     ...... initialize averages ...
c     
         avrqd=0.
         avrqn1=0.
         avrqn2=0.
         avrqn3=0.
         avrqnd = 0.0
         nodes = 0
         j=1
         call intgrt (j,yl,ylp)
         zxmax = xmax
         zymax = ymax
         j=-1
         call intgrt (j,yr,yrp )
         zth0 = th0/twopi - 1e-4
         if( ltrak ) th0 = zxmax
         if(ymax.gt.zymax .and. ltrak) th0 = xmax
         dyp = yrp*yl - ylp*yr
         djump = alfa*dyp/(yr*yl + 1.e-10)
         enorm = yr*yl + 1.e-10
         err = dyp / enorm
         almda = almda - djump/avrqd
c     
c     ...... eliminate continuum solutions
         if(almda.gt.almerc) almda=almerc
         if( (almda.eq.almerc) .and. (itry.gt.1) ) lcontm = .true.
         if( lcontm ) go to 31
         etest = abs(err)
         if ( etest .lt. error ) go to 31
 100  continue
      lconv1 = .false.
 31   continue
      alrq = (avrqn1 + avrqn2 + avrqn3)/avrqd
      al1 = avrqn1/avrqd
      al2 = avrqn2/avrqd
      al3 = avrqn3/avrqd
      ald = avrqnd/avrqd
      erl = abs(almda-alrq)
      tkham = tk/twopi
      th0 = th0/twopi
c
c     ...... store lambda  ......
c     
      if((erl/(abs(almerc)+1.e-3) .gt. 1.e-1) .and. lconv1
     .     .and. (almda.ne.almerc)) then
         ifail=1
C#         write(6,38)'a=',aesc,' q=',q,' ifail=',ifail
C# 38      format(a,1pe10.3,a,1pe10.3,a,i3," large numerical error")
      endif
      lcont = 0
      if(almda.eq.almerc) then
         lcont = 1
      endif
      if(.not.lconv1) then
         ifail=ifail+2
C#         write(6,36)'a=',aesc,' q=',q,' ifail=',ifail
C# 36      format(a,1pe10.3,a,1pe10.3,a,i3,"failure to converge")
      endif
c     
c     ...... most unstable mode ? ......
c     
      if(almda.ge.almin) go to 150
      almin = almda
      if(.not.lchek) go to 150
      do 126 ist=1,mth1
         zx = thst(ist)
         thstm(ist) = zx
         if(zx.lt.zth0) eigstm(ist) = eigst(ist)/yl
         if(zx.gt.zth0) eigstm(ist) = eigst(ist)/yr
         if(zx.lt.zth0) wstm(ist) = wst(ist)/yl/yl
         if(zx.gt.zth0) wstm(ist) = wst(ist)/yr/yr
 126  continue
 150  continue
      return
      end
c======================================================================
      subroutine intgrt(mp,yy1,yyp)
      include 'balmc.i'
      real*8 zy(7), zyp(7), work(7)
      real*8 range,thmax,h,thet1,zx,t,tts,zxsq,zj,zb0,zb1,zg0,gamad,za
     &     ,zb,zjp,zbsq,zn,zy2,zy11,zy1p,zypp,rqdi,rqni1,rqni2
     &     ,rqni3,rqnid,yy1,yyp,ynorm
      real *4 t4
      integer nord,i,ith0,ith,ith1
      integer mp
c     
c     ...... placate compiler
c     
      work(1) = 0.
      zyp(1) = 0.
      range = thpi * pi
      istep = ( range + 0.0001 ) / (two*dth)
      istep = ( (istep+1)/2 ) * 2
      thmax = istep*two*dth*mp
      h = two*dth*mp
      nord = 7
c     
c     .... initialize the integration ...
c     
      thet1 = th0 - thmax
      do  i=1,nord
         zy(i) = zero
      enddo
      zy(2) = .05
c     ...... perform integrations
      call deq(thet1,h,nord,zy,zyp,work,istep)
c     
c     ...... update derivative vector ......
c     
      zx = thet1
c     
c.... routine coeff, called by include 'coeff.i' in BALMSC.F
c     
      t = zx
      t4=t
      ith0 = (t+sign(1.0,t4)*1.0e-6 ) / dth
      ith = mod(ith0,mth)
      if ( ith .lt. 0 ) ith = mth + ith
      ith1 = ith + 1
      if(ldi) go to 10
      tts = t - tk
      zxsq = xsq(ith1)
      zj = xjacob(ith1)
      zb0 = beta0(ith1)
      zb1 = beta1(ith1) - beta1(1)
      zg0 = gama0(ith1)
      alfa = alfa0(ith1) + tts*( alfa1(ith1)+tts*alfa2(ith1) )
      beta = zb0 + zb1*tts
      gama = zg0 + zb1/zj
      gamad = gama0d(ith1)
      go to 20
c     ...... integrands for di ......
 10   continue
      zxsq = xsq(ith1)
      za = alfa2(ith1)
      zb = beta1(ith1)
      zj = xjacob(ith1)
      zyp(1) = zj*zb/za
      zyp(2) = zj/za
      zyp(3) = zj*zb*zb/za
      zyp(4) = zj*gama0(ith1) + zb
      zyp(5) = zj
c.... 
c.... zyp(4) may be written as zyp(6) ...
c     
      zjp = xjprym(ith1)
      zbsq = ( g2 + grpssq(ith1) ) / zxsq
      zyp(6) = pp * ( pp*zj/zbsq - zjp ) - gp*qp
 20   continue
c     
c.... routine der, called by include 'der.i' in BALMSC.F
c     
c     
      if(ldi) go to 30
      zyp(1) = ( beta*zy(1) + zy(2) ) / alfa * zj
c     
c     set normalization: pest-2 if lp2nrm = .t
c     
      zn = 1.0
      if ( .not. lp2nrm ) zn = alfa
      zy2 = - ( alfa*gama + beta**2 + almda*alfa*zn ) * zy(1)
      zy2 = zy2 - beta*zy(2)
      zyp(2) = zy2 / alfa * zj
c     
c.... routine igrand, called by include 'igrand.i' in BALMSC.F
c     ...... calculate integrands for rqd, rqn, averages ...
c     
      zy11 = zy(1)*zy(1)
      zy1p = zy(1)*zyp(1)
      zypp = zyp(1)*zyp(1)
c     kinetic energy for model of pest-2 if lp2nrm = .t
      zn = 1.0
      if ( .not. lp2nrm ) zn = alfa
      rqdi = zn*zy11*zj
      rqni1 = - (2.*tts*zy1p + zy11)*zb1
      rqni2 = -2.*zb0*zy1p - zg0*zy11*zj
      rqni3 = alfa*zypp/zj
      rqnid = zj * gamad * zy11
      zyp(3) = rqdi
      zyp(4) = rqni1
      zyp(5) = rqni2
      zyp(6) = rqni3
      zyp(7) = rqnid
 30   continue
      yy1 = zy(1)
c     
c     return invariant derivative
c     
      yyp = zyp(1) / zj
      ymax = ymax/( abs(yy1) + 1.e-10 )
      ynorm = (yy1*yy1 + 1.e-10) * mp
      avrqd = avrqd + zy(3) / ynorm
      avrqn1 = avrqn1 + zy(4) / ynorm
      avrqn2 = avrqn2 + zy(5) / ynorm
      avrqn3 = avrqn3 + zy(6) / ynorm
      avrqnd = avrqnd + zy(7) / ynorm
      return
      end
c======================================================================
      subroutine deq(zx,zzh,iin,zy,zyp,zzq,iins)
c     
c     fourth order runge-kutta integrator written by ( ? )
c     modified by r.l. dewar feb 21 1979 to reduce subroutine calls
c     warning: zyp does not return up-to-date values
      include 'balmc.i'
      real*8 zt,zx,zzqx,t,tts,zxsq,zj,zb0,zb1,zg0,gamad,za,zb,zjp,zbsq
     &     ,zn,zy2,zy11,zy1p,zypp,rqdi,rqni1,rqni2,rqni3,rqnid,zzd,zzh
     &     ,zzr,zyo,zya
      real *4 t4
      integer jjns,iins,jji,iin,jjj,ith0,ith,ith1,ist
      real*8 zza(4),zzb(4),zzc(4),zy(*),zyp(*),zzq(*)
      data zza/.5,.2928932188134525,1.707106781186548,.1666666666666
     1     667/ ,zzb/1.,.2928932188134525,1.707106781186548,.3333333333
     2     333333/ , zzc/.5,.2928932188134525,1.707106781186548,.5/
c     ...... clich1 contains common blocks needed by der
      zt = zx + 1.0
      do jjns=1,iins
         zzqx=0.0
         do jji=1,iin
            zzq(jji)=0.0
         enddo
         do jjj=1,4
c     ...... cliche der calculates derivative vector zyp
c     ...... cliche coeff computes coefficients depending
c     on zx, but not on zy, that are used in der
            if(zt.eq.zx) go to 750
            zt = zx
c     
c.... routine coeff, called by include 'coeff.i' in BALMSC.F
c     
            t = zx
            t4=t
            ith0 = ( t+sign(1.0,t4)*1.0e-6 ) / dth
            ith = mod(ith0,mth)
            if ( ith .lt. 0 ) ith = mth + ith
            ith1 = ith + 1
c     
            if(ldi) go to 10
            tts = t - tk
            zxsq = xsq(ith1)
            zj = xjacob(ith1)
            zb0 = beta0(ith1)
            zb1 = beta1(ith1) - beta1(1)
            zg0 = gama0(ith1)
            alfa = alfa0(ith1) + tts*( alfa1(ith1)+tts*alfa2(ith1) )
            beta = zb0 + zb1*tts
            gama = zg0 + zb1/zj
            gamad = gama0d(ith1)
            go to 20
c     ...... integrands for di ......
 10         continue
            zxsq = xsq(ith1)
            za = alfa2(ith1)
            zb = beta1(ith1)
            zj = xjacob(ith1)
            zyp(1) = zj*zb/za
            zyp(2) = zj/za
            zyp(3) = zj*zb*zb/za
            zyp(4) = zj*gama0(ith1) + zb
            zyp(5) = zj
c.... 
c.... zyp(4) may be written as zyp(6) ...
c     
            zjp = xjprym(ith1)
            zbsq = ( g2 + grpssq(ith1) ) / zxsq
            zyp(6) = pp * ( pp*zj/zbsq - zjp ) - gp*qp
 20         continue
 750        continue
c     
c.... routine der, called by include 'der.i' in BALMSC.F
c     
            if(ldi) go to 30
            zyp(1) = ( beta*zy(1) + zy(2) ) / alfa * zj
c     
c     set normalization: pest-2 if lp2nrm = .t
c     
            zn = 1.0
            if ( .not. lp2nrm ) zn = alfa
            zy2 = - ( alfa*gama + beta**2 + almda*alfa*zn ) * zy(1)
            zy2 = zy2 - beta*zy(2)
            zyp(2) = zy2 / alfa * zj
c     
c.... routine igrand, called by include 'igrand.i' in BALMSC.F
c     
c     
c     ...... calculate integrands for rqd, rqn, averages ...
c     
            zy11 = zy(1)*zy(1)
            zy1p = zy(1)*zyp(1)
            zypp = zyp(1)*zyp(1)
c     
c     kinetic energy for model of pest-2 if lp2nrm = .t
            zn = 1.0
            if ( .not. lp2nrm ) zn = alfa
            rqdi = zn*zy11*zj
            rqni1 = - (2.*tts*zy1p + zy11)*zb1
            rqni2 = -2.*zb0*zy1p - zg0*zy11*zj
            rqni3 = alfa*zypp/zj
            rqnid = zj * gamad * zy11
            zyp(3) = rqdi
            zyp(4) = rqni1
            zyp(5) = rqni2
            zyp(6) = rqni3
            zyp(7) = rqnid
 30         continue
c     
            do jji=1,iin
               zzd=zzh*zyp(jji)
               zzr=zzd*zza(jjj)-zzb(jjj)*zzq(jji)
               zy(jji)=zy(jji)+zzr
               zzq(jji)=zzq(jji)+3.0*zzr-zzc(jjj)*zzd
            enddo
            zzr=zzh*zza(jjj)-zzb(jjj)*zzqx
            zx=zx+zzr
            zzqx=zzqx+3.0*zzr-zzc(jjj)*zzh
         enddo
c     
c     ...... cliche check allows access to zy at each step
c     
c.... routine check, called by include 'check.i' in BALMSC.F
c     
c     ...... find maximum of zy(1)
c     
         if(jjns.ne.1) go to 760
         xmax = zx
         ymax = abs( zy(1) )
         zyo = zy(1)
 760     zya = abs (zy(1))
         if( ymax .gt. zya ) go to 761
         xmax = zx
         ymax = zya
 761     continue
c     
c     ...... count nodes ......
c     
         if(zyo*zy(1) .gt. zero) go to 765
         nodes = nodes + 1
 765     zyo = zy(1)
c     
c     ...... store eigenfunction ......
c     
         if(.not.lchek) go to 766
         ist = 1.0001 + ( (zx-th0)/thwin + pi )/dth
         if(ist.lt.1) go to 766
         if(ist.gt.mth1) go to 766
         thst(ist) = zx/twopi
         eigst(ist) = zy(1)
         wst (ist) = rqni1 + rqni2 + rqni3
 766     continue
c     
      enddo
      return
      end
c======================================================================
      subroutine fprime(far,xar,nls,nrs,idisp,xval,f,fpr)
c     this subroutine returns the function far & its 1st derivative
c     corresponding to the independent variable xar using a 4 point
c     lagrangian interpolation formula.(ref. abramowitz & stegun)
c     
c     it finds the nearest mesh point in the x array. idisp may be
c     used to use every idisp'th point in x. nl and nr define the
c     valid range for the search. when the nearest mesh point is found
c     then fpr is computed using the next 3 points to the right
c     
c     note: this routine assumes
c     a) x(nr) is greater than x(nl)
c     b) idsip = 1
c     these assumptions are temporary and will be eliminated.
c     author: j. manickam; modified rld apr 11 79
c**************************************************************************
c     
      implicit none
      real*8 far(*),xar(*),xval,f,fpr
      integer nls,nrs,idisp
      integer nl,nr,ntsts,ntst,n0,n1,n2,n3
      real*8 x0,x1,x2,x3,a012,a034,a135,a245,b0b1,b0b2,b0b3,b1b2,b1b3
     &     ,b2b3
      nl = nls
      nr = nrs
c     
c     first locate the nearest mesh point
c     
      if(xval .gt. xar(nl))go to 100
c     outside the range on lhs
      go to 700
 100  if(xval .lt. xar(nr))go to 200
c     outside the range on rhs
      nl = nr - 3 * idisp
      go to 700
c     
 200  continue
c     find the nearest mesh point
      ntsts = ( nr - nl )
      ntst = ( nr + nl ) / 2
c     check for value on a mesh point
      if(ntsts .lt. 2) go to 500
      if(xval-xar(ntst))250,400,300
c     value to the left, bring in nr
 250  nr = ntst
      go to 200
c     value to the right, bring in nl
 300  nl = ntst
      go to 200
 400  nl = ntst
 500  continue
      nl = nl - 1
      if(nl.lt.nls) nl = nls
c     now set n0
 600  if(nl+3*idisp-nrs)700,700,650
 650  nl = nl - idisp
      go to 600
c     now find fprime
 700  continue
c     
      n0 = nl
      n1 = nl + idisp
      n2 = nl + 2 * idisp
      n3 = nl + 3 * idisp
c     
      x0 = xar(n0)
      x1 = xar(n1)
      x2 = xar(n2)
      x3 = xar(n3)
      a012 = (x0-x1) * (x0-x2) * (x0-x3)
      a034 = (x0-x1) * (x1-x2) * (x1-x3)
      a135 = (x0-x2) * (x1-x2) * (x2-x3)
      a245 = (x0-x3) * (x1-x3) * (x2-x3)
c     
      b0b1 = (xval-x0) * (xval-x1)
      b0b2 = (xval-x0) * (xval-x2)
      b0b3 = (xval-x0) * (xval-x3)
      b1b2 = (xval-x1) * (xval-x2)
      b1b3 = (xval-x1) * (xval-x3)
      b2b3 = (xval-x2) * (xval-x3)
c     
      fpr = far(n0) * (b2b3 + b1b3 + b1b2) / a012
     .     - far(n1) * (b2b3 + b0b3 + b0b2) / a034
     .     + far(n2) * (b1b3 + b0b3 + b0b1) / a135
     .     - far(n3) * (b1b2 + b0b1 + b0b2) / a245
c     
      f = far(n0) * b1b2 * (xval-x3) / a012
     .     - far(n1) * b2b3 * (xval-x0) / a034
     .     + far(n2) * b0b3 * (xval-x1) / a135
     .     - far(n3) * b0b1 * (xval-x2) / a245
      return
      end

c============================================================================
       subroutine bootstrappp(ajavbs, iboot, ierr, nout,
     &  zeff,rm2avx,epsx,rmajx,qx,ftrapx,
     &  pvale,ppvale,tempxe,tempxpe,
     &  pvali,ppvali,tempxi,tempxpi )
         implicit none

         integer iboot,ierr,nout,igoof
         real*8 ajavbs,zeff,rm2avx,epsx,rmajx,qx,ftrapx,pvale,ppvale
     &        ,tempxe,tempxpe,pvali,ppvali,tempxi,tempxpi

         real*8 amu0,pi,echarge,emass,pmass,xeps0,zave,xf,alfi,denom
     &        ,anum1,anum2,tepote,tipoti,xmi,xclog,xdense,xdensi,xvthe
     &        ,xvthi,xtaue,xtaui,xnuse,xnusi,xl31,xl32,col1,col2
     &        ,alfhh,a13,b13,c13,a23,b23,c23,pepope,pipopi

c                                             6/27/94 scj
c
c     DESCRIPTION OF RETURNED PARAMETERS:
c     ajavbs   surface averaged J dot B due to bootstrap
c              divided by surface averaged B dot grad phi
c              (divide by mu0=4*pi*1.e-7 to get amperes/meter)
c     ierr     =0 for normal calculation
c              =1 for error exit
c
c     DESCRIPTION OF INPUT PARAMETERS:
c     iboot    =1 for collisionless bootstrap calculation
c              =2 to include collisonal corrections
c     nout     logical unit for error output
c
c     zeff     effective charge
c     rm2avx   surface average of 1/R**2
c     epsx     inverse aspect ratio of this flux surface
c     rmajx    major radius of this flux surface (m)
c     qx       safety factor
c     ftrapx   trapped particle fraction
c     pvale    electron pressure times mu0
c     ppvale   derivative of electron pressure wrt psi  (PF/radian)
c     pvali    ion pressure times mu0
c     ppvali   derivative of ion pressure wrt psi       (PF/radian)
c     tempxe   electron temperature (kev)
c     tempxpe  derivative of electron pressure wrt psi  (PF/radian)
c     tempxi   ion temperature      (kev)
c     tempxpi  derivative of ion temperature wrt psi    (PF/radian)
c
c
      data amu0/1.256637061D-6/
      data pi/3.1415926535D+0/
      data echarge/1.6022D-19/
      data emass/9.1095D-31/
      data pmass/1.6726D-27/
      data xeps0/8.854D-12/
      data igoof/0/
c
      ierr = 0
      if(tempxe.le.0 .or. tempxi.le.0) go to 101
      if(ftrapx.le.0 .or. ftrapx.ge.1) go to 102
      if(pvale .le.0 .or. pvali .le.0) go to 103
      if(rm2avx .le. 0) go to 104
      if(zeff .lt. 0.9) go to 105
      if(iboot.lt.1 .or. iboot.gt.2) go to 106
      zave = (pvale/tempxe)/(pvali/tempxi)
      xf=ftrapx/(1.0-ftrapx)
      alfi = 1.172/(1.0+0.462*xf)
      denom = 1.414*zeff + zeff**2 + xf*(0.754+2.657*zeff+2.*zeff**2)
     1                          + xf**2*(0.348+1.243*zeff+   zeff**2)
      anum1 = xf*(0.754 + 2.21*zeff + zeff**2
     1      + xf*(0.348 + 1.24*zeff + zeff**2) )
      anum2 = xf*(0.884 + 2.07*zeff)
      pepope=ppvale/pvale
      pipopi=ppvali/pvali
      tepote=tempxpe/tempxe
      tipoti=tempxpi/tempxi
      if(abs(iboot) .eq. 1) then
      ajavbs=-(pvale/rm2avx)*(anum1*(pepope+(pvali/pvale)
     1          *(pipopi - alfi*tipoti) ) - anum2*tepote) / denom
      return
      endif
c
      if(zeff .le. 2.0) then
      a13=1.02-0.13*(zeff-1.0)
      b13=1.07-0.45*(zeff-1.0)
      c13=1.07-0.38*(zeff-1.0)
      a23=0.57-0.05*(zeff-1.0)
      b23=0.61-0.27*(zeff-1.0)
      c23=0.61-0.23*(zeff-1.0)
      else
      a13=0.89-0.05*(zeff-2.0)
      b13=0.62-0.03*(zeff-2.0)
      c13=0.69-0.09*(zeff-2.0)
      a23=0.52-0.02*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      c23=0.38-0.05*(zeff-2.0)
      endif
      xmi=1.0
      xclog=17.0
      xdense=pvale/(amu0*tempxe*echarge*1.e3)
      xdensi=pvali/(amu0*tempxi*echarge*1.e3)
      xvthe=sqrt(2.0*echarge*1.e3*tempxe/emass)
      xvthi=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass))
      xtaue=(12.0*pi**1.5*xeps0**2*(emass)**2*xvthe**3)/
     +(4.0*xdense*zeff*(echarge)**4*xclog)
      xtaui=1.414*(12.0*pi**1.5*xeps0**2*(xmi*pmass)**2*
     +xvthi**3)/(4.0*xdensi*(zave*echarge)**4*xclog)
      xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*epsx**1.5)
      xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*epsx**1.5)
      xl31=anum1/denom
      xl32=anum2/denom
      col1=(1.0/(1.0+a13*sqrt(xnuse)+b13*xnuse))*
     +(1.0/(1.0+c13*xnuse*epsx**1.5))
      col2=(1.0/(1.0+a23*sqrt(xnuse)+b23*xnuse))*
     +(1.0/(1.0+c23*xnuse*epsx**1.5))
      alfhh=((alfi-0.35*sqrt(xnusi))/(1.0+0.7*sqrt(xnusi))
     +-2.1*xnusi**2*epsx**3)*(1.0/(1.0+xnusi**2*epsx**3))*
     +(1.0/(1.0+xnuse**2*epsx**3))
      ajavbs=-(pvale/rm2avx)*(xl31*col1*(pepope+(pvali/pvale)
     +*(pipopi-alfhh*tipoti))-(2.5*xl31*col1-(2.5*xl31-xl32)*col2)
     +*tepote)
c
      return
c
c.....error exit
  101 continue
      write(nout,1101) tempxe,tempxi
 1101 format(" Error in bootstrap, tempxe,tempxi=",1p2e12.4)
      go to 200
  102 continue
      write(nout,1102) ftrapx
 1102 format(" Error in bootstrap, ftrapx=",1pe12.4)
      go to 200
  103 continue
      write(nout,1103) pvale,pvali
 1103 format(" Error in bootstrap, pvale,pvali =",1p2e12.4)
      go to 200
  104 continue
      write(nout,1104) rm2avx
 1104 format(" Error in bootstrap, rm2avx=",1pe12.4)
      go to 200
  105 continue
      write(nout,1105) zeff
 1105 format(" Error in bootstrap, zeff=",1pe12.4)
      go to 200
  106 continue
      write(nout,1106) iboot
 1106 format(" Error in bootstrap, iboot=",i5)
c
  200 igoof=igoof+1
      if(igoof.gt.10) ierr=1
      return
      end
