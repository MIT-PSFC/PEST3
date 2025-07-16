      implicit none
      integer nths0,nths,nthp0,nrths
c.....nths0 is maximum number of theta points
      parameter(nths0=256)
c.....nths  is maximum array dimension for theta points
      parameter(nths=nths0+5)
c.....nthp0 is another theta dimension
      parameter(nthp0=nths0+1)
c.....nrths is another theta dimension
      parameter(nrths=nths)

      logical lconv
      logical lconv1, lchek, ldi, ltrak
      logical lmerc, lcontm, lecon,lp2nrm

      real*8
     .     grpssq(nths),xsq(nths),grpsth(nths),grthsq(nths),xsqdps(nths)
     &     ,gpsdth(nths), gptdth(nths), xjacob(nths), xjprym(nths)
     &     ,delta(nths), qdelp(nths), xsqdth(nths)

      real*8 zero,one,two,four,pi,pi2,twopi, octpi,thpi,thwin,dth,delk
      !integer*8  kdist(2)
      real*8 aesc,f,fp,p,pp,q,qp,g,gp,s
      !equivalence (kmin, kdist(1))

      common / const / zero, one, two, four, pi, pi2, twopi, octpi
      integer mth1,mth2,kmin,kmax,iotty
      common/deflti/mth1, mth2,kmin,kmax
      common/defltr/thpi,thwin,dth,delk

      common / logic /lconv, ltrak,lconv1,lchek, ldi, lp2nrm,
     .     lmerc, lcontm, lecon
      integer ntry,istep,nodes,iwords
      real*8 error,th0,thet0,almda0,tk,almda
      common / integi/ntry,istep,nodes,iwords
      common / integr /error,th0,thet0,almda0,tk,almda
      common / chnel / iotty
      common / metrc / grpssq, xsq, grpsth, grthsq, xsqdps,
     .     gpsdth, gptdth, xjacob, xjprym,
     .     delta, qdelp, xsqdth
      common/eqmap/aesc,f,fp,p,pp,q,qp,g,gp,s
      real*8  r2,r4,r6,alphap,f2,g2,q2,di,almerc,avx2,avx2p2,avp2m1
     &     ,const,r
      integer iref,mth
      common/equili/iref
      common/equilr/r2,r4,r6,alphap,f2,g2,q2
      common/mercr/ di, almerc, avx2, avx2p2, avp2m1, const
      common / equ5i /mth
      common / equ5r /r
      real*8 avc(6),avrqd,avrqn1,avrqn2,avrqn3,avrqnd
      common / avrge / avc, avrqd, avrqn1, avrqn2, avrqn3, avrqnd
      real*8 alfa0(nrths), alfa1(nrths), alfa2(nrths), beta0(nrths),
     .     beta1(nrths), gama0(nrths), gama0d(nrths)
      common / perod / alfa0, alfa1, alfa2, beta0, beta1, gama0, gama0d
      real*8 alfa,beta,gama,xmax,ymax,almin
      common / coefs / alfa, beta, gama
      common / setmx / xmax, ymax
      common / sumop /almin
      real*8 eigst(nthp0), eigstm(nthp0), thst(nthp0), thstm(nthp0),
     .     wst(nthp0), wstm(nthp0)
      common / chekst / eigst, eigstm, wst, wstm, thst, thstm

