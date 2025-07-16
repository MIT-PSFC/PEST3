!-----------------------------------------------------------------------
 SUBROUTINE pstmetric
!
! fetch profile and metric from disk ....
!
! aplet  14.2_r8 .97
!
!.......................................................................
!
 USE pstcom
 USE newmet
 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NADRES
      INTEGER LGIVUP
      INTEGER LENGTH
      INTEGER LADR2
      real*8 CONST
      real*8 GAS
      real*8 QANEW
      real*8 QPANEW
      real*8 B0SC
      real*8 TOP
      real*8 TOPS
      real*8 BOTTOM
      real*8 BOTTM2
      real*8 TOTCUR
      real*8 ZVOL
      real*8 ZPSURF
      real*8 PSQSUM
      real*8 BP2
      INTEGER NSURF
      real*8 ZHL
      real*8 ZHR
      real*8 ZDPSI
      real*8 B2SUM
      real*8 VTSUM
      real*8 SSUM
      real*8 SJPHI
      real*8 DSUBR
      real*8 DSUBEE
      real*8 BP2SUM
      real*8 RGRPS
      real*8 ZPAV
      real*8 ZB2AV
      real*8 BETA1
      real*8 BETA3
      real*8 ARAT
      real*8 XRAD
      real*8 ELLIP
      real*8 ZDELTA
      real*8 ZZMIN
      real*8 ZZMAX
      real*8 ZXMIN
      real*8 ZXMAX
      real*8 XIN
      real*8 XOUT
      real*8 ZARAD
      real*8 ZR0
      real*8 EPSILN
      real*8 ZAVERA
      real*8 ZELONG
      real*8 EBETAP
      real*8 BETAPL
      real*8 ZBETAT
      real*8 BETAST
      real*8 BETATR
      real*8 PPF
      real*8 ELI
      real*8 ZASPCT
      real*8 ZQSPEC
      real*8 ZSSPEC
      real*8 ZCUR
      real*8 ZFACT
      real*8 ZZ0
      INTEGER I
      INTEGER JNSTTP
      real*8 ZPP
      real*8 ZGGP

!
pa =  0._r8 
ppa =  0._r8 
qa =  0._r8 
qpa =  0._r8 
ga =  0._r8 
gpa =  0._r8 
fa =  0._r8 
fpa =  0._r8 
!
      nadres = 50
      CALL pstzrd(outmp1,pa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,ppa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,qa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,qpa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,ga(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,gpa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,fa(1),nosurf,nadres,lgivup,7000)
      nadres = nadres + nosurf
      CALL pstzrd(outmp1,fpa(1),nosurf,nadres,lgivup,7000)

      if ( lpstps ) then
         nadres = nadres + nosurf
         CALL pstzrd(outmp1,di(1),nosurf,nadres,lgivup,7000)
      end if

      length = nths * nsf
      ladr2 = 50 + 2*mth2 + nosurf + 8*length
      lgivup = 1
!
xa =  0._r8 
za =  0._r8 
xdth =  0._r8 
zdth =  0._r8 
xdps =  0._r8 
zdps =  0._r8 
!
      nadres =  50 + 10*nosurf
if ( lpstps )  nadres = 50 + 8*nosurf

      CALL pstzrd(outmp1,xa(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(outmp1,za(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(outmp1,xdth(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(outmp1,zdth(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(outmp1,xdps(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(outmp1,zdps(1,1),length,nadres,lgivup,7000)
!
xinf =  0._r8 
zinf =  0._r8 
grpssq =  0._r8 
psibig =  0._r8 
xsq =  0._r8 
grpsth =  0._r8 
xsqdps =  0._r8 
xjacob =  0._r8 
xjprym =  0._r8 
delta  =  0._r8 
qdelp  =  0._r8 
!
      nadres = 50
      CALL pstzrd(iomode,xinf(1),mth2,nadres,lgivup,7000)
      nadres = nadres + mth2
      CALL pstzrd(iomode,zinf(1),mth2,nadres,lgivup,7000)
      nadres = nadres + mth2
      CALL pstzrd(iomode,psibig(1),nosurf,nadres,lgivup,7000)
!
      nadres = 50 + 2*mth2 + nosurf
      CALL pstzrd(iomode,grpssq(1,1),length,nadres,lgivup,7000)
      nadres = 50 + 2*mth2 + nosurf + length
      CALL pstzrd(iomode,xsq(1,1),length,nadres,lgivup,7000)
      nadres = 50 + 2*mth2 + nosurf + 3*length
      CALL pstzrd(iomode,grpsth(1,1),length,nadres,lgivup,7000)
      nadres = 50 + 2*mth2 + nosurf + 5*length
      CALL pstzrd(iomode,xsqdps(1,1),length,nadres,lgivup,7000)

      nadres = 50 + 2*mth2 + nosurf + 8*length
      CALL pstzrd(iomode,xjacob(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(iomode,xjprym(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(iomode,delta(1,1),length,nadres,lgivup,7000)
      nadres = nadres + length
      CALL pstzrd(iomode,qdelp(1,1),length,nadres,lgivup,7000)
!
end SUBROUTINE pstmetric
