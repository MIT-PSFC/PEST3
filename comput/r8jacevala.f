      subroutine r8jacevala(DRMDROC,DYMDROC,DRMDROS,DYMDROS,
     >                    RMSPLC,YMSPLC,RMSPLS,YMSPLS,
     >                    ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
c
c  evaluate 2x2 jacobian using moments and derivatives
c    ** asymmetric formula **
c
c  evaluate moments position:  asymmetric eq. moments formula
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER nmom,imul
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zc,zs
!============
      REAL*8 rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      REAL*8 drmdroc(0:nmom),dymdroc(0:nmom)
      REAL*8 drmdros(0:nmom),dymdros(0:nmom)
      REAL*8 zcsthtk(nmom),zsnthtk(nmom)
c
      REAL*8 zjac(2,2)
c
c--------------------------------
c
      ZJAC(1,1)= DRMDROC(0)
      ZJAC(1,2)= 0.0D0
C
      ZJAC(2,1)= DYMDROC(0)
      ZJAC(2,2)= 0.0D0
C
      DO 151 IMUL=1,NMOM
C
         ZC = ZCSTHTK(IMUL)
         ZS = ZSNTHTK(IMUL)
C
         ZJAC(1,1) = ZJAC(1,1) + DRMDROC(IMUL)*ZC
     >      + DRMDROS(IMUL)*ZS
         ZJAC(1,2) = ZJAC(1,2) - RMSPLC(IMUL)*IMUL*ZS
     >      + RMSPLS(IMUL)*IMUL*ZC
C
         ZJAC(2,1) = ZJAC(2,1) + DYMDROC(IMUL)*ZC
     >      + DYMDROS(IMUL)*ZS
         ZJAC(2,2) = ZJAC(2,2) - YMSPLC(IMUL)*IMUL*ZS
     >      + YMSPLS(IMUL)*IMUL*ZC
C
 151  CONTINUE
C
      return
      end
 
! 11Jan2003 fgtok -s r8_precision.sub misc.sub "r8con.csh conversion"
