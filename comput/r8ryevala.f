      subroutine r8ryevala(rmsplc,ymsplc,rmspls,ymspls,
     >                   zsnthtk,zcsthtk,nmom,zr99,zy99)
c
c  evaluate moments position:  asymmetric eq. moments formula
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER nmom,imul
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zr99,zy99
!============
      REAL*8 rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      REAL*8 zcsthtk(nmom),zsnthtk(nmom)
c
c--------------------------------
c
      ZR99 = RMSPLC(0)
      ZY99 = YMSPLC(0)
      DO 51 IMUL=1,NMOM
         ZR99 = ZR99 + RMSPLC(IMUL)*ZCSTHTK(IMUL)
     >      + RMSPLS(IMUL)*ZSNTHTK(IMUL)
         ZY99 = ZY99 + YMSPLC(IMUL)*ZCSTHTK(IMUL)
     >      + YMSPLS(IMUL)*ZSNTHTK(IMUL)
 51   CONTINUE
c
      return
      end
! 11Jan2003 fgtok -s r8_precision.sub misc.sub "r8con.csh conversion"
