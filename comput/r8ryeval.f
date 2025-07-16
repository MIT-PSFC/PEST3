      subroutine r8ryeval(r0spl,rmspl,ymspl,zsnthtk,zcsthtk,nmom,
     >                   zr99,zy99)
c
c  evaluate moments position:  symmetric eq. moments formula
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER nmom,imul
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zr99,zy99
!============
      REAL*8 r0spl,rmspl(nmom),ymspl(nmom)
      REAL*8 zcsthtk(nmom),zsnthtk(nmom)
c
c--------------------------------
c
      ZR99=R0SPL
      ZY99 = 0.D0
      DO 50 IMUL=1,NMOM
         ZR99=ZR99+RMSPL(IMUL)*ZCSTHTK(IMUL)
         ZY99=ZY99+YMSPL(IMUL)*ZSNTHTK(IMUL)
 50   CONTINUE
c
      return
      end
! 11Jan2003 fgtok -s r8_precision.sub misc.sub "r8con.csh conversion"
