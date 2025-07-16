      subroutine r8jaceval(DR0DRO,DRMDRO,DYMDRO,RMSPL,YMSPL,
     >                    ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
c
c  evaluate 2x2 jacobian using moments and derivatives
c    ** symmetric formula **
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER nmom,imul
!============
      REAL*8 rmspl(nmom),ymspl(nmom)
      REAL*8 dr0dro,drmdro(nmom),dymdro(nmom)
      REAL*8 zcsthtk(nmom),zsnthtk(nmom)
c
      REAL*8 zjac(2,2)
c
c--------------------------------
c
      ZJAC(1,1)=DR0DRO
      ZJAC(1,2)=0.0D0
C
      ZJAC(2,1)=0.0D0
      ZJAC(2,2)=0.0D0
C
      DO 150 IMUL=1,NMOM
C
         ZJAC(1,1)=ZJAC(1,1)+DRMDRO(IMUL)*ZCSTHTK(IMUL)
         ZJAC(1,2)=ZJAC(1,2)-RMSPL(IMUL)*IMUL*ZSNTHTK(IMUL)
C
         ZJAC(2,1)=ZJAC(2,1)+DYMDRO(IMUL)*ZSNTHTK(IMUL)
         ZJAC(2,2)=ZJAC(2,2)+YMSPL(IMUL)*IMUL*ZCSTHTK(IMUL)
C
 150  CONTINUE
c
      return
      end
 
! 11Jan2003 fgtok -s r8_precision.sub misc.sub "r8con.csh conversion"
