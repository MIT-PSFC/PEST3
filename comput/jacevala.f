      subroutine jacevala(DRMDROC,DYMDROC,DRMDROS,DYMDROS,
     >                    RMSPLC,YMSPLC,RMSPLS,YMSPLS,
     >                    ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
c
c  evaluate 2x2 jacobian using moments and derivatives
c    ** asymmetric formula **
c
c  evaluate moments position:  asymmetric eq. moments formula
c
      real rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      real drmdroc(0:nmom),dymdroc(0:nmom)
      real drmdros(0:nmom),dymdros(0:nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
c
      real zjac(2,2)
c
c--------------------------------
c
      ZJAC(1,1)= DRMDROC(0)
      ZJAC(1,2)= 0.0
C
      ZJAC(2,1)= DYMDROC(0)
      ZJAC(2,2)= 0.0
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
 
