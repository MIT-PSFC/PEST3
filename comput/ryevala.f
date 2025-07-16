      subroutine ryevala(rmsplc,ymsplc,rmspls,ymspls,
     >                   zsnthtk,zcsthtk,nmom,zr99,zy99)
c
c  evaluate moments position:  asymmetric eq. moments formula
c
      real rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
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
