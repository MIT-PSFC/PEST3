      subroutine ryeval(r0spl,rmspl,ymspl,zsnthtk,zcsthtk,nmom,
     >                   zr99,zy99)
c
c  evaluate moments position:  symmetric eq. moments formula
c
      real r0spl,rmspl(nmom),ymspl(nmom)
      real zcsthtk(nmom),zsnthtk(nmom)
c
c--------------------------------
c
      ZR99=R0SPL
      ZY99 = 0.
      DO 50 IMUL=1,NMOM
         ZR99=ZR99+RMSPL(IMUL)*ZCSTHTK(IMUL)
         ZY99=ZY99+YMSPL(IMUL)*ZSNTHTK(IMUL)
 50   CONTINUE
c
      return
      end
