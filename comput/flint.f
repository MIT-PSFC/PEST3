C******************** START FILE FLINT.FOR ; GROUP SIGS2 ******************
C================================================================
C
C  FLINT
C   INTERPOLATE ON A 2D ARRAY
C
      REAL FUNCTION FLINT(NOUT,ARRAY,X,Y,IDENT,
     >   XMIN,XMAX,XLR,NX,YMIN,YMAX,YLR,NY)
C
      DIMENSION ARRAY(NX,NY)
c
c Interpolates the table ARRAY(NX,NY) to find the value at (x,y).
c If XLR>0, then the X grid is equally spaced on a logarithmic scale:
c
c   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
c
c   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
c
c For XLR<=0, the X grid is equally spaced on a linear scale from
c XMIN to XMAX.
c
c For YLR>0, Y is logarithmically spaced from YMIN to YMIN*exp(YLR)
c and for YLR<=0, Y is linearly spaced from YMIN to YMAX.
c
c NOUT is the fortran unit to write error messages to.
c IDENT is a code identifying the table ARRAY in error messages.
C
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
C
      external zcflnt                   ! block data initialization
C
C
      ZX=X
      ZY=Y
C
      if(zx .lt. xmin) then
          zzx=0.0
      else
 
      IF(XLR.GT.0.0) THEN
        ZLOGX=ALOG(ZX/XMIN)
        ZZX=1.0+ZLOGX/XLR * (NX-1)
      ELSE
        ZZX=1.0+(ZX-XMIN)/(XMAX-XMIN) * (NX-1)
      ENDIF
 
      endif
C
      if(zy .lt. ymin) then
          zzy=0.0
      else
 
      IF(YLR.GT.0.0) THEN
        ZLOGY=ALOG(ZY/YMIN)
        ZZY=1.0+ZLOGY/YLR * (NY-1)
      ELSE
        ZZY=1.0+(ZY-YMIN)/(YMAX-YMIN) * (NY-1)
      ENDIF
 
      endif
C
      IX=ZZX
      IY=ZZY
C  CHECK FOR POINT OUT OF BOUNDS OF ARRAY
      IF(IX.GE.1) GO TO 10
C
      IPXL=IPXL-1
      IF(IPXL.GE.0) WRITE(NOUT, 9001) IDENT,X,Y
 9001 FORMAT(/////' *************************************'/
     >  '  ATTEMPT TO INTERPOLATE OUT OF BOUNDS    IDENT=',I8/
     >  '  FLINT SIGMA*V TABLE LOOKUP, STANDARD FIXUP TAKEN'/
     >  '  ARGUMENTS:  X=',1PE10.3,'  Y=',1PE10.3/
     >  ' **************************************'//////)
      IX=1
      ZZX=1.0
C
 10   CONTINUE
      IF(IY.GE.1) GO TO 20
      IPYL=IPYL-1
      IF(IPYL.GE.0) WRITE(NOUT, 9001) IDENT,X,Y
      IY=1
      ZZY=1.0
C
 20   CONTINUE
      IF(IX.LT.NX) GO TO 30
      IX=NX-1
      ZZX=NX
C
      IF(ZX.GT.XMAX) THEN
        write(nout,9001) ident,x,y
        CALL BAD_EXIT  ! DEEMED FATAL - DMC 1988
      ENDIF
C
 30   CONTINUE
      IF(IY.LT.NY) GO TO 40
      IY=NY-1
      ZZY=NY
C
      IF(ZY.GT.YMAX) THEN
        write(nout,9001) ident,x,y
        CALL BAD_EXIT  ! DEEMED FATAL - DMC 1988
      ENDIF
C
C---------------------------
C  OK
C
 40   CONTINUE
      IXP1=IX+1
      IYP1=IY+1
      ZZX=ZZX-IX
      ZZY=ZZY-IY
      ZF00=(1.-ZZX)*(1.-ZZY)
      ZF01=(1.-ZZX)*ZZY
      ZF10=ZZX*(1.-ZZY)
      ZF11=ZZX*ZZY
C
C  INTERPOLATE
C
      FLINT=ZF00*ARRAY(IX,IY)+ZF01*ARRAY(IX,IYP1)+
     >        ZF10*ARRAY(IXP1,IY)+ZF11*ARRAY(IXP1,IYP1)
C
      RETURN
      	END
C
C******************** END FILE FLINT.FOR ; GROUP SIGS2 ******************
