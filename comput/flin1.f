C******************** START FILE FLIN1.FOR ; GROUP SIGS2 ******************
C================================================================
C
C  FLIN1
C   INTERPOLATE ON A 1D ARRAY
C
      REAL FUNCTION FLIN1(NOUT,ARRAY,X,IDENT,XMIN,XMAX,XLR,NX)
c
c Interpolates the table ARRAY(NX) to find the value at x.
c If XLR>0, then the X grid is equally spaced on a logarithmic scale:
c
c   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
c
c   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
c
c For XLR<=0, the X grid is equally spaced on a linear scale from
c XMIN to XMAX.
c
c NOUT is the fortran unit to write error messages to.
c IDENT is a code identifying the table ARRAY in error messages.
C
 
C
      DIMENSION ARRAY(NX)
C
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
C
      external zcflnt                   ! block data initialization
C
C  INTERPOLATE TO X  LOGARITHMIC SPACING IN X
C     LINEAR IF XLR .LE. 0.0
C
C
      ZX=X
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
      IX=ZZX
C
C  CHECK FOR POINT OUT OF BOUNDS OF ARRAY
      IF(IX.GE.1) GO TO 10
C
      IPXL=IPXL-1
      IF(IPXL.GE.0) WRITE(NOUT, 9001) IDENT,X
 9001 FORMAT(///' *************************************'/
     >  '  ATTEMPT TO INTERPOLATE OUT OF BOUNDS    IDENT=',I8/
     >  '  FLIN1 SIGMA*V TABLE LOOKUP, STANDARD FIXUP TAKEN'/
     >  '  ARGUMENT:  X=',1PE10.3/
     >  ' **************************************'///)
      IX=1
      ZZX=1.0
C
 10   CONTINUE
C
      IF(IX.LT.NX) GO TO 30
      IX=NX-1
      ZZX=NX
C
      IF(ZX.GT.XMAX) THEN
        write(nout,9001) ident,x
        call bad_exit
      ENDIF
C
C---------------------------
C  OK
C
 30   CONTINUE
C
      IXP1=IX+1
      ZZX=ZZX-IX
C
C
C  INTERPOLATE
C
      FLIN1=(1.0-ZZX)*ARRAY(IX)+ZZX*ARRAY(IXP1)
C
      RETURN
      	END
C
C******************** END FILE FLIN1.FOR ; GROUP SIGS2 ******************
C-------------------------
C  FLINT LOCAL PERM. MEMORY
C   (CF SIGS2.FOR)
C
      BLOCK DATA ZCFLNT
C
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
C
      DATA ZX/0.0/
      DATA ZY/0.0/
      DATA JDENT/0/
C
      DATA IPXL/1/
      DATA IPYL/1/
C
      DATA IPXH/10/
      DATA IPYH/10/
C
      END
