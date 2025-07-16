C******************** START FILE ISCMP0.FOR ; GROUP URDUTS ******************
C-------------
C&  ISCMP0 - ISCOMP FUNCTION WITHOUT USERIO INCLUDE COMMON BLOCKS
C
      INTEGER FUNCTION ISCMP0(DBUFF,ABLIS,N,NMAX,ND2)
C
      implicit NONE
C
      integer LMAX
      PARAMETER (LMAX=512)
C
      integer N,NMAX,ND2
      CHARACTER*(*) ABLIS(N)
      CHARACTER*1 DBUFF(LMAX)
C
C-----------------------------------
C
      CHARACTER*1 ZBUFF(LMAX),MMINUS,MPLUS,MBLANK
C
      integer imax,imaxp1,imul,i,ij,id,iz
      integer num,nd1,numz,nz1,nz2,ilena
C
      DATA MBLANK/' '/
      DATA MMINUS/'-'/
      DATA MPLUS/'+'/
C
C-----------------------------------
      ISCMP0=0
      IMAX=MIN(NMAX,LMAX)
      IMAXP1=MIN(NMAX+1,LMAX)
      ilena=min(lmax,len(ablis(1)))
      IMUL=1
C
      CALL NUMFLD(DBUFF,ND1,ND2,NUM)
      IF(ND1.LE.0) GO TO 5
      IF((DBUFF(ND1).NE.MPLUS).AND.(DBUFF(ND1).NE.MMINUS)) GO TO 5
      IF(DBUFF(ND1).EQ.MMINUS) IMUL=-1
      ND1=ND1+1
      NUM=NUM-1
 5    CONTINUE
      IF(NUM.LE.0) RETURN
C  DECODE LIST ELEMENTS AND COMPARE
      DO 20 I=1,N
          DO 10 IJ=1,IMAXP1
             ZBUFF(IJ)=MBLANK
 10       CONTINUE
C
          do ij=1,ilena
             zbuff(ij)=ablis(i)(ij:ij)
          enddo
C
12        continue
          CALL FLDLIM(ZBUFF,NZ1,NZ2,NUMZ)
C
          IF(NUM.NE.NUMZ) GO TO 20
          DO 15 IJ=1,NUM
             ID=ND1+IJ-1
             IZ=NZ1+IJ-1
             IF(ZBUFF(IZ).NE.DBUFF(ID)) GO TO 20
 15       CONTINUE
          GO TO 50
 20   CONTINUE
C
C  NOT FOUND
C
      RETURN
C
C  FOUND
C
 50   CONTINUE
      ISCMP0=IMUL*I
      RETURN
      	END
C******************** END FILE ISCMP0.FOR ; GROUP URDUTS ******************
