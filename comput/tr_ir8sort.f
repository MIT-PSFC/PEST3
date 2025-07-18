*DECK TR_IR8SORT
      SUBROUTINE TR_IR8SORT (IDX, DY, N, KFLAG)
C***BEGIN PROLOGUE  TR_IR8SORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
c
c  >> modified version (dmc) -- sort by integer key; carry REAL*8 DY array
c
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      DOUBLE PRECISION (SSORT-S, TR_DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   TR_IR8SORT sorts array IDX and optionally makes the same interchanges in
C   array DY.  The array IDX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IDX - array of values to be sorted   (usually abscissas)
C      DY - array to be (optionally) carried along
C      N  - number of values in array IDX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IDX in increasing order and carry DY along.
C            =  1  means sort IDX in increasing order (ignoring DY)
C            = -1  means sort IDX in decreasing order (ignoring DY)
C            = -2  means sort IDX in decreasing order and carry DY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to DX,DY; changed
C           code to parallel SSORT. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  TR_IR8SORT
C     .. Scalar Arguments ..
c
      implicit NONE
c
      INTEGER KFLAG, N
C     .. Array Arguments ..
      INTEGER IDX(*)
      REAL*8 DY(*)
C     .. Local Scalars ..
      REAL*8 R, TTY, TY
      integer T, TT
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  TR_IR8SORT
      NN = N
      IF (NN .LT. 1) THEN
         write(6,*) ' ?tr_ir8sort:  n.gt.0 expected, n = ',n
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         write(6,*) ' ?tr_ir8sort: kflag = ',kflag,' not in {-2,-1,1,2}'
         RETURN
      ENDIF
C
C     Alter array IDX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IDX(I) = -IDX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort IDX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375D0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IDX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IDX(I) .GT. T) THEN
         IDX(IJ) = IDX(I)
         IDX(I) = T
         T = IDX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (IDX(J) .LT. T) THEN
         IDX(IJ) = IDX(J)
         IDX(J) = T
         T = IDX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IDX(I) .GT. T) THEN
            IDX(IJ) = IDX(I)
            IDX(I) = T
            T = IDX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (IDX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (IDX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IDX(L)
         IDX(L) = IDX(K)
         IDX(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = IDX(I+1)
      IF (IDX(I) .LE. T) GO TO 80
      K = I
C
   90 IDX(K+1) = IDX(K)
      K = K-1
      IF (T .LT. IDX(K)) GO TO 90
      IDX(K+1) = T
      GO TO 80
C
C     Sort IDX and carry DY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IDX(IJ)
      TY = DY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IDX(I) .GT. T) THEN
         IDX(IJ) = IDX(I)
         IDX(I) = T
         T = IDX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IDX(J) .LT. T) THEN
         IDX(IJ) = IDX(J)
         IDX(J) = T
         T = IDX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IDX(I) .GT. T) THEN
            IDX(IJ) = IDX(I)
            IDX(I) = T
            T = IDX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IDX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IDX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IDX(L)
         IDX(L) = IDX(K)
         IDX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = IDX(I+1)
      TY = DY(I+1)
      IF (IDX(I) .LE. T) GO TO 170
      K = I
C
  180 IDX(K+1) = IDX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. IDX(K)) GO TO 180
      IDX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IDX(I) = -IDX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
