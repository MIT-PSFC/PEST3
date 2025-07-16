      INTEGER FUNCTION LKUPR(X,TABLE,N)
C
C Look up a real number X in a TABLE of length N.
C returns the value lkupr such that:
C
C	table(lkupr) <= x < table(lkupr+1)
C
C assumes that table is in increasing order and does a binary search.
C
      REAL X,TABLE(N)
C
C========================================UPPER CASE & PRETTY PRINT RTM MAY 1988
C
      IF(X .LT. TABLE(1))THEN
          LKUPR=0
          GO TO 70000		!DONE
      ENDIF
C
      IF(X .GE. TABLE(N))THEN
          LKUPR=N
          GO TO 70000		!DONE
      ENDIF
C
      N1=1
      N2=N
C
100   CONTINUE	!LOOK AGAIN LOOP
C
      N3=(N2+N1)/2
C
      IF(X .GE. TABLE(N3))THEN
          N1=N3
      ELSE
          N2=N3
      ENDIF
C
      IF(N2 .EQ. N1+1)THEN
          LKUPR=N1
          GO TO 70000	    !DONE
      ELSE
          GOTO 100	    !LOOK AGAIN
      ENDIF
 
70000 CONTINUE	!All returns from here
C
      RETURN
      END
