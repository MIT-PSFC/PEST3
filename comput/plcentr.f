      SUBROUTINE PLCENTR (R,Z,N,RCENTR,ZCENTR)
 
      IMPLICIT NONE 
 
C     USE GREEN'S THEOREM TO CALCULATE CENTROIDS
C     DICK WIELAND
C
C     RGA: Jun2010, switch to centroid of inscribed polygon
C
      REAL R(*),Z(*),RCENTR,ZCENTR
      INTEGER N
 
      REAL AREA,LINT,MR,MZ
      INTEGER I
 
C     CALCULATE AREA ENCLOSED
      AREA = 0.
      MR   = 0.
      MZ   = 0.
      DO I = 1,N-1
          LINT = R(I)*Z(I+1)-R(I+1)*Z(I)
          AREA = AREA + LINT
          MR   = MR + (R(I)+R(I+1))*LINT
          MZ   = MZ + (Z(I)+Z(I+1))*LINT
      END DO
      AREA = AREA/2.

      if (abs(AREA)>1.e-10) then
         MR = MR/(6.*AREA)
         MZ = MZ/(6.*AREA)
      else
         MR=0.
         MZ=0.
         do i = 1, N
            MR = MR+R(I)
            MZ = MZ+Z(I)
         end do
         MR = MR/max(1,N)
         MZ = MZ/max(1,N)
      end if

      RCENTR = MR
      ZCENTR = MZ
 
      RETURN
      END
 
