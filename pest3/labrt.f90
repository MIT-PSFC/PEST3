      SUBROUTINE pstlabrt(isw,lhol,inx)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER ISW
      INTEGER INX
      INTEGER NP

!
 LOGICAL	 	 	 :: ps
 LOGICAL	 	 	 :: ts
      data np/10/,ps/.true./,ts/.false./
!
 INTEGER, DIMENSION(8) 	 :: lhol
      if((isw == 0) .OR. (isw > 5))return
      go to ( 1,2,3,4,5 ), isw
    1 if ( ps  .AND.  (np  >  0) )   write (6, 27 )   lhol, inx
   27 format(1h0,9x,8a10,3x,i4)
      np=np-1
      if ( ts )   stop 'ERROR in LABRT'
      return
    2 ps=.false.
      return
    3 ps=.true.
      np=inx
      return
    4 ts=.true.
      return
    5 ts=.false.
      return
      end


