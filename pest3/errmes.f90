!......................................................................
      SUBROUTINE psterrmes ( nout,mesage,is )
!......................................................................
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NOUT
      INTEGER IS


 INTEGER, DIMENSION(2) 	 :: mesage
      write(itty,100)mesage,is
      write(nout,100)mesage,is
  100 format(" error in subroutine",1x,2a5/,1x," error # ",i3 )
      stop 'ERROR in errmes'
      end


