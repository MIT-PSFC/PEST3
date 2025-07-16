!......................................................................
      SUBROUTINE pstempty ( iarg )
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IARG

!
!      simulates empty buffer on ltss...
!
      close ( unit = iarg )
!
      return
      end


