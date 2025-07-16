!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_TIME
      SUBROUTINE R4_EQ_TIME(
     > R4_ZTIME)
 
      external EQ_TIME
 
! argument declarations
 ! floating type, input/output:
      REAL R4_ZTIME
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZTIME_act
      REAL*8 :: ZTIME_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZTIME_act=R4_ZTIME
      ZTIME_ref=R4_ZTIME
 
! call to original routine:  EQ_TIME
 
      CALL EQ_TIME(
     > ZTIME_act)
 
! copy back outputs if modified.
 
      if(ZTIME_act.ne.ZTIME_ref) then
         R4_ZTIME = ZTIME_act
      endif
 
! exit
      return
      end
