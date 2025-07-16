!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GLIMB
      SUBROUTINE R4_EQ_GLIMB(
     > R4_RHO,R4_B_MIN,R4_B_MAX,IERR)
 
      external EQ_GLIMB
 
! argument declarations
 ! floating type, input only:
      REAL R4_RHO
 
 ! floating type, output only:
      REAL R4_B_MIN
 
 ! floating type, output only:
      REAL R4_B_MAX
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: RHO_act
      REAL*8 :: B_MIN_act
      REAL*8 :: B_MAX_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      RHO_act=R4_RHO
      B_MIN_act=0
      B_MAX_act=0
! call to original routine:  EQ_GLIMB
 
      CALL EQ_GLIMB(
     > RHO_act,B_MIN_act,B_MAX_act,IERR)
 
! copy back outputs if modified.
      R4_B_MIN = B_MIN_act
      R4_B_MAX = B_MAX_act
 
! exit
      return
      end
