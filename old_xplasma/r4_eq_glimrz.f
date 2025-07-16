!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GLIMRZ
      SUBROUTINE R4_EQ_GLIMRZ(
     > R4_RHO,R4_R_MIN,R4_R_MAX,R4_Z_MIN,R4_Z_MAX,IERR)
 
      external EQ_GLIMRZ
 
! argument declarations
 ! floating type, input only:
      REAL R4_RHO
 
 ! floating type, output only:
      REAL R4_R_MIN
 
 ! floating type, output only:
      REAL R4_R_MAX
 
 ! floating type, output only:
      REAL R4_Z_MIN
 
 ! floating type, output only:
      REAL R4_Z_MAX
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: RHO_act
      REAL*8 :: R_MIN_act
      REAL*8 :: R_MAX_act
      REAL*8 :: Z_MIN_act
      REAL*8 :: Z_MAX_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      RHO_act=R4_RHO
      R_MIN_act=0
      R_MAX_act=0
      Z_MIN_act=0
      Z_MAX_act=0
! call to original routine:  EQ_GLIMRZ
 
      CALL EQ_GLIMRZ(
     > RHO_act,R_MIN_act,R_MAX_act,Z_MIN_act,Z_MAX_act,IERR)
 
! copy back outputs if modified.
      R4_R_MIN = R_MIN_act
      R4_R_MAX = R_MAX_act
      R4_Z_MIN = Z_MIN_act
      R4_Z_MAX = Z_MAX_act
 
! exit
      return
      end
