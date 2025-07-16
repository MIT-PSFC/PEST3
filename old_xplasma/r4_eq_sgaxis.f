!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_SGAXIS
      SUBROUTINE R4_EQ_SGAXIS(
     > R4_B_AXIS,R4_R_AXIS,R4_Z_AXIS)
 
      external EQ_SGAXIS
 
! argument declarations
 ! floating type, input/output:
      REAL R4_B_AXIS
 
 ! floating type, input/output:
      REAL R4_R_AXIS
 
 ! floating type, input/output:
      REAL R4_Z_AXIS
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: B_AXIS_act
      REAL*8 :: B_AXIS_ref
      REAL*8 :: R_AXIS_act
      REAL*8 :: R_AXIS_ref
      REAL*8 :: Z_AXIS_act
      REAL*8 :: Z_AXIS_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      B_AXIS_act=R4_B_AXIS
      B_AXIS_ref=R4_B_AXIS
 
      R_AXIS_act=R4_R_AXIS
      R_AXIS_ref=R4_R_AXIS
 
      Z_AXIS_act=R4_Z_AXIS
      Z_AXIS_ref=R4_Z_AXIS
 
! call to original routine:  EQ_SGAXIS
 
      CALL EQ_SGAXIS(
     > B_AXIS_act,R_AXIS_act,Z_AXIS_act)
 
! copy back outputs if modified.
 
      if(B_AXIS_act.ne.B_AXIS_ref) then
         R4_B_AXIS = B_AXIS_act
      endif
 
      if(R_AXIS_act.ne.R_AXIS_ref) then
         R4_R_AXIS = R_AXIS_act
      endif
 
      if(Z_AXIS_act.ne.Z_AXIS_ref) then
         R4_Z_AXIS = Z_AXIS_act
      endif
 
! exit
      return
      end
