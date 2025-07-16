!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RHOLIM
      SUBROUTINE R4_EQ_RHOLIM(
     > R4_ZRHO_AXIS,R4_ZRHO_BDY)
 
      external EQ_RHOLIM
 
! argument declarations
 ! floating type, output only:
      REAL R4_ZRHO_AXIS
 
 ! floating type, output only:
      REAL R4_ZRHO_BDY
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZRHO_AXIS_act
      REAL*8 :: ZRHO_BDY_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZRHO_AXIS_act=0
      ZRHO_BDY_act=0
! call to original routine:  EQ_RHOLIM
 
      CALL EQ_RHOLIM(
     > ZRHO_AXIS_act,ZRHO_BDY_act)
 
! copy back outputs if modified.
      R4_ZRHO_AXIS = ZRHO_AXIS_act
      R4_ZRHO_BDY = ZRHO_BDY_act
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RHOXLIM
      SUBROUTINE R4_EQ_RHOXLIM(
     > R4_ZRHO_BDY,R4_ZRHO_EXTEND)
 
      external EQ_RHOXLIM
 
! argument declarations
 ! floating type, input/output:
      REAL R4_ZRHO_BDY
 
 ! floating type, input/output:
      REAL R4_ZRHO_EXTEND
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZRHO_BDY_act
      REAL*8 :: ZRHO_BDY_ref
      REAL*8 :: ZRHO_EXTEND_act
      REAL*8 :: ZRHO_EXTEND_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZRHO_BDY_act=R4_ZRHO_BDY
      ZRHO_BDY_ref=R4_ZRHO_BDY
 
      ZRHO_EXTEND_act=R4_ZRHO_EXTEND
      ZRHO_EXTEND_ref=R4_ZRHO_EXTEND
 
! call to original routine:  EQ_RHOXLIM
 
      CALL EQ_RHOXLIM(
     > ZRHO_BDY_act,ZRHO_EXTEND_act)
 
! copy back outputs if modified.
 
      if(ZRHO_BDY_act.ne.ZRHO_BDY_ref) then
         R4_ZRHO_BDY = ZRHO_BDY_act
      endif
 
      if(ZRHO_EXTEND_act.ne.ZRHO_EXTEND_ref) then
         R4_ZRHO_EXTEND = ZRHO_EXTEND_act
      endif
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RHOLIMA
      SUBROUTINE R4_EQ_RHOLIMA(
     > R4_ZRHO_AXIS,R4_ZRHO_EXTEND)
 
      external EQ_RHOLIMA
 
! argument declarations
 ! floating type, input/output:
      REAL R4_ZRHO_AXIS
 
 ! floating type, input/output:
      REAL R4_ZRHO_EXTEND
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZRHO_AXIS_act
      REAL*8 :: ZRHO_AXIS_ref
      REAL*8 :: ZRHO_EXTEND_act
      REAL*8 :: ZRHO_EXTEND_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZRHO_AXIS_act=R4_ZRHO_AXIS
      ZRHO_AXIS_ref=R4_ZRHO_AXIS
 
      ZRHO_EXTEND_act=R4_ZRHO_EXTEND
      ZRHO_EXTEND_ref=R4_ZRHO_EXTEND
 
! call to original routine:  EQ_RHOLIMA
 
      CALL EQ_RHOLIMA(
     > ZRHO_AXIS_act,ZRHO_EXTEND_act)
 
! copy back outputs if modified.
 
      if(ZRHO_AXIS_act.ne.ZRHO_AXIS_ref) then
         R4_ZRHO_AXIS = ZRHO_AXIS_act
      endif
 
      if(ZRHO_EXTEND_act.ne.ZRHO_EXTEND_ref) then
         R4_ZRHO_EXTEND = ZRHO_EXTEND_act
      endif
 
! exit
      return
      end
