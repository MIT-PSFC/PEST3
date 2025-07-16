!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RZLIMS
      SUBROUTINE R4_EQ_RZLIMS(
     > R4_RMIN,R4_RMAX,R4_ZMIN,R4_ZMAX,IERR)
 
      external EQ_RZLIMS
 
! argument declarations
 ! floating type, input/output:
      REAL R4_RMIN
 
 ! floating type, input/output:
      REAL R4_RMAX
 
 ! floating type, input/output:
      REAL R4_ZMIN
 
 ! floating type, input/output:
      REAL R4_ZMAX
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: RMIN_act
      REAL*8 :: RMIN_ref
      REAL*8 :: RMAX_act
      REAL*8 :: RMAX_ref
      REAL*8 :: ZMIN_act
      REAL*8 :: ZMIN_ref
      REAL*8 :: ZMAX_act
      REAL*8 :: ZMAX_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      RMIN_act=R4_RMIN
      RMIN_ref=R4_RMIN
 
      RMAX_act=R4_RMAX
      RMAX_ref=R4_RMAX
 
      ZMIN_act=R4_ZMIN
      ZMIN_ref=R4_ZMIN
 
      ZMAX_act=R4_ZMAX
      ZMAX_ref=R4_ZMAX
 
! call to original routine:  EQ_RZLIMS
 
      CALL EQ_RZLIMS(
     > RMIN_act,RMAX_act,ZMIN_act,ZMAX_act,IERR)
 
! copy back outputs if modified.
 
      if(RMIN_act.ne.RMIN_ref) then
         R4_RMIN = RMIN_act
      endif
 
      if(RMAX_act.ne.RMAX_ref) then
         R4_RMAX = RMAX_act
      endif
 
      if(ZMIN_act.ne.ZMIN_ref) then
         R4_ZMIN = ZMIN_act
      endif
 
      if(ZMAX_act.ne.ZMAX_ref) then
         R4_ZMAX = ZMAX_act
      endif
 
! exit
      return
      end
