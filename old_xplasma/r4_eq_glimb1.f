!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GLIMB1
      SUBROUTINE R4_EQ_GLIMB1(
     > R4_RHO,R4_PHI,R4_B_MIN,R4_CHI_B_MIN,R4_B_MAX,R4_CHI_B_MAX,IERR)
 
      external EQ_GLIMB1
 
! argument declarations
 ! floating type, input/output:
      REAL R4_RHO
 
 ! floating type, input/output:
      REAL R4_PHI
 
 ! floating type, input/output:
      REAL R4_B_MIN
 
 ! floating type, input/output:
      REAL R4_B_MAX
 
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_CHI_B_MIN(1)
 
 ! floating type, input/output:
      REAL R4_CHI_B_MAX(1)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: RHO_act
      REAL*8 :: RHO_ref
      REAL*8 :: PHI_act
      REAL*8 :: PHI_ref
      REAL*8 :: B_MIN_act
      REAL*8 :: B_MIN_ref
      REAL*8, allocatable :: CHI_B_MIN_act(:)
      REAL*8, allocatable :: CHI_B_MIN_ref(:)
      REAL*8 :: B_MAX_act
      REAL*8 :: B_MAX_ref
      REAL*8, allocatable :: CHI_B_MAX_act(:)
      REAL*8, allocatable :: CHI_B_MAX_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(CHI_B_MIN_act(1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GLIMB1 -- CHI_B_MIN_act ALLOCATE error!')
      ALLOCATE(CHI_B_MIN_ref(1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GLIMB1 -- CHI_B_MIN_ref ALLOCATE error!')
      ALLOCATE(CHI_B_MAX_act(1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GLIMB1 -- CHI_B_MAX_act ALLOCATE error!')
      ALLOCATE(CHI_B_MAX_ref(1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GLIMB1 -- CHI_B_MAX_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RHO_act=R4_RHO
      RHO_ref=R4_RHO
 
      PHI_act=R4_PHI
      PHI_ref=R4_PHI
 
      B_MIN_act=R4_B_MIN
      B_MIN_ref=R4_B_MIN
 
      CHI_B_MIN_act=R4_CHI_B_MIN
      CHI_B_MIN_ref=R4_CHI_B_MIN
 
      B_MAX_act=R4_B_MAX
      B_MAX_ref=R4_B_MAX
 
      CHI_B_MAX_act=R4_CHI_B_MAX
      CHI_B_MAX_ref=R4_CHI_B_MAX
 
! call to original routine:  EQ_GLIMB1
 
      CALL EQ_GLIMB1(
     > RHO_act,PHI_act,B_MIN_act,CHI_B_MIN_act(1),B_MAX_act,
     > CHI_B_MAX_act(1),IERR)
 
! copy back outputs if modified.
 
      if(RHO_act.ne.RHO_ref) then
         R4_RHO = RHO_act
      endif
 
      if(PHI_act.ne.PHI_ref) then
         R4_PHI = PHI_act
      endif
 
      if(B_MIN_act.ne.B_MIN_ref) then
         R4_B_MIN = B_MIN_act
      endif
 
      where(CHI_B_MIN_act.ne.CHI_B_MIN_ref)
         R4_CHI_B_MIN = CHI_B_MIN_act
      end where
      DEALLOCATE(CHI_B_MIN_act)
      DEALLOCATE(CHI_B_MIN_ref)
 
      if(B_MAX_act.ne.B_MAX_ref) then
         R4_B_MAX = B_MAX_act
      endif
 
      where(CHI_B_MAX_act.ne.CHI_B_MAX_ref)
         R4_CHI_B_MAX = CHI_B_MAX_act
      end where
      DEALLOCATE(CHI_B_MAX_act)
      DEALLOCATE(CHI_B_MAX_ref)
 
! exit
      return
      end
