!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GAXIS
      SUBROUTINE R4_EQ_GAXIS(
     > IVEC,R4_PHI,R4_B_AXIS,R4_R_AXIS,R4_Z_AXIS)
 
      external EQ_GAXIS
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_PHI(IVEC)
 
 ! floating type, output only:
      REAL R4_B_AXIS(IVEC)
 
 ! floating type, output only:
      REAL R4_R_AXIS(IVEC)
 
 ! floating type, output only:
      REAL R4_Z_AXIS(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: PHI_act(:)
      REAL*8, allocatable :: B_AXIS_act(:)
      REAL*8, allocatable :: R_AXIS_act(:)
      REAL*8, allocatable :: Z_AXIS_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(PHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GAXIS -- PHI_act ALLOCATE error!')
      ALLOCATE(B_AXIS_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GAXIS -- B_AXIS_act ALLOCATE error!')
      ALLOCATE(R_AXIS_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GAXIS -- R_AXIS_act ALLOCATE error!')
      ALLOCATE(Z_AXIS_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GAXIS -- Z_AXIS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      PHI_act=R4_PHI
      B_AXIS_act=0
      R_AXIS_act=0
      Z_AXIS_act=0
! call to original routine:  EQ_GAXIS
 
      CALL EQ_GAXIS(
     > IVEC,PHI_act(1),B_AXIS_act(1),R_AXIS_act(1),Z_AXIS_act(1))
 
! copy back outputs if modified.
      DEALLOCATE(PHI_act)
      R4_B_AXIS = B_AXIS_act
      DEALLOCATE(B_AXIS_act)
      R4_R_AXIS = R_AXIS_act
      DEALLOCATE(R_AXIS_act)
      R4_Z_AXIS = Z_AXIS_act
      DEALLOCATE(Z_AXIS_act)
 
! exit
      return
      end
