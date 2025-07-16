!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_LIMCON
      SUBROUTINE R4_EQ_LIMCON(
     > IMAX,INUM,R4_RLIM,R4_ZLIM,R4_DTOL,IERR)
 
      external EQ_LIMCON
 
! argument declarations
      INTEGER IMAX
      INTEGER INUM
 ! floating type, input only:
      REAL R4_DTOL
 
      INTEGER IERR
 ! floating type, output only:
      REAL R4_RLIM(IMAX)
 
 ! floating type, output only:
      REAL R4_ZLIM(IMAX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RLIM_act(:)
      REAL*8, allocatable :: ZLIM_act(:)
      REAL*8 :: DTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(RLIM_act(IMAX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_LIMCON -- RLIM_act ALLOCATE error!')
      ALLOCATE(ZLIM_act(IMAX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_LIMCON -- ZLIM_act ALLOCATE error!')
 
! executable code:  copy for input
 
      RLIM_act=0
      ZLIM_act=0
      DTOL_act=R4_DTOL
! call to original routine:  EQ_LIMCON
 
      CALL EQ_LIMCON(
     > IMAX,INUM,RLIM_act(1),ZLIM_act(1),DTOL_act,IERR)
 
! copy back outputs if modified.
      R4_RLIM = RLIM_act
      DEALLOCATE(RLIM_act)
      R4_ZLIM = ZLIM_act
      DEALLOCATE(ZLIM_act)
 
! exit
      return
      end
