!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_CBDY_GRID
      SUBROUTINE R4_EQM_CBDY_GRID(
     > IPTS,R4_RLIM,R4_ZLIM,INUMR,INUMZ,ID_R,ID_Z,IERR)
 
      external EQM_CBDY_GRID
 
! argument declarations
      INTEGER IPTS
      INTEGER INUMR
      INTEGER INUMZ
      INTEGER ID_R
      INTEGER ID_Z
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_RLIM(IPTS)
 
 ! floating type, input/output:
      REAL R4_ZLIM(IPTS)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RLIM_act(:)
      REAL*8, allocatable :: RLIM_ref(:)
      REAL*8, allocatable :: ZLIM_act(:)
      REAL*8, allocatable :: ZLIM_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(RLIM_act(IPTS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CBDY_GRID -- RLIM_act ALLOCATE error!')
      ALLOCATE(RLIM_ref(IPTS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CBDY_GRID -- RLIM_ref ALLOCATE error!')
      ALLOCATE(ZLIM_act(IPTS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CBDY_GRID -- ZLIM_act ALLOCATE error!')
      ALLOCATE(ZLIM_ref(IPTS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CBDY_GRID -- ZLIM_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RLIM_act=R4_RLIM
      RLIM_ref=R4_RLIM
 
      ZLIM_act=R4_ZLIM
      ZLIM_ref=R4_ZLIM
 
! call to original routine:  EQM_CBDY_GRID
 
      CALL EQM_CBDY_GRID(
     > IPTS,RLIM_act(1),ZLIM_act(1),INUMR,INUMZ,ID_R,ID_Z,IERR)
 
! copy back outputs if modified.
 
      where(RLIM_act.ne.RLIM_ref)
         R4_RLIM = RLIM_act
      end where
      DEALLOCATE(RLIM_act)
      DEALLOCATE(RLIM_ref)
 
      where(ZLIM_act.ne.ZLIM_ref)
         R4_ZLIM = ZLIM_act
      end where
      DEALLOCATE(ZLIM_act)
      DEALLOCATE(ZLIM_ref)
 
! exit
      return
      end
