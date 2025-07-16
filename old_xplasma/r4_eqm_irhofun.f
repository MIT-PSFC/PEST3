!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_IRHOFUN
      SUBROUTINE R4_EQM_IRHOFUN(
     > ID_AXIS,ZLBL,INPROF,R4_ZPROF,IFLAG,ID,IERR)
 
      external EQM_IRHOFUN
 
! argument declarations
      INTEGER ID_AXIS
      CHARACTER*(*) ZLBL
      INTEGER INPROF
      INTEGER IFLAG
      INTEGER ID
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZPROF(INPROF)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZPROF_act(:)
      REAL*8, allocatable :: ZPROF_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZPROF_act(INPROF)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_IRHOFUN -- ZPROF_act ALLOCATE error!')
      ALLOCATE(ZPROF_ref(INPROF)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_IRHOFUN -- ZPROF_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZPROF_act=R4_ZPROF
      ZPROF_ref=R4_ZPROF
 
! call to original routine:  EQM_IRHOFUN
 
      CALL EQM_IRHOFUN(
     > ID_AXIS,ZLBL,INPROF,ZPROF_act(1),IFLAG,ID,IERR)
 
! copy back outputs if modified.
 
      where(ZPROF_act.ne.ZPROF_ref)
         R4_ZPROF = ZPROF_act
      end where
      DEALLOCATE(ZPROF_act)
      DEALLOCATE(ZPROF_ref)
 
! exit
      return
      end
