!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_BBIN
      SUBROUTINE R4_EQM_BBIN(
     > R4_ZBB,IAUTO,INUM,R4_ZTOL,ID,IERR)
 
      external EQM_BBIN
 
! argument declarations
      INTEGER IAUTO
      INTEGER INUM
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER ID
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZBB(INUM)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZBB_act(:)
      REAL*8, allocatable :: ZBB_ref(:)
      REAL*8 :: ZTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(ZBB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_BBIN -- ZBB_act ALLOCATE error!')
      ALLOCATE(ZBB_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_BBIN -- ZBB_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZBB_act=R4_ZBB
      ZBB_ref=R4_ZBB
 
      ZTOL_act=R4_ZTOL
! call to original routine:  EQM_BBIN
 
      CALL EQM_BBIN(
     > ZBB_act(1),IAUTO,INUM,ZTOL_act,ID,IERR)
 
! copy back outputs if modified.
 
      where(ZBB_act.ne.ZBB_ref)
         R4_ZBB = ZBB_act
      end where
      DEALLOCATE(ZBB_act)
      DEALLOCATE(ZBB_ref)
 
! exit
      return
      end
