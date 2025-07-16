!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_TBDY_GRID
      SUBROUTINE R4_EQM_TBDY_GRID(
     > ILINES,R4_RL,R4_ZL,R4_THL,ICIRCS,R4_RC,R4_ZC,R4_RAD,INUMR,INUMZ,
     > ID_R,ID_Z,IERR)
 
      external EQM_TBDY_GRID
 
! argument declarations
      INTEGER ILINES
      INTEGER ICIRCS
      INTEGER INUMR
      INTEGER INUMZ
      INTEGER ID_R
      INTEGER ID_Z
      INTEGER IERR
 ! floating type, input only:
      REAL R4_RL(ILINES)
 
 ! floating type, input only:
      REAL R4_ZL(ILINES)
 
 ! floating type, input only:
      REAL R4_THL(ILINES)
 
 ! floating type, input only:
      REAL R4_RC(ICIRCS)
 
 ! floating type, input only:
      REAL R4_ZC(ICIRCS)
 
 ! floating type, input only:
      REAL R4_RAD(ICIRCS)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RL_act(:)
      REAL*8, allocatable :: ZL_act(:)
      REAL*8, allocatable :: THL_act(:)
      REAL*8, allocatable :: RC_act(:)
      REAL*8, allocatable :: ZC_act(:)
      REAL*8, allocatable :: RAD_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(RL_act(ILINES)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- RL_act ALLOCATE error!')
      ALLOCATE(ZL_act(ILINES)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- ZL_act ALLOCATE error!')
      ALLOCATE(THL_act(ILINES)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- THL_act ALLOCATE error!')
      ALLOCATE(RC_act(ICIRCS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- RC_act ALLOCATE error!')
      ALLOCATE(ZC_act(ICIRCS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- ZC_act ALLOCATE error!')
      ALLOCATE(RAD_act(ICIRCS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_TBDY_GRID -- RAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      RL_act=R4_RL
      ZL_act=R4_ZL
      THL_act=R4_THL
      RC_act=R4_RC
      ZC_act=R4_ZC
      RAD_act=R4_RAD
! call to original routine:  EQM_TBDY_GRID
 
      CALL EQM_TBDY_GRID(
     > ILINES,RL_act(1),ZL_act(1),THL_act(1),ICIRCS,RC_act(1),ZC_act(1),
     > RAD_act(1),INUMR,INUMZ,ID_R,ID_Z,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(RL_act)
      DEALLOCATE(ZL_act)
      DEALLOCATE(THL_act)
      DEALLOCATE(RC_act)
      DEALLOCATE(ZC_act)
      DEALLOCATE(RAD_act)
 
! exit
      return
      end
