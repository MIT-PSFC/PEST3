!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FXCEP
      SUBROUTINE R4_EQ_FXCEP(
     > INUM,ITYPE,IDINIT,R4_RA,R4_ZA,R4_PHIA,R4_DA,R4_RB,R4_ZB,R4_PHIB,
     > R4_DB,R4_TOL,R4_RX,R4_ZX,R4_PHIX,ISTAT,IERR)
 
      external EQ_FXCEP
 
! argument declarations
      INTEGER INUM
      INTEGER ITYPE
      INTEGER IDINIT
 ! floating type, input only:
      REAL R4_TOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_RA(INUM)
 
 ! floating type, input only:
      REAL R4_ZA(INUM)
 
 ! floating type, input only:
      REAL R4_PHIA(INUM)
 
 ! floating type, input/output:
      REAL R4_DA(INUM)
 
 ! floating type, input only:
      REAL R4_RB(INUM)
 
 ! floating type, input only:
      REAL R4_ZB(INUM)
 
 ! floating type, input only:
      REAL R4_PHIB(INUM)
 
 ! floating type, input/output:
      REAL R4_DB(INUM)
 
 ! floating type, output only:
      REAL R4_RX(INUM)
 
 ! floating type, output only:
      REAL R4_ZX(INUM)
 
 ! floating type, output only:
      REAL R4_PHIX(INUM)
 
      INTEGER ISTAT(INUM)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RA_act(:)
      REAL*8, allocatable :: ZA_act(:)
      REAL*8, allocatable :: PHIA_act(:)
      REAL*8, allocatable :: DA_act(:)
      REAL*8, allocatable :: DA_ref(:)
      REAL*8, allocatable :: RB_act(:)
      REAL*8, allocatable :: ZB_act(:)
      REAL*8, allocatable :: PHIB_act(:)
      REAL*8, allocatable :: DB_act(:)
      REAL*8, allocatable :: DB_ref(:)
      REAL*8 :: TOL_act
      REAL*8, allocatable :: RX_act(:)
      REAL*8, allocatable :: ZX_act(:)
      REAL*8, allocatable :: PHIX_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(RA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- RA_act ALLOCATE error!')
      ALLOCATE(ZA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- ZA_act ALLOCATE error!')
      ALLOCATE(PHIA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- PHIA_act ALLOCATE error!')
      ALLOCATE(DA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- DA_act ALLOCATE error!')
      ALLOCATE(DA_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- DA_ref ALLOCATE error!')
      ALLOCATE(RB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- RB_act ALLOCATE error!')
      ALLOCATE(ZB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- ZB_act ALLOCATE error!')
      ALLOCATE(PHIB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- PHIB_act ALLOCATE error!')
      ALLOCATE(DB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- DB_act ALLOCATE error!')
      ALLOCATE(DB_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- DB_ref ALLOCATE error!')
      ALLOCATE(RX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- RX_act ALLOCATE error!')
      ALLOCATE(ZX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- ZX_act ALLOCATE error!')
      ALLOCATE(PHIX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FXCEP -- PHIX_act ALLOCATE error!')
 
! executable code:  copy for input
 
      RA_act=R4_RA
      ZA_act=R4_ZA
      PHIA_act=R4_PHIA
      DA_act=R4_DA
      DA_ref=R4_DA
 
      RB_act=R4_RB
      ZB_act=R4_ZB
      PHIB_act=R4_PHIB
      DB_act=R4_DB
      DB_ref=R4_DB
 
      TOL_act=R4_TOL
      RX_act=0
      ZX_act=0
      PHIX_act=0
! call to original routine:  EQ_FXCEP
 
      CALL EQ_FXCEP(
     > INUM,ITYPE,IDINIT,RA_act(1),ZA_act(1),PHIA_act(1),DA_act(1),
     > RB_act(1),ZB_act(1),PHIB_act(1),DB_act(1),TOL_act,RX_act(1),
     > ZX_act(1),PHIX_act(1),ISTAT,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(RA_act)
      DEALLOCATE(ZA_act)
      DEALLOCATE(PHIA_act)
 
      where(DA_act.ne.DA_ref)
         R4_DA = DA_act
      end where
      DEALLOCATE(DA_act)
      DEALLOCATE(DA_ref)
      DEALLOCATE(RB_act)
      DEALLOCATE(ZB_act)
      DEALLOCATE(PHIB_act)
 
      where(DB_act.ne.DB_ref)
         R4_DB = DB_act
      end where
      DEALLOCATE(DB_act)
      DEALLOCATE(DB_ref)
      R4_RX = RX_act
      DEALLOCATE(RX_act)
      R4_ZX = ZX_act
      DEALLOCATE(ZX_act)
      R4_PHIX = PHIX_act
      DEALLOCATE(PHIX_act)
 
! exit
      return
      end
