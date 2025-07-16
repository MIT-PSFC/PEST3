!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_FXCEP
      SUBROUTINE R4_EQXYZ_FXCEP(
     > INUM,ITYPE,IDINIT,R4_XA,R4_YA,R4_ZA,R4_DA,R4_XB,R4_YB,R4_ZB,
     > R4_DB,R4_TOL,R4_XX,R4_YX,R4_ZX,ISTAT,IERR)
 
      external EQXYZ_FXCEP
 
! argument declarations
      INTEGER INUM
      INTEGER ITYPE
      INTEGER IDINIT
 ! floating type, input only:
      REAL R4_TOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XA(INUM)
 
 ! floating type, input only:
      REAL R4_YA(INUM)
 
 ! floating type, input only:
      REAL R4_ZA(INUM)
 
 ! floating type, input/output:
      REAL R4_DA(INUM)
 
 ! floating type, input only:
      REAL R4_XB(INUM)
 
 ! floating type, input only:
      REAL R4_YB(INUM)
 
 ! floating type, input only:
      REAL R4_ZB(INUM)
 
 ! floating type, input/output:
      REAL R4_DB(INUM)
 
 ! floating type, output only:
      REAL R4_XX(INUM)
 
 ! floating type, output only:
      REAL R4_YX(INUM)
 
 ! floating type, output only:
      REAL R4_ZX(INUM)
 
      INTEGER ISTAT(INUM)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XA_act(:)
      REAL*8, allocatable :: YA_act(:)
      REAL*8, allocatable :: ZA_act(:)
      REAL*8, allocatable :: DA_act(:)
      REAL*8, allocatable :: DA_ref(:)
      REAL*8, allocatable :: XB_act(:)
      REAL*8, allocatable :: YB_act(:)
      REAL*8, allocatable :: ZB_act(:)
      REAL*8, allocatable :: DB_act(:)
      REAL*8, allocatable :: DB_ref(:)
      REAL*8 :: TOL_act
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YX_act(:)
      REAL*8, allocatable :: ZX_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- XA_act ALLOCATE error!')
      ALLOCATE(YA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- YA_act ALLOCATE error!')
      ALLOCATE(ZA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- ZA_act ALLOCATE error!')
      ALLOCATE(DA_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- DA_act ALLOCATE error!')
      ALLOCATE(DA_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- DA_ref ALLOCATE error!')
      ALLOCATE(XB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- XB_act ALLOCATE error!')
      ALLOCATE(YB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- YB_act ALLOCATE error!')
      ALLOCATE(ZB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- ZB_act ALLOCATE error!')
      ALLOCATE(DB_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- DB_act ALLOCATE error!')
      ALLOCATE(DB_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- DB_ref ALLOCATE error!')
      ALLOCATE(XX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- XX_act ALLOCATE error!')
      ALLOCATE(YX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- YX_act ALLOCATE error!')
      ALLOCATE(ZX_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FXCEP -- ZX_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XA_act=R4_XA
      YA_act=R4_YA
      ZA_act=R4_ZA
      DA_act=R4_DA
      DA_ref=R4_DA
 
      XB_act=R4_XB
      YB_act=R4_YB
      ZB_act=R4_ZB
      DB_act=R4_DB
      DB_ref=R4_DB
 
      TOL_act=R4_TOL
      XX_act=0
      YX_act=0
      ZX_act=0
! call to original routine:  EQXYZ_FXCEP
 
      CALL EQXYZ_FXCEP(
     > INUM,ITYPE,IDINIT,XA_act(1),YA_act(1),ZA_act(1),DA_act(1),
     > XB_act(1),YB_act(1),ZB_act(1),DB_act(1),TOL_act,XX_act(1),
     > YX_act(1),ZX_act(1),ISTAT,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XA_act)
      DEALLOCATE(YA_act)
      DEALLOCATE(ZA_act)
 
      where(DA_act.ne.DA_ref)
         R4_DA = DA_act
      end where
      DEALLOCATE(DA_act)
      DEALLOCATE(DA_ref)
      DEALLOCATE(XB_act)
      DEALLOCATE(YB_act)
      DEALLOCATE(ZB_act)
 
      where(DB_act.ne.DB_ref)
         R4_DB = DB_act
      end where
      DEALLOCATE(DB_act)
      DEALLOCATE(DB_ref)
      R4_XX = XX_act
      DEALLOCATE(XX_act)
      R4_YX = YX_act
      DEALLOCATE(YX_act)
      R4_ZX = ZX_act
      DEALLOCATE(ZX_act)
 
! exit
      return
      end
