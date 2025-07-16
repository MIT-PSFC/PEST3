!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FLXINT_ARR2D
      SUBROUTINE R4_EQ_FLXINT_ARR2D(
     > IORDER,IRHOAX,R4_ZDATA,ID1,ID2,IBCRHO0,IDZ0,R4_ZBCRHO0,IBCRHO1,
     > IDZ1,R4_ZBCRHO1,IWANT,R4_RESULT,INUMCHI,INUMRHO,IERR)
 
      external EQ_FLXINT_ARR2D
 
! argument declarations
      INTEGER IORDER
      INTEGER IRHOAX
      INTEGER ID1
      INTEGER ID2
      INTEGER IBCRHO0
      INTEGER IDZ0
      INTEGER IBCRHO1
      INTEGER IDZ1
      INTEGER IWANT
      INTEGER INUMCHI
      INTEGER INUMRHO
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZDATA(ID1,ID2)
 
 ! floating type, input/output:
      REAL R4_ZBCRHO0(IDZ0)
 
 ! floating type, input/output:
      REAL R4_ZBCRHO1(IDZ1)
 
 ! floating type, input/output:
      REAL R4_RESULT(INUMCHI,INUMRHO)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZDATA_act(:,:)
      REAL*8, allocatable :: ZDATA_ref(:,:)
      REAL*8, allocatable :: ZBCRHO0_act(:)
      REAL*8, allocatable :: ZBCRHO0_ref(:)
      REAL*8, allocatable :: ZBCRHO1_act(:)
      REAL*8, allocatable :: ZBCRHO1_ref(:)
      REAL*8, allocatable :: RESULT_act(:,:)
      REAL*8, allocatable :: RESULT_ref(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZDATA_act(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZDATA_act ALLOCATE error!')
      ALLOCATE(ZDATA_ref(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZDATA_ref ALLOCATE error!')
      ALLOCATE(ZBCRHO0_act(IDZ0)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZBCRHO0_act ALLOCATE error!')
      ALLOCATE(ZBCRHO0_ref(IDZ0)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZBCRHO0_ref ALLOCATE error!')
      ALLOCATE(ZBCRHO1_act(IDZ1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZBCRHO1_act ALLOCATE error!')
      ALLOCATE(ZBCRHO1_ref(IDZ1)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- ZBCRHO1_ref ALLOCATE error!')
      ALLOCATE(RESULT_act(INUMCHI,INUMRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- RESULT_act ALLOCATE error!')
      ALLOCATE(RESULT_ref(INUMCHI,INUMRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_ARR2D -- RESULT_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZDATA_act=R4_ZDATA
      ZDATA_ref=R4_ZDATA
 
      ZBCRHO0_act=R4_ZBCRHO0
      ZBCRHO0_ref=R4_ZBCRHO0
 
      ZBCRHO1_act=R4_ZBCRHO1
      ZBCRHO1_ref=R4_ZBCRHO1
 
      RESULT_act=R4_RESULT
      RESULT_ref=R4_RESULT
 
! call to original routine:  EQ_FLXINT_ARR2D
 
      CALL EQ_FLXINT_ARR2D(
     > IORDER,IRHOAX,ZDATA_act(1,1),ID1,ID2,IBCRHO0,IDZ0,ZBCRHO0_act(1),
     > IBCRHO1,IDZ1,ZBCRHO1_act(1),IWANT,RESULT_act(1,1),INUMCHI,
     > INUMRHO,IERR)
 
! copy back outputs if modified.
 
      where(ZDATA_act.ne.ZDATA_ref)
         R4_ZDATA = ZDATA_act
      end where
      DEALLOCATE(ZDATA_act)
      DEALLOCATE(ZDATA_ref)
 
      where(ZBCRHO0_act.ne.ZBCRHO0_ref)
         R4_ZBCRHO0 = ZBCRHO0_act
      end where
      DEALLOCATE(ZBCRHO0_act)
      DEALLOCATE(ZBCRHO0_ref)
 
      where(ZBCRHO1_act.ne.ZBCRHO1_ref)
         R4_ZBCRHO1 = ZBCRHO1_act
      end where
      DEALLOCATE(ZBCRHO1_act)
      DEALLOCATE(ZBCRHO1_ref)
 
      where(RESULT_act.ne.RESULT_ref)
         R4_RESULT = RESULT_act
      end where
      DEALLOCATE(RESULT_act)
      DEALLOCATE(RESULT_ref)
 
! exit
      return
      end
