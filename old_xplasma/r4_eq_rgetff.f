!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RGETFF
      SUBROUTINE R4_EQ_RGETFF(
     > IVEC,R4_ZRHO,NLIST,IFCNS,IWANT,IVECD,R4_ZVALS,IERR)
 
      external EQ_RGETFF
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IWANT
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_ZVALS(IVECD,NLIST)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGETFF -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGETFF -- ZVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZVALS_act=0
! call to original routine:  EQ_RGETFF
 
      CALL EQ_RGETFF(
     > IVEC,ZRHO_act(1),NLIST,IFCNS,IWANT,IVECD,ZVALS_act(1,1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      R4_ZVALS = ZVALS_act
      DEALLOCATE(ZVALS_act)
 
! exit
      return
      end
