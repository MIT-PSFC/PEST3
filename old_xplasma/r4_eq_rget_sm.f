!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RGET_SM
      SUBROUTINE R4_EQ_RGET_SM(
     > INRHO,R4_ZRHO,ID,R4_ZDELTA,R4_ZVAL,IERR)
 
      external EQ_RGET_SM
 
! argument declarations
      INTEGER INRHO
      INTEGER ID
 ! floating type, input only:
      REAL R4_ZDELTA
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(INRHO)
 
 ! floating type, output only:
      REAL R4_ZVAL(INRHO)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8 :: ZDELTA_act
      REAL*8, allocatable :: ZVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(INRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGET_SM -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZVAL_act(INRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGET_SM -- ZVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZDELTA_act=R4_ZDELTA
      ZVAL_act=0
! call to original routine:  EQ_RGET_SM
 
      CALL EQ_RGET_SM(
     > INRHO,ZRHO_act(1),ID,ZDELTA_act,ZVAL_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      R4_ZVAL = ZVAL_act
      DEALLOCATE(ZVAL_act)
 
! exit
      return
      end
