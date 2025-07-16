!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_NINDX
      SUBROUTINE R4_EQ_NINDX(
     > IAXIS,IVEC,R4_ZVAL,KINDX,IWARN)
 
      external EQ_NINDX
 
! argument declarations
      INTEGER IAXIS
      INTEGER IVEC
      INTEGER IWARN
 ! floating type, input only:
      REAL R4_ZVAL(IVEC)
 
      INTEGER KINDX(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_NINDX -- ZVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZVAL_act=R4_ZVAL
! call to original routine:  EQ_NINDX
 
      CALL EQ_NINDX(
     > IAXIS,IVEC,ZVAL_act(1),KINDX,IWARN)
 
! copy back outputs if modified.
      DEALLOCATE(ZVAL_act)
 
! exit
      return
      end
