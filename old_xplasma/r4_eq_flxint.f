!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FLXINT
      SUBROUTINE R4_EQ_FLXINT(
     > NAME,NOPTION,R4_RESULT,INUMCHI,INUMRHO,IERR)
 
      external EQ_FLXINT
 
! argument declarations
      CHARACTER*(*) NAME
      INTEGER NOPTION
      INTEGER INUMCHI
      INTEGER INUMRHO
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_RESULT(INUMCHI,INUMRHO)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RESULT_act(:,:)
      REAL*8, allocatable :: RESULT_ref(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(RESULT_act(INUMCHI,INUMRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT -- RESULT_act ALLOCATE error!')
      ALLOCATE(RESULT_ref(INUMCHI,INUMRHO)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT -- RESULT_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RESULT_act=R4_RESULT
      RESULT_ref=R4_RESULT
 
! call to original routine:  EQ_FLXINT
 
      CALL EQ_FLXINT(
     > NAME,NOPTION,RESULT_act(1,1),INUMCHI,INUMRHO,IERR)
 
! copy back outputs if modified.
 
      where(RESULT_act.ne.RESULT_ref)
         R4_RESULT = RESULT_act
      end where
      DEALLOCATE(RESULT_act)
      DEALLOCATE(RESULT_ref)
 
! exit
      return
      end
