!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_CHI
      SUBROUTINE R4_EQM_CHI(
     > R4_ZCHI,IAUTO,INUM,R4_ZTOL,ID,IERR)
 
      external EQM_CHI
 
! argument declarations
      INTEGER IAUTO
      INTEGER INUM
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER ID
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZCHI(INUM)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZCHI_ref(:)
      REAL*8 :: ZTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(ZCHI_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CHI -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_CHI -- ZCHI_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZCHI_act=R4_ZCHI
      ZCHI_ref=R4_ZCHI
 
      ZTOL_act=R4_ZTOL
! call to original routine:  EQM_CHI
 
      CALL EQM_CHI(
     > ZCHI_act(1),IAUTO,INUM,ZTOL_act,ID,IERR)
 
! copy back outputs if modified.
 
      where(ZCHI_act.ne.ZCHI_ref)
         R4_ZCHI = ZCHI_act
      end where
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZCHI_ref)
 
! exit
      return
      end
