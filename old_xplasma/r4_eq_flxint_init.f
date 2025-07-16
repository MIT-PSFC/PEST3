!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FLXINT_INIT
      SUBROUTINE R4_EQ_FLXINT_INIT(
     > IAUTO,INUM,R4_ZRHO,R4_ZRHOMIN,IERR)
 
      external EQ_FLXINT_INIT
 
! argument declarations
      INTEGER IAUTO
      INTEGER INUM
 ! floating type, input/output:
      REAL R4_ZRHOMIN
 
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZRHO(INUM)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZRHO_ref(:)
      REAL*8 :: ZRHOMIN_act
      REAL*8 :: ZRHOMIN_ref
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_INIT -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZRHO_ref(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_INIT -- ZRHO_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZRHO_ref=R4_ZRHO
 
      ZRHOMIN_act=R4_ZRHOMIN
      ZRHOMIN_ref=R4_ZRHOMIN
 
! call to original routine:  EQ_FLXINT_INIT
 
      CALL EQ_FLXINT_INIT(
     > IAUTO,INUM,ZRHO_act(1),ZRHOMIN_act,IERR)
 
! copy back outputs if modified.
 
      where(ZRHO_act.ne.ZRHO_ref)
         R4_ZRHO = ZRHO_act
      end where
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZRHO_ref)
 
      if(ZRHOMIN_act.ne.ZRHOMIN_ref) then
         R4_ZRHOMIN = ZRHOMIN_act
      endif
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FLXINT_CHINIT
      SUBROUTINE R4_EQ_FLXINT_CHINIT(
     > IAUTO,INUMCHIS,R4_ZCHI,IERR)
 
      external EQ_FLXINT_CHINIT
 
! argument declarations
      INTEGER IAUTO
      INTEGER INUMCHIS
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZCHI(INUMCHIS)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZCHI_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZCHI_act(INUMCHIS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_CHINIT -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_ref(INUMCHIS)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FLXINT_CHINIT -- ZCHI_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZCHI_act=R4_ZCHI
      ZCHI_ref=R4_ZCHI
 
! call to original routine:  EQ_FLXINT_CHINIT
 
      CALL EQ_FLXINT_CHINIT(
     > IAUTO,INUMCHIS,ZCHI_act(1),IERR)
 
! copy back outputs if modified.
 
      where(ZCHI_act.ne.ZCHI_ref)
         R4_ZCHI = ZCHI_act
      end where
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZCHI_ref)
 
! exit
      return
      end
