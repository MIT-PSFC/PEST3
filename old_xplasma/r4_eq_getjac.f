!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GETJAC
      SUBROUTINE R4_EQ_GETJAC(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,R4_JTENS,R4_DETJ,IERR)
 
      external EQ_GETJAC
 
! argument declarations
      INTEGER IVEC
      INTEGER IRZMODE
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_JTENS(3,3,abs(IVEC))
 
 ! floating type, output only:
      REAL R4_DETJ(abs(IVEC))
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: JTENS_act(:,:,:)
      REAL*8, allocatable :: DETJ_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETJAC -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETJAC -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETJAC -- ZPHI_act ALLOCATE error!')
      ALLOCATE(JTENS_act(3,3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETJAC -- JTENS_act ALLOCATE error!')
      ALLOCATE(DETJ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETJAC -- DETJ_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      JTENS_act=0
      DETJ_act=0
! call to original routine:  EQ_GETJAC
 
      CALL EQ_GETJAC(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,
     > JTENS_act(1,1,1),DETJ_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_JTENS = JTENS_act
      DEALLOCATE(JTENS_act)
      R4_DETJ = DETJ_act
      DEALLOCATE(DETJ_act)
 
! exit
      return
      end
