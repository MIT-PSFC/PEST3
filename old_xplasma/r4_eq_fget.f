!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FGET
      SUBROUTINE R4_EQ_FGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IFCN,R4_FVAL,NREGION,IERR)
 
      external EQ_FGET
 
! argument declarations
      INTEGER IVEC
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FVAL(abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      FVAL_act=0
! call to original routine:  EQ_FGET
 
      CALL EQ_FGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IFCN,FVAL_act(1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FGRAD
      SUBROUTINE R4_EQ_FGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,IFCN,R4_FVAL,R4_FGRAD,
     > NREGION,IERR)
 
      external EQ_FGRAD
 
! argument declarations
      INTEGER IVEC
      INTEGER IRZMODE
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FVAL(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FGRAD(3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQ_FGRAD
 
      CALL EQ_FGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,IFCN,
     > FVAL_act(1),FGRAD_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_FGET
      SUBROUTINE R4_EQMAP_FGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IFCN,R4_FVAL,NREGION,
     > IERR)
 
      external EQMAP_FGET
 
! argument declarations
      INTEGER IVEC
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FVAL(abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVAL_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      FVAL_act=0
! call to original routine:  EQMAP_FGET
 
      CALL EQMAP_FGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IFCN,FVAL_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_FGRAD
      SUBROUTINE R4_EQMAP_FGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,IFCN,R4_FVAL,
     > R4_FGRAD,NREGION,IERR)
 
      external EQMAP_FGRAD
 
! argument declarations
      INTEGER IVEC
      INTEGER IRZMODE
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FVAL(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_FGRAD(3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVAL_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQMAP_FGRAD
 
      CALL EQMAP_FGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,IFCN,FVAL_act(1),FGRAD_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
