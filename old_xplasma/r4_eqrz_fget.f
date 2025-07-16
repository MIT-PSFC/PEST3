!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_FGET
      SUBROUTINE R4_EQRZ_FGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,IFCN,R4_FVAL,NREGION,IERR)
 
      external EQRZ_FGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      FVAL_act=0
! call to original routine:  EQRZ_FGET
 
      CALL EQRZ_FGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,IFCN,FVAL_act(1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_FGRAD
      SUBROUTINE R4_EQRZ_FGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,IFCN,R4_FVAL,R4_FGRAD,NREGION,
     > IERR)
 
      external EQRZ_FGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
 ! floating type, output only:
      REAL R4_FGRAD(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQRZ_FGRAD
 
      CALL EQRZ_FGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,IFCN,FVAL_act(1),
     > FGRAD_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_FGET
      SUBROUTINE R4_EQRZMAP_FGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,IFCN,R4_FVAL,
     > NREGION,IERR)
 
      external EQRZMAP_FGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      FVAL_act=0
! call to original routine:  EQRZMAP_FGET
 
      CALL EQRZMAP_FGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,IFCN,FVAL_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_FGRAD
      SUBROUTINE R4_EQRZMAP_FGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,IFCN,R4_FVAL,
     > R4_FGRAD,NREGION,IERR)
 
      external EQRZMAP_FGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
 ! floating type, output only:
      REAL R4_FGRAD(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQRZMAP_FGRAD
 
      CALL EQRZMAP_FGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,IFCN,FVAL_act(1),FGRAD_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
