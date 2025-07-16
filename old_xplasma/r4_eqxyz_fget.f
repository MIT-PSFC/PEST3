!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_FGET
      SUBROUTINE R4_EQXYZ_FGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,IFCN,R4_FVAL,NREGION,IERR)
 
      external EQXYZ_FGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      FVAL_act=0
! call to original routine:  EQXYZ_FGET
 
      CALL EQXYZ_FGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,IFCN,FVAL_act(1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_FGRAD
      SUBROUTINE R4_EQXYZ_FGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,IFCN,R4_FVAL,R4_FGRAD,NREGION,
     > IERR)
 
      external EQXYZ_FGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
 ! floating type, output only:
      REAL R4_FGRAD(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQXYZ_FGRAD
 
      CALL EQXYZ_FGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,IFCN,FVAL_act(1),
     > FGRAD_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_FGET
      SUBROUTINE R4_EQXYZMAP_FGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,IFCN,
     > R4_FVAL,NREGION,IERR)
 
      external EQXYZMAP_FGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGET -- FVAL_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      FVAL_act=0
! call to original routine:  EQXYZMAP_FGET
 
      CALL EQXYZMAP_FGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,IFCN,FVAL_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_FGRAD
      SUBROUTINE R4_EQXYZMAP_FGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,IFCN,
     > R4_FVAL,R4_FGRAD,NREGION,IERR)
 
      external EQXYZMAP_FGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IFCN
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_FVAL(IVEC)
 
 ! floating type, output only:
      REAL R4_FGRAD(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVAL_act(:)
      REAL*8, allocatable :: FGRAD_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- FVAL_act ALLOCATE error!')
      ALLOCATE(FGRAD_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FGRAD -- FGRAD_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      FVAL_act=0
      FGRAD_act=0
! call to original routine:  EQXYZMAP_FGRAD
 
      CALL EQXYZMAP_FGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,IFCN,FVAL_act(1),FGRAD_act(1,1),NREGION,
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
      R4_FVAL = FVAL_act
      DEALLOCATE(FVAL_act)
      R4_FGRAD = FGRAD_act
      DEALLOCATE(FGRAD_act)
 
! exit
      return
      end
