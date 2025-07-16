!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQX_INV
      SUBROUTINE R4_EQX_INV(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_TOL,R4_ZRHO,R4_ZCHI,R4_ZPHI,NREGION,
     > IERR)
 
      external EQX_INV
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_TOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: TOL_act
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZRHO_ref(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZCHI_ref(:)
      REAL*8, allocatable :: ZPHI_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZRHO_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZRHO_ref ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZCHI_ref ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_INV -- ZPHI_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      TOL_act=R4_TOL
      ZRHO_act=R4_ZRHO
      ZRHO_ref=R4_ZRHO
 
      ZCHI_act=R4_ZCHI
      ZCHI_ref=R4_ZCHI
 
      ZPHI_act=0
! call to original routine:  EQX_INV
 
      CALL EQX_INV(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),TOL_act,ZRHO_act(1),
     > ZCHI_act(1),ZPHI_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
 
      where(ZRHO_act.ne.ZRHO_ref)
         R4_ZRHO = ZRHO_act
      end where
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZRHO_ref)
 
      where(ZCHI_act.ne.ZCHI_ref)
         R4_ZCHI = ZCHI_act
      end where
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZCHI_ref)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_INV
      SUBROUTINE R4_EQ_INV(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_TOL,R4_ZRHO,R4_ZCHI,NREGION,IERR)
 
      external EQ_INV
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_TOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZCHI(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: TOL_act
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZRHO_ref(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZCHI_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZRHO_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZRHO_ref ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_INV -- ZCHI_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      TOL_act=R4_TOL
      ZRHO_act=R4_ZRHO
      ZRHO_ref=R4_ZRHO
 
      ZCHI_act=R4_ZCHI
      ZCHI_ref=R4_ZCHI
 
! call to original routine:  EQ_INV
 
      CALL EQ_INV(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),TOL_act,ZRHO_act(1),
     > ZCHI_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
 
      where(ZRHO_act.ne.ZRHO_ref)
         R4_ZRHO = ZRHO_act
      end where
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZRHO_ref)
 
      where(ZCHI_act.ne.ZCHI_ref)
         R4_ZCHI = ZCHI_act
      end where
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZCHI_ref)
 
! exit
      return
      end
