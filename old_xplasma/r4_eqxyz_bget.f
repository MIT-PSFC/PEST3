!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_BGET
      SUBROUTINE R4_EQXYZ_BGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,R4_BVEC,NREGION,IERR)
 
      external EQXYZ_BGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      BVEC_act=0
! call to original routine:  EQXYZ_BGET
 
      CALL EQXYZ_BGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,BVEC_act(1,1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_BGET
      SUBROUTINE R4_EQXYZMAP_BGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,R4_BVEC,
     > NREGION,IERR)
 
      external EQXYZMAP_BGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
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
      REAL R4_BVEC(3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
! call to original routine:  EQXYZMAP_BGET
 
      CALL EQXYZMAP_BGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,BVEC_act(1,1),NREGION,IERR)
 
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
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_BGRAD
      SUBROUTINE R4_EQXYZ_BGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,R4_BVEC,R4_GBTENSR,NREGION,IERR)
 
      external EQXYZ_BGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQXYZ_BGRAD
 
      CALL EQXYZ_BGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,BVEC_act(1,1),
     > GBTENSR_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_BGRAD
      SUBROUTINE R4_EQXYZMAP_BGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,R4_BVEC,
     > R4_GBTENSR,NREGION,IERR)
 
      external EQXYZMAP_BGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
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
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQXYZMAP_BGRAD
 
      CALL EQXYZMAP_BGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,BVEC_act(1,1),GBTENSR_act(1,1,1),NREGION,
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
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
