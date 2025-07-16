!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_BGET
      SUBROUTINE R4_EQRZ_BGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,R4_BVEC,NREGION,IERR)
 
      external EQRZ_BGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      BVEC_act=0
! call to original routine:  EQRZ_BGET
 
      CALL EQRZ_BGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,BVEC_act(1,1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_BGRAD
      SUBROUTINE R4_EQRZ_BGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,R4_BVEC,R4_GBTENSR,NREGION,IERR)
 
      external EQRZ_BGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQRZ_BGRAD
 
      CALL EQRZ_BGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,BVEC_act(1,1),
     > GBTENSR_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_BGET
      SUBROUTINE R4_EQRZMAP_BGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,R4_BVEC,NREGION,
     > IERR)
 
      external EQRZMAP_BGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
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
      REAL R4_BVEC(3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
! call to original routine:  EQRZMAP_BGET
 
      CALL EQRZMAP_BGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,BVEC_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_BGRAD
      SUBROUTINE R4_EQRZMAP_BGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,R4_BVEC,
     > R4_GBTENSR,NREGION,IERR)
 
      external EQRZMAP_BGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
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
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQRZMAP_BGRAD
 
      CALL EQRZMAP_BGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,BVEC_act(1,1),GBTENSR_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
