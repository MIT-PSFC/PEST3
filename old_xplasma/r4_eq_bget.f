!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BGET
      SUBROUTINE R4_EQ_BGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,R4_BVEC,NREGION,IERR)
 
      external EQ_BGET
 
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
      REAL R4_BVEC(3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      BVEC_act=0
! call to original routine:  EQ_BGET
 
      CALL EQ_BGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,BVEC_act(1,1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BGRAD
      SUBROUTINE R4_EQ_BGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,R4_BVEC,R4_GBTENSR,NREGION,
     > IERR)
 
      external EQ_BGRAD
 
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
      REAL R4_BVEC(3,abs(IVEC))
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQ_BGRAD
 
      CALL EQ_BGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,BVEC_act(1,1),
     > GBTENSR_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_BGET
      SUBROUTINE R4_EQMAP_BGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,R4_BVEC,NREGION,
     > IERR)
 
      external EQMAP_BGET
 
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
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_BVEC(3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: BVEC_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGET -- BVEC_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      BVEC_act=0
! call to original routine:  EQMAP_BGET
 
      CALL EQMAP_BGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,BVEC_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_BGRAD
      SUBROUTINE R4_EQMAP_BGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,R4_BVEC,
     > R4_GBTENSR,NREGION,IERR)
 
      external EQMAP_BGRAD
 
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
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_BVEC(3,abs(IVEC))
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,abs(IVEC))
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BGRAD -- GBTENSR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      BVEC_act=0
      GBTENSR_act=0
! call to original routine:  EQMAP_BGRAD
 
      CALL EQMAP_BGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,BVEC_act(1,1),GBTENSR_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
 
! exit
      return
      end
